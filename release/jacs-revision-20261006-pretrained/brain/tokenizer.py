"""
CoordRep Hybrid Tokenizer

防碎策略：
1. 结构化标签作为 Special Tokens（严禁拆分）
2. 数值分桶（CShM 值 → 离散 Token）
3. 配体 SMILES 使用 BPE 切分
"""

import re
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from collections import Counter
from dataclasses import dataclass


# ============================================================
# 特殊 Token 定义
# ============================================================

# 结构化标签 - 严禁拆分
SPECIAL_TOKENS = [
    # 分隔符
    "[SEP]", "[CLS]", "[PAD]", "[MASK]", "[UNK]",

    # 金属块
    "[Fe]", "[Co]", "[Ni]", "[Cu]", "[Zn]", "[Ru]", "[Rh]", "[Pd]", "[Pt]", "[Ir]",
    "[Sc]", "[Ti]", "[V]", "[Cr]", "[Mn]", "[Y]", "[Zr]", "[Nb]", "[Mo]", "[Tc]",
    "[Ag]", "[Cd]", "[La]", "[Hf]", "[Ta]", "[W]", "[Re]", "[Os]", "[Au]", "[Hg]",

    # 氧化态
    ";ox=+1", ";ox=+2", ";ox=+3", ";ox=+4", ";ox=+5", ";ox=+6",
    ";ox=0", ";ox=-1", ";ox=-2",

    # 行/d电子
    ";row=3", ";row=4", ";row=5",
    ";d=0", ";d=1", ";d=2", ";d=3", ";d=4", ";d=5", ";d=6", ";d=7", ";d=8", ";d=9", ";d=10",

    # 配位数
    ";CN=2", ";CN=3", ";CN=4", ";CN=5", ";CN=6", ";CN=7", ";CN=8",

    # 几何类型
    "<Shape:", ">",
    "Td", "SP", "Oh", "TP", "TBP", "SPY", "TPr", "L",

    # 数值分桶 (0.00-20.00, 步长 0.2, 共 100 个)
    *[f"V_{i:03d}" for i in range(100)],

    # 约束块
    "{trans:", "{cis:", "}",
    "--",

    # 配体引用
    *[f"L{i}" for i in range(1, 21)],
    *[f":N:{i}" for i in range(1, 10)],
    *[f":O:{i}" for i in range(1, 10)],
    *[f":S:{i}" for i in range(1, 10)],
    *[f":P:{i}" for i in range(1, 10)],
    *[f":Cl:{i}" for i in range(1, 10)],

    # 配体字典
    "|", "=",
]


@dataclass
class TokenizerConfig:
    """Tokenizer 配置"""
    vocab_size: int = 10000
    min_frequency: int = 2
    special_tokens: List[str] = None
    cshm_bins: int = 100
    cshm_max: float = 20.0

    def __post_init__(self):
        if self.special_tokens is None:
            self.special_tokens = SPECIAL_TOKENS


class CoordRepTokenizer:
    """
    CoordRep 混合 Tokenizer

    策略：
    1. 先用正则提取结构化整词
    2. 对 SMILES 部分使用字符级/BPE
    3. 数值自动分桶
    """

    def __init__(self, config: TokenizerConfig = None):
        self.config = config or TokenizerConfig()

        # 词表
        self.token2id: Dict[str, int] = {}
        self.id2token: Dict[int, str] = {}

        # 特殊 Token 先占位
        for i, tok in enumerate(self.config.special_tokens):
            self.token2id[tok] = i
            self.id2token[i] = tok

        self.next_id = len(self.config.special_tokens)

        # 正则模式 - 用于提取结构化 Token
        self._build_patterns()

    def _build_patterns(self):
        """构建正则表达式模式"""
        # 金属块: [Fe;ox=+2;row=4;d=6;CN=6]
        self.metal_pattern = re.compile(
            r'\[([A-Z][a-z]?)(;ox=[+-]?\d)?(;row=\d)?(;d=\d+)?(;CN=\d+)?\]'
        )

        # 形状块: <Shape:Oh=4.99,Td=12.3>
        self.shape_pattern = re.compile(
            r'<Shape:([A-Za-z]+)=([\d.]+)(?:,([A-Za-z]+)=([\d.]+))*>'
        )

        # 约束块: {trans:L1:N:1--L2:N:1}
        self.constraint_pattern = re.compile(
            r'\{(trans|cis):(L\d+:[A-Z][a-z]?:\d+)--(L\d+:[A-Z][a-z]?:\d+)\}'
        )

        # 配体引用: (L1@1,2)
        self.ligand_ref_pattern = re.compile(r'\(L(\d+)@[\d,]+\)')

        # 配体字典: |L1=SMILES|
        self.ligand_dict_pattern = re.compile(r'\|L(\d+)=([^|]+)')

        # 数值模式
        self.number_pattern = re.compile(r'(\d+\.\d+)')

    def _quantize_cshm(self, value: float) -> str:
        """将 CShM 值量化为离散 Token"""
        # 0.00-20.00 映射到 V_000-V_099
        bin_idx = min(int(value / self.config.cshm_max * self.config.cshm_bins),
                      self.config.cshm_bins - 1)
        return f"V_{bin_idx:03d}"

    def _dequantize_cshm(self, token: str) -> float:
        """将离散 Token 还原为 CShM 值"""
        if token.startswith("V_"):
            bin_idx = int(token[2:])
            return bin_idx * self.config.cshm_max / self.config.cshm_bins
        return 0.0

    def tokenize(self, text: str) -> List[str]:
        """
        将 CoordRep 字符串切分为 Token 列表

        策略：
        1. 先提取结构化整词
        2. 对数值进行分桶
        3. 对 SMILES 字符级切分
        """
        tokens = []

        # 预处理：数值分桶
        def replace_numbers(match):
            value = float(match.group(1))
            return self._quantize_cshm(value)

        text = self.number_pattern.sub(replace_numbers, text)

        # 分段处理
        parts = self._split_into_blocks(text)

        for block_type, content in parts:
            if block_type == "metal":
                tokens.extend(self._tokenize_metal(content))
            elif block_type == "shape":
                tokens.extend(self._tokenize_shape(content))
            elif block_type == "constraint":
                tokens.extend(self._tokenize_constraint(content))
            elif block_type == "ligand_dict":
                tokens.extend(self._tokenize_ligand_dict(content))
            else:
                # 默认字符级切分
                tokens.extend(list(content))

        return tokens

    def _split_into_blocks(self, text: str) -> List[Tuple[str, str]]:
        """将文本分割为不同类型的块"""
        blocks = []
        pos = 0

        while pos < len(text):
            # 尝试匹配各种模式
            if text[pos] == '[':
                # 金属块
                end = text.find(']', pos)
                if end != -1:
                    blocks.append(("metal", text[pos:end+1]))
                    pos = end + 1
                    continue

            if text[pos] == '<':
                # 形状块
                end = text.find('>', pos)
                if end != -1:
                    blocks.append(("shape", text[pos:end+1]))
                    pos = end + 1
                    continue

            if text[pos] == '{':
                # 约束块
                end = text.find('}', pos)
                if end != -1:
                    blocks.append(("constraint", text[pos:end+1]))
                    pos = end + 1
                    continue

            if text[pos] == '|':
                # 配体字典
                end = text.find('|', pos + 1)
                if end != -1:
                    blocks.append(("ligand_dict", text[pos:end+1]))
                    pos = end + 1
                    continue

            # 默认：单字符
            blocks.append(("char", text[pos]))
            pos += 1

        return blocks

    def _tokenize_metal(self, block: str) -> List[str]:
        """切分金属块"""
        tokens = []
        # [Fe;ox=+2;row=4;d=6;CN=6]
        inner = block[1:-1]  # 去掉 []
        parts = inner.split(';')

        tokens.append(f"[{parts[0]}]")
        for part in parts[1:]:
            tokens.append(f";{part}")

        return tokens

    def _tokenize_shape(self, block: str) -> List[str]:
        """切分形状块"""
        tokens = ["<Shape:"]
        # <Shape:Oh=V_025,Td=V_062>
        inner = block[7:-1]  # 去掉 <Shape: 和 >
        for item in inner.split(','):
            if '=' in item:
                geo, val = item.split('=')
                tokens.append(geo)
                tokens.append("=")
                tokens.append(val)
                tokens.append(",")
        if tokens[-1] == ",":
            tokens.pop()
        tokens.append(">")
        return tokens

    def _tokenize_constraint(self, block: str) -> List[str]:
        """切分约束块"""
        tokens = []
        # {trans:L1:N:1--L2:N:1}
        inner = block[1:-1]  # 去掉 {}
        if inner.startswith("trans:"):
            tokens.append("{trans:")
            rest = inner[6:]
        elif inner.startswith("cis:"):
            tokens.append("{cis:")
            rest = inner[4:]
        else:
            return list(block)

        parts = rest.split("--")
        if len(parts) == 2:
            tokens.append(parts[0])
            tokens.append("--")
            tokens.append(parts[1])

        tokens.append("}")
        return tokens

    def _tokenize_ligand_dict(self, block: str) -> List[str]:
        """切分配体字典"""
        tokens = ["|"]
        # |L1=SMILES|
        inner = block[1:-1] if block.endswith('|') else block[1:]

        if '=' in inner:
            lid, smiles = inner.split('=', 1)
            tokens.append(lid)
            tokens.append("=")
            # SMILES 字符级切分
            tokens.extend(list(smiles))

        tokens.append("|")
        return tokens

    def encode(self, text: str, add_special: bool = True) -> List[int]:
        """将文本编码为 ID 序列"""
        tokens = self.tokenize(text)

        if add_special:
            tokens = ["[CLS]"] + tokens + ["[SEP]"]

        ids = []
        for tok in tokens:
            if tok in self.token2id:
                ids.append(self.token2id[tok])
            else:
                # 未知 Token - 尝试添加到词表
                if self.next_id < self.config.vocab_size:
                    self.token2id[tok] = self.next_id
                    self.id2token[self.next_id] = tok
                    ids.append(self.next_id)
                    self.next_id += 1
                else:
                    ids.append(self.token2id["[UNK]"])

        return ids

    def decode(self, ids: List[int]) -> str:
        """将 ID 序列解码为文本"""
        tokens = [self.id2token.get(i, "[UNK]") for i in ids]
        # 移除特殊 Token
        tokens = [t for t in tokens if t not in ["[CLS]", "[SEP]", "[PAD]"]]
        return "".join(tokens)

    def save(self, path: str):
        """保存 Tokenizer"""
        data = {
            "token2id": self.token2id,
            "config": {
                "vocab_size": self.config.vocab_size,
                "cshm_bins": self.config.cshm_bins,
                "cshm_max": self.config.cshm_max,
            }
        }
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load(cls, path: str) -> 'CoordRepTokenizer':
        """加载 Tokenizer"""
        with open(path, 'r') as f:
            data = json.load(f)

        config = TokenizerConfig(**data["config"])
        tokenizer = cls(config)
        tokenizer.token2id = data["token2id"]
        tokenizer.id2token = {int(k): v for k, v in
                              {v: k for k, v in data["token2id"].items()}.items()}
        tokenizer.next_id = max(tokenizer.token2id.values()) + 1
        return tokenizer

    @property
    def vocab_size(self) -> int:
        return len(self.token2id)

    @property
    def pad_token_id(self) -> int:
        return self.token2id["[PAD]"]

    @property
    def mask_token_id(self) -> int:
        return self.token2id["[MASK]"]

    @property
    def cls_token_id(self) -> int:
        return self.token2id["[CLS]"]

    @property
    def sep_token_id(self) -> int:
        return self.token2id["[SEP]"]


if __name__ == "__main__":
    # 测试
    tokenizer = CoordRepTokenizer()

    test_input = "[Fe;ox=+2;row=4;d=6;CN=6]<Shape:Oh=4.99,Td=12.30>{trans:L1:N:1--L2:N:1}|L1=NCCN|"

    print("Input:", test_input)
    tokens = tokenizer.tokenize(test_input)
    print("Tokens:", tokens)

    ids = tokenizer.encode(test_input)
    print("IDs:", ids)

    decoded = tokenizer.decode(ids)
    print("Decoded:", decoded)
