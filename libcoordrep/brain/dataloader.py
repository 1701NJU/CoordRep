"""
CoordRep

80% tmQM + 20% COD
    +
"""

import sys
import random
import warnings
import numpy as np
from pathlib import Path
from typing import List, Dict, Iterator, Optional
from dataclasses import dataclass

warnings.filterwarnings('ignore')
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import torch
from torch.utils.data import Dataset, DataLoader

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep import encode_molecule
from coordrep.io.tmqm_reader import TMQMReader
from coordrep.io.cif_reader import CIFReader
from brain.tokenizer import CoordRepTokenizer
from brain.masking import CoordinatedMasker, LigandShuffler
from brain.smart_truncation import SmartTruncator, TruncationConfig, TruncationMonitor


@dataclass
class DataConfig:
    """"""
    tmqm_dir: str = None
    cod_dir: str = None
    tmqm_ratio: float = 0.8
    max_length: int = 512
    batch_size: int = 32
    num_workers: int = 4
    shuffle_ligands: bool = True
    n_augments: int = 2


class CoordRepDataset(Dataset):
    """
    CoordRep Training data
    
    1. tmQM + COD
    2.
    3.
    """
    
    def __init__(
        self,
        config: DataConfig,
        tokenizer: CoordRepTokenizer,
        masker: CoordinatedMasker,
        split: str = "train"
    ):
        self.config = config
        self.tokenizer = tokenizer
        self.masker = masker
        self.shuffler = LigandShuffler()
        self.split = split
        
        # Load data
        self.samples = []
        self._load_data()
    
    def _load_data(self):
        """"""
        print(f"Loading {self.split} data...")
        
        if self.config.tmqm_dir:
            tmqm_samples = self._load_tmqm()
            n_tmqm = int(len(tmqm_samples) * self.config.tmqm_ratio)
            self.samples.extend(tmqm_samples[:n_tmqm])
            print(f"  Loaded {n_tmqm} tmQM samples")
        
        if self.config.cod_dir:
            cod_samples = self._load_cod()
            n_cod = int(len(cod_samples) * (1 - self.config.tmqm_ratio))
            self.samples.extend(cod_samples[:n_cod])
            print(f"  Loaded {n_cod} COD samples")
        
        print(f"Total: {len(self.samples)} samples")
    
    def _load_tmqm(self) -> List[str]:
        """ tmQM """
        reader = TMQMReader(self.config.tmqm_dir)
        reader.load()
        
        samples = []
        for mol in reader.iter_molecules(limit=10000):
            try:
                cc = encode_molecule(mol)
                if cc.metal.element != "?":
                    cc_can = cc.canonicalize()
                    string = cc_can.to_string()
                    if len(string) < self.config.max_length * 2:
                        samples.append(string)
            except:
                continue
        
        return samples
    
    def _load_cod(self) -> List[str]:
        """ COD """
        reader = CIFReader(self.config.cod_dir)
        reader.load()
        
        samples = []
        for mol in reader.iter_molecules(limit=5000):
            try:
                cc = encode_molecule(mol)
                if cc.metal.element != "?":
                    cc_can = cc.canonicalize()
                    string = cc_can.to_string()
                    if len(string) < self.config.max_length * 2:
                        samples.append(string)
            except:
                continue
        
        return samples
    
    def __len__(self) -> int:
        if self.config.shuffle_ligands:
            return len(self.samples) * (1 + self.config.n_augments)
        return len(self.samples)
    
    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        sample_idx = idx % len(self.samples)
        text = self.samples[sample_idx]
        
        if self.config.shuffle_ligands and idx >= len(self.samples):
            text = self.shuffler.shuffle(text)
        
        # Tokenize
        tokens = self.tokenizer.tokenize(text)
        
        if len(tokens) > self.config.max_length - 2: # [CLS] [SEP]
            truncator = SmartTruncator(TruncationConfig(
                max_length=self.config.max_length - 2,
                max_ligand_tokens=80
            ))
            tokens, trunc_stats = truncator.truncate(tokens)
        
        input_ids = self.tokenizer.encode(text, add_special=True)
        
        if len(input_ids) > self.config.max_length:
            input_ids = input_ids[:self.config.max_length]
            tokens = tokens[:self.config.max_length - 2]
        
        masked_ids, labels = self.masker.create_mlm_inputs(
            input_ids, 
            ["[CLS]"] + tokens + ["[SEP]"],
            self.tokenizer.mask_token_id
        )
        
        padding_length = self.config.max_length - len(masked_ids)
        attention_mask = [1] * len(masked_ids) + [0] * padding_length
        masked_ids = masked_ids + [self.tokenizer.pad_token_id] * padding_length
        labels = labels + [-100] * padding_length
        
        return {
            "input_ids": torch.tensor(masked_ids, dtype=torch.long),
            "attention_mask": torch.tensor(attention_mask, dtype=torch.long),
            "labels": torch.tensor(labels, dtype=torch.long),
        }


class CoordRepDataModule:
    """
    DataLoader
    """
    
    def __init__(
        self,
        config: DataConfig,
        tokenizer: CoordRepTokenizer = None,
        masker: CoordinatedMasker = None
    ):
        self.config = config
        self.tokenizer = tokenizer or CoordRepTokenizer()
        self.masker = masker or CoordinatedMasker()
    
    def train_dataloader(self) -> DataLoader:
        dataset = CoordRepDataset(
            self.config, self.tokenizer, self.masker, split="train"
        )
        return DataLoader(
            dataset,
            batch_size=self.config.batch_size,
            shuffle=True,
            num_workers=self.config.num_workers,
            pin_memory=True
        )
    
    def val_dataloader(self) -> DataLoader:
        dataset = CoordRepDataset(
            self.config, self.tokenizer, self.masker, split="val"
        )
        return DataLoader(
            dataset,
            batch_size=self.config.batch_size,
            shuffle=False,
            num_workers=self.config.num_workers
        )


if __name__ == "__main__":
    config = DataConfig(
        tmqm_dir="../tmQM-master/tmQM",
        batch_size=4,
        max_length=128
    )
    
    tokenizer = CoordRepTokenizer()
    masker = CoordinatedMasker()
    
    dataset = CoordRepDataset(config, tokenizer, masker, split="train")
    
    print(f"\nDataset size: {len(dataset)}")
    
    sample = dataset[0]
    print("\nSample:")
    print(f"  input_ids shape: {sample['input_ids'].shape}")
    print(f"  attention_mask shape: {sample['attention_mask'].shape}")
    print(f"  labels shape: {sample['labels'].shape}")
    print(f"  num masked: {(sample['labels'] != -100).sum()}")

