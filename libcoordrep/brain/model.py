"""
CoordRep Mini Brain: RoBERTa-Small

- Layers: 6
- Hidden: 384
- Heads: 6
- FFN: 1536

1.
2. < 1 GPU-day
3.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass
class CoordRepModelConfig:
 """"""
    vocab_size: int = 10000
    hidden_size: int = 384
    num_hidden_layers: int = 6
    num_attention_heads: int = 6
    intermediate_size: int = 1536
    hidden_dropout_prob: float = 0.1
    attention_probs_dropout_prob: float = 0.1
    max_position_embeddings: int = 512
    layer_norm_eps: float = 1e-6
    pad_token_id: int = 2
    
    @classmethod
    def small(cls, max_length: int = 768):
 """RoBERTa-Small """
        return cls(max_position_embeddings=max_length)
    
    @classmethod
    def tiny(cls, max_length: int = 512):
 """"""
        return cls(
            hidden_size=256,
            num_hidden_layers=4,
            num_attention_heads=4,
            intermediate_size=1024,
            max_position_embeddings=max_length
        )


class MultiHeadAttention(nn.Module):
 """"""
    
    def __init__(self, config: CoordRepModelConfig):
        super().__init__()
        self.num_heads = config.num_attention_heads
        self.head_dim = config.hidden_size // config.num_attention_heads
        self.scale = self.head_dim ** -0.5
        
        self.q_proj = nn.Linear(config.hidden_size, config.hidden_size)
        self.k_proj = nn.Linear(config.hidden_size, config.hidden_size)
        self.v_proj = nn.Linear(config.hidden_size, config.hidden_size)
        self.out_proj = nn.Linear(config.hidden_size, config.hidden_size)
        
        self.dropout = nn.Dropout(config.attention_probs_dropout_prob)
    
    def forward(
        self,
        hidden_states: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        batch_size, seq_len, _ = hidden_states.shape
        
        q = self.q_proj(hidden_states)
        k = self.k_proj(hidden_states)
        v = self.v_proj(hidden_states)
        
        q = q.view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        k = k.view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        v = v.view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        
        attn_weights = torch.matmul(q, k.transpose(-2, -1)) * self.scale
        
        if attention_mask is not None:
            attn_weights = attn_weights + attention_mask
        
        attn_weights = F.softmax(attn_weights, dim=-1)
        attn_weights = self.dropout(attn_weights)
        
        attn_output = torch.matmul(attn_weights, v)
        
        attn_output = attn_output.transpose(1, 2).contiguous().view(
            batch_size, seq_len, -1
        )
        
        return self.out_proj(attn_output)


class TransformerBlock(nn.Module):
 """Transformer """
    
    def __init__(self, config: CoordRepModelConfig):
        super().__init__()
        
        self.attention = MultiHeadAttention(config)
        self.ln1 = nn.LayerNorm(config.hidden_size, eps=config.layer_norm_eps)
        
        self.ffn = nn.Sequential(
            nn.Linear(config.hidden_size, config.intermediate_size),
            nn.GELU(),
            nn.Linear(config.intermediate_size, config.hidden_size),
            nn.Dropout(config.hidden_dropout_prob)
        )
        self.ln2 = nn.LayerNorm(config.hidden_size, eps=config.layer_norm_eps)
        
        self.dropout = nn.Dropout(config.hidden_dropout_prob)
    
    def forward(
        self,
        hidden_states: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        # Self-attention
        attn_output = self.attention(hidden_states, attention_mask)
        hidden_states = self.ln1(hidden_states + self.dropout(attn_output))
        
        # FFN
        ffn_output = self.ffn(hidden_states)
        hidden_states = self.ln2(hidden_states + ffn_output)
        
        return hidden_states


class CoordRepEmbeddings(nn.Module):
 """"""
    
    def __init__(self, config: CoordRepModelConfig):
        super().__init__()
        
        self.word_embeddings = nn.Embedding(
            config.vocab_size, config.hidden_size, padding_idx=config.pad_token_id
        )
        self.position_embeddings = nn.Embedding(
            config.max_position_embeddings, config.hidden_size
        )
        
        self.layer_norm = nn.LayerNorm(config.hidden_size, eps=config.layer_norm_eps)
        self.dropout = nn.Dropout(config.hidden_dropout_prob)
        
        self.register_buffer(
            "position_ids",
            torch.arange(config.max_position_embeddings).expand((1, -1))
        )
    
    def forward(self, input_ids: torch.Tensor) -> torch.Tensor:
        seq_len = input_ids.size(1)
        position_ids = self.position_ids[:, :seq_len]
        
        word_embeds = self.word_embeddings(input_ids)
        position_embeds = self.position_embeddings(position_ids)
        
        embeddings = word_embeds + position_embeds
        embeddings = self.layer_norm(embeddings)
        embeddings = self.dropout(embeddings)
        
        return embeddings


class CoordRepEncoder(nn.Module):
 """CoordRep """
    
    def __init__(self, config: CoordRepModelConfig):
        super().__init__()
        
        self.config = config
        self.embeddings = CoordRepEmbeddings(config)
        
        self.layers = nn.ModuleList([
            TransformerBlock(config) for _ in range(config.num_hidden_layers)
        ])
    
    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        hidden_states = self.embeddings(input_ids)
        
        if attention_mask is not None:
            # [batch, seq] -> [batch, 1, 1, seq]
            extended_mask = attention_mask.unsqueeze(1).unsqueeze(2)
            extended_mask = (1.0 - extended_mask) * -10000.0
        else:
            extended_mask = None
        
        for layer in self.layers:
            hidden_states = layer(hidden_states, extended_mask)
        
        return hidden_states


class CoordRepForMLM(nn.Module):
    """
 CoordRep
    
 CoordRep
    """
    
    def __init__(self, config: CoordRepModelConfig):
        super().__init__()
        
        self.config = config
        self.encoder = CoordRepEncoder(config)
        
        self.lm_head = nn.Sequential(
            nn.Linear(config.hidden_size, config.hidden_size),
            nn.GELU(),
            nn.LayerNorm(config.hidden_size, eps=config.layer_norm_eps),
            nn.Linear(config.hidden_size, config.vocab_size)
        )
        
        self.lm_head[-1].weight = self.encoder.embeddings.word_embeddings.weight
    
    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
        labels: Optional[torch.Tensor] = None
    ) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        hidden_states = self.encoder(input_ids, attention_mask)
        
        # Predict
        logits = self.lm_head(hidden_states)
        
        loss = None
        if labels is not None:
            loss_fct = nn.CrossEntropyLoss(ignore_index=-100)
            loss = loss_fct(logits.view(-1, self.config.vocab_size), labels.view(-1))
        
        return logits, loss
    
    def get_embeddings(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
        pooling: str = "mean"
    ) -> torch.Tensor:
        """
        """
        hidden_states = self.encoder(input_ids, attention_mask)
        
        if pooling == "cls":
            return hidden_states[:, 0]  # [CLS] token
        elif pooling == "mean":
            if attention_mask is not None:
                mask = attention_mask.unsqueeze(-1).float()
                return (hidden_states * mask).sum(1) / mask.sum(1)
            return hidden_states.mean(1)
        else:
            raise ValueError(f"Unknown pooling: {pooling}")


class CoordRepForRegression(nn.Module):
    """
 CoordRep Predict
    
 +
    """
    
    def __init__(
        self,
        encoder: CoordRepEncoder,
        hidden_size: int = 384,
        freeze_encoder: bool = True
    ):
        super().__init__()
        
        self.encoder = encoder
        self.freeze_encoder = freeze_encoder
        
        if freeze_encoder:
            for param in self.encoder.parameters():
                param.requires_grad = False
        
        self.regressor = nn.Sequential(
            nn.Linear(hidden_size, hidden_size // 2),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_size // 2, 1)
        )
    
    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
        labels: Optional[torch.Tensor] = None
    ) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        if self.freeze_encoder:
            with torch.no_grad():
                hidden_states = self.encoder(input_ids, attention_mask)
        else:
            hidden_states = self.encoder(input_ids, attention_mask)
        
        if attention_mask is not None:
            mask = attention_mask.unsqueeze(-1).float()
            pooled = (hidden_states * mask).sum(1) / mask.sum(1)
        else:
            pooled = hidden_states.mean(1)
        
        # Predict
        predictions = self.regressor(pooled).squeeze(-1)
        
        loss = None
        if labels is not None:
            loss = F.mse_loss(predictions, labels.float())
        
        return predictions, loss


def count_parameters(model: nn.Module) -> int:
 """Statistics"""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


if __name__ == "__main__":
    config = CoordRepModelConfig.small()
    model = CoordRepForMLM(config)
    
    print("=== CoordRep Mini Brain ===")
    print(f"Config: {config}")
    print(f"Parameters: {count_parameters(model):,}")
    
    batch_size = 4
    seq_len = 128
    
    input_ids = torch.randint(0, config.vocab_size, (batch_size, seq_len))
    attention_mask = torch.ones(batch_size, seq_len)
    labels = input_ids.clone()
 labels[labels != 4] = -100 # label
    
    logits, loss = model(input_ids, attention_mask, labels)
    
    print(f"\nOutput shape: {logits.shape}")
    print(f"Loss: {loss.item():.4f}")
    
    embeddings = model.get_embeddings(input_ids, attention_mask)
    print(f"Embeddings shape: {embeddings.shape}")

