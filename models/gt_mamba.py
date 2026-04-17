# models/gt_mamba.py
import torch
import torch.nn as nn
from mambapy.mamba import Mamba, MambaConfig
from configs.config import CONFIG

class GraphMaskedAttention(nn.Module):
    def __init__(self, input_dim, num_heads=4, dropout=0.1):
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = input_dim // num_heads
        self.scale = self.head_dim ** -0.5
        self.q_proj = nn.Linear(input_dim, input_dim)
        self.k_proj = nn.Linear(input_dim, input_dim)
        self.v_proj = nn.Linear(input_dim, input_dim)
        self.out = nn.Linear(input_dim, input_dim)
        self.dropout = nn.Dropout(dropout)
        self.norm = nn.LayerNorm(input_dim)

    def forward(self, x, adj_mask):
        B, N, C = x.shape
        residual = x
        q = self.q_proj(x).view(B, N, self.num_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(B, N, self.num_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(B, N, self.num_heads, self.head_dim).transpose(1, 2)

        scores = (q @ k.transpose(-2, -1)) * self.scale
        if adj_mask is not None:
            scores = scores + adj_mask.unsqueeze(0).unsqueeze(0)

        attn = torch.softmax(scores, dim=-1)
        attn = self.dropout(attn)

        out = (attn @ v).transpose(1, 2).reshape(B, N, C)
        out = self.out(out)
        out = self.dropout(out)
        return self.norm(residual + out)

class GraphMambaRegressor(nn.Module):
    def __init__(self, num_nodes=CONFIG['NUM_NODES']):
        super().__init__()
        d_model = CONFIG['D_MODEL']

        mask = torch.full((num_nodes, num_nodes), float('-inf'))
        mask.fill_diagonal_(0.0)
        self.register_buffer('adj_mask', mask)

        self.node_embedding = nn.Linear(1, d_model)
        self.norm_in = nn.LayerNorm(d_model)

        self.gat = GraphMaskedAttention(d_model, num_heads=CONFIG['NUM_HEADS'], dropout=CONFIG['DROPOUT'])

        config = MambaConfig(
            d_model=d_model,
            n_layers=2,
            d_state=CONFIG['D_STATE'],
            d_conv=CONFIG['D_CONV'],
            pscan=True
        )
        self.mamba = Mamba(config)

        self.head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Dropout(0.1),
            nn.Linear(d_model, 32),
            nn.GELU(),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        x = x.unsqueeze(-1)
        x = self.node_embedding(x)
        x = self.norm_in(x)

        x = self.gat(x, self.adj_mask)
        x = self.mamba(x)

        x = x.mean(dim=1)
        return self.head(x).squeeze(-1)