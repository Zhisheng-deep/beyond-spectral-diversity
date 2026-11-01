"""DK-CAT feature extraction and CatBoost species classification (Section 3.1.2 of the paper).

Stage 1: DK-CAT, a transformer backbone with self-attention paired with a
knowledge-guided cross-attention module. The parallel knowledge stream encodes
the spectral regions of pigment absorption in the visible (450-500 nm), the red
edge (700-750 nm) and the near-infrared (> 750 nm), and spectral descriptors
(derivatives and statistical moments).
Stage 2: CatBoost classifier, trained with classes weighted by the square root
of their inverse frequency, assigns each pixel to the species with the highest
predicted probability (Section 3.1.3).
"""

import json
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from catboost import CatBoostClassifier
from scipy.stats import kurtosis, skew


def knowledge_ids(wavelength_nm):
    """Spectral region of each band: 0 other bands, 1 visible pigment absorption
    (450-500 nm), 2 red edge (700-750 nm), 3 near-infrared (> 750 nm)."""
    wl = np.asarray(wavelength_nm, dtype=np.float32)
    ids = np.zeros(len(wl), dtype=np.int64)
    ids[(wl >= 450.0) & (wl <= 500.0)] = 1
    ids[(wl >= 700.0) & (wl <= 750.0)] = 2
    ids[wl > 750.0] = 3
    return torch.from_numpy(ids)


def brightness_normalize(spectra):
    """Divide each spectrum by its Euclidean norm (Section 2.3)."""
    values = np.asarray(spectra, dtype=np.float32)
    return (values / np.linalg.norm(values, axis=1, keepdims=True)).astype(np.float32)


def spectral_descriptors(spectra):
    """Standardized spectrum, its first and second derivatives, and the standard
    deviation, skewness and kurtosis of each spectrum."""
    mean = spectra.mean(axis=1, keepdims=True)
    sd = spectra.std(axis=1, keepdims=True)
    z = ((spectra - mean) / (sd + 1e-8)).astype(np.float32)
    d1 = np.gradient(z, axis=1).astype(np.float32)
    d2 = np.gradient(d1, axis=1).astype(np.float32)
    sk = np.nan_to_num(skew(z, axis=1, bias=False), nan=0.0).reshape(-1, 1)
    ku = np.nan_to_num(kurtosis(z, axis=1, bias=False), nan=0.0).reshape(-1, 1)
    return np.column_stack([z, d1, d2, sd, sk, ku]).astype(np.float32)


class DKCAT(nn.Module):
    """Domain Knowledge-infused Cross-Attention Transformer."""

    def __init__(self, num_bands, num_classes, knowledge_ids, embed_dim=32):
        super().__init__()
        self.register_buffer("knowledge_ids", knowledge_ids)
        self.spectral_embedding = nn.Linear(1, embed_dim)
        self.knowledge_embedding = nn.Embedding(4, embed_dim)
        self.cross_attn = nn.MultiheadAttention(embed_dim, num_heads=4, batch_first=True)
        self.norm1 = nn.LayerNorm(embed_dim)
        layer = nn.TransformerEncoderLayer(d_model=embed_dim, nhead=4, dim_feedforward=128,
                                           dropout=0.1, batch_first=True, norm_first=True)
        self.transformer_encoder = nn.TransformerEncoder(layer, num_layers=2, enable_nested_tensor=False)
        self.gap = nn.AdaptiveAvgPool1d(1)
        self.descriptor_stream = nn.Sequential(
            nn.Linear(num_bands * 3 + 3, 1024), nn.GELU(), nn.LayerNorm(1024),
            nn.Linear(1024, 1024), nn.GELU(), nn.LayerNorm(1024),
            nn.Linear(1024, embed_dim),
        )
        self.classifier = nn.Linear(embed_dim, num_classes)

    def forward(self, spectra, descriptors, extract_features=False):
        spectral = self.spectral_embedding(spectra.unsqueeze(-1))
        ids = self.knowledge_ids.unsqueeze(0).expand(spectra.shape[0], -1)
        knowledge = self.knowledge_embedding(ids)
        attended, _ = self.cross_attn(query=knowledge, key=spectral, value=spectral)
        hidden = self.transformer_encoder(self.norm1(spectral + attended))
        features = self.gap(hidden.transpose(1, 2)).squeeze(-1) + self.descriptor_stream(descriptors)
        return features if extract_features else self.classifier(features)


class SpeciesModel:
    """Trained DK-CAT encoder and CatBoost classifier read from `weights/`."""

    def __init__(self, weights_dir, device="cpu"):
        weights_dir = Path(weights_dir)
        meta = json.loads((weights_dir / "model.json").read_text(encoding="utf-8"))
        self.num_bands = int(meta["band_count"])
        self.wavelength_nm = np.asarray(meta["wavelength_nm"], dtype=np.float32)
        self.species = list(meta["species"])
        self.codes = np.asarray(meta["output_codes"], dtype=np.float32)
        self.batch_size = int(meta["batch_size"])
        self.device = torch.device(device)
        self.encoder = DKCAT(self.num_bands, len(self.species), knowledge_ids(self.wavelength_nm),
                             int(meta["embed_dim"])).to(self.device)
        state = torch.load(weights_dir / "encoder.pt", map_location=self.device, weights_only=True)
        self.encoder.load_state_dict(state)
        self.encoder.eval()
        self.classifier = CatBoostClassifier()
        self.classifier.load_model(str(weights_dir / "classifier.cbm"))

    def features(self, spectra):
        """DK-CAT features of brightness-normalized spectra (N x bands)."""
        bn = brightness_normalize(spectra)
        desc = spectral_descriptors(bn)
        out = []
        with torch.inference_mode():
            for i in range(0, len(bn), self.batch_size):
                x = torch.from_numpy(bn[i:i + self.batch_size]).to(self.device)
                d = torch.from_numpy(desc[i:i + self.batch_size]).to(self.device)
                out.append(self.encoder(x, d, extract_features=True).cpu().numpy())
        return np.vstack(out)

    def predict(self, spectra):
        """Species code of each spectrum (the species with the highest probability)."""
        prob = self.classifier.predict_proba(self.features(spectra))
        return self.codes[np.argmax(prob, axis=1)]
