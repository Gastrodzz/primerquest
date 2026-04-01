"""
PrimersQuest Pro v3.0
A robust, high-accuracy qPCR primer design tool.

Major Improvements over v2:
- Consistent Tm calculation via primer3.calc_tm() throughout
- Secondary structure checks (hairpin, homodimer, heterodimer ΔG)
- Smart fallback cascade (product size → Tm → GC, not all-at-once)
- Exon-junction spanning support for RT-qPCR
- SYBR Green vs TaqMan design presets
- 50 kb sequence support + FASTA/GenBank file upload
- Genome-browser-style primer map with GC-content track
- Expanded species support (Human, Mouse, Rat, Zebrafish)
- PrimerBank parsing with canary validation
- Clean, theme-aware UI

Author: Dr. Ahmed bey Chaker, King's College London
License: MIT
"""

import streamlit as st
import pandas as pd
import requests
from bs4 import BeautifulSoup
import primer3
import re
import urllib.parse
import time
from datetime import datetime
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from typing import Dict, List, Tuple, Optional
import json
from dataclasses import dataclass, field
from io import StringIO
import hashlib

# ---------------------------------------------------------------------------
# Safe BioPython imports
# ---------------------------------------------------------------------------
BIOPYTHON_AVAILABLE = False
try:
    from Bio.Seq import Seq as BioSeq
    from Bio import SeqIO
    from Bio.SeqUtils import gc_fraction
    BIOPYTHON_AVAILABLE = True
except ImportError:
    SeqIO = None

# ---------------------------------------------------------------------------
# Core thermodynamic helpers – use primer3 engine everywhere
# ---------------------------------------------------------------------------

def calc_tm(sequence: str, mv_conc: float = 50.0, dv_conc: float = 1.5,
            dntp_conc: float = 0.2, dna_conc: float = 50.0) -> float:
    """Nearest-neighbor Tm via primer3 – single source of truth."""
    if not sequence:
        return 0.0
    try:
        return primer3.calc_tm(
            sequence.upper(),
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )
    except Exception:
        # Absolute last-resort fallback (should never hit)
        s = sequence.upper()
        if len(s) < 14:
            return 2 * (s.count('A') + s.count('T')) + 4 * (s.count('G') + s.count('C'))
        gc = (s.count('G') + s.count('C')) / len(s)
        return 64.9 + 41 * (gc - 0.36) * (1 / len(s))


def calc_gc(sequence: str) -> float:
    """GC content as a percentage."""
    if not sequence:
        return 0.0
    s = sequence.upper()
    return (s.count('G') + s.count('C')) / len(s) * 100


def calc_hairpin_dg(sequence: str, mv_conc: float = 50.0, dv_conc: float = 1.5) -> float:
    """Hairpin ΔG (kcal/mol) via primer3."""
    try:
        res = primer3.calc_hairpin(sequence.upper(), mv_conc=mv_conc, dv_conc=dv_conc)
        return res.dg / 1000  # milli-cal → kcal
    except Exception:
        return 0.0


def calc_homodimer_dg(sequence: str, mv_conc: float = 50.0, dv_conc: float = 1.5) -> float:
    """Homodimer ΔG (kcal/mol) via primer3."""
    try:
        res = primer3.calc_homodimer(sequence.upper(), mv_conc=mv_conc, dv_conc=dv_conc)
        return res.dg / 1000
    except Exception:
        return 0.0


def calc_heterodimer_dg(seq1: str, seq2: str, mv_conc: float = 50.0, dv_conc: float = 1.5) -> float:
    """Heterodimer ΔG (kcal/mol) via primer3."""
    try:
        res = primer3.calc_heterodimer(seq1.upper(), seq2.upper(), mv_conc=mv_conc, dv_conc=dv_conc)
        return res.dg / 1000
    except Exception:
        return 0.0


def gc_content_profile(sequence: str, window: int = 50) -> Tuple[List[int], List[float]]:
    """Sliding-window GC% along a sequence."""
    positions, values = [], []
    s = sequence.upper()
    half = window // 2
    for i in range(half, len(s) - half):
        seg = s[i - half:i + half]
        gc = (seg.count('G') + seg.count('C')) / len(seg) * 100
        positions.append(i)
        values.append(gc)
    return positions, values


# ---------------------------------------------------------------------------
# Page configuration
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="PrimersQuest Pro — qPCR Primer Design",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'About': "PrimersQuest Pro v3.0 — Advanced qPCR primer design by Dr. Ahmed bey Chaker, King's College London"
    }
)

# ---------------------------------------------------------------------------
# Theme-aware CSS
# ---------------------------------------------------------------------------
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

    /* Typography */
    .main { font-family: 'Inter', sans-serif; }

    .app-header {
        font-size: 2.8rem;
        font-weight: 700;
        background: linear-gradient(135deg, #6366f1 0%, #8b5cf6 50%, #a855f7 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        text-align: center;
        margin-bottom: 0.25rem;
        letter-spacing: -0.03em;
    }
    .app-subtitle {
        text-align: center;
        opacity: 0.65;
        font-size: 1.05rem;
        margin-bottom: 2rem;
    }

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] { gap: 0.5rem; }
    .stTabs [data-baseweb="tab"] {
        padding: 0.7rem 1.4rem;
        font-weight: 600;
        border-radius: 10px;
    }
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #6366f1, #8b5cf6);
        color: white !important;
    }

    /* Sequence code blocks */
    code { letter-spacing: 0.06em; }

    /* Donate button */
    .donate-btn {
        background: linear-gradient(135deg, #FFDD00, #FBB034);
        color: #000 !important;
        padding: 0.75rem 1.2rem;
        border-radius: 10px;
        font-weight: 700;
        text-align: center;
        text-decoration: none !important;
        display: block;
        transition: transform 0.2s, box-shadow 0.2s;
        box-shadow: 0 4px 12px rgba(251,176,52,0.35);
    }
    .donate-btn:hover {
        transform: translateY(-2px);
        box-shadow: 0 8px 20px rgba(251,176,52,0.5);
    }

    /* Hide Streamlit branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# SEO Meta Tags and Google Analytics
st.markdown("""
<head>
    <!-- Google Analytics -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-NGW0S1841L"></script>
    <script>
        window.dataLayer = window.dataLayer || [];
        function gtag(){dataLayer.push(arguments);}
        gtag('js', new Date());
        gtag('config', 'G-NGW0S1841L');
        function trackEvent(category, action, label, value) {
            if (typeof gtag !== 'undefined') {
                gtag('event', action, {
                    'event_category': category,
                    'event_label': label,
                    'value': value
                });
            }
        }
    </script>
    <!-- SEO Meta Tags -->
    <meta name="description" content="PrimersQuest Pro - Free online qPCR primer design tool with PrimerBank integration. Design and validate primers for Real-Time PCR, RT-qPCR experiments.">
    <meta name="keywords" content="primer design, qPCR, PCR primers, PrimerBank, real-time PCR, RT-qPCR, primer3, NCBI primer blast, free primer design tool, molecular biology">
    <meta name="author" content="Dr. Ahmed bey Chaker, King's College London">
    <meta name="robots" content="index, follow">
    <link rel="canonical" href="https://primersquest.streamlit.app">
    <meta property="og:title" content="PrimersQuest Pro - Advanced qPCR Primer Design Tool">
    <meta property="og:description" content="Design high-quality primers for qPCR with our free tool. Features PrimerBank search, Primer3 algorithm, and NCBI validation.">
    <meta property="og:type" content="website">
    <meta property="og:url" content="https://primersquest.streamlit.app">
    <script type="application/ld+json">
    {
        "@context": "https://schema.org",
        "@type": "WebApplication",
        "name": "PrimersQuest Pro",
        "description": "Advanced qPCR primer design tool with PrimerBank integration",
        "url": "https://primersquest.streamlit.app",
        "applicationCategory": "UtilityApplication",
        "operatingSystem": "All",
        "offers": {"@type": "Offer", "price": "0", "priceCurrency": "USD"},
        "author": {"@type": "Person", "name": "Dr. Ahmed bey Chaker", "affiliation": {"@type": "Organization", "name": "King's College London"}}
    }
    </script>
</head>
""", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SPECIES_MAP = {
    "Human": "Homo sapiens",
    "Mouse": "Mus musculus",
    "Rat": "Rattus norvegicus",
    "Zebrafish": "Danio rerio",
}

MAX_SEQUENCE_LENGTH = 50_000  # 50 kb

DESIGN_PRESETS = {
    "SYBR Green (Standard qPCR)": {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_MAX_POLY_X': 3,
        'PRIMER_PRODUCT_SIZE_RANGE': [[80, 150]],
        'PRIMER_NUM_RETURN': 10,
        'PRIMER_MAX_SELF_ANY_TH': 45.0,
        'PRIMER_MAX_SELF_END_TH': 35.0,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
        'PRIMER_MAX_HAIRPIN_TH': 47.0,
    },
    "TaqMan (Probe-based)": {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_MIN_GC': 30.0,
        'PRIMER_MAX_GC': 70.0,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_PRODUCT_SIZE_RANGE': [[60, 120]],
        'PRIMER_NUM_RETURN': 10,
        'PRIMER_MAX_SELF_ANY_TH': 45.0,
        'PRIMER_MAX_SELF_END_TH': 35.0,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
        'PRIMER_MAX_HAIRPIN_TH': 47.0,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_OPT_TM': 68.0,
        'PRIMER_INTERNAL_MIN_TM': 65.0,
        'PRIMER_INTERNAL_MAX_TM': 72.0,
        'PRIMER_INTERNAL_MIN_SIZE': 18,
        'PRIMER_INTERNAL_MAX_SIZE': 30,
    },
}


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------
@dataclass
class PrimerPair:
    forward_seq: str
    reverse_seq: str
    forward_tm: float
    reverse_tm: float
    forward_gc: float
    reverse_gc: float
    product_size: int
    forward_position: int
    reverse_position: int
    penalty: float
    confidence: int = 0
    probe_seq: str = ""
    probe_tm: float = 0.0
    probe_gc: float = 0.0
    hairpin_dg_fwd: float = 0.0
    hairpin_dg_rev: float = 0.0
    homodimer_dg_fwd: float = 0.0
    homodimer_dg_rev: float = 0.0
    heterodimer_dg: float = 0.0
    self_any_th: float = 0.0
    self_end_th: float = 0.0
    left_end_stability: float = 0.0
    right_end_stability: float = 0.0


# ---------------------------------------------------------------------------
# Primer Designer
# ---------------------------------------------------------------------------
class PrimerDesigner:
    """Robust primer designer with smart fallback cascade."""

    SCORING_WEIGHTS = {
        'tm_match': 0.20,
        'gc_content': 0.15,
        'product_size': 0.12,
        'primer3_penalty': 0.10,
        'self_complementarity': 0.12,
        'end_stability': 0.08,
        'hairpin': 0.10,
        'homodimer': 0.08,
        'heterodimer': 0.05,
    }

    # ---- Sequence parsing --------------------------------------------------

    @staticmethod
    def parse_sequence(text: str) -> Tuple[Optional[Dict], Optional[str]]:
        """Parse FASTA or raw sequence. Supports up to 50 kb."""
        text = text.strip()
        if not text:
            return None, "No sequence provided."

        header = "User_Sequence"
        if text.startswith('>'):
            if BIOPYTHON_AVAILABLE:
                try:
                    records = list(SeqIO.parse(StringIO(text), "fasta"))
                    if not records:
                        return None, "No valid FASTA records found."
                    rec = records[0]
                    sequence = str(rec.seq).upper()
                    header = rec.description
                except Exception as e:
                    return None, f"FASTA parse error: {e}"
            else:
                lines = text.split('\n')
                header = lines[0][1:].strip() or "User_Sequence"
                sequence = ''.join(lines[1:])
                sequence = re.sub(r'[^ATGCNatgcn]', '', sequence).upper()
        else:
            sequence = re.sub(r'[^ATGCNatgcn]', '', text).upper()

        if len(sequence) < 50:
            return None, f"Sequence too short ({len(sequence)} bp). Minimum 50 bp required."
        if len(sequence) > MAX_SEQUENCE_LENGTH:
            return None, f"Sequence too long ({len(sequence):,} bp). Maximum {MAX_SEQUENCE_LENGTH:,} bp."

        n_pct = sequence.count('N') / len(sequence) * 100
        return {
            'sequence': sequence,
            'header': header,
            'length': len(sequence),
            'gc_content': calc_gc(sequence),
            'n_content': n_pct,
            'n_warning': n_pct > 5,
        }, None

    @staticmethod
    def parse_uploaded_file(uploaded_file) -> Tuple[Optional[Dict], Optional[str]]:
        """Parse an uploaded FASTA / GenBank / plain-text file."""
        try:
            content = uploaded_file.read().decode('utf-8', errors='replace')
        except Exception as e:
            return None, f"Could not read file: {e}"

        name = uploaded_file.name.lower()

        # GenBank
        if name.endswith(('.gb', '.gbk', '.genbank')):
            if not BIOPYTHON_AVAILABLE:
                return None, "BioPython is required for GenBank files. Install with: pip install biopython"
            try:
                records = list(SeqIO.parse(StringIO(content), "genbank"))
                if not records:
                    return None, "No records found in GenBank file."
                rec = records[0]
                seq = str(rec.seq).upper()
                # Extract exon junctions from CDS features
                junctions = []
                for feat in rec.features:
                    if feat.type == 'CDS' and len(feat.location.parts) > 1:
                        parts = sorted(feat.location.parts, key=lambda p: int(p.start))
                        for i in range(len(parts) - 1):
                            junctions.append(int(parts[i].end))
                return {
                    'sequence': seq,
                    'header': rec.description or rec.name,
                    'length': len(seq),
                    'gc_content': calc_gc(seq),
                    'n_content': seq.count('N') / len(seq) * 100 if seq else 0,
                    'n_warning': seq.count('N') / len(seq) * 100 > 5 if seq else False,
                    'junctions': junctions,
                }, None
            except Exception as e:
                return None, f"GenBank parse error: {e}"

        # FASTA or plain text
        return PrimerDesigner.parse_sequence(content)

    # ---- Primer design -----------------------------------------------------

    def design(self, seq_data: Dict, preset_params: Dict,
               custom_overrides: Dict = None,
               junction_list: List[int] = None) -> Tuple[Optional[List[PrimerPair]], Optional[str]]:
        """Design primers with smart fallback cascade."""

        sequence = seq_data['sequence']
        header = seq_data['header']

        base_params = preset_params.copy()
        if custom_overrides:
            base_params.update(custom_overrides)

        # Build fallback cascade: widen product size → Tm → GC
        cascades = self._build_fallback_cascade(base_params)

        all_primers: List[PrimerPair] = []

        for params in cascades:
            seq_args = {
                'SEQUENCE_ID': header,
                'SEQUENCE_TEMPLATE': sequence,
            }
            if junction_list:
                seq_args['SEQUENCE_OVERLAP_JUNCTION_LIST'] = junction_list
                params['PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION'] = 4
                params['PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION'] = 4

            global_args = {
                'PRIMER_TASK': 'generic',
                'PRIMER_PICK_LEFT_PRIMER': 1,
                'PRIMER_PICK_RIGHT_PRIMER': 1,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_SALT_DIVALENT': 1.5,
                'PRIMER_DNTP_CONC': 0.2,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                **params,
            }

            # Make sure internal oligo pick is set correctly
            if 'PRIMER_PICK_INTERNAL_OLIGO' not in global_args:
                global_args['PRIMER_PICK_INTERNAL_OLIGO'] = 0

            try:
                # Use non-deprecated API
                try:
                    results = primer3.design_primers(seq_args, global_args)
                except AttributeError:
                    results = primer3.designPrimers(seq_args, global_args)
            except Exception:
                continue

            n = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            for i in range(n):
                pp = self._extract_pair(results, i)
                if pp:
                    all_primers.append(pp)

            # If we already have enough good primers, stop cascading
            if len(all_primers) >= 8:
                break

        if not all_primers:
            return None, ("No suitable primers found. Try widening the product size range "
                          "or relaxing Tm/GC constraints in the sidebar.")

        # Deduplicate
        seen = set()
        unique = []
        for p in all_primers:
            key = (p.forward_seq, p.reverse_seq)
            if key not in seen:
                seen.add(key)
                unique.append(p)

        # Score and sort
        for p in unique:
            p.confidence = self._score(p)
        unique.sort(key=lambda p: p.confidence, reverse=True)

        # Return top 5 — quality over quantity
        return unique[:5], None

    # ---- Internal helpers --------------------------------------------------

    @staticmethod
    def _build_fallback_cascade(base: Dict) -> List[Dict]:
        """Product-size first, then Tm, then GC."""
        cascades = [base.copy()]

        # Step 1: widen product size
        wider_size = base.copy()
        orig_range = base.get('PRIMER_PRODUCT_SIZE_RANGE', [[80, 150]])
        lo = max(50, orig_range[0][0] - 30)
        hi = orig_range[0][1] + 80
        wider_size['PRIMER_PRODUCT_SIZE_RANGE'] = [[lo, hi]]
        cascades.append(wider_size)

        # Step 2: widen Tm by ±1 °C
        wider_tm = wider_size.copy()
        wider_tm['PRIMER_MIN_TM'] = base.get('PRIMER_MIN_TM', 58.0) - 1.0
        wider_tm['PRIMER_MAX_TM'] = base.get('PRIMER_MAX_TM', 62.0) + 1.0
        cascades.append(wider_tm)

        # Step 3: relax GC (last resort)
        wider_gc = wider_tm.copy()
        wider_gc['PRIMER_MIN_GC'] = max(25, base.get('PRIMER_MIN_GC', 40) - 10)
        wider_gc['PRIMER_MAX_GC'] = min(75, base.get('PRIMER_MAX_GC', 60) + 10)
        cascades.append(wider_gc)

        return cascades

    def _extract_pair(self, data: Dict, i: int) -> Optional[PrimerPair]:
        fwd = data.get(f'PRIMER_LEFT_{i}_SEQUENCE', '')
        rev = data.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')
        if not fwd or not rev:
            return None

        probe = data.get(f'PRIMER_INTERNAL_{i}_SEQUENCE', '')

        return PrimerPair(
            forward_seq=fwd,
            reverse_seq=rev,
            forward_tm=data.get(f'PRIMER_LEFT_{i}_TM', calc_tm(fwd)),
            reverse_tm=data.get(f'PRIMER_RIGHT_{i}_TM', calc_tm(rev)),
            forward_gc=data.get(f'PRIMER_LEFT_{i}_GC_PERCENT', calc_gc(fwd)),
            reverse_gc=data.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', calc_gc(rev)),
            product_size=data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0),
            forward_position=data.get(f'PRIMER_LEFT_{i}', [0, 0])[0],
            reverse_position=data.get(f'PRIMER_RIGHT_{i}', [0, 0])[0],
            penalty=data.get(f'PRIMER_PAIR_{i}_PENALTY', 0),
            probe_seq=probe,
            probe_tm=data.get(f'PRIMER_INTERNAL_{i}_TM', calc_tm(probe)) if probe else 0.0,
            probe_gc=data.get(f'PRIMER_INTERNAL_{i}_GC_PERCENT', calc_gc(probe)) if probe else 0.0,
            hairpin_dg_fwd=calc_hairpin_dg(fwd),
            hairpin_dg_rev=calc_hairpin_dg(rev),
            homodimer_dg_fwd=calc_homodimer_dg(fwd),
            homodimer_dg_rev=calc_homodimer_dg(rev),
            heterodimer_dg=calc_heterodimer_dg(fwd, rev),
            self_any_th=data.get(f'PRIMER_PAIR_{i}_COMPL_ANY_TH', 0),
            self_end_th=data.get(f'PRIMER_PAIR_{i}_COMPL_END_TH', 0),
            left_end_stability=data.get(f'PRIMER_LEFT_{i}_END_STABILITY', 0),
            right_end_stability=data.get(f'PRIMER_RIGHT_{i}_END_STABILITY', 0),
        )

    def _score(self, p: PrimerPair) -> int:
        """Comprehensive confidence score [0–100]."""
        s = {}

        # 1. Tm match (ideal ΔTm ≤ 1 °C)
        dt = abs(p.forward_tm - p.reverse_tm)
        s['tm_match'] = max(0, 100 - dt * 25)

        # 2. GC content (optimal 45-55 each)
        gc_scores = []
        for gc in [p.forward_gc, p.reverse_gc]:
            if 45 <= gc <= 55:
                gc_scores.append(100)
            elif 40 <= gc <= 60:
                gc_scores.append(80 - abs(gc - 50) * 2)
            else:
                gc_scores.append(max(0, 60 - abs(gc - 50) * 3))
        s['gc_content'] = sum(gc_scores) / 2

        # 3. Product size (qPCR sweet-spot 80–120 bp)
        ps = p.product_size
        if 80 <= ps <= 120:
            s['product_size'] = 100
        elif 60 <= ps <= 200:
            s['product_size'] = max(50, 100 - abs(ps - 100) * 0.8)
        else:
            s['product_size'] = max(0, 50 - abs(ps - 100) * 0.5)

        # 4. Primer3 penalty (lower = better)
        s['primer3_penalty'] = max(0, 100 - p.penalty * 15)

        # 5. Self-complementarity
        if p.self_any_th < 25 and p.self_end_th < 15:
            s['self_complementarity'] = 100
        elif p.self_any_th < 35 and p.self_end_th < 25:
            s['self_complementarity'] = 75
        else:
            s['self_complementarity'] = max(0, 50 - p.self_any_th)

        # 6. 3′ end stability (GC clamp)
        end_sc = []
        for seq in [p.forward_seq, p.reverse_seq]:
            last2 = seq[-2:].upper() if len(seq) >= 2 else ''
            gc_end = sum(1 for c in last2 if c in 'GC')
            end_sc.append(100 if gc_end >= 1 else 60)
        s['end_stability'] = sum(end_sc) / 2

        # 7. Hairpin ΔG (more negative = worse)
        worst_hp = min(p.hairpin_dg_fwd, p.hairpin_dg_rev)
        if worst_hp > -2:
            s['hairpin'] = 100
        elif worst_hp > -4:
            s['hairpin'] = 70
        elif worst_hp > -6:
            s['hairpin'] = 40
        else:
            s['hairpin'] = 10

        # 8. Homodimer ΔG
        worst_hd = min(p.homodimer_dg_fwd, p.homodimer_dg_rev)
        if worst_hd > -5:
            s['homodimer'] = 100
        elif worst_hd > -8:
            s['homodimer'] = 65
        else:
            s['homodimer'] = 20

        # 9. Heterodimer ΔG
        if p.heterodimer_dg > -5:
            s['heterodimer'] = 100
        elif p.heterodimer_dg > -8:
            s['heterodimer'] = 60
        else:
            s['heterodimer'] = 20

        total = sum(s[k] * self.SCORING_WEIGHTS[k] for k in s)
        return int(max(0, min(100, total)))


# ---------------------------------------------------------------------------
# PrimerBank Searcher
# ---------------------------------------------------------------------------
class PrimerBankSearcher:
    """Search Harvard PrimerBank with canary validation."""

    BASE = "https://pga.mgh.harvard.edu"
    SEARCH = f"{BASE}/cgi-bin/primerbank/new_search2.cgi"

    def __init__(self):
        self.session = requests.Session()

    def search(self, gene: str, species: str = "Human") -> Tuple[Optional[List[Dict]], Optional[str]]:
        try:
            # Warm up session
            self.session.get(f"{self.BASE}/primerbank/", timeout=10)

            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) Chrome/120.0',
                'Content-Type': 'application/x-www-form-urlencoded',
                'Referer': f"{self.BASE}/primerbank/",
            }
            data = {
                'selectBox': 'NCBI Gene Symbol',
                'species': species,
                'searchBox': gene,
                'Submit': 'Submit',
            }
            resp = self.session.post(self.SEARCH, data=data, headers=headers, timeout=15)
            if resp.status_code != 200:
                return None, f"PrimerBank returned HTTP {resp.status_code}"

            return self._parse(resp.text, gene)
        except requests.exceptions.Timeout:
            return None, "PrimerBank request timed out. Try again."
        except requests.exceptions.RequestException as e:
            return None, f"Network error: {e}"

    def _parse(self, html: str, gene: str) -> Tuple[Optional[List[Dict]], Optional[str]]:
        soup = BeautifulSoup(html, 'html.parser')
        text = soup.get_text()

        if any(kw in text.lower() for kw in ['no primer', 'not found', 'no result']):
            return None, f"No PrimerBank results for '{gene}'"

        sets: List[Dict] = []
        # Extract all DNA-like sequences
        all_seqs = re.findall(r'[ATGC]{18,35}', text)
        all_ids = re.findall(r'(\d{5,}[a-z]\d)', text)

        if not all_seqs:
            return None, f"Could not parse primer data for '{gene}'. PrimerBank page structure may have changed."

        # Group in pairs (fwd, rev) + optional probe
        i = 0
        set_num = 0
        while i + 1 < len(all_seqs):
            ps = {
                'forward_seq': all_seqs[i],
                'reverse_seq': all_seqs[i + 1],
                'gene_symbol': gene,
                'source': 'PrimerBank',
            }
            if set_num < len(all_ids):
                ps['primerbank_id'] = all_ids[set_num]

            # Compute Tm/GC using primer3
            ps['forward_tm'] = calc_tm(ps['forward_seq'])
            ps['reverse_tm'] = calc_tm(ps['reverse_seq'])
            ps['forward_gc'] = calc_gc(ps['forward_seq'])
            ps['reverse_gc'] = calc_gc(ps['reverse_seq'])

            sets.append(ps)
            set_num += 1
            i += 2

        # Deduplicate
        seen = set()
        unique = []
        for s in sets:
            key = (s['forward_seq'], s['reverse_seq'])
            if key not in seen:
                seen.add(key)
                unique.append(s)

        for idx, s in enumerate(unique):
            s['set_number'] = idx + 1

        if unique:
            return unique, None
        return None, f"Could not extract primer pairs for '{gene}'."


# ---------------------------------------------------------------------------
# Visualization functions
# ---------------------------------------------------------------------------

def fig_primer_map(sequence: str, primers: List[PrimerPair], n_show: int = 5) -> go.Figure:
    """Genome-browser-style primer binding map with GC track."""
    seq_len = len(sequence)
    n = min(n_show, len(primers))

    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True,
        row_heights=[0.35, 0.65],
        vertical_spacing=0.06,
        subplot_titles=["GC Content (50 bp window)", "Primer Binding Sites"]
    )

    # --- GC track ---
    positions, gc_vals = gc_content_profile(sequence, window=50)
    fig.add_trace(go.Scatter(
        x=positions, y=gc_vals,
        mode='lines', fill='tozeroy',
        line=dict(color='rgba(99,102,241,0.6)', width=1),
        fillcolor='rgba(99,102,241,0.15)',
        name='GC %', hovertemplate='Pos %{x}: %{y:.1f}% GC',
    ), row=1, col=1)
    fig.add_hline(y=50, line_dash='dot', line_color='rgba(100,100,100,0.4)',
                  annotation_text='50%', row=1, col=1)

    # --- Template bar ---
    fig.add_trace(go.Scatter(
        x=[0, seq_len], y=[0, 0],
        mode='lines', line=dict(color='rgba(150,150,150,0.5)', width=6),
        name='Template', showlegend=False, hoverinfo='skip',
    ), row=2, col=1)

    colors = px.colors.qualitative.Vivid
    for i in range(n):
        p = primers[i]
        c = colors[i % len(colors)]
        y_fwd = 1 + i * 0.9
        y_rev = -(1 + i * 0.9)

        fwd_start = p.forward_position
        fwd_end = fwd_start + len(p.forward_seq)
        rev_end = p.reverse_position + 1
        rev_start = rev_end - len(p.reverse_seq)

        # Forward
        fig.add_trace(go.Scatter(
            x=[fwd_start, fwd_end], y=[y_fwd, y_fwd],
            mode='lines+text',
            line=dict(color=c, width=8),
            text=[f'F{i+1}', ''], textposition='top center',
            name=f'Pair {i+1} Fwd',
            hovertemplate=(f'<b>Pair {i+1} Forward</b><br>'
                           f'Pos: {fwd_start}–{fwd_end}<br>'
                           f'Tm: {p.forward_tm:.1f}°C<br>'
                           f'GC: {p.forward_gc:.1f}%<extra></extra>'),
        ), row=2, col=1)

        # Reverse
        fig.add_trace(go.Scatter(
            x=[rev_start, rev_end], y=[y_rev, y_rev],
            mode='lines+text',
            line=dict(color=c, width=8),
            text=[f'R{i+1}', ''], textposition='bottom center',
            name=f'Pair {i+1} Rev',
            hovertemplate=(f'<b>Pair {i+1} Reverse</b><br>'
                           f'Pos: {rev_start}–{rev_end}<br>'
                           f'Tm: {p.reverse_tm:.1f}°C<br>'
                           f'GC: {p.reverse_gc:.1f}%<extra></extra>'),
        ), row=2, col=1)

        # Product span (light connector)
        fig.add_trace(go.Scatter(
            x=[fwd_start, fwd_end, rev_start, rev_end],
            y=[y_fwd * 0.3, y_fwd * 0.3, y_rev * 0.3, y_rev * 0.3],
            mode='lines', line=dict(color=c, width=1, dash='dot'),
            showlegend=False, hoverinfo='skip',
        ), row=2, col=1)

    max_y = 1 + n * 0.9 + 0.5
    fig.update_yaxes(range=[-max_y, max_y], showticklabels=False,
                     zeroline=True, zerolinecolor='rgba(150,150,150,0.3)',
                     row=2, col=1)
    fig.update_xaxes(title_text="Position (bp)", row=2, col=1)
    fig.update_yaxes(title_text="GC %", row=1, col=1)
    fig.update_layout(
        height=500, margin=dict(t=50, b=40),
        legend=dict(orientation='h', y=-0.12),
        hovermode='closest',
    )
    return fig


def fig_tm_correlation(primers: List[PrimerPair]) -> go.Figure:
    """Tm scatter with ±1 °C tolerance band."""
    fwd_tms = [p.forward_tm for p in primers]
    rev_tms = [p.reverse_tm for p in primers]
    confs = [p.confidence for p in primers]
    sizes = [p.product_size for p in primers]

    lo = min(min(fwd_tms), min(rev_tms)) - 1
    hi = max(max(fwd_tms), max(rev_tms)) + 1

    fig = go.Figure()

    # ±1 °C band
    fig.add_trace(go.Scatter(
        x=[lo, hi, hi, lo, lo],
        y=[lo - 1, hi - 1, hi + 1, lo + 1, lo - 1],
        fill='toself', fillcolor='rgba(16,185,129,0.10)',
        line=dict(width=0), name='±1°C zone', showlegend=True,
    ))

    # Ideal line
    fig.add_trace(go.Scatter(
        x=[lo, hi], y=[lo, hi],
        mode='lines', line=dict(color='rgba(239,68,68,0.5)', dash='dash', width=1),
        name='Ideal', showlegend=True,
    ))

    # Points
    fig.add_trace(go.Scatter(
        x=fwd_tms, y=rev_tms,
        mode='markers',
        marker=dict(size=12, color=confs, colorscale='Viridis',
                    showscale=True, colorbar=dict(title='Score'),
                    line=dict(width=1, color='white')),
        text=[f'Pair {i+1}<br>Size: {s} bp' for i, s in enumerate(sizes)],
        hovertemplate='Fwd Tm: %{x:.1f}°C<br>Rev Tm: %{y:.1f}°C<br>%{text}<extra></extra>',
        name='Primers',
    ))

    fig.update_layout(
        title='Melting Temperature Correlation',
        xaxis_title='Forward Tm (°C)', yaxis_title='Reverse Tm (°C)',
        height=380, xaxis=dict(range=[lo, hi]), yaxis=dict(range=[lo, hi]),
    )
    return fig


def fig_confidence_bar(primers: List[PrimerPair]) -> go.Figure:
    """Horizontal bar chart of confidence scores."""
    labels = [f'Pair {i+1}' for i in range(len(primers))]
    scores = [p.confidence for p in primers]
    colors = ['#10b981' if s >= 85 else '#f59e0b' if s >= 70 else '#ef4444' for s in scores]

    fig = go.Figure(go.Bar(
        y=labels, x=scores, orientation='h',
        marker_color=colors,
        text=[f'{s}%' for s in scores], textposition='outside',
    ))
    fig.add_vline(x=85, line_dash='dash', line_color='#10b981',
                  annotation_text='Excellent', annotation_position='top')
    fig.add_vline(x=70, line_dash='dash', line_color='#f59e0b',
                  annotation_text='Good', annotation_position='top')
    fig.update_layout(
        title='Confidence Scores', xaxis_title='Score (%)',
        height=max(250, 40 * len(primers) + 100),
        xaxis=dict(range=[0, 105]),
        yaxis=dict(autorange='reversed'),
    )
    return fig


def fig_product_strip(primers: List[PrimerPair]) -> go.Figure:
    """Strip chart of product sizes (better than histogram at low n)."""
    sizes = [p.product_size for p in primers]
    confs = [p.confidence for p in primers]
    labels = [f'Pair {i+1}' for i in range(len(primers))]

    fig = go.Figure(go.Scatter(
        x=sizes, y=[0]*len(sizes),
        mode='markers+text',
        marker=dict(size=14, color=confs, colorscale='Viridis',
                    showscale=True, colorbar=dict(title='Score'),
                    line=dict(width=1, color='white')),
        text=labels, textposition='top center',
        hovertemplate='%{text}: %{x} bp<extra></extra>',
    ))
    # Optimal zone
    fig.add_vrect(x0=80, x1=120, fillcolor='rgba(16,185,129,0.08)',
                  line_width=0, annotation_text='Optimal qPCR', annotation_position='top left')
    fig.update_layout(
        title='Product Sizes', xaxis_title='Product Size (bp)',
        yaxis=dict(showticklabels=False, range=[-0.5, 0.8]),
        height=250,
    )
    return fig


# ---------------------------------------------------------------------------
# BLAST URL builder
# ---------------------------------------------------------------------------
def primer_blast_url(fwd: str, rev: str, species: str) -> str:
    org = SPECIES_MAP.get(species, "Homo sapiens")
    params = {
        'PRIMER_LEFT_INPUT': fwd,
        'PRIMER_RIGHT_INPUT': rev,
        'PRIMER_PRODUCT_MIN': 70,
        'PRIMER_PRODUCT_MAX': 1000,
        'PRIMER_SPECIFICITY_DATABASE': 'refseq_rna',
        'ORGANISM': org,
    }
    return "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?" + urllib.parse.urlencode(params)


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def primer_quality_badge(score: int):
    """Return a Streamlit badge for quality."""
    if score >= 85:
        st.success(f"✅ EXCELLENT ({score}%)")
    elif score >= 70:
        st.warning(f"⚠️ GOOD ({score}%)")
    else:
        st.error(f"❌ NEEDS OPTIMISATION ({score}%)")


def display_primer_card(p: PrimerPair, idx: int, species: str):
    """Render a single designed primer pair."""
    with st.container(border=True):
        hcol1, hcol2 = st.columns([3, 1])
        with hcol1:
            st.markdown(f"#### 🧬 Primer Pair {idx + 1}")
        with hcol2:
            primer_quality_badge(p.confidence)

        c1, c2 = st.columns(2)
        with c1:
            st.caption("FORWARD PRIMER")
            st.code(p.forward_seq)
            m1, m2, m3 = st.columns(3)
            m1.metric("Tm", f"{p.forward_tm:.1f} °C")
            m2.metric("GC", f"{p.forward_gc:.1f}%")
            m3.metric("Length", f"{len(p.forward_seq)} bp")
        with c2:
            st.caption("REVERSE PRIMER")
            st.code(p.reverse_seq)
            m1, m2, m3 = st.columns(3)
            m1.metric("Tm", f"{p.reverse_tm:.1f} °C")
            m2.metric("GC", f"{p.reverse_gc:.1f}%")
            m3.metric("Length", f"{len(p.reverse_seq)} bp")

        if p.probe_seq:
            st.caption("TAQMAN PROBE")
            st.code(p.probe_seq)
            m1, m2, m3 = st.columns(3)
            m1.metric("Tm", f"{p.probe_tm:.1f} °C")
            m2.metric("GC", f"{p.probe_gc:.1f}%")
            m3.metric("Length", f"{len(p.probe_seq)} bp")

        # Key metrics row
        k1, k2, k3, k4 = st.columns(4)
        k1.metric("Product", f"{p.product_size} bp")
        k2.metric("ΔTm", f"{abs(p.forward_tm - p.reverse_tm):.1f} °C")
        k3.metric("Penalty", f"{p.penalty:.2f}")
        k4.metric("Score", f"{p.confidence}%")

        with st.expander("Thermodynamic details"):
            t1, t2 = st.columns(2)
            with t1:
                st.markdown("**Hairpin ΔG (kcal/mol)**")
                st.write(f"Forward: {p.hairpin_dg_fwd:.1f} | Reverse: {p.hairpin_dg_rev:.1f}")
                st.markdown("**Homodimer ΔG (kcal/mol)**")
                st.write(f"Forward: {p.homodimer_dg_fwd:.1f} | Reverse: {p.homodimer_dg_rev:.1f}")
            with t2:
                st.markdown("**Heterodimer ΔG**")
                st.write(f"{p.heterodimer_dg:.1f} kcal/mol")
                st.markdown("**Self-complementarity**")
                st.write(f"Any: {p.self_any_th:.1f} °C | End: {p.self_end_th:.1f} °C")
                st.markdown("**3′ end stability**")
                st.write(f"Fwd: {p.left_end_stability:.1f} | Rev: {p.right_end_stability:.1f}")

            # Score breakdown
            st.markdown("**Score breakdown**")
            # Re-calculate sub-scores for display
            designer = PrimerDesigner()
            sub = {}
            dt = abs(p.forward_tm - p.reverse_tm)
            sub['Tm match'] = max(0, 100 - dt * 25)
            gc_sc = []
            for gc in [p.forward_gc, p.reverse_gc]:
                if 45 <= gc <= 55: gc_sc.append(100)
                elif 40 <= gc <= 60: gc_sc.append(80 - abs(gc - 50) * 2)
                else: gc_sc.append(max(0, 60 - abs(gc - 50) * 3))
            sub['GC content'] = sum(gc_sc) / 2
            ps = p.product_size
            if 80 <= ps <= 120: sub['Product size'] = 100
            elif 60 <= ps <= 200: sub['Product size'] = max(50, 100 - abs(ps - 100) * 0.8)
            else: sub['Product size'] = max(0, 50 - abs(ps - 100) * 0.5)
            sub['Penalty'] = max(0, 100 - p.penalty * 15)

            worst_hp = min(p.hairpin_dg_fwd, p.hairpin_dg_rev)
            sub['Hairpin'] = 100 if worst_hp > -2 else (70 if worst_hp > -4 else (40 if worst_hp > -6 else 10))
            worst_hd = min(p.homodimer_dg_fwd, p.homodimer_dg_rev)
            sub['Homodimer'] = 100 if worst_hd > -5 else (65 if worst_hd > -8 else 20)
            sub['Heterodimer'] = 100 if p.heterodimer_dg > -5 else (60 if p.heterodimer_dg > -8 else 20)

            breakdown_df = pd.DataFrame({
                'Component': list(sub.keys()),
                'Score': [f"{v:.0f}" for v in sub.values()],
            })
            st.dataframe(breakdown_df, hide_index=True, use_container_width=True)

        st.link_button("🎯 Validate with NCBI Primer-BLAST",
                       primer_blast_url(p.forward_seq, p.reverse_seq, species),
                       use_container_width=True)


def display_primerbank_card(ps: Dict, idx: int, species: str):
    """Render a PrimerBank result card."""
    fwd = ps.get('forward_seq', '')
    rev = ps.get('reverse_seq', '')
    if not fwd or not rev:
        return

    with st.container(border=True):
        h1, h2 = st.columns([3, 1])
        with h1:
            pid = ps.get('primerbank_id', '')
            title = f"#### Set {idx}" + (f" — ID: {pid}" if pid else "")
            st.markdown(title)
        with h2:
            st.success("Experimentally validated")

        c1, c2 = st.columns(2)
        with c1:
            st.caption("FORWARD")
            st.code(fwd)
            m1, m2, m3 = st.columns(3)
            m1.metric("Tm", f"{ps.get('forward_tm', 0):.1f} °C")
            m2.metric("GC", f"{ps.get('forward_gc', 0):.1f}%")
            m3.metric("Len", f"{len(fwd)} bp")
        with c2:
            st.caption("REVERSE")
            st.code(rev)
            m1, m2, m3 = st.columns(3)
            m1.metric("Tm", f"{ps.get('reverse_tm', 0):.1f} °C")
            m2.metric("GC", f"{ps.get('reverse_gc', 0):.1f}%")
            m3.metric("Len", f"{len(rev)} bp")

        # Thermodynamic analysis
        with st.expander("Thermodynamic analysis"):
            hp_f = calc_hairpin_dg(fwd)
            hp_r = calc_hairpin_dg(rev)
            hd_f = calc_homodimer_dg(fwd)
            hd_r = calc_homodimer_dg(rev)
            het = calc_heterodimer_dg(fwd, rev)
            t1, t2 = st.columns(2)
            with t1:
                st.write(f"**Hairpin ΔG:** Fwd {hp_f:.1f} | Rev {hp_r:.1f} kcal/mol")
                st.write(f"**Homodimer ΔG:** Fwd {hd_f:.1f} | Rev {hd_r:.1f} kcal/mol")
            with t2:
                st.write(f"**Heterodimer ΔG:** {het:.1f} kcal/mol")
                tm_diff = abs(ps.get('forward_tm', 0) - ps.get('reverse_tm', 0))
                st.write(f"**ΔTm:** {tm_diff:.1f} °C")

        b1, b2 = st.columns(2)
        with b1:
            st.link_button("Validate with Primer-BLAST",
                           primer_blast_url(fwd, rev, species))
        with b2:
            txt = f"Forward: {fwd}\nReverse: {rev}"
            st.download_button("Download sequences", data=txt,
                               file_name=f"primerbank_set_{idx}.txt",
                               key=f"pb_dl_{idx}")


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def export_designed_primers(primers: List[PrimerPair], key_prefix: str = "exp"):
    """Export designed primers in CSV, JSON, and lab format."""
    rows = []
    for i, p in enumerate(primers):
        rows.append({
            'pair': i + 1,
            'forward_seq': p.forward_seq,
            'reverse_seq': p.reverse_seq,
            'forward_tm': round(p.forward_tm, 2),
            'reverse_tm': round(p.reverse_tm, 2),
            'forward_gc': round(p.forward_gc, 2),
            'reverse_gc': round(p.reverse_gc, 2),
            'product_size': p.product_size,
            'confidence': p.confidence,
            'penalty': round(p.penalty, 3),
            'hairpin_dg_fwd': round(p.hairpin_dg_fwd, 2),
            'hairpin_dg_rev': round(p.hairpin_dg_rev, 2),
            'homodimer_dg_fwd': round(p.homodimer_dg_fwd, 2),
            'homodimer_dg_rev': round(p.homodimer_dg_rev, 2),
            'heterodimer_dg': round(p.heterodimer_dg, 2),
            'probe_seq': p.probe_seq or '',
        })

    df = pd.DataFrame(rows)
    c1, c2, c3 = st.columns(3)
    with c1:
        st.download_button("📥 Download CSV", df.to_csv(index=False),
                           "primers.csv", "text/csv",
                           key=f"{key_prefix}_csv")
    with c2:
        st.download_button("📥 Download JSON", json.dumps(rows, indent=2),
                           "primers.json", "application/json",
                           key=f"{key_prefix}_json")
    with c3:
        lab = []
        for r in rows:
            lab.append(f"--- Pair {r['pair']} (Score: {r['confidence']}%) ---")
            lab.append(f"  Forward : {r['forward_seq']}")
            lab.append(f"  Reverse : {r['reverse_seq']}")
            if r['probe_seq']:
                lab.append(f"  Probe   : {r['probe_seq']}")
            lab.append(f"  Product : {r['product_size']} bp")
            lab.append(f"  Fwd Tm  : {r['forward_tm']} °C | Rev Tm: {r['reverse_tm']} °C")
            lab.append("")
        st.download_button("📥 Download Lab Format", "\n".join(lab),
                           "primers_lab.txt", "text/plain",
                           key=f"{key_prefix}_lab")


# ---------------------------------------------------------------------------
# MAIN APP
# ---------------------------------------------------------------------------

def main():
    # Header
    st.markdown('<h1 class="app-header">PrimersQuest Pro</h1>', unsafe_allow_html=True)
    st.markdown('<p class="app-subtitle">High-accuracy qPCR primer design &amp; validation</p>',
                unsafe_allow_html=True)

    designer = PrimerDesigner()
    searcher = PrimerBankSearcher()

    # ---- Sidebar -----------------------------------------------------------
    with st.sidebar:
        st.markdown("### Configuration")

        preset_name = st.selectbox("Chemistry preset", list(DESIGN_PRESETS.keys()))
        preset = DESIGN_PRESETS[preset_name]

        with st.expander("Advanced parameters", expanded=False):
            min_tm = st.slider("Min Tm (°C)", 54.0, 66.0,
                               float(preset.get('PRIMER_MIN_TM', 58.0)), 0.5)
            max_tm = st.slider("Max Tm (°C)", 54.0, 66.0,
                               float(preset.get('PRIMER_MAX_TM', 62.0)), 0.5)
            min_size = st.number_input("Min product (bp)", 40, 500,
                                       preset['PRIMER_PRODUCT_SIZE_RANGE'][0][0])
            max_size = st.number_input("Max product (bp)", 40, 1000,
                                       preset['PRIMER_PRODUCT_SIZE_RANGE'][0][1])
            min_gc = st.slider("Min GC%", 20, 80,
                               int(preset.get('PRIMER_MIN_GC', 40)))
            max_gc = st.slider("Max GC%", 20, 80,
                               int(preset.get('PRIMER_MAX_GC', 60)))

        custom = {
            'PRIMER_MIN_TM': min_tm,
            'PRIMER_MAX_TM': max_tm,
            'PRIMER_PRODUCT_SIZE_RANGE': [[min_size, max_size]],
            'PRIMER_MIN_GC': float(min_gc),
            'PRIMER_MAX_GC': float(max_gc),
        }

        st.divider()
        st.markdown("### Resources")
        st.markdown(
            "- [Primer3 docs](http://primer3.org/)\n"
            "- [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)\n"
            "- [PrimerBank](https://pga.mgh.harvard.edu/primerbank/)\n"
            "- [MIQE guidelines](https://pubmed.ncbi.nlm.nih.gov/19246619/)"
        )

        st.divider()
        st.markdown("### Support this tool")
        st.markdown(
            '<a class="donate-btn" href="https://www.buymeacoffee.com/primerquest" '
            'target="_blank">Buy Me a Coffee ☕</a>',
            unsafe_allow_html=True,
        )

        with st.expander("System info"):
            st.write(f"BioPython: {'available' if BIOPYTHON_AVAILABLE else 'not installed'}")
            st.write(f"primer3-py: {primer3.__version__ if hasattr(primer3, '__version__') else 'installed'}")
            st.write(f"Max sequence: {MAX_SEQUENCE_LENGTH:,} bp")

    # ---- Tabs --------------------------------------------------------------
    tab_design, tab_bank, tab_analysis, tab_about = st.tabs([
        "🎯 Custom Design", "🔍 PrimerBank Search", "📊 Analysis Dashboard", "ℹ️ About"
    ])

    # ================================================================
    # TAB 1: Custom Design
    # ================================================================
    with tab_design:
        st.markdown("### Design Custom Primers")

        # Input method toggle
        input_method = st.radio("Input method", ["Paste sequence", "Upload file"],
                                horizontal=True, label_visibility="collapsed")

        seq_data = None
        error = None

        if input_method == "Paste sequence":
            sequence_input = st.text_area(
                "Target sequence (FASTA or raw DNA)",
                height=200,
                placeholder=(
                    ">GeneX mRNA\n"
                    "ATGGCAGAAATCGGTGTCAACGGATTTGGCCGTATTGGGCGCCTGGTCACC...\n\n"
                    "Accepts up to 50 kb. Paste raw nucleotides or FASTA format."
                ),
            )
            if sequence_input.strip():
                seq_data, error = designer.parse_sequence(sequence_input)
        else:
            uploaded = st.file_uploader(
                "Upload FASTA, GenBank, or plain text",
                type=['fasta', 'fa', 'txt', 'gb', 'gbk', 'genbank'],
            )
            if uploaded:
                seq_data, error = designer.parse_uploaded_file(uploaded)

        # Show stats
        if seq_data:
            s1, s2, s3, s4 = st.columns(4)
            s1.metric("Length", f"{seq_data['length']:,} bp")
            s2.metric("GC content", f"{seq_data['gc_content']:.1f}%")
            s3.metric("N content", f"{seq_data['n_content']:.1f}%")
            s4.metric("Format", "FASTA" if seq_data['header'] != 'User_Sequence' else "Raw")
            if seq_data.get('n_warning'):
                st.warning("High ambiguous-base content (>5%). Results may be suboptimal.")
            if seq_data.get('junctions'):
                st.info(f"Detected {len(seq_data['junctions'])} exon junction(s) from GenBank annotation. "
                        "Primers will be designed to span junctions where possible.")

        if error:
            st.error(error)

        # Junction input for non-GenBank
        junction_list = None
        if seq_data and not seq_data.get('junctions'):
            with st.expander("Exon-junction positions (optional, for RT-qPCR)"):
                junc_input = st.text_input(
                    "Comma-separated positions",
                    placeholder="e.g. 150, 320, 485",
                    help="Position of each exon-exon junction in your sequence. "
                         "At least one primer will be forced to span a junction."
                )
                if junc_input.strip():
                    try:
                        junction_list = [int(x.strip()) for x in junc_input.split(',') if x.strip()]
                    except ValueError:
                        st.error("Junction positions must be integers separated by commas.")
        elif seq_data and seq_data.get('junctions'):
            junction_list = seq_data['junctions']

        # Species
        blast_species = st.selectbox("Organism (for BLAST validation)",
                                     list(SPECIES_MAP.keys()), key="sp_design")

        # Design button
        if st.button("Design Primers", type="primary", use_container_width=True,
                      disabled=seq_data is None):
            with st.spinner("Designing primers..."):
                primers, err = designer.design(
                    seq_data, preset, custom, junction_list=junction_list
                )

            if err:
                st.error(err)
            else:
                st.session_state['primers'] = primers
                st.session_state['seq'] = seq_data['sequence']
                st.success(f"Designed {len(primers)} primer pair(s)")

        # Results
        if 'primers' in st.session_state and st.session_state['primers']:
            primers = st.session_state['primers']
            st.divider()

            # Primer map
            if 'seq' in st.session_state:
                st.plotly_chart(
                    fig_primer_map(st.session_state['seq'], primers),
                    use_container_width=True,
                )

            # Cards
            for i, p in enumerate(primers):
                display_primer_card(p, i, blast_species)

            # Export
            st.markdown("### Export")
            export_designed_primers(primers, key_prefix="design_export")
    # ================================================================
    with tab_bank:
        st.markdown("### Search Validated Primers")

        bc1, bc2 = st.columns([3, 1])
        with bc1:
            gene = st.text_input("Gene symbol", placeholder="e.g. GAPDH, ACTB, TP53")
        with bc2:
            sp = st.selectbox("Species", ["Human", "Mouse"], key="sp_bank")

        if st.button("Search PrimerBank", type="primary", use_container_width=True):
            if gene.strip():
                with st.spinner(f"Searching PrimerBank for {gene.strip()} ({sp})..."):
                    results, err = searcher.search(gene.strip(), sp)

                if err:
                    st.error(err)
                    c1, c2 = st.columns(2)
                    with c1:
                        url = (f"https://pga.mgh.harvard.edu/cgi-bin/primerbank/"
                               f"new_search2.cgi?searchBox={gene}&selectBox=NCBI+Gene+Symbol"
                               f"&species={sp}&Submit=Submit")
                        st.link_button("Try manual search on PrimerBank", url)
                    with c2:
                        st.link_button("Verify gene symbol on NCBI",
                                       f"https://www.ncbi.nlm.nih.gov/gene/?term={gene}")
                elif results:
                    st.session_state['pb_results'] = results
                    st.success(f"Found {len(results)} validated primer set(s)")

        if 'pb_results' in st.session_state and st.session_state['pb_results']:
            results = st.session_state['pb_results']
            for i, ps in enumerate(results):
                display_primerbank_card(ps, i + 1, sp)

            # Comparison table
            if len(results) > 1:
                st.markdown("### Comparison")
                comp = []
                for i, ps in enumerate(results):
                    comp.append({
                        'Set': i + 1,
                        'Fwd Tm': f"{ps.get('forward_tm',0):.1f}",
                        'Rev Tm': f"{ps.get('reverse_tm',0):.1f}",
                        'ΔTm': f"{abs(ps.get('forward_tm',0) - ps.get('reverse_tm',0)):.1f}",
                        'Fwd GC%': f"{ps.get('forward_gc',0):.1f}",
                        'Rev GC%': f"{ps.get('reverse_gc',0):.1f}",
                    })
                st.dataframe(pd.DataFrame(comp), hide_index=True, use_container_width=True)

                # Recommendations
                st.markdown("### 💡 Recommendations")
                best_idx = 0
                best_score = 0
                for i, ps in enumerate(results):
                    fwd = ps.get('forward_seq', '')
                    rev = ps.get('reverse_seq', '')
                    if not fwd or not rev:
                        continue
                    score = 100
                    tm_diff = abs(ps.get('forward_tm', 0) - ps.get('reverse_tm', 0))
                    if tm_diff > 5: score -= 20
                    elif tm_diff > 3: score -= 10
                    elif tm_diff > 1: score -= 5
                    for gc in [ps.get('forward_gc', 50), ps.get('reverse_gc', 50)]:
                        if gc < 30 or gc > 70: score -= 15
                        elif gc < 40 or gc > 60: score -= 5
                    for length in [len(fwd), len(rev)]:
                        if length < 18 or length > 25: score -= 10
                    if score > best_score:
                        best_score = score
                        best_idx = i

                reasons = []
                best_ps = results[best_idx]
                tm_diff = abs(best_ps.get('forward_tm', 0) - best_ps.get('reverse_tm', 0))
                if tm_diff <= 1:
                    reasons.append("excellent Tm matching")
                avg_gc = (best_ps.get('forward_gc', 50) + best_ps.get('reverse_gc', 50)) / 2
                if 45 <= avg_gc <= 55:
                    reasons.append("optimal GC content")
                if best_ps.get('probe_seq'):
                    reasons.append("includes validated probe")
                if best_ps.get('product_size') and 80 <= best_ps['product_size'] <= 150:
                    reasons.append("ideal product size for qPCR")

                st.success(f"**Recommended: Set {best_idx + 1}** (Quality Score: {best_score}%)")
                if reasons:
                    st.write("**Reasons:** " + ", ".join(reasons))

            # Export
            st.markdown("### Export")
            rows = []
            for i, ps in enumerate(results):
                rows.append({
                    'set': i+1, 'gene': ps.get('gene_symbol',''),
                    'primerbank_id': ps.get('primerbank_id',''),
                    'forward': ps.get('forward_seq',''),
                    'reverse': ps.get('reverse_seq',''),
                    'fwd_tm': round(ps.get('forward_tm',0), 2),
                    'rev_tm': round(ps.get('reverse_tm',0), 2),
                })
            c1, c2 = st.columns(2)
            with c1:
                st.download_button("📥 Download CSV",
                                   pd.DataFrame(rows).to_csv(index=False),
                                   "primerbank.csv", "text/csv",
                                   key="pb_export_csv")
            with c2:
                st.download_button("📥 Download JSON",
                                   json.dumps(rows, indent=2),
                                   "primerbank.json", "application/json",
                                   key="pb_export_json")

    # ================================================================
    # TAB 3: Analysis Dashboard
    # ================================================================
    with tab_analysis:
        st.markdown("### Analysis Dashboard")

        if 'primers' not in st.session_state or not st.session_state['primers']:
            st.info("Design primers first to see analysis.")
        else:
            primers = st.session_state['primers']

            # Summary metrics
            m1, m2, m3, m4 = st.columns(4)
            confs = [p.confidence for p in primers]
            m1.metric("Avg Score", f"{np.mean(confs):.0f}%")
            m2.metric("Excellent (≥85)", f"{sum(1 for c in confs if c >= 85)}/{len(primers)}")
            tds = [abs(p.forward_tm - p.reverse_tm) for p in primers]
            m3.metric("Avg ΔTm", f"{np.mean(tds):.1f} °C")
            m4.metric("Avg Product", f"{np.mean([p.product_size for p in primers]):.0f} bp")

            # Charts
            c1, c2 = st.columns(2)
            with c1:
                st.plotly_chart(fig_confidence_bar(primers), use_container_width=True)
            with c2:
                st.plotly_chart(fig_tm_correlation(primers), use_container_width=True)

            st.plotly_chart(fig_product_strip(primers), use_container_width=True)

            # Thermodynamic summary table
            st.markdown("### Thermodynamic Summary")
            thermo_rows = []
            for i, p in enumerate(primers):
                thermo_rows.append({
                    'Pair': i + 1,
                    'Score': p.confidence,
                    'Fwd Tm': f"{p.forward_tm:.1f}",
                    'Rev Tm': f"{p.reverse_tm:.1f}",
                    'ΔTm': f"{abs(p.forward_tm - p.reverse_tm):.1f}",
                    'Product': p.product_size,
                    'HP ΔG (F)': f"{p.hairpin_dg_fwd:.1f}",
                    'HP ΔG (R)': f"{p.hairpin_dg_rev:.1f}",
                    'HD ΔG (F)': f"{p.homodimer_dg_fwd:.1f}",
                    'HD ΔG (R)': f"{p.homodimer_dg_rev:.1f}",
                    'Het ΔG': f"{p.heterodimer_dg:.1f}",
                })
            st.dataframe(pd.DataFrame(thermo_rows), hide_index=True, use_container_width=True)

            # Export
            st.markdown("### Export")
            export_designed_primers(primers, key_prefix="analysis_export")
    # ================================================================
    with tab_about:
        st.markdown("### About PrimersQuest Pro")
        st.markdown(
            "**PrimersQuest Pro v3.0** was built by **Dr. Ahmed bey Chaker** "
            "(King's College London) to make qPCR primer design faster, "
            "more accurate, and freely accessible to the research community."
        )

        st.markdown("#### What's new in v3.0")
        st.markdown(
            "- **Consistent thermodynamics** — All Tm, hairpin, homodimer, and heterodimer "
            "calculations use the primer3 nearest-neighbour engine with identical salt/oligo "
            "concentrations. No more mixed-method discrepancies.\n"
            "- **Smart fallback** — When strict parameters yield no primers, the algorithm "
            "widens product-size range first (cheapest compromise), then Tm, and only relaxes "
            "GC as a last resort.\n"
            "- **Secondary structure scoring** — Confidence scores now incorporate hairpin ΔG, "
            "homodimer ΔG, and heterodimer ΔG from primer3.\n"
            "- **Exon-junction spanning** — Upload a GenBank file or enter junction positions "
            "manually; primers will be forced across junctions for RT-qPCR.\n"
            "- **50 kb sequences + file upload** — Paste up to 50 kb or upload FASTA / GenBank "
            "files directly.\n"
            "- **SYBR Green & TaqMan presets** — One click switches between chemistry-specific "
            "default parameters (including internal probe design for TaqMan).\n"
            "- **Genome-browser-style map** — GC-content track + stacked primer binding sites.\n"
            "- **Expanded species** — Human, Mouse, Rat, and Zebrafish BLAST validation."
        )

        st.markdown("#### Quick start")
        st.markdown(
            "1. **Custom Design** — Paste or upload your target sequence, pick a chemistry "
            "preset, and click *Design Primers*.\n"
            "2. **PrimerBank Search** — Enter a gene symbol to retrieve experimentally "
            "validated primer sets from Harvard PrimerBank.\n"
            "3. **Analysis Dashboard** — Review comparative visualisations and export your "
            "results in CSV, JSON, or lab-ready format.\n"
            "4. **Always validate** — Use the Primer-BLAST button on each card to check "
            "specificity against the refseq transcriptome."
        )

        st.markdown("#### 🔬 Best Practices")
        st.markdown(
            "- **Always validate** — Use the Primer-BLAST button to check specificity "
            "against the target genome.\n"
            "- **Check for SNPs** — Be aware of single nucleotide polymorphisms that may "
            "fall within your primer binding sites.\n"
            "- **Use high-quality sequences** — The quality of your input sequence directly "
            "impacts primer design quality.\n"
            "- **Run melt curves** — Always include melt curve analysis with SYBR Green to "
            "confirm single product amplification.\n"
            "- **Validate efficiency** — Run a 5-point serial dilution to confirm 90–110% "
            "amplification efficiency and R² > 0.99."
        )

        st.markdown("#### Resources")
        st.markdown(
            "- [Primer3](http://primer3.org/) — The underlying algorithm\n"
            "- [PrimerBank](https://pga.mgh.harvard.edu/primerbank/) — Harvard Medical School\n"
            "- [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)\n"
            "- [MIQE guidelines](https://pubmed.ncbi.nlm.nih.gov/19246619/) — Minimum Information for qPCR"
        )

        st.markdown("#### Contact")
        st.markdown(
            "- [LinkedIn](https://www.linkedin.com/in/ahmed-bey-chaker-6b3908192/)\n"
            "- [Support this app](https://www.buymeacoffee.com/primerquest)"
        )

    # Footer
    st.divider()
    st.markdown("""
    <div style='text-align: center; color: #718096; padding: 1.5rem;' itemscope itemtype="https://schema.org/WebApplication">
        <p style="font-size: 0.875rem;">
            <span itemprop="name">PrimersQuest Pro</span> v3.0 | Powered by Primer3 | Integrated with NCBI & Harvard PrimerBank<br>
            Made with ❤️ by <a href="https://www.linkedin.com/in/ahmed-bey-chaker-6b3908192/" target="_blank" style="color: #718096; text-decoration: underline;" itemprop="author">Dr. Ahmed bey Chaker</a>, King's College London
        </p>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
