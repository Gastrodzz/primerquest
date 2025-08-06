"""
PrimersQuest Pro v3.0
Enhanced with modern UI, dimer analysis, and targeted position features
Author: Dr. Ahmed bey Chaker, King's College London
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
from typing import Dict, List, Tuple, Optional, Set
import json
from dataclasses import dataclass, field
from io import StringIO
import hashlib
from collections import defaultdict

# Configure page
st.set_page_config(
    page_title="PrimersQuest Pro v3.0 - Advanced qPCR Primer Design",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Enhanced CSS for modern social media-like design
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap');
    
    /* Main Theme */
    .main {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
        background: linear-gradient(135deg, #667eea15 0%, #764ba215 100%);
    }
    
    /* Header Styles */
    .main-header {
        font-size: 3.5rem;
        font-weight: 800;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        text-align: center;
        margin-bottom: 1rem;
        letter-spacing: -0.03em;
    }
    
    .version-badge {
        display: inline-block;
        background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        color: white;
        padding: 0.25rem 0.75rem;
        border-radius: 20px;
        font-size: 0.85rem;
        font-weight: 600;
        margin-left: 1rem;
    }
    
    /* Modern Card Design */
    .primer-card {
        background: white;
        border-radius: 16px;
        padding: 1.5rem;
        margin-bottom: 1rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.12);
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        border: 1px solid rgba(0,0,0,0.06);
    }
    
    .primer-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 12px 24px rgba(0,0,0,0.1), 0 4px 8px rgba(0,0,0,0.08);
        border-color: #667eea;
    }
    
    /* Thread-like Container */
    .thread-container {
        position: relative;
        padding-left: 2rem;
    }
    
    .thread-container::before {
        content: '';
        position: absolute;
        left: 0.75rem;
        top: 2rem;
        bottom: 0;
        width: 2px;
        background: linear-gradient(180deg, #667eea 0%, #764ba2 100%);
        opacity: 0.3;
    }
    
    .thread-node {
        position: absolute;
        left: 0.25rem;
        width: 1rem;
        height: 1rem;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        border-radius: 50%;
        border: 3px solid white;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    
    /* Status Badges */
    .status-badge {
        display: inline-flex;
        align-items: center;
        padding: 0.375rem 0.75rem;
        border-radius: 20px;
        font-size: 0.875rem;
        font-weight: 600;
        letter-spacing: 0.025em;
    }
    
    .status-excellent {
        background: linear-gradient(135deg, #10b981 0%, #059669 100%);
        color: white;
    }
    
    .status-good {
        background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%);
        color: white;
    }
    
    .status-warning {
        background: linear-gradient(135deg, #ef4444 0%, #dc2626 100%);
        color: white;
    }
    
    /* Sequence Display */
    .sequence-display {
        font-family: 'SF Mono', 'Monaco', 'Courier New', monospace;
        background: linear-gradient(135deg, #f3f4f6 0%, #e5e7eb 100%);
        padding: 1rem;
        border-radius: 8px;
        font-size: 0.95rem;
        letter-spacing: 0.05em;
        word-break: break-all;
        border: 1px solid #d1d5db;
    }
    
    /* Dimer Analysis Styles */
    .dimer-warning {
        background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%);
        border-left: 4px solid #f59e0b;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
    }
    
    .dimer-critical {
        background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%);
        border-left: 4px solid #ef4444;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
    }
    
    /* Metrics Grid */
    .metrics-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
        gap: 1rem;
        margin: 1rem 0;
    }
    
    .metric-item {
        background: linear-gradient(135deg, #f9fafb 0%, #f3f4f6 100%);
        padding: 0.75rem;
        border-radius: 8px;
        text-align: center;
        border: 1px solid #e5e7eb;
    }
    
    .metric-label {
        font-size: 0.75rem;
        color: #6b7280;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        margin-bottom: 0.25rem;
    }
    
    .metric-value {
        font-size: 1.25rem;
        font-weight: 700;
        color: #111827;
    }
    
    /* Interactive Elements */
    .expand-button {
        background: transparent;
        border: none;
        color: #667eea;
        font-weight: 600;
        cursor: pointer;
        padding: 0.5rem;
        border-radius: 8px;
        transition: all 0.2s;
    }
    
    .expand-button:hover {
        background: rgba(102, 126, 234, 0.1);
    }
    
    /* Tab Pills */
    .tab-pills {
        display: flex;
        gap: 0.5rem;
        padding: 0.25rem;
        background: #f3f4f6;
        border-radius: 12px;
        margin-bottom: 1rem;
    }
    
    .tab-pill {
        padding: 0.5rem 1rem;
        border-radius: 8px;
        font-weight: 600;
        transition: all 0.2s;
        cursor: pointer;
    }
    
    .tab-pill.active {
        background: white;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }
    
    /* Animation */
    @keyframes slideIn {
        from {
            opacity: 0;
            transform: translateY(10px);
        }
        to {
            opacity: 1;
            transform: translateY(0);
        }
    }
    
    .animate-in {
        animation: slideIn 0.3s ease-out;
    }
</style>
""", unsafe_allow_html=True)

# BioPython fallbacks
try:
    from Bio.Seq import Seq
    from Bio import SeqIO
except ImportError:
    class Seq:
        def __init__(self, data):
            self.data = str(data).upper()
        def __str__(self):
            return self.data
        def reverse_complement(self):
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
            return Seq(''.join(complement.get(base, 'N') for base in self.data[::-1]))
    SeqIO = None

@dataclass
class DimerAnalysis:
    """Results from dimer analysis"""
    dimer_type: str  # 'self', 'hetero'
    score: float
    alignment: str
    delta_g: float
    tm: float
    position: Tuple[int, int]
    severity: str  # 'low', 'medium', 'high', 'critical'

@dataclass
class PrimerPair:
    """Enhanced primer pair with dimer analysis"""
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
    confidence: int
    target_position: Optional[Tuple[int, int]] = None
    self_dimer_forward: Optional[DimerAnalysis] = None
    self_dimer_reverse: Optional[DimerAnalysis] = None
    hetero_dimer: Optional[DimerAnalysis] = None
    hairpin_forward: Optional[Dict] = None
    hairpin_reverse: Optional[Dict] = None
    additional_info: Dict = field(default_factory=dict)

class DimerAnalyzer:
    """Advanced dimer and secondary structure analyzer"""
    
    def __init__(self):
        self.complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        self.dg_values = {
            'AA': 1.2, 'TT': 1.2, 'AT': 0.9, 'TA': 0.9,
            'CA': 1.7, 'TG': 1.7, 'GT': 1.5, 'AC': 1.5,
            'CT': 1.3, 'AG': 1.3, 'GA': 1.3, 'TC': 1.3,
            'CG': 2.8, 'GC': 2.3, 'GG': 2.1, 'CC': 2.1
        }
    
    def analyze_primer_dimers(self, primer1: str, primer2: str = None) -> DimerAnalysis:
        """Comprehensive dimer analysis between primers"""
        if primer2 is None:
            primer2 = primer1
            dimer_type = 'self'
        else:
            dimer_type = 'hetero'
        
        # Reverse complement of primer2 for alignment
        primer2_rc = self.reverse_complement(primer2)
        
        # Find best alignment
        best_score = 0
        best_alignment = None
        best_position = (0, 0)
        
        for i in range(len(primer1)):
            for j in range(len(primer2_rc)):
                score, alignment = self.calculate_alignment_score(
                    primer1[i:], primer2_rc[j:]
                )
                if score > best_score:
                    best_score = score
                    best_alignment = alignment
                    best_position = (i, j)
        
        # Calculate thermodynamics
        delta_g = self.calculate_delta_g(best_alignment)
        tm = self.calculate_dimer_tm(best_alignment, delta_g)
        
        # Determine severity
        severity = self.assess_dimer_severity(best_score, delta_g, best_position)
        
        return DimerAnalysis(
            dimer_type=dimer_type,
            score=best_score,
            alignment=best_alignment,
            delta_g=delta_g,
            tm=tm,
            position=best_position,
            severity=severity
        )
    
    def analyze_hairpin(self, sequence: str) -> Dict:
        """Analyze potential hairpin structures"""
        min_stem = 4
        min_loop = 3
        best_hairpin = {'score': 0, 'structure': None, 'tm': 0}
        
        for i in range(len(sequence) - min_stem * 2 - min_loop):
            for j in range(i + min_stem + min_loop, len(sequence)):
                stem_length = min(j - i - min_loop, len(sequence) - j)
                
                if stem_length >= min_stem:
                    stem1 = sequence[i:i + stem_length]
                    stem2 = sequence[j:j + stem_length]
                    stem2_rc = self.reverse_complement(stem2)
                    
                    score = sum(1 for a, b in zip(stem1, stem2_rc) if a == b)
                    
                    if score > best_hairpin['score']:
                        loop = sequence[i + stem_length:j]
                        structure = f"{stem1}-{loop}-{stem2}"
                        tm = self.calculate_hairpin_tm(stem1, stem2_rc)
                        
                        best_hairpin = {
                            'score': score,
                            'structure': structure,
                            'stem_length': stem_length,
                            'loop_length': len(loop),
                            'position': (i, j),
                            'tm': tm
                        }
        
        return best_hairpin
    
    def reverse_complement(self, seq: str) -> str:
        """Get reverse complement of sequence"""
        return ''.join(self.complement.get(base, 'N') for base in seq[::-1])
    
    def calculate_alignment_score(self, seq1: str, seq2: str) -> Tuple[float, str]:
        """Calculate alignment score for two sequences"""
        min_len = min(len(seq1), len(seq2))
        matches = 0
        alignment = []
        
        for i in range(min_len):
            if seq1[i] == seq2[i]:
                matches += 1
                alignment.append('|')
            else:
                alignment.append('.')
        
        score = matches / min_len if min_len > 0 else 0
        alignment_str = f"{seq1[:min_len]}\n{''.join(alignment)}\n{seq2[:min_len]}"
        
        return score * 100, alignment_str
    
    def calculate_delta_g(self, alignment: str) -> float:
        """Calculate ΔG for dimer formation"""
        if not alignment:
            return 0
        
        lines = alignment.split('\n')
        if len(lines) < 3:
            return 0
        
        seq1 = lines[0]
        seq2 = lines[2]
        
        delta_g = 0
        for i in range(len(seq1) - 1):
            if i < len(seq2) - 1:
                dinuc = seq1[i:i+2] + seq2[i:i+2]
                if dinuc[:2] in self.dg_values:
                    delta_g -= self.dg_values[dinuc[:2]]
        
        return delta_g
    
    def calculate_dimer_tm(self, alignment: str, delta_g: float) -> float:
        """Calculate melting temperature of dimer"""
        if delta_g >= 0:
            return 0
        
        # Simplified Tm calculation for dimers
        # Tm = ΔH / (ΔS + R*ln(C))
        # Using approximation: Tm ≈ -ΔG * 1000 / (10.8)
        tm = abs(delta_g) * 1000 / 10.8
        
        return min(tm, 85)  # Cap at 85°C
    
    def calculate_hairpin_tm(self, stem1: str, stem2: str) -> float:
        """Calculate hairpin melting temperature"""
        matches = sum(1 for a, b in zip(stem1, stem2) if a == b)
        gc_count = sum(1 for a, b in zip(stem1, stem2) if a == b and a in 'GC')
        
        # Wallace rule approximation
        tm = 4 * gc_count + 2 * (matches - gc_count)
        
        return tm
    
    def assess_dimer_severity(self, score: float, delta_g: float, position: Tuple[int, int]) -> str:
        """Assess the severity of dimer formation"""
        # Check 3' end involvement (most critical)
        is_3_prime = position[0] > 15 or position[1] > 15
        
        if score > 80 and abs(delta_g) > 8 and is_3_prime:
            return 'critical'
        elif score > 60 and abs(delta_g) > 6:
            return 'high'
        elif score > 40 and abs(delta_g) > 4:
            return 'medium'
        else:
            return 'low'

class TargetedPrimerDesigner:
    """Design primers for specific target regions"""
    
    def __init__(self):
        self.designer = EnhancedPrimerDesigner()
        self.analyzer = DimerAnalyzer()
    
    def design_for_region(self, sequence: str, target_start: int, target_end: int, 
                          params: Dict = None) -> List[PrimerPair]:
        """Design primers targeting a specific region"""
        
        # Adjust design parameters for targeted region
        targeted_params = {
            'SEQUENCE_TARGET': f"{target_start},{target_end - target_start}",
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_PRODUCT_SIZE_RANGE': [[50, target_end - target_start + 100]]
        }
        
        if params:
            targeted_params.update(params)
        
        # Design primers
        primers = self.designer.design_primers_enhanced(
            {'sequence': sequence, 'header': 'Targeted_Region'},
            targeted_params
        )
        
        # Add target position info
        if primers[0]:
            for primer in primers[0]:
                primer.target_position = (target_start, target_end)
                
                # Perform dimer analysis
                primer.self_dimer_forward = self.analyzer.analyze_primer_dimers(primer.forward_seq)
                primer.self_dimer_reverse = self.analyzer.analyze_primer_dimers(primer.reverse_seq)
                primer.hetero_dimer = self.analyzer.analyze_primer_dimers(
                    primer.forward_seq, primer.reverse_seq
                )
                
                # Hairpin analysis
                primer.hairpin_forward = self.analyzer.analyze_hairpin(primer.forward_seq)
                primer.hairpin_reverse = self.analyzer.analyze_hairpin(primer.reverse_seq)
        
        return primers[0] if primers[0] else []

class ModernPrimerCard:
    """Generate modern social media-style primer cards"""
    
    @staticmethod
    def render_primer_thread(primers: List[PrimerPair], show_details: bool = False):
        """Render primers as a social media thread"""
        
        for idx, primer in enumerate(primers[:10]):  # Show top 10
            with st.container():
                # Thread node
                if idx > 0:
                    st.markdown(f'<div class="thread-node" style="top: {idx * 12}rem;"></div>', 
                               unsafe_allow_html=True)
                
                # Main card
                st.markdown('<div class="primer-card animate-in">', unsafe_allow_html=True)
                
                # Card header
                col1, col2, col3 = st.columns([2, 3, 2])
                
                with col1:
                    st.markdown(f"### 🧬 Primer Set {idx + 1}")
                    ModernPrimerCard._render_confidence_badge(primer.confidence)
                
                with col2:
                    if primer.target_position:
                        st.caption(f"📍 Target: {primer.target_position[0]}-{primer.target_position[1]} bp")
                    st.caption(f"📏 Product: {primer.product_size} bp")
                
                with col3:
                    # Quick metrics
                    st.markdown('<div class="metrics-grid">', unsafe_allow_html=True)
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.metric("ΔTm", f"{abs(primer.forward_tm - primer.reverse_tm):.1f}°C", 
                                 delta_color="inverse")
                    with col_b:
                        st.metric("Avg GC", f"{(primer.forward_gc + primer.reverse_gc)/2:.0f}%")
                    st.markdown('</div>', unsafe_allow_html=True)
                
                # Expandable details
                with st.expander("🔍 View Details", expanded=show_details):
                    tabs = st.tabs(["📝 Sequences", "🔬 Dimer Analysis", "📊 Properties", "🧪 Validation"])
                    
                    with tabs[0]:
                        ModernPrimerCard._render_sequences(primer)
                    
                    with tabs[1]:
                        ModernPrimerCard._render_dimer_analysis(primer)
                    
                    with tabs[2]:
                        ModernPrimerCard._render_properties(primer)
                    
                    with tabs[3]:
                        ModernPrimerCard._render_validation_tools(primer, idx)
                
                st.markdown('</div>', unsafe_allow_html=True)
                
                if idx < len(primers) - 1:
                    st.markdown("---")
    
    @staticmethod
    def _render_confidence_badge(confidence: int):
        """Render confidence score as a badge"""
        if confidence >= 85:
            badge_class = "status-excellent"
            icon = "✅"
            text = "EXCELLENT"
        elif confidence >= 70:
            badge_class = "status-good"
            icon = "⚠️"
            text = "GOOD"
        else:
            badge_class = "status-warning"
            icon = "⚡"
            text = "CHECK"
        
        st.markdown(
            f'<span class="status-badge {badge_class}">{icon} {text} ({confidence}%)</span>',
            unsafe_allow_html=True
        )
    
    @staticmethod
    def _render_sequences(primer: PrimerPair):
        """Render primer sequences with visual enhancements"""
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Forward Primer (5' → 3')**")
            st.markdown(f'<div class="sequence-display">{primer.forward_seq}</div>', 
                       unsafe_allow_html=True)
            
            # Metrics
            metrics_html = f"""
            <div class="metrics-grid">
                <div class="metric-item">
                    <div class="metric-label">Length</div>
                    <div class="metric-value">{len(primer.forward_seq)}</div>
                </div>
                <div class="metric-item">
                    <div class="metric-label">Tm</div>
                    <div class="metric-value">{primer.forward_tm:.1f}°C</div>
                </div>
                <div class="metric-item">
                    <div class="metric-label">GC%</div>
                    <div class="metric-value">{primer.forward_gc:.0f}%</div>
                </div>
            </div>
            """
            st.markdown(metrics_html, unsafe_allow_html=True)
        
        with col2:
            st.markdown("**Reverse Primer (5' → 3')**")
            st.markdown(f'<div class="sequence-display">{primer.reverse_seq}</div>', 
                       unsafe_allow_html=True)
            
            # Metrics
            metrics_html = f"""
            <div class="metrics-grid">
                <div class="metric-item">
                    <div class="metric-label">Length</div>
                    <div class="metric-value">{len(primer.reverse_seq)}</div>
                </div>
                <div class="metric-item">
                    <div class="metric-label">Tm</div>
                    <div class="metric-value">{primer.reverse_tm:.1f}°C</div>
                </div>
                <div class="metric-item">
                    <div class="metric-label">GC%</div>
                    <div class="metric-value">{primer.reverse_gc:.0f}%</div>
                </div>
            </div>
            """
            st.markdown(metrics_html, unsafe_allow_html=True)
    
    @staticmethod
    def _render_dimer_analysis(primer: PrimerPair):
        """Render comprehensive dimer analysis"""
        
        # Self-dimers
        st.markdown("#### Self-Dimer Analysis")
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Forward Primer Self-Dimer**")
            if primer.self_dimer_forward:
                ModernPrimerCard._display_dimer_result(primer.self_dimer_forward)
        
        with col2:
            st.markdown("**Reverse Primer Self-Dimer**")
            if primer.self_dimer_reverse:
                ModernPrimerCard._display_dimer_result(primer.self_dimer_reverse)
        
        # Heterodimer
        st.markdown("#### Heterodimer Analysis")
        if primer.hetero_dimer:
            ModernPrimerCard._display_dimer_result(primer.hetero_dimer)
        
        # Hairpin structures
        st.markdown("#### Hairpin Structure Analysis")
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Forward Primer Hairpin**")
            if primer.hairpin_forward and primer.hairpin_forward['score'] > 0:
                st.info(f"Stem: {primer.hairpin_forward['stem_length']} bp | "
                       f"Loop: {primer.hairpin_forward['loop_length']} bp | "
                       f"Tm: {primer.hairpin_forward['tm']:.1f}°C")
        
        with col2:
            st.markdown("**Reverse Primer Hairpin**")
            if primer.hairpin_reverse and primer.hairpin_reverse['score'] > 0:
                st.info(f"Stem: {primer.hairpin_reverse['stem_length']} bp | "
                       f"Loop: {primer.hairpin_reverse['loop_length']} bp | "
                       f"Tm: {primer.hairpin_reverse['tm']:.1f}°C")
    
    @staticmethod
    def _display_dimer_result(dimer: DimerAnalysis):
        """Display dimer analysis result"""
        severity_colors = {
            'low': 'info',
            'medium': 'warning',
            'high': 'warning',
            'critical': 'error'
        }
        
        severity_icons = {
            'low': '✓',
            'medium': '⚠️',
            'high': '⚠️',
            'critical': '⛔'
        }
        
        # Create severity message
        severity_msg = f"{severity_icons[dimer.severity]} {dimer.severity.upper()}"
        
        if dimer.severity in ['high', 'critical']:
            st.error(f"{severity_msg} | Score: {dimer.score:.1f}% | ΔG: {dimer.delta_g:.1f} kcal/mol")
        elif dimer.severity == 'medium':
            st.warning(f"{severity_msg} | Score: {dimer.score:.1f}% | ΔG: {dimer.delta_g:.1f} kcal/mol")
        else:
            st.success(f"{severity_msg} | Score: {dimer.score:.1f}% | ΔG: {dimer.delta_g:.1f} kcal/mol")
        
        # Show alignment in code block
        if dimer.alignment:
            st.code(dimer.alignment, language=None)
    
    @staticmethod
    def _render_properties(primer: PrimerPair):
        """Render detailed primer properties"""
        
        # Create properties dataframe
        properties = {
            'Property': [
                'Product Size',
                'Forward Position',
                'Reverse Position',
                'Penalty Score',
                'Tm Difference',
                'GC Difference',
                '3\' End Stability (F)',
                '3\' End Stability (R)'
            ],
            'Value': [
                f"{primer.product_size} bp",
                f"{primer.forward_position}",
                f"{primer.reverse_position}",
                f"{primer.penalty:.2f}",
                f"{abs(primer.forward_tm - primer.reverse_tm):.1f}°C",
                f"{abs(primer.forward_gc - primer.reverse_gc):.1f}%",
                primer.additional_info.get('left_end_stability', 'N/A'),
                primer.additional_info.get('right_end_stability', 'N/A')
            ]
        }
        
        df = pd.DataFrame(properties)
        st.dataframe(df, hide_index=True, use_container_width=True)
    
    @staticmethod
    def _render_validation_tools(primer: PrimerPair, idx: int):
        """Render validation and export tools"""
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            # BLAST validation
            species = st.selectbox(
                "Organism",
                ["Homo sapiens", "Mus musculus"],
                key=f"blast_species_{idx}"
            )
            
            blast_url = ModernPrimerCard._create_blast_url(
                primer.forward_seq, primer.reverse_seq, species
            )
            st.link_button("🎯 Validate with BLAST", blast_url, use_container_width=True)
        
        with col2:
            # Copy sequences
            sequences = f"Forward: {primer.forward_seq}\nReverse: {primer.reverse_seq}"
            st.download_button(
                "📋 Copy Sequences",
                data=sequences,
                file_name=f"primer_set_{idx + 1}.txt",
                key=f"copy_{idx}",
                use_container_width=True
            )
        
        with col3:
            # Export JSON
            primer_data = {
                'forward_seq': primer.forward_seq,
                'reverse_seq': primer.reverse_seq,
                'forward_tm': primer.forward_tm,
                'reverse_tm': primer.reverse_tm,
                'forward_gc': primer.forward_gc,
                'reverse_gc': primer.reverse_gc,
                'product_size': primer.product_size,
                'confidence': primer.confidence
            }
            
            st.download_button(
                "💾 Export JSON",
                data=json.dumps(primer_data, indent=2),
                file_name=f"primer_{idx + 1}.json",
                key=f"export_{idx}",
                use_container_width=True
            )
    
    @staticmethod
    def _create_blast_url(forward_seq: str, reverse_seq: str, organism: str) -> str:
        """Create NCBI Primer-BLAST URL"""
        base_url = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi"
        params = {
            'PRIMER_LEFT_INPUT': forward_seq,
            'PRIMER_RIGHT_INPUT': reverse_seq,
            'PRIMER_PRODUCT_MIN': 70,
            'PRIMER_PRODUCT_MAX': 500,
            'ORGANISM': organism,
            'PRIMER_SPECIFICITY_DATABASE': 'refseq_rna'
        }
        return f"{base_url}?" + urllib.parse.urlencode(params)

class EnhancedAnalysisDashboard:
    """Create beautiful analysis dashboards"""
    
    @staticmethod
    def render_comprehensive_analysis(primers: List[PrimerPair]):
        """Render comprehensive analysis with modern visualizations"""
        
        st.markdown("### 📊 Comprehensive Primer Analysis")
        
        # Summary metrics
        col1, col2, col3, col4 = st.columns(4)
        
        excellent_count = sum(1 for p in primers if p.confidence >= 85)
        avg_confidence = np.mean([p.confidence for p in primers])
        avg_tm_diff = np.mean([abs(p.forward_tm - p.reverse_tm) for p in primers])
        dimer_issues = sum(1 for p in primers 
                          if p.hetero_dimer and p.hetero_dimer.severity in ['high', 'critical'])
        
        with col1:
            st.metric("Excellent Primers", f"{excellent_count}/{len(primers)}", 
                     f"{excellent_count/len(primers)*100:.0f}%")
        
        with col2:
            st.metric("Avg Confidence", f"{avg_confidence:.1f}%",
                     "Good" if avg_confidence >= 70 else "Needs Review")
        
        with col3:
            st.metric("Avg ΔTm", f"{avg_tm_diff:.2f}°C",
                     "Optimal" if avg_tm_diff <= 1 else "Check")
        
        with col4:
            st.metric("Dimer Issues", dimer_issues,
                     "Clear" if dimer_issues == 0 else f"-{dimer_issues}")
        
        # Visualizations
        tabs = st.tabs(["📈 Overview", "🔬 Dimer Analysis", "🎯 Target Coverage", "📊 Distributions"])
        
        with tabs[0]:
            EnhancedAnalysisDashboard._render_overview_charts(primers)
        
        with tabs[1]:
            EnhancedAnalysisDashboard._render_dimer_heatmap(primers)
        
        with tabs[2]:
            EnhancedAnalysisDashboard._render_target_coverage(primers)
        
        with tabs[3]:
            EnhancedAnalysisDashboard._render_distributions(primers)
    
    @staticmethod
    def _render_overview_charts(primers: List[PrimerPair]):
        """Render overview charts"""
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Confidence radar chart
            fig = go.Figure()
            
            categories = ['Confidence', 'Tm Match', 'GC Balance', 'Size Optimal', 'No Dimers']
            
            for i, primer in enumerate(primers[:5]):  # Top 5
                values = [
                    primer.confidence,
                    100 - min(abs(primer.forward_tm - primer.reverse_tm) * 10, 100),
                    100 - abs(50 - (primer.forward_gc + primer.reverse_gc) / 2) * 2,
                    100 if 80 <= primer.product_size <= 150 else 50,
                    100 if primer.hetero_dimer and primer.hetero_dimer.severity == 'low' else 50
                ]
                
                fig.add_trace(go.Scatterpolar(
                    r=values,
                    theta=categories,
                    fill='toself',
                    name=f'Set {i+1}',
                    opacity=0.6
                ))
            
            fig.update_layout(
                polar=dict(
                    radialaxis=dict(
                        visible=True,
                        range=[0, 100]
                    )),
                showlegend=True,
                title="Primer Quality Radar",
                height=400
            )
            
            st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            # Tm correlation scatter
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x=[p.forward_tm for p in primers],
                y=[p.reverse_tm for p in primers],
                mode='markers',
                marker=dict(
                    size=[p.confidence/5 for p in primers],
                    color=[p.confidence for p in primers],
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title="Confidence")
                ),
                text=[f"Set {i+1}<br>Product: {p.product_size} bp" 
                      for i, p in enumerate(primers)],
                hovertemplate='%{text}<br>Forward Tm: %{x:.1f}°C<br>Reverse Tm: %{y:.1f}°C'
            ))
            
            # Add ideal line
            min_tm = min([p.forward_tm for p in primers] + [p.reverse_tm for p in primers])
            max_tm = max([p.forward_tm for p in primers] + [p.reverse_tm for p in primers])
            
            fig.add_trace(go.Scatter(
                x=[min_tm, max_tm],
                y=[min_tm, max_tm],
                mode='lines',
                line=dict(dash='dash', color='red'),
                name='Ideal Match',
                showlegend=True
            ))
            
            fig.update_layout(
                title="Melting Temperature Correlation",
                xaxis_title="Forward Tm (°C)",
                yaxis_title="Reverse Tm (°C)",
                height=400
            )
            
            st.plotly_chart(fig, use_container_width=True)
    
    @staticmethod
    def _render_dimer_heatmap(primers: List[PrimerPair]):
        """Render dimer analysis heatmap"""
        
        # Create dimer matrix
        n = min(len(primers), 10)
        dimer_matrix = np.zeros((n * 2, n * 2))
        labels = []
        
        for i in range(n):
            labels.append(f"F{i+1}")
            labels.append(f"R{i+1}")
        
        analyzer = DimerAnalyzer()
        
        # Calculate all pairwise dimers
        all_sequences = []
        for p in primers[:n]:
            all_sequences.append(p.forward_seq)
            all_sequences.append(p.reverse_seq)
        
        for i, seq1 in enumerate(all_sequences):
            for j, seq2 in enumerate(all_sequences):
                if i != j:
                    dimer = analyzer.analyze_primer_dimers(seq1, seq2)
                    dimer_matrix[i, j] = dimer.score
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=dimer_matrix,
            x=labels,
            y=labels,
            colorscale='RdYlGn_r',
            text=np.round(dimer_matrix, 1),
            texttemplate='%{text}',
            textfont={"size": 8},
            colorbar=dict(title="Dimer Score %")
        ))
        
        fig.update_layout(
            title="Primer-Primer Dimer Matrix",
            xaxis_title="Primer",
            yaxis_title="Primer",
            height=500
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Dimer warnings
        critical_dimers = []
        for i in range(len(dimer_matrix)):
            for j in range(i+1, len(dimer_matrix)):
                if dimer_matrix[i, j] > 70:
                    critical_dimers.append(f"{labels[i]} × {labels[j]}: {dimer_matrix[i, j]:.1f}%")
        
        if critical_dimers:
            st.warning(f"⚠️ Critical dimer pairs detected: {', '.join(critical_dimers[:3])}")
    
    @staticmethod
    def _render_target_coverage(primers: List[PrimerPair]):
        """Render target coverage visualization"""
        
        if not any(p.target_position for p in primers):
            st.info("No target position information available")
            return
        
        # Create Gantt-like chart for primer coverage
        fig = go.Figure()
        
        for i, primer in enumerate(primers[:10]):
            if primer.target_position:
                # Forward primer
                fig.add_trace(go.Scatter(
                    x=[primer.forward_position, 
                       primer.forward_position + len(primer.forward_seq)],
                    y=[i, i],
                    mode='lines',
                    line=dict(color='blue', width=10),
                    name=f'F{i+1}',
                    showlegend=i == 0,
                    legendgroup='forward',
                    hovertemplate=f'Forward {i+1}<br>Pos: %{{x}}'
                ))
                
                # Reverse primer
                fig.add_trace(go.Scatter(
                    x=[primer.reverse_position - len(primer.reverse_seq),
                       primer.reverse_position],
                    y=[i, i],
                    mode='lines',
                    line=dict(color='red', width=10),
                    name=f'R{i+1}',
                    showlegend=i == 0,
                    legendgroup='reverse',
                    hovertemplate=f'Reverse {i+1}<br>Pos: %{{x}}'
                ))
                
                # Target region
                if primer.target_position:
                    fig.add_shape(
                        type="rect",
                        x0=primer.target_position[0],
                        x1=primer.target_position[1],
                        y0=i - 0.2,
                        y1=i + 0.2,
                        fillcolor="lightgreen",
                        opacity=0.3,
                        line_width=0,
                    )
        
        fig.update_layout(
            title="Primer Binding Sites & Target Coverage",
            xaxis_title="Position (bp)",
            yaxis_title="Primer Set",
            height=400,
            yaxis=dict(
                tickmode='array',
                tickvals=list(range(min(10, len(primers)))),
                ticktext=[f'Set {i+1}' for i in range(min(10, len(primers)))]
            )
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    @staticmethod
    def _render_distributions(primers: List[PrimerPair]):
        """Render distribution plots"""
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Product Size Distribution', 'GC Content Distribution',
                           'Tm Distribution', 'Confidence Score Distribution')
        )
        
        # Product size
        fig.add_trace(
            go.Histogram(
                x=[p.product_size for p in primers],
                nbinsx=20,
                marker_color='lightblue',
                name='Product Size'
            ),
            row=1, col=1
        )
        
        # GC content
        gc_values = []
        for p in primers:
            gc_values.extend([p.forward_gc, p.reverse_gc])
        
        fig.add_trace(
            go.Histogram(
                x=gc_values,
                nbinsx=20,
                marker_color='lightgreen',
                name='GC Content'
            ),
            row=1, col=2
        )
        
        # Tm distribution
        tm_values = []
        for p in primers:
            tm_values.extend([p.forward_tm, p.reverse_tm])
        
        fig.add_trace(
            go.Histogram(
                x=tm_values,
                nbinsx=20,
                marker_color='coral',
                name='Tm'
            ),
            row=2, col=1
        )
        
        # Confidence scores
        fig.add_trace(
            go.Histogram(
                x=[p.confidence for p in primers],
                nbinsx=20,
                marker_color='purple',
                name='Confidence'
            ),
            row=2, col=2
        )
        
        fig.update_layout(height=600, showlegend=False)
        fig.update_xaxes(title_text="Product Size (bp)", row=1, col=1)
        fig.update_xaxes(title_text="GC %", row=1, col=2)
        fig.update_xaxes(title_text="Tm (°C)", row=2, col=1)
        fig.update_xaxes(title_text="Confidence %", row=2, col=2)
        
        st.plotly_chart(fig, use_container_width=True)

# Keep the original helper classes
class EnhancedPrimerDesigner:
    """Enhanced primer designer with improved algorithms"""
    
    def __init__(self):
        self.optimal_params = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 23,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 58.0,
            'PRIMER_MAX_TM': 62.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_OPT_GC_PERCENT': 50.0,
            'PRIMER_MAX_POLY_X': 3,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY_TH': 40.0,
            'PRIMER_MAX_SELF_END_TH': 35.0,
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': 40.0,
            'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[80, 150]],
            'PRIMER_NUM_RETURN': 15,
            'PRIMER_MAX_HAIRPIN_TH': 40.0,
            'PRIMER_INTERNAL_MAX_HAIRPIN_TH': 40.0,
            'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 7,
            'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION': 7
        }
        
        # qPCR-specific scoring weights
        self.scoring_weights = {
            'tm_difference': 0.25,
            'gc_content': 0.20,
            'product_size': 0.15,
            'self_complementarity': 0.15,
            'end_stability': 0.10,
            'complexity': 0.10,
            'penalty': 0.05
        }
        
        self.dimer_analyzer = DimerAnalyzer()
    
    def parse_sequence_input(self, sequence_input: str) -> Tuple[Optional[Dict], Optional[str]]:
        """Enhanced sequence parser supporting multiple formats"""
        try:
            sequence_input = sequence_input.strip()
            
            # Try to parse as FASTA
            if sequence_input.startswith('>'):
                if SeqIO:
                    # Use BioPython if available
                    fasta_io = StringIO(sequence_input)
                    records = list(SeqIO.parse(fasta_io, "fasta"))
                    
                    if not records:
                        return None, "No valid sequences found in FASTA input"
                    
                    # Use the first record
                    record = records[0]
                    sequence = str(record.seq).upper()
                    header = record.description
                else:
                    # Manual FASTA parsing if BioPython not available
                    lines = sequence_input.strip().split('\n')
                    header = lines[0][1:].strip() if lines[0].startswith('>') else "User_Sequence"
                    sequence = ''.join(lines[1:])
                    sequence = re.sub(r'[^ATGCNatgcn]', '', sequence).upper()
            else:
                # Assume raw sequence
                sequence = re.sub(r'[^ATGCNatgcn]', '', sequence_input).upper()
                header = "User_Sequence"
            
            # Validate sequence
            if len(sequence) < 100:
                return None, f"Sequence too short ({len(sequence)} bp). Minimum 100 bp required."
            
            if len(sequence) > 10000:
                return None, f"Sequence too long ({len(sequence)} bp). Maximum 10,000 bp allowed."
            
            # Check for excessive N content
            n_content = sequence.count('N') / len(sequence) * 100
            if n_content > 10:
                return None, f"Sequence contains too many ambiguous bases (N content: {n_content:.1f}%)"
            
            return {
                'sequence': sequence,
                'header': header,
                'length': len(sequence),
                'gc_content': self.calculate_gc(sequence),
                'n_content': n_content
            }, None
            
        except Exception as e:
            return None, f"Error parsing sequence: {str(e)}"
    
    def calculate_gc(self, sequence: str) -> float:
        """Calculate GC content"""
        if not sequence:
            return 0.0
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    def calculate_tm(self, sequence: str, dnac: float = 50, saltc: float = 50) -> float:
        """Calculate melting temperature using nearest neighbor method"""
        sequence = str(sequence).upper()
        if not sequence:
            return 0.0
        
        # Simple Tm calculation for DNA
        if len(sequence) < 14:
            return 2 * (sequence.count('A') + sequence.count('T')) + 4 * (sequence.count('G') + sequence.count('C'))
        
        # Using GC content method for longer sequences
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        
        # Basic formula: Tm = 81.5 + 0.41(%GC) - 675/length + salt correction
        basic_tm = 81.5 + (0.41 * gc_content * 100) - (675 / len(sequence))
        
        # Salt correction (simplified)
        salt_correction = 12.5 * np.log10(saltc / 1000)
        
        return basic_tm + salt_correction
    
    def design_primers_enhanced(self, sequence_data: Dict, custom_params: Dict = None) -> Tuple[Optional[List[PrimerPair]], Optional[str]]:
        """Enhanced primer design with dimer analysis"""
        try:
            sequence = sequence_data['sequence']
            seq_id = sequence_data['header']
            
            # Merge custom parameters if provided
            params = self.optimal_params.copy()
            if custom_params:
                params.update(custom_params)
            
            # Multiple design attempts with different parameters
            all_primers = []
            param_sets = [
                params,  # Original parameters
                {**params, 'PRIMER_PRODUCT_SIZE_RANGE': [[70, 200]]},  # Wider range
                {**params, 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0},  # Wider Tm range
                {**params, 'PRIMER_MIN_GC': 35.0, 'PRIMER_MAX_GC': 65.0}  # Wider GC range
            ]
            
            for param_set in param_sets:
                seq_args = {
                    'SEQUENCE_ID': seq_id,
                    'SEQUENCE_TEMPLATE': sequence,
                }
                
                global_args = {
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_PICK_INTERNAL_OLIGO': 0,
                    **param_set
                }
                
                try:
                    # Try different Primer3 API calls
                    try:
                        results = primer3.designPrimers(seq_args, global_args)
                    except AttributeError:
                        try:
                            results = primer3.bindings.design_primers(seq_args, global_args)
                        except AttributeError:
                            try:
                                results = primer3.bindings.designPrimers(seq_args, global_args)
                            except:
                                combined_args = {**seq_args, **global_args}
                                results = primer3.bindings.designPrimers(combined_args)
                    
                    num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
                    
                    # Process primers if found
                    if num_returned > 0:
                        for i in range(num_returned):
                            primer_pair = self._extract_primer_pair(results, i)
                            if primer_pair:
                                # Add dimer analysis
                                primer_pair.self_dimer_forward = self.dimer_analyzer.analyze_primer_dimers(
                                    primer_pair.forward_seq
                                )
                                primer_pair.self_dimer_reverse = self.dimer_analyzer.analyze_primer_dimers(
                                    primer_pair.reverse_seq
                                )
                                primer_pair.hetero_dimer = self.dimer_analyzer.analyze_primer_dimers(
                                    primer_pair.forward_seq, primer_pair.reverse_seq
                                )
                                primer_pair.hairpin_forward = self.dimer_analyzer.analyze_hairpin(
                                    primer_pair.forward_seq
                                )
                                primer_pair.hairpin_reverse = self.dimer_analyzer.analyze_hairpin(
                                    primer_pair.reverse_seq
                                )
                                
                                all_primers.append(primer_pair)
                
                except Exception as e:
                    continue
            
            if not all_primers:
                return None, "No suitable primers found. Try adjusting sequence or parameters."
            
            # Remove duplicates and sort by confidence
            unique_primers = self._remove_duplicate_primers(all_primers)
            unique_primers.sort(key=lambda x: x.confidence, reverse=True)
            
            # Return top primers
            return unique_primers[:10], None
            
        except Exception as e:
            return None, f"Error in primer design: {str(e)}"
    
    def _extract_primer_pair(self, primer_data: Dict, pair_num: int) -> Optional[PrimerPair]:
        """Extract and create PrimerPair object from Primer3 results"""
        try:
            forward_seq = primer_data.get(f'PRIMER_LEFT_{pair_num}_SEQUENCE', '')
            reverse_seq = primer_data.get(f'PRIMER_RIGHT_{pair_num}_SEQUENCE', '')
            
            if not forward_seq or not reverse_seq:
                return None
            
            # Calculate advanced confidence score
            confidence = self._calculate_advanced_confidence(primer_data, pair_num)
            
            return PrimerPair(
                forward_seq=forward_seq,
                reverse_seq=reverse_seq,
                forward_tm=primer_data.get(f'PRIMER_LEFT_{pair_num}_TM', 60),
                reverse_tm=primer_data.get(f'PRIMER_RIGHT_{pair_num}_TM', 60),
                forward_gc=primer_data.get(f'PRIMER_LEFT_{pair_num}_GC_PERCENT', 50),
                reverse_gc=primer_data.get(f'PRIMER_RIGHT_{pair_num}_GC_PERCENT', 50),
                product_size=primer_data.get(f'PRIMER_PAIR_{pair_num}_PRODUCT_SIZE', 100),
                forward_position=primer_data.get(f'PRIMER_LEFT_{pair_num}', [0])[0],
                reverse_position=primer_data.get(f'PRIMER_RIGHT_{pair_num}', [0])[0],
                penalty=primer_data.get(f'PRIMER_PAIR_{pair_num}_PENALTY', 0),
                confidence=confidence,
                additional_info={
                    'self_any_th': primer_data.get(f'PRIMER_PAIR_{pair_num}_COMPL_ANY_TH', 0),
                    'self_end_th': primer_data.get(f'PRIMER_PAIR_{pair_num}_COMPL_END_TH', 0),
                    'left_end_stability': primer_data.get(f'PRIMER_LEFT_{pair_num}_END_STABILITY', 0),
                    'right_end_stability': primer_data.get(f'PRIMER_RIGHT_{pair_num}_END_STABILITY', 0)
                }
            )
        except:
            return None
    
    def _calculate_advanced_confidence(self, primer_data: Dict, pair_num: int) -> int:
        """Calculate comprehensive confidence score for qPCR primers"""
        scores = {}
        
        # Get primer properties
        left_tm = primer_data.get(f'PRIMER_LEFT_{pair_num}_TM', 60)
        right_tm = primer_data.get(f'PRIMER_RIGHT_{pair_num}_TM', 60)
        left_gc = primer_data.get(f'PRIMER_LEFT_{pair_num}_GC_PERCENT', 50)
        right_gc = primer_data.get(f'PRIMER_RIGHT_{pair_num}_GC_PERCENT', 50)
        product_size = primer_data.get(f'PRIMER_PAIR_{pair_num}_PRODUCT_SIZE', 100)
        penalty = primer_data.get(f'PRIMER_PAIR_{pair_num}_PENALTY', 0)
        left_seq = primer_data.get(f'PRIMER_LEFT_{pair_num}_SEQUENCE', '')
        right_seq = primer_data.get(f'PRIMER_RIGHT_{pair_num}_SEQUENCE', '')
        
        # 1. Tm difference score (optimal: ≤1°C)
        tm_diff = abs(left_tm - right_tm)
        if tm_diff <= 0.5:
            scores['tm_difference'] = 100
        elif tm_diff <= 1.0:
            scores['tm_difference'] = 90
        elif tm_diff <= 2.0:
            scores['tm_difference'] = 70
        elif tm_diff <= 3.0:
            scores['tm_difference'] = 50
        else:
            scores['tm_difference'] = max(0, 30 - tm_diff * 5)
        
        # 2. GC content score (optimal: 45-55%)
        gc_scores = []
        for gc in [left_gc, right_gc]:
            if 45 <= gc <= 55:
                gc_scores.append(100)
            elif 40 <= gc <= 60:
                gc_scores.append(80)
            elif 35 <= gc <= 65:
                gc_scores.append(60)
            else:
                gc_scores.append(40)
        scores['gc_content'] = sum(gc_scores) / 2
        
        # 3. Product size score (optimal: 80-120 bp for qPCR)
        if 80 <= product_size <= 120:
            scores['product_size'] = 100
        elif 70 <= product_size <= 150:
            scores['product_size'] = 80
        elif 60 <= product_size <= 200:
            scores['product_size'] = 60
        else:
            scores['product_size'] = 40
        
        # 4. Self-complementarity score
        self_any_th = primer_data.get(f'PRIMER_PAIR_{pair_num}_COMPL_ANY_TH', 0)
        self_end_th = primer_data.get(f'PRIMER_PAIR_{pair_num}_COMPL_END_TH', 0)
        
        if self_any_th < 30 and self_end_th < 20:
            scores['self_complementarity'] = 100
        elif self_any_th < 35 and self_end_th < 25:
            scores['self_complementarity'] = 80
        elif self_any_th < 40 and self_end_th < 30:
            scores['self_complementarity'] = 60
        else:
            scores['self_complementarity'] = 40
        
        # 5. 3' end stability score
        end_scores = []
        for seq in [left_seq, right_seq]:
            if seq and seq[-1] in ['G', 'C']:
                if seq[-2:] == 'GC' or seq[-2:] == 'CG':
                    end_scores.append(100)
                else:
                    end_scores.append(80)
            else:
                end_scores.append(60)
        scores['end_stability'] = sum(end_scores) / 2
        
        # 6. Sequence complexity score
        complexity_scores = []
        for seq in [left_seq, right_seq]:
            if seq:
                # Check for problematic patterns
                problems = 0
                if 'AAAA' in seq or 'TTTT' in seq:
                    problems += 1
                if 'CCCC' in seq or 'GGGG' in seq:
                    problems += 1
                if re.search(r'(.)\1{4,}', seq):  # 5+ repeats
                    problems += 2
                
                complexity_scores.append(max(40, 100 - problems * 20))
        scores['complexity'] = sum(complexity_scores) / 2 if complexity_scores else 50
        
        # 7. Penalty score
        if penalty <= 0.5:
            scores['penalty'] = 100
        elif penalty <= 1.0:
            scores['penalty'] = 80
        elif penalty <= 2.0:
            scores['penalty'] = 60
        else:
            scores['penalty'] = max(20, 60 - penalty * 10)
        
        # Calculate weighted final score
        final_score = sum(scores[key] * self.scoring_weights[key] for key in scores)
        
        return int(max(0, min(100, final_score)))
    
    def _remove_duplicate_primers(self, primers: List[PrimerPair]) -> List[PrimerPair]:
        """Remove duplicate primer pairs"""
        unique_primers = []
        seen = set()
        
        for primer in primers:
            key = (primer.forward_seq, primer.reverse_seq)
            if key not in seen:
                seen.add(key)
                unique_primers.append(primer)
        
        return unique_primers

def main():
    # Header with version badge
    st.markdown(
        '<h1 class="main-header">PrimersQuest Pro <span class="version-badge">v3.0</span></h1>',
        unsafe_allow_html=True
    )
    st.markdown("""
    <p style="text-align: center; color: #718096; font-size: 1.125rem;">
    Advanced qPCR primer design with dimer analysis and modern visualization
    </p>
    """, unsafe_allow_html=True)
    
    # Initialize tools
    designer = EnhancedPrimerDesigner()
    targeted_designer = TargetedPrimerDesigner()
    analyzer = DimerAnalyzer()
    
    # Sidebar
    with st.sidebar:
        st.markdown("### ⚙️ Configuration")
        
        # Design mode
        design_mode = st.radio(
            "Design Mode",
            ["Standard", "Targeted Region", "Dimer-Optimized"],
            help="Choose primer design strategy"
        )
        
        # Parameters
        with st.expander("🔧 Design Parameters", expanded=False):
            custom_params = {}
            
            st.markdown("#### Temperature")
            min_tm = st.slider("Min Tm (°C)", 55.0, 65.0, 58.0, 0.5)
            max_tm = st.slider("Max Tm (°C)", 55.0, 65.0, 62.0, 0.5)
            custom_params['PRIMER_MIN_TM'] = min_tm
            custom_params['PRIMER_MAX_TM'] = max_tm
            
            st.markdown("#### Product Size")
            min_size = st.number_input("Min Size (bp)", 50, 300, 80)
            max_size = st.number_input("Max Size (bp)", 50, 500, 150)
            custom_params['PRIMER_PRODUCT_SIZE_RANGE'] = [[min_size, max_size]]
            
            st.markdown("#### Dimer Thresholds")
            max_self_any = st.slider("Max Self Any (°C)", 30.0, 50.0, 40.0, 1.0)
            max_self_end = st.slider("Max Self End (°C)", 20.0, 40.0, 35.0, 1.0)
            custom_params['PRIMER_MAX_SELF_ANY_TH'] = max_self_any
            custom_params['PRIMER_MAX_SELF_END_TH'] = max_self_end
        
        st.markdown("### 📊 Display Options")
        show_details = st.checkbox("Show detailed cards", value=False)
        max_primers = st.slider("Max primers to display", 5, 20, 10)
        
        st.markdown("### 🔗 Resources")
        st.markdown("""
        - [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
        - [Primer3 Docs](http://primer3.org/)
        - [qPCR Guidelines](https://pubmed.ncbi.nlm.nih.gov/19246619/)
        """)
    
    # Main tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "🎯 Design Primers",
        "🔬 Dimer Analysis",
        "📊 Dashboard",
        "📚 Help"
    ])
    
    with tab1:
        st.markdown("### Design Custom Primers")
        
        # Sequence input
        sequence_input = st.text_area(
            "Enter target sequence (FASTA or raw):",
            height=200,
            placeholder=">Gene_Name\nATGCGATCG..."
        )
        
        # Targeted region options
        if design_mode == "Targeted Region":
            col1, col2 = st.columns(2)
            with col1:
                target_start = st.number_input("Target Start (bp)", 1, 10000, 100)
            with col2:
                target_end = st.number_input("Target End (bp)", 1, 10000, 300)
        
        if st.button("🚀 Design Primers", type="primary", use_container_width=True):
            if sequence_input.strip():
                with st.spinner("🧬 Analyzing sequence and designing primers..."):
                    # Parse sequence
                    seq_data, error = designer.parse_sequence_input(sequence_input)
                    
                    if error:
                        st.error(f"❌ {error}")
                    else:
                        # Design primers based on mode
                        if design_mode == "Targeted Region":
                            primers = targeted_designer.design_for_region(
                                seq_data['sequence'],
                                target_start,
                                target_end,
                                custom_params
                            )
                            primers = (primers, None) if primers else (None, "No primers found")
                        else:
                            primers = designer.design_primers_enhanced(seq_data, custom_params)
                        
                        if primers[1]:  # Error
                            st.error(f"❌ {primers[1]}")
                        elif primers[0]:  # Success
                            st.success(f"✅ Successfully designed {len(primers[0])} primer pairs!")
                            st.session_state['designed_primers'] = primers[0]
                            st.session_state['target_sequence'] = seq_data['sequence']
                            
                            # Display results as modern thread
                            st.markdown("---")
                            st.markdown("### 🧬 Primer Design Results")
                            ModernPrimerCard.render_primer_thread(
                                primers[0][:max_primers], 
                                show_details
                            )
            else:
                st.error("❌ Please enter a sequence")
    
    with tab2:
        st.markdown("### Comprehensive Dimer Analysis")
        
        # Manual primer input for dimer checking
        col1, col2 = st.columns(2)
        
        with col1:
            primer1 = st.text_input("Primer 1 (5' → 3')", placeholder="ATGCGATCGATCGATCG")
        
        with col2:
            primer2 = st.text_input("Primer 2 (5' → 3')", placeholder="CGTAGCTAGCTAGCTAG")
        
        if st.button("Analyze Dimers", type="primary"):
            if primer1:
                with st.spinner("Analyzing dimer formation..."):
                    # Self-dimer for primer 1
                    if primer1:
                        self_dimer1 = analyzer.analyze_primer_dimers(primer1)
                        
                        st.markdown("#### Primer 1 Self-Dimer")
                        ModernPrimerCard._display_dimer_result(self_dimer1)
                    
                    # Self-dimer for primer 2
                    if primer2:
                        self_dimer2 = analyzer.analyze_primer_dimers(primer2)
                        
                        st.markdown("#### Primer 2 Self-Dimer")
                        ModernPrimerCard._display_dimer_result(self_dimer2)
                    
                    # Heterodimer
                    if primer1 and primer2:
                        hetero_dimer = analyzer.analyze_primer_dimers(primer1, primer2)
                        
                        st.markdown("#### Heterodimer")
                        ModernPrimerCard._display_dimer_result(hetero_dimer)
                    
                    # Hairpin analysis
                    st.markdown("#### Hairpin Structures")
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        if primer1:
                            hairpin1 = analyzer.analyze_hairpin(primer1)
                            st.markdown("**Primer 1 Hairpin**")
                            if hairpin1['score'] > 0:
                                st.info(f"Stem: {hairpin1['stem_length']} bp | "
                                       f"Loop: {hairpin1['loop_length']} bp | "
                                       f"Tm: {hairpin1['tm']:.1f}°C")
                            else:
                                st.success("No significant hairpin detected")
                    
                    with col2:
                        if primer2:
                            hairpin2 = analyzer.analyze_hairpin(primer2)
                            st.markdown("**Primer 2 Hairpin**")
                            if hairpin2['score'] > 0:
                                st.info(f"Stem: {hairpin2['stem_length']} bp | "
                                       f"Loop: {hairpin2['loop_length']} bp | "
                                       f"Tm: {hairpin2['tm']:.1f}°C")
                            else:
                                st.success("No significant hairpin detected")
    
    with tab3:
        st.markdown("### Analysis Dashboard")
        
        if 'designed_primers' in st.session_state and st.session_state['designed_primers']:
            EnhancedAnalysisDashboard.render_comprehensive_analysis(
                st.session_state['designed_primers']
            )
            
            # Export options
            st.markdown("### 💾 Export Results")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                # CSV export
                export_data = []
                for i, p in enumerate(st.session_state['designed_primers']):
                    export_data.append({
                        'Set': i + 1,
                        'Forward_Seq': p.forward_seq,
                        'Reverse_Seq': p.reverse_seq,
                        'Product_Size': p.product_size,
                        'Confidence': p.confidence,
                        'Forward_Tm': p.forward_tm,
                        'Reverse_Tm': p.reverse_tm,
                        'Forward_GC': p.forward_gc,
                        'Reverse_GC': p.reverse_gc,
                        'Dimer_Severity': p.hetero_dimer.severity if p.hetero_dimer else 'N/A'
                    })
                
                df = pd.DataFrame(export_data)
                csv = df.to_csv(index=False)
                
                st.download_button(
                    "📥 Download CSV",
                    data=csv,
                    file_name="primer_results.csv",
                    mime="text/csv"
                )
            
            with col2:
                # JSON export
                json_data = json.dumps(export_data, indent=2)
                st.download_button(
                    "📥 Download JSON",
                    data=json_data,
                    file_name="primer_results.json",
                    mime="application/json"
                )
            
            with col3:
                # Lab format
                lab_text = []
                for i, p in enumerate(st.session_state['designed_primers'][:10]):
                    lab_text.append(f"Set {i+1}:")
                    lab_text.append(f"  Forward: {p.forward_seq}")
                    lab_text.append(f"  Reverse: {p.reverse_seq}")
                    lab_text.append(f"  Product: {p.product_size} bp")
                    lab_text.append(f"  Confidence: {p.confidence}%")
                    lab_text.append("")
                
                st.download_button(
                    "📥 Download Lab Format",
                    data="\n".join(lab_text),
                    file_name="primer_lab_format.txt",
                    mime="text/plain"
                )
        else:
            st.info("🔬 No primer data available. Design primers first!")
    
    with tab4:
        st.markdown("""
        ### About PrimersQuest Pro v3.0
        
        **New Features:**
        - 🔬 **Comprehensive Dimer Analysis**: Detect self-dimers, heterodimers, and hairpin structures
        - 📍 **Targeted Position Design**: Design primers for specific genomic regions
        - 🎨 **Modern Thread UI**: Social media-inspired interface for better visualization
        - 📊 **Enhanced Dashboard**: Advanced analytics and visualizations
        - 🧬 **Improved Scoring**: Updated confidence algorithm with dimer consideration
        
        **Dimer Analysis Interpretation:**
        - **Low Severity**: Minimal impact on PCR efficiency
        - **Medium Severity**: May reduce efficiency, consider alternatives
        - **High Severity**: Likely to cause issues, avoid if possible
        - **Critical**: Will significantly impact PCR, must use different primers
        
        **Best Practices:**
        1. Aim for confidence scores > 85%
        2. Keep ΔTm < 1°C between primer pairs
        3. Avoid dimers with scores > 60%
        4. Check for hairpins with Tm > 50°C
        5. Validate with NCBI Primer-BLAST
        
        **Developer:** Dr. Ahmed bey Chaker, King's College London
        """)
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #718096; padding: 1rem;'>
        PrimersQuest Pro v3.0 | Enhanced with Dimer Analysis & Modern UI<br>
        Made with ❤️ by Dr. Ahmed bey Chaker
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()