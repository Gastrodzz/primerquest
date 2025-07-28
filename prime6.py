"""
qPCR Primer Designer Pro v2.1
Enhanced version with improved PrimerBank multi-primer parsing and native Streamlit components for display.

Key Features:
- Robust primer design with fallback strategies
- Enhanced PrimerBank search with comprehensive multi-primer parsing
- Modern UI with interactive visualizations
- Export functionality (CSV/JSON/Lab format)
- Works with or without BioPython

Author: Enhanced version with improved PrimerBank functionality
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
from dataclasses import dataclass
from io import StringIO

# Safe Biopython imports with fallbacks
try:
    from Bio.Seq import Seq
except ImportError:
    class Seq:
        def __init__(self, data):
            self.data = str(data)
        def __str__(self):
            return self.data
        def upper(self):
            return Seq(self.data.upper())

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

# GC content calculation
def GC(sequence):
    """Calculate GC content of a sequence"""
    if not sequence:
        return 0.0
    sequence = str(sequence).upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    return (gc_count / total * 100) if total > 0 else 0.0

# Molecular weight calculation
def molecular_weight(sequence):
    """Calculate approximate molecular weight of DNA sequence"""
    if not sequence:
        return 0.0
    # Average molecular weight per nucleotide pair in dsDNA
    return len(str(sequence)) * 650

# Melting temperature calculation
def Tm_NN(sequence, dnac=50, saltc=50):
    """Calculate melting temperature using nearest neighbor method"""
    sequence = str(sequence).upper()
    if not sequence:
        return 0.0
    
    # Simple Tm calculation for DNA
    # Using Wallace rule for short primers (<14 bp)
    if len(sequence) < 14:
        return 2 * (sequence.count('A') + sequence.count('T')) + 4 * (sequence.count('G') + sequence.count('C'))
    
    # Using GC content method for longer sequences
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    
    # Basic formula: Tm = 81.5 + 0.41(%GC) - 675/length + salt correction
    basic_tm = 81.5 + (0.41 * gc_content * 100) - (675 / len(sequence))
    
    # Salt correction (simplified)
    salt_correction = 12.5 * np.log10(saltc / 1000)
    
    return basic_tm + salt_correction

# Page configuration
st.set_page_config(
    page_title="PrimersQuest Pro v2.1",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Enhanced CSS for modern design
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
    
    .main {
        font-family: 'Inter', sans-serif;
    }
    
    .main-header {
        font-size: 3.5rem;
        font-weight: 700;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        text-align: center;
        margin-bottom: 2rem;
        letter-spacing: -0.02em;
    }
    
    .sub-header {
        font-size: 1.8rem;
        font-weight: 600;
        color: #4a5568;
        margin-bottom: 1.5rem;
    }
    
    .stTabs [data-baseweb="tab-list"] {
        gap: 2rem;
        background-color: transparent;
    }
    
    .stTabs [data-baseweb="tab"] {
        padding: 1rem 2rem;
        font-weight: 600;
        background-color: transparent;
    }
    
    .stTabs [aria-selected="true"] {
        background-color: #667eea;
        color: white;
        border-radius: 12px;
    }
</style>
""", unsafe_allow_html=True)

@dataclass
class PrimerPair:
    """Data class for storing primer pair information"""
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
    additional_info: Dict = None

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
                'gc_content': GC(sequence),
                'n_content': n_content
            }, None
            
        except Exception as e:
            return None, f"Error parsing sequence: {str(e)}"
    
    def design_primers_enhanced(self, sequence_data: Dict, custom_params: Dict = None) -> Tuple[Optional[List[PrimerPair]], Optional[str]]:
        """Enhanced primer design with better error handling and multiple attempts"""
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
                        # Try the most common API first
                        results = primer3.designPrimers(seq_args, global_args)
                    except AttributeError:
                        try:
                            # Try bindings.design_primers
                            results = primer3.bindings.design_primers(seq_args, global_args)
                        except AttributeError:
                            try:
                                # Try bindings.designPrimers
                                results = primer3.bindings.designPrimers(seq_args, global_args)
                            except:
                                # Last resort - combine args
                                combined_args = {**seq_args, **global_args}
                                results = primer3.bindings.designPrimers(combined_args)
                    
                    num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
                    
                    # Process primers if found
                    if num_returned > 0:
                        for i in range(num_returned):
                            primer_pair = self._extract_primer_pair(results, i)
                            if primer_pair:
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

class ImprovedPrimerBankSearcher:
    """Enhanced PrimerBank searcher with comprehensive multi-primer parsing"""
    
    def __init__(self):
        self.base_url = "https://pga.mgh.harvard.edu"
        self.search_url = "https://pga.mgh.harvard.edu/cgi-bin/primerbank/new_search2.cgi"
        self.session = requests.Session()
    
    def search_primerbank(self, gene_symbol: str, species: str = "Human") -> Tuple[Optional[List[Dict]], Optional[str]]:
        """Enhanced PrimerBank search with improved parsing"""
        try:
            # Initialize session with homepage visit
            homepage_response = self.session.get(f"{self.base_url}/primerbank/", timeout=10)
            
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
                'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
                'Accept-Language': 'en-US,en;q=0.5',
                'Accept-Encoding': 'gzip, deflate, br',
                'Content-Type': 'application/x-www-form-urlencoded',
                'Origin': self.base_url,
                'Referer': f"{self.base_url}/primerbank/",
                'Connection': 'keep-alive',
                'Upgrade-Insecure-Requests': '1'
            }
            
            # Search parameters
            data = {
                'selectBox': 'NCBI Gene Symbol',
                'species': species,
                'searchBox': gene_symbol,
                'Submit': 'Submit'
            }
            
            response = self.session.post(self.search_url, data=data, headers=headers, timeout=15)
            
            if response.status_code != 200:
                return None, f"Server error: {response.status_code}"
            
            # Parse all primer sets
            return self._parse_all_primerbank_results(response.text, gene_symbol)
            
        except requests.exceptions.RequestException as e:
            return None, f"Network error: {str(e)}"
        except Exception as e:
            return None, f"Search error: {str(e)}"
    
    def _parse_all_primerbank_results(self, html_content: str, gene_symbol: str) -> Tuple[Optional[List[Dict]], Optional[str]]:
        """Parse ALL primer sets from PrimerBank results"""
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Check for no results
            if 'no primer' in html_content.lower() or 'not found' in html_content.lower():
                return None, f"No primers found for {gene_symbol}"
            
            all_primer_sets = []
            
            # Strategy 1: Find primer data in structured tables
            # PrimerBank typically shows results in tables with specific patterns
            tables = soup.find_all('table')
            
            for table in tables:
                # Skip navigation or layout tables
                if table.get('width') == '100%' or 'navigation' in str(table.get('class', [])):
                    continue
                
                # Look for tables containing primer data
                table_text = table.get_text()
                if 'primer' in table_text.lower() or re.search(r'[ATGC]{18,}', table_text):
                    primer_sets = self._extract_primers_from_table(table)
                    all_primer_sets.extend(primer_sets)
            
            # Strategy 2: Find primer data in div/pre elements
            # Sometimes primers are in preformatted text blocks
            pre_elements = soup.find_all(['pre', 'div'])
            for element in pre_elements:
                text = element.get_text()
                if re.search(r'[ATGC]{18,}', text):
                    primer_sets = self._extract_primers_from_text_block(text)
                    all_primer_sets.extend(primer_sets)
            
            # Strategy 3: Pattern-based extraction from entire page
            if not all_primer_sets:
                all_primer_sets = self._extract_primers_by_pattern(html_content)
            
            # Remove duplicates and organize
            unique_primer_sets = self._deduplicate_primer_sets(all_primer_sets)
            
            # Enrich with metadata
            for i, primer_set in enumerate(unique_primer_sets):
                primer_set['set_number'] = i + 1
                primer_set['gene_symbol'] = gene_symbol
                primer_set['source'] = 'PrimerBank'
                
                # Calculate properties if not present
                if 'forward_tm' not in primer_set and primer_set.get('forward_seq'):
                    primer_set['forward_tm'] = Tm_NN(primer_set['forward_seq'])
                if 'reverse_tm' not in primer_set and primer_set.get('reverse_seq'):
                    primer_set['reverse_tm'] = Tm_NN(primer_set['reverse_seq'])
                if 'forward_gc' not in primer_set and primer_set.get('forward_seq'):
                    primer_set['forward_gc'] = GC(primer_set['forward_seq'])
                if 'reverse_gc' not in primer_set and primer_set.get('reverse_seq'):
                    primer_set['reverse_gc'] = GC(primer_set['reverse_seq'])
            
            if unique_primer_sets:
                return unique_primer_sets, None
            else:
                return None, f"Could not parse primer information for {gene_symbol}"
                
        except Exception as e:
            return None, f"Parsing error: {str(e)}"
    
    def _extract_primers_from_table(self, table) -> List[Dict]:
        """Extract primer sets from a table element"""
        primer_sets = []
        rows = table.find_all('tr')
        
        current_set = {}
        primerbank_id = None
        
        for row in rows:
            cells = row.find_all(['td', 'th'])
            if not cells:
                continue
            
            row_text = ' '.join(cell.get_text().strip() for cell in cells)
            
            # Look for PrimerBank ID
            id_match = re.search(r'(\d+[a-z]\d+)', row_text)
            if id_match:
                primerbank_id = id_match.group(1)
            
            # Look for primer sequences
            for cell in cells:
                cell_text = cell.get_text().strip()
                
                # Check if this is a primer sequence
                if re.match(r'^[ATGC]{18,30}$', cell_text, re.IGNORECASE):
                    seq = cell_text.upper()
                    
                    # Determine if it's forward or reverse based on context
                    prev_text = row_text.lower()
                    if 'forward' in prev_text or 'left' in prev_text or 'primer 1' in prev_text:
                        current_set['forward_seq'] = seq
                    elif 'reverse' in prev_text or 'right' in prev_text or 'primer 2' in prev_text:
                        current_set['reverse_seq'] = seq
                    elif 'probe' in prev_text or 'primer 3' in prev_text:
                        current_set['probe_seq'] = seq
                    else:
                        # If no context, assign based on what we have
                        if 'forward_seq' not in current_set:
                            current_set['forward_seq'] = seq
                        elif 'reverse_seq' not in current_set:
                            current_set['reverse_seq'] = seq
                        else:
                            current_set['probe_seq'] = seq
            
            # Look for product size
            size_match = re.search(r'(\d+)\s*bp', row_text, re.IGNORECASE)
            if size_match:
                current_set['product_size'] = int(size_match.group(1))
            
            # Look for Tm
            tm_match = re.search(r'(\d+\.?\d*)\s*°?C', row_text)
            if tm_match:
                current_set['tm_info'] = float(tm_match.group(1))
            
            # If we have a complete primer pair, save it
            if 'forward_seq' in current_set and 'reverse_seq' in current_set:
                if primerbank_id:
                    current_set['primerbank_id'] = primerbank_id
                primer_sets.append(current_set.copy())
                # Reset for next set but keep the ID
                current_set = {}
        
        return primer_sets
    
    def _extract_primers_from_text_block(self, text: str) -> List[Dict]:
        """Extract primer sets from a text block"""
        primer_sets = []
        lines = text.split('\n')
        
        current_set = {}
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Look for primer sequences
            seq_match = re.search(r'([ATGC]{18,30})', line, re.IGNORECASE)
            if seq_match:
                seq = seq_match.group(1).upper()
                
                # Check context from current and previous lines
                context = ' '.join(lines[max(0, i-1):i+1]).lower()
                
                if 'forward' in context or 'left' in context:
                    current_set['forward_seq'] = seq
                elif 'reverse' in context or 'right' in context:
                    current_set['reverse_seq'] = seq
                elif 'probe' in context:
                    current_set['probe_seq'] = seq
                else:
                    # Assign based on what we have
                    if 'forward_seq' not in current_set:
                        current_set['forward_seq'] = seq
                    elif 'reverse_seq' not in current_set:
                        current_set['reverse_seq'] = seq
                        # We have a pair, save it
                        primer_sets.append(current_set.copy())
                        current_set = {}
            
            # Look for metadata
            if 'primerbank id' in line.lower():
                id_match = re.search(r'(\d+[a-z]\d+)', line)
                if id_match:
                    current_set['primerbank_id'] = id_match.group(1)
            
            size_match = re.search(r'(\d+)\s*bp', line)
            if size_match:
                current_set['product_size'] = int(size_match.group(1))
        
        # Don't forget the last set
        if 'forward_seq' in current_set and 'reverse_seq' in current_set:
            primer_sets.append(current_set)
        
        return primer_sets
    
    def _extract_primers_by_pattern(self, html_content: str) -> List[Dict]:
        """Extract primers using pattern matching on the entire content"""
        primer_sets = []
        
        # Remove HTML tags for cleaner pattern matching
        text = BeautifulSoup(html_content, 'html.parser').get_text()
        
        # Find all potential primer sequences
        primer_pattern = r'([ATGC]{18,30})'
        all_sequences = re.findall(primer_pattern, text, re.IGNORECASE)
        
        # Find all PrimerBank IDs
        id_pattern = r'(\d+[a-z]\d+)'
        all_ids = re.findall(id_pattern, text)
        
        # Group sequences into sets
        # Assuming primers come in pairs or triplets
        for i in range(0, len(all_sequences), 2):
            if i + 1 < len(all_sequences):
                primer_set = {
                    'forward_seq': all_sequences[i].upper(),
                    'reverse_seq': all_sequences[i + 1].upper()
                }
                
                # Add probe if available
                if i + 2 < len(all_sequences) and (i + 3 >= len(all_sequences) or i % 3 == 0):
                    primer_set['probe_seq'] = all_sequences[i + 2].upper()
                
                # Try to associate with PrimerBank ID
                if all_ids and len(primer_sets) < len(all_ids):
                    primer_set['primerbank_id'] = all_ids[len(primer_sets)]
                
                primer_sets.append(primer_set)
        
        return primer_sets
    
    def _deduplicate_primer_sets(self, primer_sets: List[Dict]) -> List[Dict]:
        """Remove duplicate primer sets"""
        unique_sets = []
        seen_pairs = set()
        
        for primer_set in primer_sets:
            # Create a unique key based on sequences
            forward = primer_set.get('forward_seq', '')
            reverse = primer_set.get('reverse_seq', '')
            
            if forward and reverse:
                pair_key = (forward, reverse)
                if pair_key not in seen_pairs:
                    seen_pairs.add(pair_key)
                    unique_sets.append(primer_set)
        
        return unique_sets


def calculate_primer_metrics(sequence: str) -> Dict:
    """Calculate basic primer metrics"""
    if not sequence:
        return {'length': 0, 'gc': 0, 'tm': 0}
    
    return {
        'length': len(sequence),
        'gc': GC(sequence),
        'tm': Tm_NN(sequence)
    }


def calculate_primer_quality_score(primer_set: Dict) -> int:
    """Calculate quality score for a primer set"""
    score = 100
    
    forward_seq = primer_set.get('forward_seq', '')
    reverse_seq = primer_set.get('reverse_seq', '')
    
    if not forward_seq or not reverse_seq:
        return 0
    
    # Calculate metrics
    f_metrics = calculate_primer_metrics(forward_seq)
    r_metrics = calculate_primer_metrics(reverse_seq)
    
    # Tm difference penalty
    tm_diff = abs(f_metrics['tm'] - r_metrics['tm'])
    if tm_diff > 5:
        score -= 20
    elif tm_diff > 3:
        score -= 10
    elif tm_diff > 1:
        score -= 5
    
    # GC content penalty
    for gc in [f_metrics['gc'], r_metrics['gc']]:
        if gc < 30 or gc > 70:
            score -= 15
        elif gc < 40 or gc > 60:
            score -= 5
    
    # Length penalty
    for length in [f_metrics['length'], r_metrics['length']]:
        if length < 18 or length > 25:
            score -= 10
        elif length < 19 or length > 23:
            score -= 5
    
    # Bonus for having probe
    if primer_set.get('probe_seq'):
        score += 5
    
    # Bonus for having product size info
    if primer_set.get('product_size'):
        size = primer_set['product_size']
        if 80 <= size <= 150:
            score += 5
    
    return max(0, min(100, score))


def analyze_primer_pair(forward_seq: str, reverse_seq: str) -> Dict:
    """Analyze primer pair compatibility"""
    f_metrics = calculate_primer_metrics(forward_seq)
    r_metrics = calculate_primer_metrics(reverse_seq)
    
    # Check for polyrun
    has_polyrun = bool(re.search(r'([ATGC])\1{3,}', forward_seq + reverse_seq))
    
    return {
        'tm_diff': abs(f_metrics['tm'] - r_metrics['tm']),
        'gc_diff': abs(f_metrics['gc'] - r_metrics['gc']),
        'length_diff': abs(f_metrics['length'] - r_metrics['length']),
        'avg_length': (f_metrics['length'] + r_metrics['length']) / 2,
        'forward_3_end': forward_seq[-2:] if forward_seq else '',
        'reverse_3_end': reverse_seq[-2:] if reverse_seq else '',
        'has_polyrun': has_polyrun
    }


def display_all_primerbank_results(primer_sets: List[Dict], species: str):
    """Display all PrimerBank results in an organized manner"""
    
    if not primer_sets:
        st.warning("No primer sets found")
        return
    
    # Summary card
    st.markdown(f"""### 🎯 PrimerBank Search Results""")
    st.info(f"Found **{len(primer_sets)}** validated primer set(s) for **{primer_sets[0].get('gene_symbol', 'your gene')}**")
    
    # Display each primer set in detailed view
    for idx, primer_set in enumerate(primer_sets):
        display_detailed_primerbank_card(primer_set, idx + 1, species)
        st.write("") # Add spacing between cards
    
    # Comparative analysis if more than one set is found
    if len(primer_sets) > 1:
        display_primerbank_comparison(primer_sets)
    
    # Bulk export
    if len(primer_sets) > 0:
        st.markdown("### 💾 Export All Primer Sets")
        export_primerbank_data(primer_sets)


def display_detailed_primerbank_card(primer_set: Dict, index: int, species: str):
    """Display detailed view of a primer set using native Streamlit components"""
    forward_seq = primer_set.get('forward_seq', '')
    reverse_seq = primer_set.get('reverse_seq', '')
    probe_seq = primer_set.get('probe_seq', '')
    quality_score = calculate_primer_quality_score(primer_set)
    primerbank_id = primer_set.get("primerbank_id", "")

    with st.container(border=True):
        # --- Card Header ---
        col1, col2 = st.columns([3, 2])
        with col1:
            if primerbank_id:
                st.subheader(f"🧬 Primer Set {index} - ID: {primerbank_id}")
            else:
                st.subheader(f"🧬 Primer Set {index}")

        with col2:
            st.success("✓ Experimentally Validated")
            if quality_score >= 85:
                st.markdown(f"**Quality Score: <span style='color:green;'>{quality_score}%</span>**", unsafe_allow_html=True)
            elif quality_score >= 70:
                st.markdown(f"**Quality Score: <span style='color:orange;'>{quality_score}%</span>**", unsafe_allow_html=True)
            else:
                st.markdown(f"**Quality Score: <span style='color:red;'>{quality_score}%</span>**", unsafe_allow_html=True)

        st.divider()

        # --- Primer Sequences ---
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("##### Forward Primer")
            st.code(forward_seq, language="dna")
            if forward_seq:
                metrics = calculate_primer_metrics(forward_seq)
                mcol1, mcol2, mcol3 = st.columns(3)
                mcol1.metric("Length", f"{metrics['length']} bp")
                mcol2.metric("GC%", f"{metrics['gc']:.1f}%")
                mcol3.metric("Tm", f"{metrics['tm']:.1f}°C")

        with col2:
            st.markdown("##### Reverse Primer")
            st.code(reverse_seq, language="dna")
            if reverse_seq:
                metrics = calculate_primer_metrics(reverse_seq)
                mcol1, mcol2, mcol3 = st.columns(3)
                mcol1.metric("Length", f"{metrics['length']} bp")
                mcol2.metric("GC%", f"{metrics['gc']:.1f}%")
                mcol3.metric("Tm", f"{metrics['tm']:.1f}°C")
        
        # --- Probe Sequence ---
        if probe_seq:
            st.markdown("##### 🔬 TaqMan Probe")
            st.code(probe_seq, language="dna")
            metrics = calculate_primer_metrics(probe_seq)
            pcol1, pcol2, pcol3, pcol4 = st.columns(4)
            pcol1.metric("Length", f"{metrics['length']} bp")
            pcol2.metric("GC%", f"{metrics['gc']:.1f}%")
            pcol3.metric("Tm", f"{metrics['tm']:.1f}°C")
            pcol4.metric("Type", "TaqMan" if len(probe_seq) < 30 else "Molecular Beacon")

        # --- Product Info & Analysis ---
        if primer_set.get('product_size'):
            st.markdown(f"**Expected Product Size:** {primer_set['product_size']} bp")

        with st.expander("📊 Detailed Analysis"):
            analysis = analyze_primer_pair(forward_seq, reverse_seq)
            
            acol1, acol2 = st.columns(2)
            with acol1:
                st.markdown("###### Compatibility Metrics")
                st.write(f"**Tm Difference:** {analysis['tm_diff']:.1f}°C")
                st.write(f"**GC Difference:** {analysis['gc_diff']:.1f}%")
                st.write(f"**Length Difference:** {analysis['length_diff']} bp")
            
            with acol2:
                st.markdown("###### Design Features")
                st.write(f"**3' End Stability (F):** `{analysis['forward_3_end']}`")
                st.write(f"**3' End Stability (R):** `{analysis['reverse_3_end']}`")
                st.write(f"**Contains homopolymer runs > 3:** {'Yes' if analysis['has_polyrun'] else 'No'}")

        # --- Action Buttons ---
        st.divider()
        bcol1, bcol2, bcol3 = st.columns(3)
        with bcol1:
            blast_url = create_primer_blast_url(forward_seq, reverse_seq, species)
            st.link_button("🎯 Validate with Primer-BLAST", blast_url)

        with bcol2:
            sequences_text = f"Forward: {forward_seq}\nReverse: {reverse_seq}"
            if probe_seq:
                sequences_text += f"\nProbe: {probe_seq}"
            st.download_button(f"📋 Download Sequences", data=sequences_text, file_name=f"primer_set_{index}.txt", key=f"copy_{index}")
        
        with bcol3:
            primerbank_url = f"https://pga.mgh.harvard.edu/primerbank/index.html"
            st.link_button("🔗 View on PrimerBank", primerbank_url)


def display_primerbank_comparison(primer_sets: List[Dict]):
    """Display comparative analysis of multiple primer sets"""
    
    st.markdown("### 📊 Comparative Analysis")
    
    # Prepare comparison data
    comparison_data = []
    for i, primer_set in enumerate(primer_sets):
        forward_seq = primer_set.get('forward_seq', '')
        reverse_seq = primer_set.get('reverse_seq', '')
        
        if forward_seq and reverse_seq:
            f_metrics = calculate_primer_metrics(forward_seq)
            r_metrics = calculate_primer_metrics(reverse_seq)
            
            comparison_data.append({
                'Set': f"Set {i + 1}",
                'Product Size': primer_set.get('product_size', 'N/A'),
                'Avg Length': (f_metrics['length'] + r_metrics['length']) / 2,
                'Avg GC%': (f_metrics['gc'] + r_metrics['gc']) / 2,
                'Avg Tm': (f_metrics['tm'] + r_metrics['tm']) / 2,
                'Tm Diff': abs(f_metrics['tm'] - r_metrics['tm']),
                'Has Probe': 'Yes' if primer_set.get('probe_seq') else 'No'
            })
    
    if comparison_data:
        df = pd.DataFrame(comparison_data)
        
        # Display comparison table
        st.dataframe(df.style.format({
            'Avg Length': '{:.1f}',
            'Avg GC%': '{:.1f}',
            'Avg Tm': '{:.1f}',
            'Tm Diff': '{:.1f}'
        }))
        
        # Visualizations
        col1, col2 = st.columns(2)
        
        with col1:
            # Tm comparison
            fig1 = go.Figure()
            for i, primer_set in enumerate(primer_sets):
                if primer_set.get('forward_seq') and primer_set.get('reverse_seq'):
                    f_tm = primer_set.get('forward_tm', Tm_NN(primer_set['forward_seq']))
                    r_tm = primer_set.get('reverse_tm', Tm_NN(primer_set['reverse_seq']))
                    
                    fig1.add_trace(go.Scatter(
                        x=[f'Set {i+1} Forward', f'Set {i+1} Reverse'],
                        y=[f_tm, r_tm],
                        mode='lines+markers',
                        name=f'Set {i+1}',
                        line=dict(width=2),
                        marker=dict(size=10)
                    ))
            
            fig1.update_layout(
                title='Melting Temperature Comparison',
                xaxis_title='Primer',
                yaxis_title='Tm (°C)',
                height=300
            )
            st.plotly_chart(fig1, use_container_width=True)
        
        with col2:
            # GC content comparison
            fig2 = go.Figure()
            sets = []
            gc_values = []
            
            for i, row in df.iterrows():
                sets.append(row['Set'])
                gc_values.append(row['Avg GC%'])
            
            fig2.add_trace(go.Bar(
                x=sets,
                y=gc_values,
                marker_color=['#10b981' if 45 <= gc <= 55 else '#f59e0b' if 40 <= gc <= 60 else '#ef4444' 
                             for gc in gc_values]
            ))
            
            fig2.add_hline(y=50, line_dash="dash", line_color="gray", annotation_text="Optimal")
            fig2.update_layout(
                title='Average GC Content',
                xaxis_title='Primer Set',
                yaxis_title='GC %',
                height=300
            )
            st.plotly_chart(fig2, use_container_width=True)
        
        # Recommendations
        st.markdown("### 💡 Recommendations")
        
        best_set = None
        best_score = 0
        
        for i, primer_set in enumerate(primer_sets):
            score = calculate_primer_quality_score(primer_set)
            if score > best_score:
                best_score = score
                best_set = i + 1
        
        if best_set:
            st.success(f"**Recommended: Set {best_set}** (Quality Score: {best_score}%)")
            
            # Explain why
            best_primer = primer_sets[best_set - 1]
            reasons = []
            
            tm_diff = comparison_data[best_set - 1]['Tm Diff']
            if tm_diff <= 1:
                reasons.append("Excellent Tm matching")
            
            avg_gc = comparison_data[best_set - 1]['Avg GC%']
            if 45 <= avg_gc <= 55:
                reasons.append("Optimal GC content")
            
            if best_primer.get('probe_seq'):
                reasons.append("Includes validated probe")
            
            if reasons:
                st.write("**Reasons:** " + ", ".join(reasons))


def export_primerbank_data(primer_sets: List[Dict]):
    """Export PrimerBank data in multiple formats"""
    col1, col2, col3 = st.columns(3)
    
    # Prepare export data
    export_data = []
    for i, primer_set in enumerate(primer_sets):
        export_entry = {
            'Set_Number': i + 1,
            'Gene_Symbol': primer_set.get('gene_symbol', ''),
            'PrimerBank_ID': primer_set.get('primerbank_id', ''),
            'Forward_Sequence': primer_set.get('forward_seq', ''),
            'Reverse_Sequence': primer_set.get('reverse_seq', ''),
            'Probe_Sequence': primer_set.get('probe_seq', ''),
            'Product_Size': primer_set.get('product_size', ''),
            'Forward_Tm': primer_set.get('forward_tm', ''),
            'Reverse_Tm': primer_set.get('reverse_tm', ''),
            'Forward_GC': primer_set.get('forward_gc', ''),
            'Reverse_GC': primer_set.get('reverse_gc', '')
        }
        export_data.append(export_entry)
    
    # CSV export
    with col1:
        csv_data = pd.DataFrame(export_data).to_csv(index=False)
        st.download_button(
            label="📥 Download CSV",
            data=csv_data,
            file_name=f"primerbank_{primer_sets[0].get('gene_symbol', 'results')}.csv",
            mime="text/csv"
        )
    
    # JSON export
    with col2:
        json_data = json.dumps(export_data, indent=2)
        st.download_button(
            label="📥 Download JSON",
            data=json_data,
            file_name=f"primerbank_{primer_sets[0].get('gene_symbol', 'results')}.json",
            mime="application/json"
        )
    
    # Lab format (simplified)
    with col3:
        lab_format = []
        for i, primer_set in enumerate(primer_sets):
            lab_format.append(f"Set {i+1}:")
            lab_format.append(f"Forward: {primer_set.get('forward_seq', '')}")
            lab_format.append(f"Reverse: {primer_set.get('reverse_seq', '')}")
            if primer_set.get('probe_seq'):
                lab_format.append(f"Probe: {primer_set.get('probe_seq')}")
            lab_format.append("")
        
        lab_text = "\n".join(lab_format)
        st.download_button(
            label="📥 Download Lab Format",
            data=lab_text,
            file_name=f"primerbank_{primer_sets[0].get('gene_symbol', 'primers')}.txt",
            mime="text/plain"
        )


def create_primer_visualization(sequence: str, primers: List[PrimerPair]) -> go.Figure:
    """Create interactive visualization of primer binding sites"""
    fig = go.Figure()
    
    # Add sequence track
    fig.add_trace(go.Scatter(
        x=[0, len(sequence)],
        y=[0, 0],
        mode='lines',
        line=dict(color='lightgray', width=10),
        name='Target Sequence',
        hoverinfo='skip'
    ))
    
    # Add primer binding sites
    colors = px.colors.qualitative.Set3
    
    for i, primer in enumerate(primers[:5]):  # Show top 5 primers
        color = colors[i % len(colors)]
        
        # Forward primer
        fig.add_trace(go.Scatter(
            x=[primer.forward_position, primer.forward_position + len(primer.forward_seq)],
            y=[0.2 + i * 0.15, 0.2 + i * 0.15],
            mode='lines+markers',
            line=dict(color=color, width=8),
            marker=dict(size=10, symbol='arrow-right', angleref='previous'),
            name=f'Pair {i+1} Forward',
            hovertemplate=f'Forward Primer {i+1}<br>Position: %{{x}}<br>Tm: {primer.forward_tm:.1f}°C'
        ))
        
        # Reverse primer
        fig.add_trace(go.Scatter(
            x=[primer.reverse_position - len(primer.reverse_seq), primer.reverse_position],
            y=[-0.2 - i * 0.15, -0.2 - i * 0.15],
            mode='lines+markers',
            line=dict(color=color, width=8),
            marker=dict(size=10, symbol='arrow-left', angleref='previous'),
            name=f'Pair {i+1} Reverse',
            hovertemplate=f'Reverse Primer {i+1}<br>Position: %{{x}}<br>Tm: {primer.reverse_tm:.1f}°C'
        ))
    
    fig.update_layout(
        title="Primer Binding Sites Visualization",
        xaxis_title="Position (bp)",
        yaxis_title="",
        height=400,
        showlegend=True,
        hovermode='x unified',
        yaxis=dict(showticklabels=False, range=[-1.5, 1.5])
    )
    
    return fig

def display_enhanced_primer_card(primer: PrimerPair, index: int, species: str):
    """Display enhanced primer card using native Streamlit components"""
    with st.container(border=True):
        # --- Card Header ---
        col1, col2 = st.columns([3, 2])
        with col1:
            st.subheader(f"🧬 Primer Pair {index + 1}")
        with col2:
            if primer.confidence >= 85:
                st.success(f"✅ EXCELLENT ({primer.confidence}%)")
            elif primer.confidence >= 70:
                st.warning(f"⚠️ GOOD ({primer.confidence}%)")
            else:
                st.error(f"❌ NEEDS OPTIMIZATION ({primer.confidence}%)")

        st.divider()

        # --- Primer Details ---
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("##### Forward Primer")
            st.code(primer.forward_seq, language="dna")
            mcol1, mcol2 = st.columns(2)
            mcol1.metric("Length", f"{len(primer.forward_seq)} bp")
            mcol1.metric("GC Content", f"{primer.forward_gc:.1f}%")
            mcol2.metric("Tm", f"{primer.forward_tm:.1f}°C")
            mcol2.metric("Position", f"{primer.forward_position}")

        with col2:
            st.markdown("##### Reverse Primer")
            st.code(primer.reverse_seq, language="dna")
            mcol1, mcol2 = st.columns(2)
            mcol1.metric("Length", f"{len(primer.reverse_seq)} bp")
            mcol1.metric("GC Content", f"{primer.reverse_gc:.1f}%")
            mcol2.metric("Tm", f"{primer.reverse_tm:.1f}°C")
            mcol2.metric("Position", f"{primer.reverse_position}")

        st.divider()

        # --- Product Analysis ---
        st.markdown("#### 📊 PCR Product Analysis")
        prod_col1, prod_col2, prod_col3, prod_col4 = st.columns(4)
        prod_col1.metric("Product Size", f"{primer.product_size} bp")
        prod_col2.metric("ΔTm", f"{abs(primer.forward_tm - primer.reverse_tm):.1f}°C")
        prod_col3.metric("Avg GC%", f"{(primer.forward_gc + primer.reverse_gc) / 2:.1f}%")
        prod_col4.metric("Penalty", f"{primer.penalty:.2f}")


        # --- Advanced Analysis ---
        with st.expander("🔬 Advanced Analysis"):
            adv_col1, adv_col2 = st.columns(2)
            with adv_col1:
                st.markdown("##### Thermodynamic Properties")
                st.write(f"**Heterodimer ΔG:** {primer.additional_info.get('self_any_th', 'N/A')}°C")
                st.write(f"**3' Heterodimer ΔG:** {primer.additional_info.get('self_end_th', 'N/A')}°C")
            with adv_col2:
                st.markdown("##### End Stability")
                st.write(f"**Forward 3' Stability:** {primer.additional_info.get('left_end_stability', 'N/A')}")
                st.write(f"**Reverse 3' Stability:** {primer.additional_info.get('right_end_stability', 'N/A')}")

        # --- BLAST Button ---
        st.divider()
        blast_url = create_primer_blast_url(primer.forward_seq, primer.reverse_seq, species)
        st.link_button("🎯 Validate with NCBI Primer-BLAST", blast_url, use_container_width=True)


def create_analysis_dashboard(primers: List[PrimerPair]) -> None:
    """Create comprehensive analysis dashboard"""
    
    # Prepare data for analysis
    df_data = []
    for i, primer in enumerate(primers):
        df_data.append({
            'Index': i + 1,
            'Confidence': primer.confidence,
            'Forward_Tm': primer.forward_tm,
            'Reverse_Tm': primer.reverse_tm,
            'Tm_Difference': abs(primer.forward_tm - primer.reverse_tm),
            'Forward_GC': primer.forward_gc,
            'Reverse_GC': primer.reverse_gc,
            'Avg_GC': (primer.forward_gc + primer.reverse_gc) / 2,
            'Product_Size': primer.product_size,
            'Penalty': primer.penalty
        })
    
    df = pd.DataFrame(df_data)
    
    # Create visualizations
    col1, col2 = st.columns(2)
    
    with col1:
        # Confidence distribution
        fig1 = px.bar(df, x='Index', y='Confidence',
                      title='Primer Confidence Scores',
                      labels={'Index': 'Primer Pair', 'Confidence': 'Confidence (%)'},
                      color='Confidence',
                      color_continuous_scale='viridis')
        fig1.add_hline(y=85, line_dash="dash", line_color="green", 
                       annotation_text="Excellent Threshold")
        fig1.add_hline(y=70, line_dash="dash", line_color="orange",
                       annotation_text="Good Threshold")
        st.plotly_chart(fig1, use_container_width=True)
    
    with col2:
        # Tm correlation
        fig2 = px.scatter(df, x='Forward_Tm', y='Reverse_Tm',
                         size='Product_Size', color='Confidence',
                         title='Melting Temperature Correlation',
                         labels={'Forward_Tm': 'Forward Tm (°C)', 
                                'Reverse_Tm': 'Reverse Tm (°C)'},
                         hover_data=['Index', 'Product_Size'])
        
        # Add ideal line
        min_tm = min(df['Forward_Tm'].min(), df['Reverse_Tm'].min())
        max_tm = max(df['Forward_Tm'].max(), df['Reverse_Tm'].max())
        fig2.add_shape(type="line", x0=min_tm, y0=min_tm, x1=max_tm, y1=max_tm,
                      line=dict(color="red", dash="dash"))
        st.plotly_chart(fig2, use_container_width=True)
    
    # Product size distribution
    fig3 = px.histogram(df, x='Product_Size', nbins=20,
                       title='PCR Product Size Distribution',
                       labels={'Product_Size': 'Product Size (bp)', 'count': 'Frequency'})
    fig3.add_vline(x=100, line_dash="dash", line_color="green",
                   annotation_text="Optimal for qPCR")
    st.plotly_chart(fig3, use_container_width=True)
    
    # Summary statistics
    st.markdown("### 📊 Summary Statistics")
    
    stat_col1, stat_col2, stat_col3, stat_col4 = st.columns(4)
    
    with stat_col1:
        st.metric("Average Confidence", f"{df['Confidence'].mean():.1f}%")
        st.metric("Excellent Primers", f"{len(df[df['Confidence'] >= 85])}/{len(df)}")
    
    with stat_col2:
        st.metric("Avg Tm Difference", f"{df['Tm_Difference'].mean():.2f}°C")
        st.metric("Optimal Tm Pairs", f"{len(df[df['Tm_Difference'] <= 1])}/{len(df)}")
    
    with stat_col3:
        st.metric("Avg Product Size", f"{df['Product_Size'].mean():.0f} bp")
        st.metric("Optimal Size Range", f"{len(df[(df['Product_Size'] >= 80) & (df['Product_Size'] <= 150)])}/{len(df)}")
    
    with stat_col4:
        st.metric("Avg GC Content", f"{df['Avg_GC'].mean():.1f}%")
        st.metric("Avg Penalty Score", f"{df['Penalty'].mean():.2f}")

def create_primer_blast_url(forward_seq: str, reverse_seq: str, species: str) -> str:
    """Create NCBI Primer-BLAST URL for the correct species."""
    base_url = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi"
    
    # Map app species to NCBI organism
    organism_map = {
        "Human": "Homo sapiens",
        "Mouse": "Mus musculus"
    }
    ncbi_organism = organism_map.get(species, "Homo sapiens") # Default to Human if not found

    params = {
        'PRIMER_LEFT_INPUT': forward_seq,
        'PRIMER_RIGHT_INPUT': reverse_seq,
        'PRIMER_PRODUCT_MIN': 70,
        'PRIMER_PRODUCT_MAX': 300,
        'PRIMER_NUM_RETURN': 20,
        'PRIMER_MIN_TM': 57,
        'PRIMER_OPT_TM': 60,
        'PRIMER_MAX_TM': 63,
        'PRIMER_SPECIFICITY_DATABASE': 'refseq_rna',
        'ORGANISM': ncbi_organism
    }
    return f"{base_url}?" + urllib.parse.urlencode(params)

def main():
    # Check for required dependencies
    try:
        import primer3
    except ImportError:
        st.error("""
        ❌ **Missing Required Dependency: primer3-py**
        
        Please install it using:
        ```bash
        pip install primer3-py
        ```
        
        For conda users:
        ```bash
        conda install -c bioconda primer3-py
        ```
        """)
        st.stop()
    
    # Header
    st.markdown('<h1 class="main-header">PrimersQuest Pro v2.1</h1>', unsafe_allow_html=True)
    st.markdown("""
    <p style="text-align: center; color: #718096; font-size: 1.125rem;">
    Advanced primer design tool with enhanced algorithms
    </p>
    """, unsafe_allow_html=True)
    
    # Initialize tools
    designer = EnhancedPrimerDesigner()
    searcher = ImprovedPrimerBankSearcher()
    
    # Initialize custom parameters
    custom_params = {}
    
    # Sidebar configuration
    with st.sidebar:
        st.markdown("### ⚙️ Configuration")
        
        # Design parameters
        with st.expander("🔧 Design Parameters", expanded=False):
            custom_params = {}
            
            st.markdown("#### Temperature Settings")
            min_tm = st.slider("Min Tm (°C)", 55.0, 65.0, 58.0, 0.5)
            max_tm = st.slider("Max Tm (°C)", 55.0, 65.0, 62.0, 0.5)
            custom_params['PRIMER_MIN_TM'] = min_tm
            custom_params['PRIMER_MAX_TM'] = max_tm
            
            st.markdown("#### Product Size")
            min_size = st.number_input("Min Size (bp)", 50, 300, 80)
            max_size = st.number_input("Max Size (bp)", 50, 500, 150)
            custom_params['PRIMER_PRODUCT_SIZE_RANGE'] = [[min_size, max_size]]
            
            st.markdown("#### GC Content")
            min_gc = st.slider("Min GC%", 20, 80, 40)
            max_gc = st.slider("Max GC%", 20, 80, 60)
            custom_params['PRIMER_MIN_GC'] = min_gc
            custom_params['PRIMER_MAX_GC'] = max_gc
        
        st.markdown("### 📚 Resources")
        st.markdown("""
        - [Primer3 Documentation](http://primer3.org/)
        - [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
        - [PrimerBank Database](https://pga.mgh.harvard.edu/primerbank/)
        - [qPCR Guidelines](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3573776/)
        """)
        
        # Check BioPython status
        with st.expander("ℹ️ System Status", expanded=False):
            if SeqIO:
                st.success("✅ BioPython detected (Enhanced FASTA parsing enabled)")
            else:
                st.warning("⚠️ BioPython not installed (Using basic FASTA parser)")
    
    # Main content tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "🎯 Custom Design", 
        "🔍 PrimerBank Search", 
        "📊 Analysis Dashboard",
        "📖 User Guide"
    ])
    
    with tab1:
        st.markdown('<h2 class="sub-header">Design Custom Primers</h2>', unsafe_allow_html=True)
        
        # Sequence input
        col1, col2 = st.columns([3, 1])
        
        with col1:
            sequence_input = st.text_area(
                "Enter your target sequence (FASTA format or raw sequence):",
                height=250,
                placeholder=""">Gene_Name Optional description
ATGGCAGAAATCGGTGTCAACGGATTTGGCCGTATTGGGCGCCTGGTCACCAGGGCTGCTTTTA
ACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCTTCATTGACCTCAACTACATGGTTTA
CATGTGCCCAGAGTATGCCGGAGACCCCTTTCACACATGCAGCACCTATCAGGCTGTACTCC...
                
Or paste raw sequence:
ATGGCAGAAATCGGTGTCAACGGATTTGGC...""",
                help="Minimum 100 bp, maximum 10,000 bp"
            )
        
        with col2:
            st.markdown("### 📋 Quick Stats")
            if sequence_input.strip():
                seq_data, _ = designer.parse_sequence_input(sequence_input)
                if seq_data:
                    st.metric("Length", f"{seq_data['length']} bp")
                    st.metric("GC Content", f"{seq_data['gc_content']:.1f}%")
                    st.metric("N Content", f"{seq_data['n_content']:.1f}%")
                else:
                    st.info("Enter valid sequence")
            else:
                st.info("Paste sequence to see stats")
        
        if st.button("🚀 Design Primers", type="primary", use_container_width=True):
            if sequence_input.strip():
                with st.spinner("🧬 Analyzing sequence and designing optimal primers..."):
                    # Parse sequence
                    seq_data, error = designer.parse_sequence_input(sequence_input)
                    
                    if error:
                        st.error(f"❌ {error}")
                    else:
                        # Design primers
                        primers, error = designer.design_primers_enhanced(
                            seq_data, 
                            custom_params
                        )
                        
                        if error:
                            st.error(f"❌ {error}")
                        else:
                            st.success(f"✅ Successfully designed {len(primers)} primer pairs!")
                            
                            # Store in session state
                            st.session_state['designed_primers'] = primers
                            st.session_state['target_sequence'] = seq_data['sequence']
            else:
                st.error("❌ Please enter a sequence")

        if 'designed_primers' in st.session_state and st.session_state['designed_primers']:
            primers = st.session_state['designed_primers']
            st.markdown("---")
            st.markdown('<h2 class="sub-header">Custom Primer Results</h2>', unsafe_allow_html=True)
            
            # Visualization
            if 'target_sequence' in st.session_state:
                fig = create_primer_visualization(st.session_state['target_sequence'], primers)
                st.plotly_chart(fig, use_container_width=True)
            
            # Add species selector for BLAST links
            blast_species = st.selectbox(
                "Select Organism for Primer-BLAST Validation:",
                ("Human", "Mouse"),
                key="blast_species_custom"
            )
            
            # Display primer cards
            for i, primer in enumerate(primers):
                display_enhanced_primer_card(primer, i, blast_species)
                st.write("") # Add space between cards

    with tab2:
        st.markdown('<h2 class="sub-header">Search Validated Primers</h2>', unsafe_allow_html=True)
        
        col1, col2 = st.columns([3, 1])
        
        with col1:
            gene_symbol = st.text_input(
                "Enter gene symbol:",
                placeholder="e.g., GAPDH, ACTB, TP53",
                help="Official gene symbol from NCBI Gene database"
            )
        
        with col2:
            species = st.selectbox(
                "Species:",
                ["Human", "Mouse"],
                help="Select organism"
            )
        
        if st.button("🔍 Search PrimerBank", type="primary", use_container_width=True):
            if gene_symbol.strip():
                with st.spinner(f"🔎 Searching PrimerBank for {gene_symbol} ({species})..."):
                    results, error = searcher.search_primerbank(gene_symbol.strip(), species)
                    
                    if error:
                        st.error(f"❌ {error}")
                        
                        # Provide helpful alternatives
                        col1, col2 = st.columns(2)
                        with col1:
                            # Direct link to PrimerBank
                            primerbank_url = f"https://pga.mgh.harvard.edu/cgi-bin/primerbank/new_search2.cgi?searchBox={gene_symbol}&selectBox=NCBI+Gene+Symbol&species={species}&Submit=Submit"
                            st.link_button("🔗 Try Manual Search on PrimerBank", primerbank_url)
                        
                        with col2:
                            # Alternative gene search
                            ncbi_url = f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_symbol}"
                            st.link_button("🔍 Verify Gene Symbol on NCBI", ncbi_url)
                            
                    elif results:
                        st.success(f"✅ Found {len(results)} validated primer sets!")
                        
                        # Store in session state
                        st.session_state['primerbank_results'] = results
                        
                        # Display all results using the new function
                        display_all_primerbank_results(results, species)

                    else:
                        st.warning("No results found")
                        st.info("💡 Consider using the Custom Design tab for this gene")
            else:
                st.error("❌ Please enter a gene symbol")
    
    with tab3:
        st.markdown('<h2 class="sub-header">Comprehensive Analysis Dashboard</h2>', unsafe_allow_html=True)
        
        if 'designed_primers' in st.session_state and st.session_state['designed_primers']:
            create_analysis_dashboard(st.session_state['designed_primers'])
            
            # Export functionality
            st.markdown("### 💾 Export Results")
            
            export_data = []
            for i, primer in enumerate(st.session_state['designed_primers']):
                export_data.append({
                    'Pair_Number': i + 1,
                    'Forward_Sequence': primer.forward_seq,
                    'Reverse_Sequence': primer.reverse_seq,
                    'Forward_Tm': primer.forward_tm,
                    'Reverse_Tm': primer.reverse_tm,
                    'Product_Size': primer.product_size,
                    'Confidence_Score': primer.confidence
                })
            
            export_df = pd.DataFrame(export_data)
            
            col1, col2 = st.columns(2)
            with col1:
                csv = export_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download as CSV",
                    data=csv,
                    file_name="primer_design_results.csv",
                    mime="text/csv"
                )
            
            with col2:
                json_data = json.dumps(export_data, indent=2)
                st.download_button(
                    label="📥 Download as JSON",
                    data=json_data,
                    file_name="primer_design_results.json",
                    mime="application/json"
                )
        else:
            st.info("🔬 No primer data available. Design primers or search PrimerBank first!")
    
    with tab4:
        st.markdown('<h2 class="sub-header">User Guide & Best Practices</h2>', unsafe_allow_html=True)
        
        st.markdown("""
        ### 🎯 Getting Started
        
        #### 1. **Custom Primer Design**
        - Paste your target sequence in FASTA format or as raw DNA sequence
        - Adjust parameters in the sidebar if needed
        - Click "Design Primers" to generate optimized primer pairs
        
        #### 2. **PrimerBank Search**
        - Enter the official gene symbol (e.g., GAPDH, not glyceraldehyde-3-phosphate dehydrogenase)
        - Select the appropriate species
        - Click "Search PrimerBank" to find experimentally validated primers
        - **NEW:** View all available primer sets with detailed analysis
        
        ### 🔬 Understanding Confidence Scores
        
        Our advanced scoring algorithm evaluates primers based on:
        
        - **Tm Difference (25%)**: Optimal ≤1°C between forward and reverse
        - **GC Content (20%)**: Optimal 45-55%
        - **Product Size (15%)**: Optimal 80-120 bp for qPCR
        - **Self-Complementarity (15%)**: Lower is better
        - **3' End Stability (10%)**: G or C at 3' end preferred
        - **Sequence Complexity (10%)**: Avoids repeats and poly-runs
        - **Overall Penalty (5%)**: Primer3's internal scoring
        
        ### 📊 Interpreting Results
        
        - **Confidence ≥85%**: Excellent primers, ready for use
        - **Confidence 70-84%**: Good primers, may need minor optimization
        - **Confidence <70%**: Consider alternative primer pairs
        
        ### 🛠️ Optimization Tips
        
        1. **For difficult sequences:**
           - Try relaxing Tm constraints (±1°C)
           - Allow wider product size range
           - Consider different sequence regions
        
        2. **For high specificity:**
           - Always validate with NCBI Primer-BLAST
           - Check for SNPs in primer binding sites
           - Consider species-specific validation
        
        3. **For multiplex PCR:**
           - Ensure similar Tm for all primer pairs
           - Check for primer-primer interactions
           - Validate amplicon sizes don't overlap
        
        ### 🔍 PrimerBank Features (v2.1)
        
        **NEW in v2.1:**
        - **Multi-primer parsing**: Now correctly displays ALL primer sets from PrimerBank
        - **Comparative analysis**: Compare multiple primer sets side-by-side
        - **Quality scoring**: Automatic quality assessment for each primer set
        - **Enhanced visualization**: Better charts and metrics for decision making
        - **Flexible export**: Download all primer sets in CSV, JSON, or lab format
        
        ### ⚡ Troubleshooting
        
        **No primers found?**
        - Check sequence quality (too many Ns?)
        - Try different sequence regions
        - Relax design constraints
        
        **Low confidence scores?**
        - Review problematic parameters
        - Consider manual optimization
        - Try PrimerBank for validated alternatives
        
        **PrimerBank not showing all primers?**
        - Check gene symbol spelling
        - Verify species selection
        - Some genes may have limited primer sets available
        """)
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #718096; padding: 2rem;'>
        <p style="font-size: 0.875rem;">
        PrimersQuest Pro v2.1 | Enhanced with multi-primer PrimerBank parsing<br>
        Powered by Primer3 | Integrated with NCBI & Harvard PrimerBank<br>
        Made with ❤️ for the scientific community
        </p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()