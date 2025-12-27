#!/usr/bin/env python3
"""
VarConfirm
==========
Primer design tool for Sanger sequencing validation of NGS-detected variants.

Workflow: NGS → Variant Detection → VarConfirm → Sanger Validation

Features:
- Multiple input formats (rsID, HGVS, VCF, genomic coordinates)
- BLAST-based specificity validation (local or web)
- Off-target amplicon analysis
- Gel electrophoresis simulation

Author: Burak Özkan
License: MIT
"""

import os
import sys
import re
import json
import time
import logging
from io import StringIO
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any
from datetime import datetime
from enum import Enum
import xml.etree.ElementTree as ET

import requests
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW, NCBIXML
import primer3
import myvariant


# ============== ANSI Color Codes ==============
class Colors:
    """Terminal color codes for formatted output"""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
    # Foreground colors
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    GRAY = '\033[90m'
    
    # Background colors
    BG_RED = '\033[41m'
    BG_GREEN = '\033[42m'
    BG_YELLOW = '\033[43m'
    BG_BLUE = '\033[44m'
    BG_MAGENTA = '\033[45m'
    BG_CYAN = '\033[46m'
    
    @staticmethod
    def primer_fwd(text: str) -> str:
        """Format forward primer (green background)"""
        return f"{Colors.BG_GREEN}{Colors.BOLD}{text}{Colors.RESET}"
    
    @staticmethod
    def primer_rev(text: str) -> str:
        """Format reverse primer (blue background)"""
        return f"{Colors.BG_BLUE}{Colors.BOLD}{text}{Colors.RESET}"
    
    @staticmethod
    def mutation(text: str) -> str:
        """Format mutation site (red background, white text)"""
        return f"{Colors.BG_RED}{Colors.BOLD}{Colors.WHITE}{text}{Colors.RESET}"
    
    @staticmethod
    def success(text: str) -> str:
        """Format success message (green)"""
        return f"{Colors.GREEN}{Colors.BOLD}{text}{Colors.RESET}"
    
    @staticmethod
    def error(text: str) -> str:
        """Format error message (red)"""
        return f"{Colors.RED}{Colors.BOLD}{text}{Colors.RESET}"
    
    @staticmethod
    def warning(text: str) -> str:
        """Format warning message (yellow)"""
        return f"{Colors.YELLOW}{text}{Colors.RESET}"
    
    @staticmethod
    def info(text: str) -> str:
        """Format info message (cyan)"""
        return f"{Colors.CYAN}{text}{Colors.RESET}"
    
    @staticmethod
    def header(text: str) -> str:
        """Format header (bold cyan)"""
        return f"{Colors.BOLD}{Colors.CYAN}{text}{Colors.RESET}"


# ============== 20 Known Test Variants ==============
KNOWN_TEST_VARIANTS = [
    # Cancer Hotspot Mutations (COSMIC)
    "BRAF:V600E",           # Melanoma, colorectal - most common
    "KRAS:G12D",            # Pancreatic, colorectal
    "KRAS:G12V",            # Lung, colorectal
    "KRAS:G12C",            # Lung cancer - Sotorasib target
    "EGFR:T790M",           # Lung cancer resistance mutation
    "EGFR:L858R",           # Lung cancer - TKI sensitive
    "PIK3CA:H1047R",        # Breast cancer hotspot
    "PIK3CA:E545K",         # Breast, colorectal
    "TP53:R175H",           # Pan-cancer hotspot
    "TP53:R248Q",           # Pan-cancer hotspot
    
    # rsID Variants (dbSNP/ClinVar)
    "rs121913529",          # BRAF V600E rsID
    "rs121917757",          # TP53 R175H rsID
    "rs28934576",           # BRCA1 pathogenic
    "rs80357906",           # BRCA2 pathogenic
    "rs113488022",          # BRAF V600E alternative
    
    # Genomic Coordinates (hg38)
    "chr7:140753336:A:T",   # BRAF V600E genomic
    "chr12:25245350:C:T",   # KRAS G12D genomic
    "chr17:7675088:C:T",    # TP53 R175H genomic
    
    # HGVS Coding Notation
    "NM_004333.6(BRAF):c.1799T>A",  # BRAF V600E HGVS
    "NM_033360.4(KRAS):c.35G>A",    # KRAS G12D HGVS
]


# Common variants database with coordinates (hg38)
COMMON_VARIANTS = {
    'BRAF:V600E': {'gene': 'BRAF', 'chrom': 'chr7', 'pos': 140753336, 'ref': 'A', 'alt': 'T', 'hgvs_p': 'p.Val600Glu'},
    'BRAF:V600K': {'gene': 'BRAF', 'chrom': 'chr7', 'pos': 140753335, 'ref': 'CA', 'alt': 'TT', 'hgvs_p': 'p.Val600Lys'},
    'KRAS:G12D': {'gene': 'KRAS', 'chrom': 'chr12', 'pos': 25245350, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Gly12Asp'},
    'KRAS:G12V': {'gene': 'KRAS', 'chrom': 'chr12', 'pos': 25245350, 'ref': 'C', 'alt': 'A', 'hgvs_p': 'p.Gly12Val'},
    'KRAS:G12C': {'gene': 'KRAS', 'chrom': 'chr12', 'pos': 25245350, 'ref': 'C', 'alt': 'G', 'hgvs_p': 'p.Gly12Cys'},
    'KRAS:G12R': {'gene': 'KRAS', 'chrom': 'chr12', 'pos': 25245349, 'ref': 'CC', 'alt': 'TT', 'hgvs_p': 'p.Gly12Arg'},
    'KRAS:G13D': {'gene': 'KRAS', 'chrom': 'chr12', 'pos': 25245347, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Gly13Asp'},
    'NRAS:Q61K': {'gene': 'NRAS', 'chrom': 'chr1', 'pos': 114713909, 'ref': 'G', 'alt': 'T', 'hgvs_p': 'p.Gln61Lys'},
    'NRAS:Q61R': {'gene': 'NRAS', 'chrom': 'chr1', 'pos': 114713908, 'ref': 'T', 'alt': 'C', 'hgvs_p': 'p.Gln61Arg'},
    'EGFR:T790M': {'gene': 'EGFR', 'chrom': 'chr7', 'pos': 55181378, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Thr790Met'},
    'EGFR:L858R': {'gene': 'EGFR', 'chrom': 'chr7', 'pos': 55191822, 'ref': 'T', 'alt': 'G', 'hgvs_p': 'p.Leu858Arg'},
    'EGFR:C797S': {'gene': 'EGFR', 'chrom': 'chr7', 'pos': 55181392, 'ref': 'A', 'alt': 'T', 'hgvs_p': 'p.Cys797Ser'},
    'PIK3CA:H1047R': {'gene': 'PIK3CA', 'chrom': 'chr3', 'pos': 179234297, 'ref': 'A', 'alt': 'G', 'hgvs_p': 'p.His1047Arg'},
    'PIK3CA:E545K': {'gene': 'PIK3CA', 'chrom': 'chr3', 'pos': 179218303, 'ref': 'G', 'alt': 'A', 'hgvs_p': 'p.Glu545Lys'},
    'PIK3CA:E542K': {'gene': 'PIK3CA', 'chrom': 'chr3', 'pos': 179218294, 'ref': 'G', 'alt': 'A', 'hgvs_p': 'p.Glu542Lys'},
    'TP53:R175H': {'gene': 'TP53', 'chrom': 'chr17', 'pos': 7675088, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Arg175His'},
    'TP53:R248Q': {'gene': 'TP53', 'chrom': 'chr17', 'pos': 7674220, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Arg248Gln'},
    'TP53:R248W': {'gene': 'TP53', 'chrom': 'chr17', 'pos': 7674221, 'ref': 'G', 'alt': 'T', 'hgvs_p': 'p.Arg248Trp'},
    'TP53:R273H': {'gene': 'TP53', 'chrom': 'chr17', 'pos': 7673802, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Arg273His'},
    'TP53:R273C': {'gene': 'TP53', 'chrom': 'chr17', 'pos': 7673803, 'ref': 'G', 'alt': 'A', 'hgvs_p': 'p.Arg273Cys'},
    'TP53:Y220C': {'gene': 'TP53', 'chrom': 'chr17', 'pos': 7674872, 'ref': 'A', 'alt': 'G', 'hgvs_p': 'p.Tyr220Cys'},
    'IDH1:R132H': {'gene': 'IDH1', 'chrom': 'chr2', 'pos': 208248389, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Arg132His'},
    'IDH1:R132C': {'gene': 'IDH1', 'chrom': 'chr2', 'pos': 208248388, 'ref': 'G', 'alt': 'A', 'hgvs_p': 'p.Arg132Cys'},
    'IDH2:R140Q': {'gene': 'IDH2', 'chrom': 'chr15', 'pos': 90088702, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Arg140Gln'},
    'IDH2:R172K': {'gene': 'IDH2', 'chrom': 'chr15', 'pos': 90088606, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Arg172Lys'},
    'JAK2:V617F': {'gene': 'JAK2', 'chrom': 'chr9', 'pos': 5073770, 'ref': 'G', 'alt': 'T', 'hgvs_p': 'p.Val617Phe'},
    'MPL:W515L': {'gene': 'MPL', 'chrom': 'chr1', 'pos': 43349338, 'ref': 'G', 'alt': 'T', 'hgvs_p': 'p.Trp515Leu'},
    'FLT3:D835Y': {'gene': 'FLT3', 'chrom': 'chr13', 'pos': 28018505, 'ref': 'C', 'alt': 'A', 'hgvs_p': 'p.Asp835Tyr'},
    'KIT:D816V': {'gene': 'KIT', 'chrom': 'chr4', 'pos': 54727461, 'ref': 'T', 'alt': 'A', 'hgvs_p': 'p.Asp816Val'},
    'ABL1:T315I': {'gene': 'ABL1', 'chrom': 'chr9', 'pos': 130872896, 'ref': 'G', 'alt': 'A', 'hgvs_p': 'p.Thr315Ile'},
    'MET:Y1003F': {'gene': 'MET', 'chrom': 'chr7', 'pos': 116771961, 'ref': 'A', 'alt': 'T', 'hgvs_p': 'p.Tyr1003Phe'},
    'ALK:F1174L': {'gene': 'ALK', 'chrom': 'chr2', 'pos': 29220830, 'ref': 'G', 'alt': 'A', 'hgvs_p': 'p.Phe1174Leu'},
    'RET:M918T': {'gene': 'RET', 'chrom': 'chr10', 'pos': 43114500, 'ref': 'A', 'alt': 'G', 'hgvs_p': 'p.Met918Thr'},
    'AKT1:E17K': {'gene': 'AKT1', 'chrom': 'chr14', 'pos': 104780214, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Glu17Lys'},
    'PTEN:R130Q': {'gene': 'PTEN', 'chrom': 'chr10', 'pos': 87933147, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Arg130Gln'},
    'DNMT3A:R882H': {'gene': 'DNMT3A', 'chrom': 'chr2', 'pos': 25234373, 'ref': 'C', 'alt': 'T', 'hgvs_p': 'p.Arg882His'},
    'SF3B1:K700E': {'gene': 'SF3B1', 'chrom': 'chr2', 'pos': 197402634, 'ref': 'T', 'alt': 'C', 'hgvs_p': 'p.Lys700Glu'},
    'GNAS:R201C': {'gene': 'GNAS', 'chrom': 'chr20', 'pos': 58909365, 'ref': 'G', 'alt': 'A', 'hgvs_p': 'p.Arg201Cys'},
    'HRAS:G12V': {'gene': 'HRAS', 'chrom': 'chr11', 'pos': 534288, 'ref': 'C', 'alt': 'A', 'hgvs_p': 'p.Gly12Val'},
}


# Logging setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler('primer_pipeline.log')]
)
logger = logging.getLogger(__name__)


class VariantFormat(Enum):
    """Supported variant input formats"""
    RSID = "rsid"
    HGVS_CODING = "hgvs_coding"
    HGVS_PROTEIN = "hgvs_protein"
    GENOMIC_COORD = "genomic_coord"
    STRUCTURAL = "structural"
    SEQUENCE = "sequence"
    UNKNOWN = "unknown"


class ValidationStatus(Enum):
    """Primer validation status levels"""
    PASS = "PASS"
    WARN = "WARNING"
    FAIL = "FAIL"
    PENDING = "PENDING"


@dataclass
class PipelineConfig:
    """Configuration for primer design pipeline"""
    ncbi_email: str = "user@example.com"
    ncbi_api_key: str = ""  # Optional: increases rate limit
    min_primer_size: int = 18
    opt_primer_size: int = 20
    max_primer_size: int = 25
    min_primer_tm: float = 58.0
    opt_primer_tm: float = 60.0
    max_primer_tm: float = 62.0
    min_primer_gc: float = 40.0
    max_primer_gc: float = 60.0
    product_size_min: int = 150
    product_size_max: int = 400
    flank_size: int = 300
    output_dir: str = "./results"
    cache_dir: str = "./cache"
    blast_enabled: bool = True
    blast_timeout: int = 120  # seconds
    min_specificity_score: float = 85.0
    # Local BLAST settings
    use_local_blast: bool = True  # Use local BLAST instead of web
    blast_db_path: str = os.path.expanduser("~/blast_db/human_genome")
    blast_threads: int = 4  # Number of threads for local BLAST


@dataclass
class Variant:
    """Parsed variant information"""
    id: str
    original_input: str
    format_type: VariantFormat
    chrom: str = ""
    pos: int = 0
    end_pos: int = 0
    ref: str = ""
    alt: str = ""
    variant_type: str = ""
    gene: Optional[str] = None
    transcript: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    sequence: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    is_valid: bool = True


@dataclass
class Primer:
    """Single primer information"""
    sequence: str
    start: int
    end: int
    tm: float
    gc: float
    length: int
    direction: str


@dataclass
class BlastHit:
    """BLAST alignment hit"""
    accession: str
    description: str
    chromosome: str
    start: int
    end: int
    identity: float
    alignment_length: int
    mismatches: int
    gaps: int
    e_value: float
    bit_score: float
    is_on_target: bool = False


@dataclass
class SpecificityResult:
    """Primer specificity validation result"""
    status: ValidationStatus
    score: float
    on_target_hits: int
    off_target_hits: int
    blast_hits: List[BlastHit] = field(default_factory=list)
    potential_amplicons: List[Dict] = field(default_factory=list)  # Detailed amplicon info
    fwd_binding_sites: int = 0  # Number of forward primer binding sites
    rev_binding_sites: int = 0  # Number of reverse primer binding sites
    three_prime_unique: bool = True
    details: List[str] = field(default_factory=list)
    blast_time: float = 0.0


@dataclass
class PrimerPair:
    """Primer pair with validation results"""
    variant_id: str
    forward: Primer
    reverse: Primer
    product_size: int
    product_sequence: str
    variant_position_in_product: int
    amplicon_validated: bool = False
    specificity: Optional[SpecificityResult] = None
    validation_details: Dict = field(default_factory=dict)


class VariantParser:
    """Parse multiple variant input formats"""
    
    PATTERNS = {
        'rsid': re.compile(r'^rs\d+$', re.IGNORECASE),
        'hgvs_coding_full': re.compile(
            r'^(NM_\d+(?:\.\d+)?)\(([A-Z0-9]+)\):c\.(.+)$', re.IGNORECASE
        ),
        'hgvs_coding_gene': re.compile(
            r'^([A-Z][A-Z0-9]+):c\.(.+)$', re.IGNORECASE
        ),
        'hgvs_protein': re.compile(
            r'^([A-Z][A-Z0-9]+):([A-Z])(\d+)([A-Z])$', re.IGNORECASE
        ),
        'hgvs_protein_long': re.compile(
            r'^([A-Z][A-Z0-9]+):p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|fs)$', re.IGNORECASE
        ),
        'genomic_vcf': re.compile(
            r'^(?:chr)?(\d+|[XYM]):(\d+):([ACGT]+):([ACGT]+)$', re.IGNORECASE
        ),
        'genomic_dash': re.compile(
            r'^(?:chr)?(\d+|[XYM])-(\d+)-([ACGT]*)-([ACGT]*)$', re.IGNORECASE
        ),
        'genomic_simple': re.compile(r'^(?:chr)?(\d+|[XYM]):(\d+)$'),
        'structural_del': re.compile(
            r'^(?:chr)?(\d+|[XYM]):(\d+):(\d+):DEL$', re.IGNORECASE
        ),
        'dna_sequence': re.compile(r'^[ACGTN]{30,}$', re.IGNORECASE),
    }
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.mv = myvariant.MyVariantInfo()
    
    def parse(self, variant_str: str) -> Optional[Variant]:
        """Parse variant string into Variant object"""
        variant_str = variant_str.strip()
        if not variant_str or variant_str.startswith('#'):
            return None
        
        # Check common variants first (fastest)
        upper_var = variant_str.upper()
        if upper_var in COMMON_VARIANTS:
            return self._parse_common_variant(variant_str)
        
        # Try each pattern
        if self.PATTERNS['rsid'].match(variant_str):
            return self._parse_rsid(variant_str)
        
        if self.PATTERNS['hgvs_coding_full'].match(variant_str):
            return self._parse_hgvs_coding(variant_str)
        
        if self.PATTERNS['hgvs_coding_gene'].match(variant_str):
            return self._parse_hgvs_coding(variant_str)
        
        if self.PATTERNS['hgvs_protein'].match(variant_str):
            return self._parse_protein_variant(variant_str)
        
        if self.PATTERNS['hgvs_protein_long'].match(variant_str):
            return self._parse_protein_variant(variant_str)
        
        for pattern_name in ['genomic_vcf', 'genomic_dash']:
            if self.PATTERNS[pattern_name].match(variant_str):
                return self._parse_genomic(variant_str, pattern_name)
        
        if self.PATTERNS['genomic_simple'].match(variant_str):
            return self._parse_genomic_simple(variant_str)
        
        if self.PATTERNS['structural_del'].match(variant_str):
            return self._parse_structural(variant_str)
        
        if self.PATTERNS['dna_sequence'].match(variant_str):
            return self._parse_sequence(variant_str)
        
        return Variant(
            id=variant_str, original_input=variant_str,
            format_type=VariantFormat.UNKNOWN, is_valid=False,
            warnings=["Unrecognized variant format"]
        )
    
    def _parse_common_variant(self, variant_str: str) -> Variant:
        """Parse from common variants database"""
        info = COMMON_VARIANTS[variant_str.upper()]
        return Variant(
            id=variant_str,
            original_input=variant_str,
            format_type=VariantFormat.HGVS_PROTEIN,
            chrom=info['chrom'],
            pos=info['pos'],
            ref=info['ref'],
            alt=info['alt'],
            gene=info['gene'],
            hgvs_p=info.get('hgvs_p', ''),
            variant_type='SNV' if len(info['ref']) == len(info['alt']) == 1 else 'INDEL',
            is_valid=True
        )
    
    def _parse_rsid(self, rsid: str) -> Variant:
        """Parse rsID using MyVariant.info API"""
        try:
            result = self.mv.getvariant(rsid, assembly='hg38')
            if result and isinstance(result, dict):
                chrom = str(result.get('chrom', ''))
                if chrom and not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                
                hg38 = result.get('hg38', {})
                pos = hg38.get('start', 0) if isinstance(hg38, dict) else 0
                
                gene_data = result.get('gene', {})
                gene = gene_data.get('symbol') if isinstance(gene_data, dict) else None
                
                return Variant(
                    id=rsid,
                    original_input=rsid,
                    format_type=VariantFormat.RSID,
                    chrom=chrom,
                    pos=pos,
                    ref=result.get('ref', ''),
                    alt=result.get('alt', ''),
                    gene=gene,
                    is_valid=bool(chrom and pos > 0)
                )
        except Exception as e:
            logger.warning(f"rsID lookup failed: {rsid}: {e}")
        
        return Variant(
            id=rsid, original_input=rsid,
            format_type=VariantFormat.RSID, is_valid=False,
            warnings=["Could not resolve rsID"]
        )
    
    def _parse_hgvs_coding(self, hgvs: str) -> Variant:
        """Parse HGVS coding notation"""
        match = self.PATTERNS['hgvs_coding_full'].match(hgvs)
        if match:
            transcript, gene, change = match.groups()
        else:
            match = self.PATTERNS['hgvs_coding_gene'].match(hgvs)
            if match:
                gene, change = match.groups()
                transcript = None
            else:
                return Variant(
                    id=hgvs, original_input=hgvs,
                    format_type=VariantFormat.HGVS_CODING,
                    is_valid=False, warnings=["HGVS parsing error"]
                )
        
        coords = self._hgvs_to_coords(hgvs)
        
        return Variant(
            id=hgvs,
            original_input=hgvs,
            format_type=VariantFormat.HGVS_CODING,
            chrom=coords.get('chrom', '') if coords else '',
            pos=coords.get('pos', 0) if coords else 0,
            ref=coords.get('ref', '') if coords else '',
            alt=coords.get('alt', '') if coords else '',
            gene=gene,
            transcript=transcript,
            hgvs_c=f"c.{change}",
            is_valid=bool(coords and coords.get('pos', 0) > 0),
            warnings=[] if coords else ["Could not resolve coordinates"]
        )
    
    def _parse_protein_variant(self, variant_str: str) -> Variant:
        """Parse protein variant notation"""
        upper_var = variant_str.upper()
        if upper_var in COMMON_VARIANTS:
            return self._parse_common_variant(variant_str)
        
        return Variant(
            id=variant_str, original_input=variant_str,
            format_type=VariantFormat.HGVS_PROTEIN,
            is_valid=False,
            warnings=["Protein variant not in database"]
        )
    
    def _parse_genomic(self, coord_str: str, pattern_name: str) -> Variant:
        """Parse genomic coordinates"""
        match = self.PATTERNS[pattern_name].match(coord_str)
        if match:
            groups = match.groups()
            chrom = f"chr{groups[0]}"
            pos = int(groups[1])
            ref = groups[2].upper() if groups[2] else ''
            alt = groups[3].upper() if groups[3] else ''
            
            return Variant(
                id=coord_str,
                original_input=coord_str,
                format_type=VariantFormat.GENOMIC_COORD,
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                variant_type=self._get_variant_type(ref, alt),
                is_valid=True,
                warnings=[] if (ref and alt) else ["Missing ref/alt alleles"]
            )
        return None
    
    def _parse_genomic_simple(self, coord_str: str) -> Variant:
        """Parse simple position notation"""
        match = self.PATTERNS['genomic_simple'].match(coord_str)
        if match:
            return Variant(
                id=coord_str,
                original_input=coord_str,
                format_type=VariantFormat.GENOMIC_COORD,
                chrom=f"chr{match.group(1)}",
                pos=int(match.group(2)),
                is_valid=True,
                warnings=["Position only - no ref/alt specified"]
            )
        return None
    
    def _parse_structural(self, coord_str: str) -> Variant:
        """Parse structural variant"""
        match = self.PATTERNS['structural_del'].match(coord_str)
        if match:
            return Variant(
                id=coord_str,
                original_input=coord_str,
                format_type=VariantFormat.STRUCTURAL,
                chrom=f"chr{match.group(1)}",
                pos=int(match.group(2)),
                end_pos=int(match.group(3)),
                variant_type='DEL',
                is_valid=True,
                warnings=[f"Large deletion ({int(match.group(3))-int(match.group(2))} bp)"]
            )
        return None
    
    def _parse_sequence(self, seq: str) -> Variant:
        """Parse raw DNA sequence"""
        return Variant(
            id=f"SEQ_{len(seq)}bp",
            original_input=seq,
            format_type=VariantFormat.SEQUENCE,
            sequence=seq.upper(),
            is_valid=True,
            warnings=[f"Raw sequence input ({len(seq)} bp)"]
        )
    
    def _hgvs_to_coords(self, hgvs: str) -> Optional[Dict]:
        """Convert HGVS to genomic coordinates using Ensembl VEP"""
        try:
            url = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs}"
            r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=15)
            if r.status_code == 200:
                data = r.json()
                if data:
                    item = data[0]
                    alleles = item.get('allele_string', '').split('/')
                    return {
                        'chrom': f"chr{item.get('seq_region_name', '')}",
                        'pos': item.get('start', 0),
                        'ref': alleles[0] if alleles else '',
                        'alt': alleles[1] if len(alleles) > 1 else ''
                    }
        except Exception as e:
            logger.debug(f"Ensembl VEP error: {e}")
        return None
    
    def _get_variant_type(self, ref: str, alt: str) -> str:
        """Determine variant type from ref/alt"""
        if not ref or not alt:
            return "UNKNOWN"
        if len(ref) == len(alt) == 1:
            return "SNV"
        elif len(ref) > len(alt):
            return "DEL"
        elif len(ref) < len(alt):
            return "INS"
        return "COMPLEX"


class SequenceFetcher:
    """Fetch sequences from Ensembl REST API"""
    
    def __init__(self):
        self.cache = {}
    
    def get_sequence(self, chrom: str, start: int, end: int) -> Optional[str]:
        """Fetch genomic sequence from Ensembl"""
        cache_key = f"{chrom}:{start}-{end}"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        chrom_clean = chrom.replace('chr', '')
        url = f"https://rest.ensembl.org/sequence/region/human/{chrom_clean}:{start}..{end}:1"
        
        try:
            r = requests.get(url, headers={"Content-Type": "text/plain"}, timeout=15)
            if r.status_code == 200:
                seq = r.text.upper()
                self.cache[cache_key] = seq
                return seq
        except Exception as e:
            logger.error(f"Sequence fetch error: {e}")
        return None


class BlastValidator:
    """
    BLAST-based primer specificity validation using Primer-BLAST approach.
    
    Instead of BLASTing the amplicon (which gives many false positives),
    this class BLASTs each primer separately and checks if both primers
    can bind within amplifiable distance (~3kb) at off-target sites.
    
    This is the same approach used by NCBI Primer-BLAST.
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.last_blast_time = 0
        self.min_interval = 10  # NCBI requires 10s between requests (web only)
        self.max_amplicon_size = 3000  # Max distance between primers for amplification
        
        # Check local BLAST availability
        if self.config.use_local_blast:
            self._check_local_blast()
    
    def _check_local_blast(self):
        """Verify local BLAST installation and database"""
        import subprocess
        
        # Check blastn command
        try:
            result = subprocess.run(['blastn', '-version'], capture_output=True, text=True)
            if result.returncode == 0:
                version = result.stdout.split('\n')[0]
                logger.info(f"Local BLAST found: {version}")
            else:
                logger.warning("blastn command failed, falling back to web BLAST")
                self.config.use_local_blast = False
                return
        except FileNotFoundError:
            logger.warning("blastn not found in PATH, falling back to web BLAST")
            self.config.use_local_blast = False
            return
        
        # Check database
        db_path = self.config.blast_db_path
        if not (os.path.exists(f"{db_path}.nin") or os.path.exists(f"{db_path}.nhr") or 
                os.path.exists(f"{db_path}.00.nin") or os.path.exists(f"{db_path}.00.nhr")):
            logger.warning(f"BLAST database not found at {db_path}, falling back to web BLAST")
            self.config.use_local_blast = False
    
    def validate_primers(self, primer_pair: PrimerPair, 
                        expected_chrom: str, expected_pos: int) -> SpecificityResult:
        """
        Validate primer specificity using Primer-BLAST approach.
        
        1. BLAST forward primer → get all binding sites
        2. BLAST reverse primer → get all binding sites  
        3. Find locations where BOTH primers bind within 3kb
        4. Only those are potential off-target amplicons
        
        Returns SpecificityResult with score and details
        """
        if not self.config.blast_enabled:
            return SpecificityResult(
                status=ValidationStatus.PENDING,
                score=0.0,
                on_target_hits=0,
                off_target_hits=0,
                details=["BLAST validation disabled"]
            )
        
        if self.config.use_local_blast:
            return self._run_primer_blast_local(primer_pair, expected_chrom, expected_pos)
        else:
            return self._run_primer_blast_web(primer_pair, expected_chrom, expected_pos)
    
    def _run_primer_blast_local(self, primer_pair: PrimerPair,
                                expected_chrom: str, expected_pos: int) -> SpecificityResult:
        """Run Primer-BLAST approach using local database"""
        import subprocess
        import tempfile
        
        start_time = time.time()
        
        try:
            print(Colors.info(f"    Running local Primer-BLAST validation..."))
            
            # BLAST forward primer
            print(Colors.gray(f"      Checking forward primer binding sites..."))
            fwd_hits = self._blast_primer_local(primer_pair.forward.sequence, "forward")
            
            # BLAST reverse primer (reverse complement for binding)
            print(Colors.gray(f"      Checking reverse primer binding sites..."))
            rev_hits = self._blast_primer_local(primer_pair.reverse.sequence, "reverse")
            
            blast_time = time.time() - start_time
            
            # Analyze primer pair specificity
            result = self._analyze_primer_pair_specificity(
                fwd_hits, rev_hits, expected_chrom, expected_pos, primer_pair
            )
            result.blast_time = blast_time
            
            # Add detailed report
            result.details = [
                f"Primer-BLAST completed in {blast_time:.1f}s",
                f"Forward primer binding sites: {len(fwd_hits)}",
                f"Reverse primer binding sites: {len(rev_hits)}",
                f"Potential amplicons: {result.on_target_hits + result.off_target_hits}",
                f"  - On-target: {result.on_target_hits}",
                f"  - Off-target: {result.off_target_hits}"
            ]
            
            return result
            
        except Exception as e:
            logger.error(f"Primer-BLAST error: {e}")
            return SpecificityResult(
                status=ValidationStatus.PENDING,
                score=0.0,
                on_target_hits=0,
                off_target_hits=0,
                details=[f"Primer-BLAST error: {str(e)}"]
            )
    
    def _blast_primer_local(self, primer_seq: str, primer_name: str) -> List[Dict]:
        """BLAST a single primer against local database"""
        import subprocess
        import tempfile
        
        hits = []
        
        # Create temp query file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">{primer_name}\n{primer_seq}\n")
            query_path = f.name
        
        output_path = query_path + ".tsv"
        
        try:
            # BLAST with tabular output for easier parsing
            # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand
            cmd = [
                'blastn',
                '-query', query_path,
                '-db', self.config.blast_db_path,
                '-out', output_path,
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand sseq',
                '-word_size', '7',
                '-evalue', '1000',
                '-max_target_seqs', '100',
                '-num_threads', str(self.config.blast_threads),
                '-task', 'blastn-short',
                '-dust', 'no'  # Don't mask low complexity for short primers
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            if result.returncode != 0 and result.stderr:
                logger.warning(f"BLAST warning: {result.stderr}")
            
            # Parse tabular output
            with open(output_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 13:
                        chrom = self._extract_chromosome(parts[1])
                        hit = {
                            'chromosome': chrom,
                            'identity': float(parts[2]),
                            'alignment_length': int(parts[3]),
                            'mismatches': int(parts[4]),
                            'gaps': int(parts[5]),
                            'query_start': int(parts[6]),
                            'query_end': int(parts[7]),
                            'subject_start': int(parts[8]),
                            'subject_end': int(parts[9]),
                            'e_value': float(parts[10]),
                            'bit_score': float(parts[11]),
                            'strand': parts[12],  # 'plus' or 'minus'
                            'subject_seq': parts[13] if len(parts) > 13 else '',
                            'accession': parts[1]
                        }
                        hits.append(hit)
        finally:
            # Cleanup
            try:
                os.unlink(query_path)
                os.unlink(output_path)
            except:
                pass
        
        return hits
    
    def _analyze_primer_pair_specificity(self, fwd_hits: List[Dict], rev_hits: List[Dict],
                                         expected_chrom: str, expected_pos: int,
                                         primer_pair: PrimerPair) -> SpecificityResult:
        """
        Analyze if primer pairs can form off-target amplicons.
        
        For amplification to occur:
        1. Forward primer binds on plus strand
        2. Reverse primer binds on minus strand (or vice versa)
        3. They are within max_amplicon_size bp of each other
        4. They are oriented toward each other
        """
        
        potential_amplicons = []
        on_target_count = 0
        off_target_count = 0
        scaffold_filtered = 0  # Track filtered scaffolds
        
        # Target product size for comparison
        target_size = primer_pair.product_size
        
        # Group hits by chromosome
        fwd_by_chrom = {}
        for hit in fwd_hits:
            chrom = hit['chromosome']
            if chrom not in fwd_by_chrom:
                fwd_by_chrom[chrom] = []
            fwd_by_chrom[chrom].append(hit)
        
        rev_by_chrom = {}
        for hit in rev_hits:
            chrom = hit['chromosome']
            if chrom not in rev_by_chrom:
                rev_by_chrom[chrom] = []
            rev_by_chrom[chrom].append(hit)
        
        # Find potential amplicons on each chromosome
        for chrom in set(fwd_by_chrom.keys()) & set(rev_by_chrom.keys()):
            # Check if this is a primary chromosome or scaffold
            is_primary = self._is_primary_chromosome(chrom)
            
            for fwd in fwd_by_chrom[chrom]:
                for rev in rev_by_chrom[chrom]:
                    # Check if primers are on opposite strands
                    if fwd['strand'] == rev['strand']:
                        continue  # Same strand - can't amplify
                    
                    # Calculate distance and orientation
                    if fwd['strand'] == 'plus':
                        # Forward on plus strand should be upstream of reverse
                        fwd_pos = fwd['subject_start']
                        rev_pos = rev['subject_end']
                        distance = rev_pos - fwd_pos
                    else:
                        # Forward on minus strand 
                        fwd_pos = fwd['subject_end']
                        rev_pos = rev['subject_start']
                        distance = fwd_pos - rev_pos
                    
                    # Check if within amplifiable distance
                    if 0 < distance <= self.max_amplicon_size:
                        # Calculate size ratio for risk assessment
                        size_ratio = abs(distance) / target_size if target_size > 0 else 1
                        
                        amplicon_info = {
                            'chromosome': chrom,
                            'fwd_pos': fwd['subject_start'],
                            'rev_pos': rev['subject_start'],
                            'size': abs(distance),
                            'size_ratio': size_ratio,  # How much larger than target
                            'fwd_identity': fwd['identity'],
                            'rev_identity': rev['identity'],
                            'fwd_mismatches': fwd['mismatches'],
                            'rev_mismatches': rev['mismatches'],
                            'fwd_3prime_match': self._check_3prime_match(fwd, primer_pair.forward.sequence),
                            'rev_3prime_match': self._check_3prime_match(rev, primer_pair.reverse.sequence),
                            'is_primary_chrom': is_primary,
                            'is_scaffold': not is_primary,
                            # Add alignment data for visualization
                            'fwd_subject_seq': fwd.get('subject_seq', ''),
                            'rev_subject_seq': rev.get('subject_seq', ''),
                            'fwd_query_start': fwd.get('query_start', 1),
                            'fwd_query_end': fwd.get('query_end', len(primer_pair.forward.sequence)),
                            'rev_query_start': rev.get('query_start', 1),
                            'rev_query_end': rev.get('query_end', len(primer_pair.reverse.sequence)),
                            'fwd_alignment_len': fwd.get('alignment_length', len(primer_pair.forward.sequence)),
                            'rev_alignment_len': rev.get('alignment_length', len(primer_pair.reverse.sequence)),
                            'fwd_gaps': fwd.get('gaps', 0),
                            'rev_gaps': rev.get('gaps', 0),
                            'fwd_strand': fwd.get('strand', 'plus'),
                            'rev_strand': rev.get('strand', 'minus')
                        }
                        
                        # Check if on-target
                        is_on_target = self._is_on_target(
                            chrom, fwd['subject_start'], rev['subject_start'],
                            expected_chrom, expected_pos
                        )
                        amplicon_info['is_on_target'] = is_on_target
                        
                        # Filter scaffolds for off-targets but keep for reporting
                        if not is_on_target and not is_primary:
                            scaffold_filtered += 1
                            amplicon_info['filtered'] = True
                        else:
                            amplicon_info['filtered'] = False
                        
                        potential_amplicons.append(amplicon_info)
                        
                        if is_on_target:
                            on_target_count += 1
                        elif is_primary:  # Only count primary chromosome off-targets
                            off_target_count += 1
        
        # Store scaffold count for reporting
        self._last_scaffold_filtered = scaffold_filtered
        
        # Calculate specificity score based on potential off-target amplicons
        score = self._calculate_primer_pair_score(potential_amplicons, on_target_count, off_target_count)
        
        # Count actual high-risk amplicons (will truly amplify)
        high_risk_count = sum(1 for amp in potential_amplicons 
                              if not amp['is_on_target'] 
                              and amp['fwd_identity'] >= 95 
                              and amp['rev_identity'] >= 95
                              and amp.get('fwd_3prime_match', 0) >= 5
                              and amp.get('rev_3prime_match', 0) >= 5)
        
        # Determine status based on ACTUAL risk, not just presence of off-targets
        if high_risk_count == 0 and on_target_count >= 1:
            status = ValidationStatus.PASS
        elif high_risk_count <= 2:
            status = ValidationStatus.WARN
        else:
            status = ValidationStatus.FAIL
        
        # Create detailed blast hits for reporting
        blast_hits = []
        for amp in potential_amplicons[:10]:  # Top 10 for report
            hit = BlastHit(
                accession=amp['chromosome'],
                description=f"Potential amplicon: {amp['size']}bp",
                chromosome=amp['chromosome'],
                start=amp['fwd_pos'],
                end=amp['rev_pos'],
                identity=(amp['fwd_identity'] + amp['rev_identity']) / 2,
                alignment_length=amp['size'],
                mismatches=amp['fwd_mismatches'] + amp['rev_mismatches'],
                gaps=0,
                e_value=0,
                bit_score=0,
                is_on_target=amp['is_on_target']
            )
            blast_hits.append(hit)
        
        # Check 3' uniqueness
        three_prime_unique = all(
            amp['fwd_3prime_match'] < 6 or amp['rev_3prime_match'] < 6 
            for amp in potential_amplicons if not amp['is_on_target']
        ) if potential_amplicons else True
        
        # Count binding sites
        fwd_binding_count = sum(len(hits) for hits in fwd_by_chrom.values())
        rev_binding_count = sum(len(hits) for hits in rev_by_chrom.values())
        
        return SpecificityResult(
            status=status,
            score=score,
            on_target_hits=on_target_count,
            off_target_hits=off_target_count,
            blast_hits=blast_hits,
            potential_amplicons=potential_amplicons,  # Full list for visualization
            fwd_binding_sites=fwd_binding_count,
            rev_binding_sites=rev_binding_count,
            three_prime_unique=three_prime_unique,
            details=[],  # Will be filled by caller
            blast_time=0  # Will be filled by caller
        )
    
    def _check_3prime_match(self, hit: Dict, primer_seq: str) -> int:
        """Check how many 3' bases match perfectly"""
        if 'subject_seq' not in hit or not hit['subject_seq']:
            return hit['alignment_length'] - hit['mismatches']
        
        # For primer binding, 3' end is critical
        # Count matching bases from 3' end
        subject = hit['subject_seq'].upper()
        primer = primer_seq.upper()
        
        # Align from 3' end
        matches = 0
        for i in range(1, min(len(subject), len(primer), 8) + 1):
            if primer[-i] == subject[-i]:
                matches += 1
            else:
                break
        return matches
    
    def _calculate_primer_pair_score(self, amplicons: List[Dict], 
                                     on_target: int, off_target: int) -> float:
        """
        Calculate specificity score based on potential amplicons.
        
        Scoring focuses on ACTUAL amplification risk:
        - Only high-identity sites with good 3' match are real risks
        - Low-identity sites won't amplify and shouldn't penalize heavily
        - Score represents confidence level, not strict pass/fail
        - Larger off-target products are less competitive with smaller target
        - Scaffolds/unknown sequences are not penalized (not in reference genome)
        """
        score = 100.0
        high_risk_count = 0
        
        for amp in amplicons:
            if amp['is_on_target']:
                continue
            
            # Skip scaffolds - they're not relevant for analysis
            if amp.get('is_scaffold', False) or amp.get('filtered', False):
                continue
            
            fwd_id = amp['fwd_identity']
            rev_id = amp['rev_identity']
            fwd_3p = amp.get('fwd_3prime_match', 0)
            rev_3p = amp.get('rev_3prime_match', 0)
            size_ratio = amp.get('size_ratio', 1.0)
            
            # Size-based penalty reduction
            # If off-target is much larger than target, it won't compete well
            # Target always wins in PCR when it's shorter
            size_penalty_factor = 1.0
            if size_ratio > 3.0:
                size_penalty_factor = 0.1  # 10x larger = negligible risk
            elif size_ratio > 2.0:
                size_penalty_factor = 0.3  # 2-3x larger = low risk
            elif size_ratio > 1.5:
                size_penalty_factor = 0.6  # 1.5-2x larger = moderate reduction
            
            # HIGH RISK: Both primers >95% AND good 3' match - WILL amplify
            if fwd_id >= 95 and rev_id >= 95 and fwd_3p >= 5 and rev_3p >= 5:
                penalty = 15 * size_penalty_factor
                score -= penalty
                if size_penalty_factor >= 0.6:  # Only count if similar size
                    high_risk_count += 1
            # MEDIUM RISK: Both primers >90% - MAY amplify
            elif fwd_id >= 90 and rev_id >= 90:
                if fwd_3p >= 5 and rev_3p >= 5:
                    score -= 8 * size_penalty_factor
                else:
                    score -= 4 * size_penalty_factor
            # LOW RISK: 85-90% identity - unlikely to amplify
            elif fwd_id >= 85 and rev_id >= 85:
                score -= 1 * size_penalty_factor
            # NEGLIGIBLE: <85% - won't amplify, no penalty
        
        # Bonus for clean target
        if on_target >= 1 and high_risk_count == 0:
            score = min(100, score + 10)
        
        return max(0, score)
    
    def _run_primer_blast_web(self, primer_pair: PrimerPair,
                              expected_chrom: str, expected_pos: int) -> SpecificityResult:
        """Run Primer-BLAST approach using NCBI web service"""
        # Rate limiting
        self._wait_for_rate_limit()
        
        start_time = time.time()
        
        try:
            print(Colors.info(f"    Running web Primer-BLAST validation (this may take 1-2 minutes)..."))
            
            # BLAST forward primer
            print(Colors.gray(f"      Checking forward primer..."))
            fwd_hits = self._blast_primer_web(primer_pair.forward.sequence)
            self._wait_for_rate_limit()
            
            # BLAST reverse primer
            print(Colors.gray(f"      Checking reverse primer..."))
            rev_hits = self._blast_primer_web(primer_pair.reverse.sequence)
            
            blast_time = time.time() - start_time
            
            # Convert web hits to local format
            fwd_hits_converted = self._convert_web_hits(fwd_hits)
            rev_hits_converted = self._convert_web_hits(rev_hits)
            
            # Analyze
            result = self._analyze_primer_pair_specificity(
                fwd_hits_converted, rev_hits_converted, 
                expected_chrom, expected_pos, primer_pair
            )
            result.blast_time = blast_time
            result.details = [
                f"Web Primer-BLAST completed in {blast_time:.1f}s",
                f"Forward primer hits: {len(fwd_hits_converted)}",
                f"Reverse primer hits: {len(rev_hits_converted)}",
                f"Potential off-target amplicons: {result.off_target_hits}"
            ]
            
            return result
            
        except Exception as e:
            logger.error(f"Web Primer-BLAST error: {e}")
            return SpecificityResult(
                status=ValidationStatus.PENDING,
                score=0.0,
                on_target_hits=0,
                off_target_hits=0,
                details=[f"Web Primer-BLAST error: {str(e)}"]
            )
    
    def _blast_primer_web(self, primer_seq: str) -> List:
        """BLAST primer using NCBI web service"""
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="refseq_genomes",
            sequence=primer_seq,
            entrez_query="Homo sapiens[organism]",
            word_size=7,
            expect=1000,
            hitlist_size=100,
            format_type="XML"
        )
        self.last_blast_time = time.time()
        
        blast_records = list(NCBIXML.parse(result_handle))
        return blast_records
    
    def _convert_web_hits(self, blast_records: List) -> List[Dict]:
        """Convert NCBIXML hits to local format"""
        hits = []
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    chrom = self._extract_chromosome(alignment.title)
                    strand = 'plus' if hsp.sbjct_start < hsp.sbjct_end else 'minus'
                    hit = {
                        'chromosome': chrom,
                        'identity': (hsp.identities / hsp.align_length) * 100,
                        'alignment_length': hsp.align_length,
                        'mismatches': hsp.align_length - hsp.identities,
                        'gaps': hsp.gaps,
                        'query_start': hsp.query_start,
                        'query_end': hsp.query_end,
                        'subject_start': min(hsp.sbjct_start, hsp.sbjct_end),
                        'subject_end': max(hsp.sbjct_start, hsp.sbjct_end),
                        'e_value': hsp.expect,
                        'bit_score': hsp.bits,
                        'strand': strand,
                        'accession': alignment.accession
                    }
                    hits.append(hit)
        return hits
    
    def _wait_for_rate_limit(self):
        """Enforce NCBI rate limiting (web BLAST only)"""
        elapsed = time.time() - self.last_blast_time
        if elapsed < self.min_interval:
            wait_time = self.min_interval - elapsed
            print(Colors.gray(f"    Waiting {wait_time:.1f}s for NCBI rate limit..."))
            time.sleep(wait_time)
    
    def _extract_chromosome(self, description: str) -> str:
        """Extract chromosome from BLAST hit description"""
        # Look for chromosome patterns - handle RefSeq format
        patterns = [
            r'chromosome (\d+|X|Y|M|MT)',
            r'chr(\d+|X|Y|M|MT)',
            r'NC_0000(\d+)\.',  # NC_000001.11 = chr1
            r'NC_0000(\d+)',
        ]
        
        for pattern in patterns:
            match = re.search(pattern, description, re.IGNORECASE)
            if match:
                chrom = match.group(1)
                # Handle NC accession numbers
                if chrom.isdigit():
                    num = int(chrom)
                    if num == 23:
                        return "chrX"
                    elif num == 24:
                        return "chrY"
                    elif num == 12920:  # MT
                        return "chrM"
                    elif num <= 22:
                        return f"chr{num}"
                return f"chr{chrom}"
        
        return "unknown"
    
    def _is_on_target(self, hit_chrom: str, hit_start: int, hit_end: int,
                     expected_chrom: str, expected_pos: int) -> bool:
        """Check if hit is on expected target"""
        # Normalize chromosome names
        hit_chrom = hit_chrom.lower().replace('chr', '')
        exp_chrom = expected_chrom.lower().replace('chr', '')
        
        if hit_chrom != exp_chrom:
            return False
        
        # Check if position is within expected range (±1000 bp)
        tolerance = 1000
        hit_min = min(hit_start, hit_end)
        hit_max = max(hit_start, hit_end)
        return (hit_min <= expected_pos + tolerance) and (hit_max >= expected_pos - tolerance)
    
    def _is_primary_chromosome(self, chrom: str) -> bool:
        """
        Check if chromosome is a primary reference chromosome.
        
        Filters out:
        - Unplaced scaffolds (chrUn_*)
        - Alt haplotypes (*_alt, *_hap*)
        - Random chromosomes (*_random)
        - Unknown/unlocalized contigs
        - Patch sequences (*_fix, *_novel)
        """
        chrom_lower = chrom.lower()
        
        # List of valid primary chromosomes
        primary_chroms = {
            'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
            'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
            'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
            'chrx', 'chry', 'chrm', 'chrmt',
            # Also accept without 'chr' prefix
            '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
            '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
            'x', 'y', 'm', 'mt'
        }
        
        # Check if it's a primary chromosome
        if chrom_lower in primary_chroms:
            return True
        
        # Filter patterns for non-primary sequences
        filter_patterns = [
            'unknown', 'un_', 'random', '_alt', '_hap', 
            '_fix', '_novel', 'scaffold', 'contig', 'patch'
        ]
        
        for pattern in filter_patterns:
            if pattern in chrom_lower:
                return False
        
        # If it starts with 'chr' followed by just a number or X/Y/M, it's primary
        if chrom_lower.startswith('chr'):
            suffix = chrom_lower[3:]
            if suffix.isdigit() or suffix in ['x', 'y', 'm', 'mt']:
                return True
        
        return False


class InSilicoPCR:
    """In silico PCR validation"""
    
    @staticmethod
    def validate_primers(forward_seq: str, reverse_seq: str,
                        template: str, expected_size: int) -> Dict:
        """Validate primer binding and product formation"""
        
        result = {
            'valid': False,
            'forward_binds': False,
            'reverse_binds': False,
            'forward_pos': -1,
            'reverse_pos': -1,
            'actual_size': 0,
            'size_match': False,
            'details': []
        }
        
        template = template.upper()
        forward_seq = forward_seq.upper()
        reverse_seq = reverse_seq.upper()
        
        # Check forward primer binding
        fwd_pos = template.find(forward_seq)
        if fwd_pos >= 0:
            result['forward_binds'] = True
            result['forward_pos'] = fwd_pos
        else:
            # Try fuzzy match
            fwd_pos = InSilicoPCR._fuzzy_find(template, forward_seq, max_mismatch=2)
            if fwd_pos >= 0:
                result['forward_binds'] = True
                result['forward_pos'] = fwd_pos
                result['details'].append("Forward primer: 1-2 mismatches allowed")
        
        # Check reverse primer binding (reverse complement)
        rev_comp = str(Seq(reverse_seq).reverse_complement())
        rev_pos = template.find(rev_comp)
        if rev_pos >= 0:
            result['reverse_binds'] = True
            result['reverse_pos'] = rev_pos
        else:
            rev_pos = InSilicoPCR._fuzzy_find(template, rev_comp, max_mismatch=2)
            if rev_pos >= 0:
                result['reverse_binds'] = True
                result['reverse_pos'] = rev_pos
                result['details'].append("Reverse primer: 1-2 mismatches allowed")
        
        # Calculate actual product size
        if result['forward_binds'] and result['reverse_binds']:
            result['actual_size'] = result['reverse_pos'] + len(reverse_seq) - result['forward_pos']
            
            # Check size tolerance (10%)
            size_diff = abs(result['actual_size'] - expected_size)
            result['size_match'] = size_diff <= expected_size * 0.1
            result['valid'] = result['size_match']
            
            if result['size_match']:
                result['details'].append(f"Product size: {result['actual_size']} bp (expected: {expected_size})")
            else:
                result['details'].append(f"Size mismatch: {result['actual_size']} vs {expected_size}")
        
        return result
    
    @staticmethod
    def _fuzzy_find(text: str, pattern: str, max_mismatch: int = 2) -> int:
        """Find pattern allowing mismatches"""
        for i in range(len(text) - len(pattern) + 1):
            mismatches = sum(1 for a, b in zip(text[i:i+len(pattern)], pattern) if a != b)
            if mismatches <= max_mismatch:
                return i
        return -1


class GelSimulator:
    """Simulate agarose gel electrophoresis results"""
    
    @staticmethod
    def simulate(target_size: int, off_target_amplicons: List[Dict], 
                 gel_length: int = 25, max_bands: int = 10) -> str:
        """
        Generate ASCII art gel simulation.
        
        Shows expected bands based on PCR product sizes.
        Smaller products migrate faster (appear lower on gel).
        """
        
        output = []
        output.append("")
        output.append(Colors.header("=" * 50))
        output.append(Colors.header("  GEL ELECTROPHORESIS SIMULATION"))
        output.append(Colors.header("=" * 50))
        output.append("")
        
        # Collect all products with their risk levels
        products = []
        
        # Add target product
        products.append({
            'size': target_size,
            'label': 'TARGET',
            'type': 'target',
            'intensity': 'strong'
        })
        
        # Add off-target products (only primary chromosomes, sorted by risk)
        for amp in off_target_amplicons:
            if amp.get('filtered', False) or amp.get('is_scaffold', False):
                continue
            if amp.get('is_on_target', False):
                continue
                
            # Determine intensity based on identity and 3' match
            fwd_id = amp.get('fwd_identity', 0)
            rev_id = amp.get('rev_identity', 0)
            fwd_3p = amp.get('fwd_3prime_match', 0)
            rev_3p = amp.get('rev_3prime_match', 0)
            
            if fwd_id >= 95 and rev_id >= 95 and fwd_3p >= 5 and rev_3p >= 5:
                intensity = 'strong'
                risk = 'HIGH'
            elif fwd_id >= 90 and rev_id >= 90:
                intensity = 'medium'
                risk = 'MED'
            else:
                intensity = 'weak'
                risk = 'LOW'
            
            products.append({
                'size': amp['size'],
                'label': f"{amp['chromosome']}:{amp['size']}bp",
                'type': 'off-target',
                'intensity': intensity,
                'risk': risk
            })
        
        # Limit to max_bands
        products = products[:max_bands]
        
        if not products:
            output.append("  No products to display")
            return "\n".join(output)
        
        # Define size range for gel
        all_sizes = [p['size'] for p in products]
        min_size = min(all_sizes)
        max_size = max(all_sizes)
        
        # Standard ladder positions (100bp ladder)
        ladder = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000]
        
        # Gel header
        output.append("         Ladder    Sample")
        output.append("         ------    ------")
        output.append("    (+)  │    │    │    │")
        
        # Calculate band positions using log scale (more realistic)
        import math
        
        def size_to_position(size, gel_len=gel_length):
            # Log scale: smaller = further traveled (higher position number)
            if size <= 0:
                return gel_len - 1
            log_size = math.log10(size)
            log_min = math.log10(50)  # Minimum visible size
            log_max = math.log10(5000)  # Maximum visible size
            # Invert: smaller size = larger position (migrated further)
            pos = int((log_max - log_size) / (log_max - log_min) * (gel_len - 1))
            return max(0, min(gel_len - 1, pos))
        
        # Create gel lanes
        ladder_lane = [' '] * gel_length
        sample_lane = [' '] * gel_length
        band_labels = [''] * gel_length
        
        # Add ladder bands
        for lsize in ladder:
            if 50 <= lsize <= 5000:
                pos = size_to_position(lsize)
                ladder_lane[pos] = '═'
        
        # Add sample bands
        for prod in products:
            pos = size_to_position(prod['size'])
            
            if prod['type'] == 'target':
                sample_lane[pos] = '█'  # Strong target band
                band_labels[pos] = f" ← {Colors.success(str(prod['size']) + 'bp TARGET')}"
            else:
                if prod['intensity'] == 'strong':
                    sample_lane[pos] = '█'
                    band_labels[pos] = f" ← {Colors.error(str(prod['size']) + 'bp')} ({prod.get('risk', '')})"
                elif prod['intensity'] == 'medium':
                    sample_lane[pos] = '▓'
                    band_labels[pos] = f" ← {Colors.warning(str(prod['size']) + 'bp')} ({prod.get('risk', '')})"
                else:
                    sample_lane[pos] = '░'
                    band_labels[pos] = f" ← {prod['size']}bp (weak)"
        
        # Draw gel
        for i in range(gel_length):
            # Find ladder size at this position
            ladder_label = ""
            for lsize in ladder:
                if size_to_position(lsize) == i:
                    ladder_label = f"{lsize:>4}"
                    break
            
            ladder_char = ladder_lane[i] if ladder_lane[i] != ' ' else '│'
            sample_char = sample_lane[i] if sample_lane[i] != ' ' else '│'
            
            line = f"  {ladder_label:>4}  │ {ladder_char}  │    │ {sample_char}  │{band_labels[i]}"
            output.append(line)
        
        output.append("    (-)  └────┘    └────┘")
        output.append("")
        
        # Legend
        output.append("  Legend:")
        output.append(f"    █ = Strong band (will see on gel)")
        output.append(f"    ▓ = Medium band (may see under optimal conditions)")
        output.append(f"    ░ = Weak band (unlikely to see)")
        output.append("")
        
        # Interpretation
        output.append(Colors.bold("  INTERPRETATION:"))
        
        target_only = all(p['type'] == 'target' for p in products if p['intensity'] == 'strong')
        
        if target_only:
            output.append(f"    {Colors.success('✓')} Clean gel expected - only target band visible")
        else:
            strong_offtargets = [p for p in products if p['type'] == 'off-target' and p['intensity'] == 'strong']
            if strong_offtargets:
                output.append(f"    {Colors.warning('⚠')} Multiple bands expected!")
                for p in strong_offtargets[:3]:
                    size_diff = abs(p['size'] - target_size)
                    if size_diff < 50:
                        output.append(f"       {Colors.error('!')} {p['size']}bp band very close to target - hard to distinguish!")
                    else:
                        output.append(f"       • {p['size']}bp off-target band")
        
        output.append("")
        output.append(Colors.header("=" * 50))
        output.append("")
        
        return "\n".join(output)


class PrimerBindingAnalyzer:
    """Analyze and visualize primer binding details for amplicons"""
    
    @staticmethod
    def analyze_binding(primer_seq: str, template_seq: str, is_reverse: bool = False) -> Dict:
        """
        Analyze how a primer binds to a template sequence.
        
        Returns detailed binding information including:
        - Match/mismatch positions
        - 3' end match count
        - Binding strength estimate
        """
        primer = primer_seq.upper()
        template = template_seq.upper()
        
        if is_reverse:
            # For reverse primer, we compare to reverse complement
            primer = str(Seq(primer).reverse_complement())
        
        # Find best alignment position
        best_pos = -1
        best_matches = 0
        best_alignment = None
        
        for i in range(len(template) - len(primer) + 1):
            region = template[i:i+len(primer)]
            matches = sum(1 for a, b in zip(primer, region) if a == b)
            if matches > best_matches:
                best_matches = matches
                best_pos = i
                best_alignment = region
        
        if best_pos == -1:
            return {
                'bound': False,
                'identity': 0,
                'matches': 0,
                'mismatches': len(primer),
                'primer_len': len(primer),
                '3prime_matches': 0,
                'alignment': None
            }
        
        # Analyze the alignment
        alignment_details = []
        mismatches = []
        for j, (p, t) in enumerate(zip(primer, best_alignment)):
            if p == t:
                alignment_details.append({'pos': j, 'primer': p, 'template': t, 'match': True})
            else:
                alignment_details.append({'pos': j, 'primer': p, 'template': t, 'match': False})
                mismatches.append(j)
        
        # Count 3' end matches (critical for extension)
        three_prime_matches = 0
        for j in range(len(primer) - 1, -1, -1):
            if primer[j] == best_alignment[j]:
                three_prime_matches += 1
            else:
                break
        
        # Calculate binding strength
        identity = (best_matches / len(primer)) * 100
        
        # Binding strength factors:
        # - Overall identity
        # - 3' end matches (most critical)
        # - Mismatch positions (3' mismatches are worse)
        binding_strength = identity * 0.6  # Base on identity
        binding_strength += (three_prime_matches / len(primer)) * 40  # 3' end bonus
        
        # Penalty for 3' mismatches
        for mm_pos in mismatches:
            distance_from_3prime = len(primer) - 1 - mm_pos
            if distance_from_3prime < 3:
                binding_strength -= 15  # Strong penalty for last 3 bases
            elif distance_from_3prime < 6:
                binding_strength -= 5   # Moderate penalty for bases 4-6
        
        return {
            'bound': True,
            'identity': identity,
            'matches': best_matches,
            'mismatches': len(primer) - best_matches,
            'mismatch_positions': mismatches,
            'primer_len': len(primer),
            '3prime_matches': three_prime_matches,
            'binding_strength': max(0, min(100, binding_strength)),
            'template_pos': best_pos,
            'alignment': alignment_details,
            'primer_seq': primer_seq,
            'aligned_template': best_alignment
        }
    
    @staticmethod
    def visualize_binding(fwd_analysis: Dict, rev_analysis: Dict, 
                         amplicon_info: Dict, target_size: int,
                         forward_primer: str, reverse_primer: str) -> str:
        """Generate BLAST-style visualization of primer binding for an amplicon"""
        
        lines = []
        
        chrom = amplicon_info.get('chromosome', 'unknown')
        size = amplicon_info.get('size', 0)
        is_on_target = amplicon_info.get('is_on_target', False)
        
        # Header
        if is_on_target:
            lines.append(Colors.success(f"  ┌{'─'*64}┐"))
            lines.append(Colors.success(f"  │  ON-TARGET AMPLICON: {chrom}:{amplicon_info.get('fwd_pos', '?'):<30}  │"))
            lines.append(Colors.success(f"  └{'─'*64}┘"))
        else:
            lines.append(Colors.warning(f"  ┌{'─'*64}┐"))
            lines.append(Colors.warning(f"  │  OFF-TARGET AMPLICON: {chrom}:{amplicon_info.get('fwd_pos', '?'):<29}  │"))
            lines.append(Colors.warning(f"  └{'─'*64}┘"))
        
        # Amplicon size comparison
        size_ratio = size / target_size if target_size > 0 else 1
        lines.append("")
        lines.append(f"  Product Size: {size} bp  (Target: {target_size} bp, Ratio: {size_ratio:.1f}x)")
        
        if size_ratio > 3:
            lines.append(f"  {Colors.success('✓')} Much larger than target - will be out-competed in PCR")
        elif size_ratio > 2:
            lines.append(f"  {Colors.info('○')} Larger product - reduced amplification efficiency")
        elif size_ratio > 1.5:
            lines.append(f"  {Colors.warning('⚠')} Moderately larger - may still amplify")
        elif 0.7 <= size_ratio <= 1.3:
            lines.append(f"  {Colors.error('!')} Similar size to target - will compete!")
        else:
            lines.append(f"  {Colors.error('!')} Smaller than target - may out-compete!")
        
        lines.append("")
        lines.append("  " + "─" * 64)
        
        # Forward primer BLAST-style alignment
        lines.append("")
        lines.append(Colors.bold("  FORWARD PRIMER:"))
        fwd_subject = amplicon_info.get('fwd_subject_seq', '')
        fwd_alignment = PrimerBindingAnalyzer._create_blast_alignment(
            forward_primer, 
            fwd_subject,
            amplicon_info.get('fwd_query_start', 1),
            amplicon_info.get('fwd_query_end', len(forward_primer)),
            amplicon_info.get('fwd_identity', 0),
            amplicon_info.get('fwd_mismatches', 0),
            amplicon_info.get('fwd_gaps', 0),
            amplicon_info.get('fwd_3prime_match', 0)
        )
        lines.append(fwd_alignment)
        
        lines.append("")
        lines.append("  " + "─" * 64)
        
        # Reverse primer BLAST-style alignment  
        lines.append("")
        lines.append(Colors.bold("  REVERSE PRIMER:"))
        rev_subject = amplicon_info.get('rev_subject_seq', '')
        rev_alignment = PrimerBindingAnalyzer._create_blast_alignment(
            reverse_primer,
            rev_subject,
            amplicon_info.get('rev_query_start', 1),
            amplicon_info.get('rev_query_end', len(reverse_primer)),
            amplicon_info.get('rev_identity', 0),
            amplicon_info.get('rev_mismatches', 0),
            amplicon_info.get('rev_gaps', 0),
            amplicon_info.get('rev_3prime_match', 0)
        )
        lines.append(rev_alignment)
        
        lines.append("")
        lines.append("  " + "─" * 64)
        
        # Overall amplification prediction
        lines.append("")
        lines.append(Colors.bold("  AMPLIFICATION PREDICTION:"))
        
        fwd_id = amplicon_info.get('fwd_identity', 0)
        rev_id = amplicon_info.get('rev_identity', 0)
        fwd_3p = amplicon_info.get('fwd_3prime_match', 0)
        rev_3p = amplicon_info.get('rev_3prime_match', 0)
        
        # Determine risk level
        if fwd_id >= 95 and rev_id >= 95 and fwd_3p >= 5 and rev_3p >= 5 and size_ratio <= 1.5:
            lines.append(f"  {Colors.error('█████')} HIGH RISK - Will likely amplify!")
            lines.append(f"       • Both primers: >95% identity")
            lines.append(f"       • 3' ends: well-matched ({fwd_3p}bp, {rev_3p}bp)")
            lines.append(f"       • Size: competitive with target")
        elif fwd_id >= 90 and rev_id >= 90:
            if size_ratio > 2:
                lines.append(f"  {Colors.warning('███░░')} MEDIUM RISK - May amplify weakly")
                lines.append(f"       • Good primer binding but larger product")
                lines.append(f"       • Target ({target_size}bp) will out-compete")
            else:
                lines.append(f"  {Colors.warning('████░')} MEDIUM-HIGH RISK")
                lines.append(f"       • Both primers bind well")
                lines.append(f"       • Similar size - may see on gel")
        elif fwd_id >= 85 or rev_id >= 85:
            lines.append(f"  {Colors.info('██░░░')} LOW RISK - Unlikely to amplify")
            lines.append(f"       • Primer binding is suboptimal")
        else:
            lines.append(f"  {Colors.gray('█░░░░')} NEGLIGIBLE - Will not amplify")
            lines.append(f"       • Poor primer binding (<85%)")
        
        lines.append("")
        
        return "\n".join(lines)
    
    @staticmethod
    def _create_blast_alignment(primer_seq: str, subject_seq: str,
                                 query_start: int, query_end: int,
                                 identity: float, mismatches: int, gaps: int,
                                 three_prime_match: int) -> str:
        """Create BLAST-style alignment visualization"""
        
        lines = []
        primer_len = len(primer_seq)
        
        # Stats line
        lines.append(f"  Identity: {identity:.1f}% ({primer_len - mismatches}/{primer_len}), "
                    f"Mismatches: {mismatches}, Gaps: {gaps}")
        lines.append(f"  3' end match: {three_prime_match} bp " + 
                    (Colors.success("(good)") if three_prime_match >= 5 else 
                     Colors.warning("(moderate)") if three_prime_match >= 3 else 
                     Colors.error("(poor - extension blocked!)")))
        lines.append("")
        
        # If we have subject sequence, show alignment
        if subject_seq and len(subject_seq) >= 10:
            # Build alignment
            query_line = f"  Query  {query_start:>3}  "
            match_line = f"             "
            sbjct_line = f"  Sbjct       "
            
            # Pad subject if needed
            subject_padded = subject_seq.upper()
            primer_upper = primer_seq.upper()
            
            # Ensure same length
            min_len = min(len(primer_upper), len(subject_padded))
            
            mismatch_positions = []
            for i in range(min_len):
                p = primer_upper[i]
                s = subject_padded[i] if i < len(subject_padded) else '-'
                
                if p == s:
                    query_line += Colors.success(p)
                    match_line += "|"
                    sbjct_line += Colors.success(s)
                else:
                    query_line += Colors.error(p)
                    match_line += " "
                    sbjct_line += Colors.error(s)
                    mismatch_positions.append(i)
            
            # Add remaining primer bases if any
            for i in range(min_len, len(primer_upper)):
                query_line += Colors.gray(primer_upper[i])
                match_line += " "
                sbjct_line += Colors.gray("-")
            
            query_line += f"  {query_end}"
            
            lines.append(query_line)
            lines.append(match_line)
            lines.append(sbjct_line)
            
            # Show 3' end indicator
            lines.append("")
            three_prime_indicator = "  " + " " * 13
            for i in range(min_len):
                if i >= min_len - three_prime_match:
                    three_prime_indicator += Colors.success("↑")
                elif i in mismatch_positions:
                    three_prime_indicator += Colors.error("×")
                else:
                    three_prime_indicator += " "
            
            lines.append(three_prime_indicator)
            
            # 3' region box
            if three_prime_match >= 5:
                lines.append(f"  {'':>13}{' ' * (min_len - 6)}└─3' end─┘ {Colors.success('OK')}")
            elif three_prime_match >= 3:
                lines.append(f"  {'':>13}{' ' * (min_len - 6)}└─3' end─┘ {Colors.warning('Weak')}")
            else:
                lines.append(f"  {'':>13}{' ' * (min_len - 6)}└─3' end─┘ {Colors.error('BLOCKED')}")
            
            # Explain mismatches
            if mismatch_positions:
                lines.append("")
                mm_desc = []
                for pos in mismatch_positions:
                    dist_from_3p = primer_len - 1 - pos
                    if dist_from_3p < 3:
                        mm_desc.append(f"pos {pos+1} " + Colors.error("(3' critical!)"))
                    elif dist_from_3p < 6:
                        mm_desc.append(f"pos {pos+1} " + Colors.warning("(3' region)"))
                    else:
                        mm_desc.append(f"pos {pos+1}")
                lines.append(f"  Mismatch positions: {', '.join(mm_desc)}")
        
        else:
            # No subject sequence - show stats only
            lines.append(f"  Query  1    {primer_seq}  {primer_len}")
            lines.append(f"              {'|' * (primer_len - mismatches)}{'X' * mismatches}")
            lines.append(f"  Sbjct       {'?' * primer_len}")
            lines.append("")
            lines.append(f"  (Subject sequence not available - showing estimate)")
        
        return "\n".join(lines)
    
    @staticmethod
    def generate_detailed_report(amplicons: List[Dict], 
                                  forward_primer: str, 
                                  reverse_primer: str,
                                  target_size: int,
                                  max_display: int = 5) -> str:
        """Generate detailed binding analysis for top amplicons with BLAST-style alignments"""
        
        lines = []
        lines.append("")
        lines.append(Colors.header("=" * 70))
        lines.append(Colors.header("  DETAILED PRIMER BINDING ANALYSIS"))
        lines.append(Colors.header("=" * 70))
        lines.append("")
        
        # Primer info
        lines.append(Colors.bold("PRIMERS:"))
        lines.append(f"  Forward (5'→3'): {Colors.primer_fwd(forward_primer)} ({len(forward_primer)}bp)")
        lines.append(f"  Reverse (5'→3'): {Colors.primer_rev(reverse_primer)} ({len(reverse_primer)}bp)")
        lines.append(f"  Target Product:  {target_size} bp")
        lines.append("")
        
        # Analyze on-target first
        on_targets = [a for a in amplicons if a.get('is_on_target', False)]
        off_targets = [a for a in amplicons if not a.get('is_on_target', False) 
                       and a.get('is_primary_chrom', True)]
        
        # Sort off-targets by risk (identity)
        off_targets_sorted = sorted(off_targets, 
                                    key=lambda x: -(x.get('fwd_identity', 0) + x.get('rev_identity', 0)))
        
        total_to_show = len(on_targets) + min(max_display, len(off_targets_sorted))
        lines.append(Colors.bold(f"ALIGNMENT ANALYSIS ({total_to_show} sites):"))
        lines.append("")
        
        # Analyze on-target
        for amp in on_targets:
            lines.append(PrimerBindingAnalyzer.visualize_binding(
                {}, {},  # Old analysis dicts not needed anymore
                amp, 
                target_size,
                forward_primer,
                reverse_primer
            ))
            lines.append("")
        
        # Analyze top off-targets
        for i, amp in enumerate(off_targets_sorted[:max_display], 1):
            lines.append(f"  Off-target #{i}:")
            lines.append(PrimerBindingAnalyzer.visualize_binding(
                {}, {},
                amp,
                target_size,
                forward_primer,
                reverse_primer
            ))
            lines.append("")
        
        if len(off_targets_sorted) > max_display:
            remaining = len(off_targets_sorted) - max_display
            lines.append(f"  ... and {remaining} more off-target sites (not shown)")
            lines.append("")
        
        lines.append(Colors.header("=" * 70))
        lines.append("")
        
        return "\n".join(lines)


class ProductVisualizer:
    """Visualize PCR product with colored primers and mutation"""
    
    @staticmethod
    def visualize(product_seq: str, forward_primer: str, reverse_primer: str,
                 mutation_pos: int, ref: str = "", alt: str = "",
                 line_width: int = 60) -> str:
        """
        Create colored visualization of PCR product
        
        Colors:
        - Green background: Forward primer binding site
        - Blue background: Reverse primer binding site
        - Red background: Mutation position
        """
        
        output_lines = []
        output_lines.append("")
        output_lines.append(Colors.header("=" * 70))
        output_lines.append(Colors.header("  PCR PRODUCT VISUALIZATION"))
        output_lines.append(Colors.header("=" * 70))
        output_lines.append("")
        
        # Legend
        output_lines.append("  Color Legend:")
        output_lines.append(f"    {Colors.primer_fwd('XXXXX')} = Forward Primer Binding Site")
        output_lines.append(f"    {Colors.primer_rev('XXXXX')} = Reverse Primer Binding Site")
        output_lines.append(f"    {Colors.mutation('X')} = Mutation Position")
        output_lines.append("")
        
        # Find primer positions
        fwd_start = 0
        fwd_end = len(forward_primer)
        
        rev_comp = str(Seq(reverse_primer).reverse_complement())
        rev_start = product_seq.upper().find(rev_comp.upper())
        if rev_start == -1:
            rev_start = len(product_seq) - len(reverse_primer)
        rev_end = rev_start + len(reverse_primer)
        
        # Product info
        output_lines.append(f"  Product Length: {len(product_seq)} bp")
        output_lines.append(f"  Mutation Position: base {mutation_pos + 1}")
        if ref and alt:
            output_lines.append(f"  Variant: {ref} > {alt}")
        output_lines.append("")
        output_lines.append("  Sequence:")
        output_lines.append("  " + "-" * 66)
        
        # Position ruler
        ruler = "  "
        for i in range(0, len(product_seq), 10):
            ruler += f"{i+1:<10}"
        output_lines.append(Colors.warning(ruler[:68]))
        
        # Build colored sequence line by line
        current_line = "  "
        char_count = 0
        
        for i, base in enumerate(product_seq):
            # Determine color
            if i >= fwd_start and i < fwd_end:
                if i == mutation_pos:
                    current_line += Colors.mutation(base)
                else:
                    current_line += Colors.primer_fwd(base)
            elif i >= rev_start and i < rev_end:
                if i == mutation_pos:
                    current_line += Colors.mutation(base)
                else:
                    current_line += Colors.primer_rev(base)
            elif i == mutation_pos:
                current_line += Colors.mutation(base)
            else:
                current_line += base
            
            char_count += 1
            
            if char_count >= line_width:
                output_lines.append(current_line)
                current_line = "  "
                char_count = 0
        
        if char_count > 0:
            output_lines.append(current_line)
        
        output_lines.append("  " + "-" * 66)
        output_lines.append("")
        
        # Primer details
        output_lines.append("  Primer Sequences:")
        output_lines.append(f"    Forward (5'→3'): {Colors.primer_fwd(forward_primer)}")
        output_lines.append(f"    Reverse (5'→3'): {Colors.primer_rev(reverse_primer)}")
        output_lines.append("")
        
        return "\n".join(output_lines)
    
    @staticmethod
    def text_report(product_seq: str, forward_primer: str, reverse_primer: str,
                   mutation_pos: int, ref: str = "", alt: str = "") -> str:
        """Plain text report for file output"""
        
        lines = []
        lines.append("=" * 70)
        lines.append("PCR PRODUCT DETAILS")
        lines.append("=" * 70)
        lines.append(f"Product Length: {len(product_seq)} bp")
        lines.append(f"Mutation Position: base {mutation_pos + 1}")
        if ref and alt:
            lines.append(f"Variant: {ref} > {alt}")
        lines.append("")
        lines.append(f"Forward Primer: {forward_primer}")
        lines.append(f"Reverse Primer: {reverse_primer}")
        lines.append("")
        lines.append("Sequence (mutation marked with [X]):")
        lines.append("-" * 70)
        
        # Mark mutation
        marked_seq = list(product_seq)
        if 0 <= mutation_pos < len(marked_seq):
            marked_seq[mutation_pos] = f"[{marked_seq[mutation_pos]}]"
        
        seq_str = "".join(marked_seq)
        
        # Wrap sequence
        for i in range(0, len(seq_str), 60):
            lines.append(seq_str[i:i+60])
        
        lines.append("-" * 70)
        
        return "\n".join(lines)


class PrimerDesigner:
    """Design primers using Primer3"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.seq_fetcher = SequenceFetcher()
    
    def design(self, variant: Variant) -> Optional[PrimerPair]:
        """Design primers for a variant"""
        if not variant.is_valid:
            return None
        
        if variant.format_type == VariantFormat.SEQUENCE:
            return self._design_from_sequence(variant)
        
        if not variant.chrom or variant.pos <= 0:
            return None
        
        # Fetch flanking sequence
        start = max(1, variant.pos - self.config.flank_size)
        end = variant.pos + self.config.flank_size
        
        sequence = self.seq_fetcher.get_sequence(variant.chrom, start, end)
        if not sequence:
            return None
        
        return self._run_primer3(variant, sequence, start)
    
    def _design_from_sequence(self, variant: Variant) -> Optional[PrimerPair]:
        """Design primers from raw sequence"""
        sequence = variant.sequence
        if len(sequence) < 100:
            return None
        target_pos = len(sequence) // 2
        return self._run_primer3_on_seq(variant.id, sequence, target_pos, 0)
    
    def _run_primer3(self, variant: Variant, sequence: str, seq_start: int) -> Optional[PrimerPair]:
        """Run Primer3 with variant information"""
        target_pos = variant.pos - seq_start
        return self._run_primer3_on_seq(variant.id, sequence, target_pos, seq_start)
    
    def _run_primer3_on_seq(self, var_id: str, sequence: str,
                           target_pos: int, seq_start: int) -> Optional[PrimerPair]:
        """Execute Primer3 design"""
        
        seq_args = {
            'SEQUENCE_ID': var_id,
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_TARGET': [max(0, target_pos - 30), 60],
        }
        
        global_args = {
            'PRIMER_OPT_SIZE': self.config.opt_primer_size,
            'PRIMER_MIN_SIZE': self.config.min_primer_size,
            'PRIMER_MAX_SIZE': self.config.max_primer_size,
            'PRIMER_OPT_TM': self.config.opt_primer_tm,
            'PRIMER_MIN_TM': self.config.min_primer_tm,
            'PRIMER_MAX_TM': self.config.max_primer_tm,
            'PRIMER_MIN_GC': self.config.min_primer_gc,
            'PRIMER_MAX_GC': self.config.max_primer_gc,
            'PRIMER_PRODUCT_SIZE_RANGE': [[self.config.product_size_min, self.config.product_size_max]],
            'PRIMER_NUM_RETURN': 5,
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_MAX_NS_ACCEPTED': 0,
        }
        
        try:
            results = primer3.design_primers(seq_args, global_args)
            
            if results.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                left = results.get('PRIMER_LEFT_0')
                right = results.get('PRIMER_RIGHT_0')
                
                if left and right:
                    fwd = Primer(
                        sequence=results['PRIMER_LEFT_0_SEQUENCE'],
                        start=left[0],
                        end=left[0] + left[1],
                        tm=results['PRIMER_LEFT_0_TM'],
                        gc=results['PRIMER_LEFT_0_GC_PERCENT'],
                        length=left[1],
                        direction='forward'
                    )
                    
                    rev = Primer(
                        sequence=results['PRIMER_RIGHT_0_SEQUENCE'],
                        start=right[0] - right[1] + 1,
                        end=right[0] + 1,
                        tm=results['PRIMER_RIGHT_0_TM'],
                        gc=results['PRIMER_RIGHT_0_GC_PERCENT'],
                        length=right[1],
                        direction='reverse'
                    )
                    
                    product_size = results['PRIMER_PAIR_0_PRODUCT_SIZE']
                    product_seq = sequence[fwd.start:rev.end]
                    mutation_in_product = target_pos - fwd.start
                    
                    # In silico PCR validation
                    validation = InSilicoPCR.validate_primers(
                        fwd.sequence, rev.sequence,
                        sequence, product_size
                    )
                    
                    return PrimerPair(
                        variant_id=var_id,
                        forward=fwd,
                        reverse=rev,
                        product_size=product_size,
                        product_sequence=product_seq,
                        variant_position_in_product=mutation_in_product,
                        amplicon_validated=validation['valid'],
                        validation_details=validation
                    )
                    
        except Exception as e:
            logger.error(f"Primer3 error: {e}")
        
        return None


class ValidationReport:
    """Generate validation reports"""
    
    @staticmethod
    def generate_report(variant: Variant, primer_pair: PrimerPair,
                       specificity: Optional[SpecificityResult] = None) -> str:
        """Generate formatted validation report"""
        
        lines = []
        lines.append("")
        lines.append(Colors.header("═" * 70))
        lines.append(Colors.header("  PRIMER SPECIFICITY VALIDATION REPORT"))
        lines.append(Colors.header("  Specificity Assessment"))
        lines.append(Colors.header("═" * 70))
        lines.append("")
        
        # Variant Information
        lines.append(Colors.bold("VARIANT INFORMATION:"))
        lines.append(f"  ID: {variant.id}")
        lines.append(f"  Gene: {variant.gene or 'N/A'}")
        lines.append(f"  Location: {variant.chrom}:{variant.pos}")
        if variant.ref and variant.alt:
            lines.append(f"  Change: {variant.ref} > {variant.alt}")
        if variant.hgvs_p:
            lines.append(f"  Protein: {variant.hgvs_p}")
        lines.append("")
        
        # Primer Information
        lines.append(Colors.bold("DESIGNED PRIMERS:"))
        lines.append(f"  Forward: {Colors.primer_fwd(primer_pair.forward.sequence)}")
        lines.append(f"           Tm={primer_pair.forward.tm:.1f}°C, GC={primer_pair.forward.gc:.0f}%, Length={primer_pair.forward.length}bp")
        lines.append(f"  Reverse: {Colors.primer_rev(primer_pair.reverse.sequence)}")
        lines.append(f"           Tm={primer_pair.reverse.tm:.1f}°C, GC={primer_pair.reverse.gc:.0f}%, Length={primer_pair.reverse.length}bp")
        lines.append(f"  Product Size: {primer_pair.product_size} bp")
        lines.append("")
        
        # Validation Results
        lines.append(Colors.bold("VALIDATION RESULTS:"))
        
        # In silico PCR
        if primer_pair.amplicon_validated:
            lines.append(f"  {Colors.success('✓')} In Silico PCR: PASS")
        else:
            lines.append(f"  {Colors.warning('⚠')} In Silico PCR: Check Required")
        
        # BLAST Specificity
        if specificity:
            if specificity.status == ValidationStatus.PASS:
                status_icon = Colors.success('✓')
                status_text = "PASS"
            elif specificity.status == ValidationStatus.WARN:
                status_icon = Colors.warning('⚠')
                status_text = "WARNING"
            elif specificity.status == ValidationStatus.FAIL:
                status_icon = Colors.error('✗')
                status_text = "FAIL"
            else:
                status_icon = Colors.info('○')
                status_text = "PENDING"
            
            lines.append(f"  {status_icon} BLAST Specificity: {status_text}")
            lines.append(f"      - On-target hits: {specificity.on_target_hits}")
            lines.append(f"      - Off-target hits: {specificity.off_target_hits}")
            if specificity.blast_time > 0:
                lines.append(f"      - Analysis time: {specificity.blast_time:.1f}s")
            
            # 3' End Check
            if specificity.three_prime_unique:
                lines.append(f"  {Colors.success('✓')} 3' End Uniqueness: PASS")
            else:
                lines.append(f"  {Colors.warning('⚠')} 3' End Uniqueness: Multiple matches")
            
            lines.append("")
            
            # Final Score
            lines.append(Colors.header("═" * 70))
            lines.append(f"  SPECIFICITY SCORE: {Colors.bold(f'{specificity.score:.1f}%')}")
            
            if specificity.status == ValidationStatus.PASS:
                lines.append(f"  STATUS: {Colors.success('✓ VALIDATED - HIGH SPECIFICITY')}")
            elif specificity.status == ValidationStatus.WARN:
                lines.append(f"  STATUS: {Colors.warning('⚠ REVIEW RECOMMENDED')}")
            else:
                lines.append(f"  STATUS: {Colors.error('✗ LOW SPECIFICITY - REVIEW REQUIRED')}")
            
            lines.append(Colors.header("═" * 70))
        else:
            lines.append(f"  {Colors.info('○')} BLAST Specificity: Not performed")
        
        lines.append("")
        
        return "\n".join(lines)
    
    @staticmethod
    def generate_text_report(variant: Variant, primer_pair: PrimerPair,
                            specificity: Optional[SpecificityResult] = None) -> str:
        """Generate plain text report for file output"""
        
        lines = []
        lines.append("=" * 70)
        lines.append("PRIMER SPECIFICITY VALIDATION REPORT")
        lines.append("Specificity Assessment")
        lines.append("=" * 70)
        lines.append("")
        
        lines.append("VARIANT INFORMATION:")
        lines.append(f"  ID: {variant.id}")
        lines.append(f"  Gene: {variant.gene or 'N/A'}")
        lines.append(f"  Location: {variant.chrom}:{variant.pos}")
        if variant.ref and variant.alt:
            lines.append(f"  Change: {variant.ref} > {variant.alt}")
        lines.append("")
        
        lines.append("PRIMERS:")
        lines.append(f"  Forward: {primer_pair.forward.sequence}")
        lines.append(f"           Tm={primer_pair.forward.tm:.1f}C, GC={primer_pair.forward.gc:.0f}%")
        lines.append(f"  Reverse: {primer_pair.reverse.sequence}")
        lines.append(f"           Tm={primer_pair.reverse.tm:.1f}C, GC={primer_pair.reverse.gc:.0f}%")
        lines.append(f"  Product: {primer_pair.product_size} bp")
        lines.append("")
        
        lines.append("VALIDATION:")
        lines.append(f"  In Silico PCR: {'PASS' if primer_pair.amplicon_validated else 'CHECK'}")
        
        if specificity:
            lines.append(f"  BLAST Status: {specificity.status.value}")
            lines.append(f"  Specificity Score: {specificity.score:.1f}%")
            lines.append(f"  On-target hits: {specificity.on_target_hits}")
            lines.append(f"  Off-target hits: {specificity.off_target_hits}")
        
        lines.append("")
        lines.append("=" * 70)
        
        return "\n".join(lines)
    
    @staticmethod
    def visualize_blast_results(specificity: SpecificityResult, 
                                 variant: Variant,
                                 primer_pair: PrimerPair) -> str:
        """Generate detailed BLAST results visualization"""
        
        lines = []
        lines.append("")
        lines.append(Colors.header("=" * 70))
        lines.append(Colors.header("  BLAST SPECIFICITY ANALYSIS - DETAILED RESULTS"))
        lines.append(Colors.header("=" * 70))
        lines.append("")
        
        # Summary Statistics
        lines.append(Colors.bold("PRIMER BINDING SITE ANALYSIS:"))
        lines.append(f"  Forward primer: {Colors.primer_fwd(primer_pair.forward.sequence)}")
        lines.append(f"    → Found {specificity.fwd_binding_sites} binding sites in genome")
        lines.append(f"  Reverse primer: {Colors.primer_rev(primer_pair.reverse.sequence)}")
        lines.append(f"    → Found {specificity.rev_binding_sites} binding sites in genome")
        lines.append("")
        
        # Amplification Analysis
        lines.append(Colors.bold("POTENTIAL AMPLIFICATION PRODUCTS:"))
        total_amplicons = specificity.on_target_hits + specificity.off_target_hits
        lines.append(f"  Total potential amplicons: {total_amplicons}")
        lines.append(f"    {Colors.success('✓')} On-target (expected): {specificity.on_target_hits}")
        if specificity.off_target_hits > 0:
            lines.append(f"    {Colors.error('✗')} Off-target (unintended): {specificity.off_target_hits}")
        else:
            lines.append(f"    {Colors.success('✓')} Off-target (unintended): 0")
        lines.append("")
        
        # Chromosome distribution
        if specificity.potential_amplicons:
            chrom_counts = {}
            for amp in specificity.potential_amplicons:
                chrom = amp['chromosome']
                if chrom not in chrom_counts:
                    chrom_counts[chrom] = {'on_target': 0, 'off_target': 0}
                if amp['is_on_target']:
                    chrom_counts[chrom]['on_target'] += 1
                else:
                    chrom_counts[chrom]['off_target'] += 1
            
            lines.append(Colors.bold("CHROMOSOME DISTRIBUTION:"))
            lines.append("  " + "-" * 50)
            lines.append(f"  {'Chromosome':<12} {'On-Target':<12} {'Off-Target':<12} {'Risk':<10}")
            lines.append("  " + "-" * 50)
            
            for chrom in sorted(chrom_counts.keys(), key=lambda x: (0 if x == variant.chrom else 1, x)):
                counts = chrom_counts[chrom]
                on_t = counts['on_target']
                off_t = counts['off_target']
                
                if chrom == variant.chrom:
                    chrom_display = f"{chrom} ←"
                else:
                    chrom_display = chrom
                
                if off_t == 0:
                    risk = Colors.success("None")
                elif off_t <= 2:
                    risk = Colors.warning("Low")
                else:
                    risk = Colors.error("High")
                
                lines.append(f"  {chrom_display:<12} {on_t:<12} {off_t:<12} {risk}")
            
            lines.append("  " + "-" * 50)
            lines.append("")
        
        # Count scaffolds
        all_off_targets = [a for a in specificity.potential_amplicons if not a['is_on_target']]
        scaffold_count = sum(1 for a in all_off_targets if a.get('is_scaffold', False) or not a.get('is_primary_chrom', True))
        
        # Detailed off-target list (top 10) - PRIMARY CHROMOSOMES ONLY
        off_targets = [a for a in specificity.potential_amplicons 
                       if not a['is_on_target'] and a.get('is_primary_chrom', True)]
        
        if scaffold_count > 0:
            lines.append(Colors.gray(f"  Note: {scaffold_count} off-targets on scaffolds/alt-haplotypes filtered out"))
            lines.append("")
        
        if off_targets:
            lines.append(Colors.bold(f"OFF-TARGET AMPLICONS - PRIMARY CHROMOSOMES (Top 10):"))
            lines.append("  " + "-" * 72)
            lines.append(f"  {'#':<3} {'Location':<22} {'Size':<8} {'Ratio':<6} {'Fwd %':<7} {'Rev %':<7} {'3p Risk':<8}")
            lines.append("  " + "-" * 72)
            
            # Sort by risk (high identity first, then by size ratio)
            off_targets_sorted = sorted(off_targets, 
                                        key=lambda x: (-(x['fwd_identity'] + x['rev_identity']), x.get('size_ratio', 1)))
            
            target_size = primer_pair.product_size
            
            for i, amp in enumerate(off_targets_sorted[:10], 1):
                loc = f"{amp['chromosome']}:{amp['fwd_pos']}"
                if len(loc) > 21:
                    loc = loc[:18] + "..."
                
                fwd_id = f"{amp['fwd_identity']:.0f}%"
                rev_id = f"{amp['rev_identity']:.0f}%"
                size_ratio = amp.get('size_ratio', amp['size'] / target_size if target_size > 0 else 1)
                ratio_str = f"{size_ratio:.1f}x"
                
                # 3' risk assessment
                fwd_3p = amp.get('fwd_3prime_match', 0)
                rev_3p = amp.get('rev_3prime_match', 0)
                if fwd_3p >= 6 and rev_3p >= 6:
                    risk_3p = Colors.error("HIGH")
                elif fwd_3p >= 4 or rev_3p >= 4:
                    risk_3p = Colors.warning("MED")
                else:
                    risk_3p = Colors.success("LOW")
                
                # Size-based note
                if size_ratio > 3:
                    size_note = Colors.gray(" (too large)")
                elif size_ratio > 2:
                    size_note = Colors.info(" (large)")
                else:
                    size_note = ""
                
                lines.append(f"  {i:<3} {loc:<22} {amp['size']:<8} {ratio_str:<6} {fwd_id:<7} {rev_id:<7} {risk_3p}{size_note}")
            
            lines.append("  " + "-" * 72)
            
            if len(off_targets) > 10:
                lines.append(f"  ... and {len(off_targets) - 10} more primary chromosome off-targets")
            lines.append("")
        else:
            lines.append(Colors.success("  ✓ No off-targets on primary chromosomes!"))
            lines.append("")
        
        # On-target details
        on_targets = [a for a in specificity.potential_amplicons if a['is_on_target']]
        if on_targets:
            lines.append(Colors.bold("ON-TARGET AMPLICON:"))
            for amp in on_targets:
                lines.append(f"  Location: {amp['chromosome']}:{amp['fwd_pos']}-{amp['rev_pos']}")
                lines.append(f"  Size: {amp['size']} bp")
                lines.append(f"  Forward primer identity: {amp['fwd_identity']:.1f}%")
                lines.append(f"  Reverse primer identity: {amp['rev_identity']:.1f}%")
            lines.append("")
        
        # Detailed Risk Analysis
        lines.append(Colors.bold("DETAILED RISK ANALYSIS:"))
        lines.append("")
        
        # Use only primary chromosome off-targets for risk analysis
        primary_off_targets = [a for a in specificity.potential_amplicons 
                               if not a['is_on_target'] and a.get('is_primary_chrom', True)]
        
        # Initialize risk categories
        high_risk = []
        medium_risk = []
        low_risk = []
        negligible = []
        target_size = primer_pair.product_size
        
        if len(primary_off_targets) == 0:
            lines.append(f"  {Colors.success('✓')} No off-target amplification predicted on primary chromosomes")
            lines.append(f"  {Colors.success('✓')} Primers are highly specific for target region")
        else:
            # Categorize off-targets by actual risk (considering size ratio)
            for amp in primary_off_targets:
                fwd_id = amp['fwd_identity']
                rev_id = amp['rev_identity']
                fwd_3p = amp.get('fwd_3prime_match', 0)
                rev_3p = amp.get('rev_3prime_match', 0)
                size_ratio = amp.get('size_ratio', amp['size'] / target_size if target_size > 0 else 1)
                
                # Very large products won't compete - move to negligible
                if size_ratio > 3.0:
                    negligible.append(amp)
                # Categorize based on actual amplification likelihood
                elif fwd_id >= 95 and rev_id >= 95 and fwd_3p >= 5 and rev_3p >= 5:
                    if size_ratio <= 1.5:
                        high_risk.append(amp)
                    else:
                        medium_risk.append(amp)  # Large but high identity
                elif fwd_id >= 90 and rev_id >= 90:
                    medium_risk.append(amp)
                elif fwd_id >= 85 and rev_id >= 85:
                    low_risk.append(amp)
                else:
                    negligible.append(amp)
            
            # Summary Table
            lines.append("  ┌─────────────────────────────────────────────────────────────┐")
            lines.append("  │        OFF-TARGET RISK SUMMARY (Primary Chromosomes)        │")
            lines.append("  ├─────────────────────────────────────────────────────────────┤")
            lines.append(f"  │  {Colors.error('■')} HIGH RISK (will likely amplify):     {len(high_risk):>3} sites          │")
            lines.append(f"  │  {Colors.warning('■')} MEDIUM RISK (may amplify):          {len(medium_risk):>3} sites          │")
            lines.append(f"  │  {Colors.info('■')} LOW RISK (unlikely to amplify):     {len(low_risk):>3} sites          │")
            lines.append(f"  │  {Colors.gray('■')} NEGLIGIBLE (won't amplify):         {len(negligible):>3} sites          │")
            lines.append("  └─────────────────────────────────────────────────────────────┘")
            lines.append("")
            
            # Detailed explanation for each category
            if high_risk:
                lines.append(f"  {Colors.error('HIGH RISK')} - {len(high_risk)} site(s):")
                lines.append(f"      Both primers bind with >95% identity, good 3' match, similar size")
                lines.append(f"      → {Colors.error('These WILL produce non-specific products!')}")
                for i, amp in enumerate(high_risk[:3], 1):
                    ratio = amp.get('size_ratio', 1)
                    lines.append(f"      {i}. {amp['chromosome']}:{amp['fwd_pos']} ({amp['size']}bp, {ratio:.1f}x target)")
                if len(high_risk) > 3:
                    lines.append(f"      ... and {len(high_risk)-3} more")
                lines.append("")
            
            if medium_risk:
                lines.append(f"  {Colors.warning('MEDIUM RISK')} - {len(medium_risk)} site(s):")
                lines.append(f"      Both primers bind with >90% identity OR larger products")
                lines.append(f"      → May amplify under permissive conditions (low Tm, high Mg2+)")
                lines.append(f"      → Target product ({target_size}bp) will out-compete in most cases")
                lines.append("")
            
            if low_risk:
                lines.append(f"  {Colors.info('LOW RISK')} - {len(low_risk)} site(s):")
                lines.append(f"      85-90% primer identity")
                lines.append(f"      → Unlikely to amplify under standard conditions")
                lines.append("")
            
            if negligible:
                lines.append(f"  {Colors.gray('NEGLIGIBLE')} - {len(negligible)} site(s):")
                lines.append(f"      Poor primer binding OR product >3x larger than target")
                lines.append(f"      → Will NOT produce visible bands")
                lines.append("")
        
        # Filter Information
        lines.append(Colors.bold("FILTERING APPLIED:"))
        lines.append(f"  • Scaffolds/Alt-haplotypes: {scaffold_count} sites excluded (not relevant for analysis)")
        lines.append(f"  • Max amplicon size: 3000 bp (larger products filtered out)")
        lines.append(f"  • Size-based risk: Products >3x target size penalized less")
        lines.append(f"  • Both primers must bind on opposite strands")
        lines.append("")
        
        # Get counts for decision
        high_risk_count = len(high_risk)
        medium_risk_count = len(medium_risk)
        
        # User Decision Section
        lines.append(Colors.bold("YOUR DECISION:"))
        lines.append("  ┌─────────────────────────────────────────────────────────────┐")
        
        if len(primary_off_targets) == 0 or high_risk_count == 0:
            lines.append(f"  │  {Colors.success('✓ RECOMMENDED')}: Safe for diagnostic use                    │")
            lines.append("  │     No high-risk off-target sites on primary chromosomes    │")
        elif high_risk_count <= 2:
            lines.append(f"  │  {Colors.warning('⚠ CONDITIONAL')}: May be acceptable with validation         │")
            lines.append(f"  │     {high_risk_count} high-risk site(s) - verify by sequencing              │")
        else:
            lines.append(f"  │  {Colors.error('⚠ CAUTION')}: Review carefully before use                  │")
            lines.append(f"  │     {high_risk_count} high-risk sites detected                              │")
        
        lines.append("  └─────────────────────────────────────────────────────────────┘")
        lines.append("")
        
        # Actionable Recommendations
        lines.append(Colors.bold("SUGGESTED ACTIONS:"))
        if high_risk_count == 0 and medium_risk_count == 0:
            lines.append(f"  {Colors.success('✓')} Primers are highly specific - proceed with confidence")
        elif high_risk_count == 0:
            lines.append(f"  {Colors.success('✓')} Use recommended annealing temperature ({primer_pair.forward.tm:.0f}°C)")
            lines.append(f"  {Colors.success('✓')} Standard Mg2+ concentration (1.5-2.0 mM)")
            lines.append(f"  {Colors.info('○')} Target product is smaller and will out-compete medium-risk products")
        else:
            lines.append(f"  {Colors.warning('1.')} Sequence the PCR product to confirm target amplification")
            lines.append(f"  {Colors.warning('2.')} Use touchdown PCR to increase specificity")
            lines.append(f"  {Colors.warning('3.')} Consider redesigning primers with:")
            lines.append(f"       - Mutation site closer to 3' end")
            lines.append(f"       - Longer primers (23-25 bp)")
            lines.append(f"       - Higher Tm primers (62-65°C)")
            if high_risk_count > 5:
                lines.append(f"  {Colors.info('4.')} Alternative: Use nested PCR or probe-based detection")
        
        lines.append("")
        lines.append(Colors.header("=" * 70))
        
        # Add gel simulation
        gel_sim = GelSimulator.simulate(
            primer_pair.product_size,
            specificity.potential_amplicons
        )
        lines.append(gel_sim)
        
        return "\n".join(lines)


# Add Colors.bold and Colors.gray
Colors.bold = lambda text: f"{Colors.BOLD}{text}{Colors.RESET}"
Colors.gray = lambda text: f"{Colors.GRAY}{text}{Colors.RESET}"


class PrimerPipeline:
    """Main primer design and validation pipeline"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.parser = VariantParser(config)
        self.designer = PrimerDesigner(config)
        self.blast_validator = BlastValidator(config)
        os.makedirs(config.output_dir, exist_ok=True)
    
    def run(self, input_file: str, 
            run_blast: bool = True,
            show_visualization: bool = True) -> Dict:
        """Run the complete pipeline"""
        
        logger.info(f"Processing: {input_file}")
        
        with open(input_file, 'r') as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
        
        print(Colors.header(f"\n{'='*70}"))
        print(Colors.header(f"  VARCONFIRM v6.0"))
        print(Colors.header(f"  Processing {len(lines)} variants..."))
        print(Colors.header(f"{'='*70}\n"))
        
        if run_blast:
            print(Colors.info("  BLAST validation enabled - each primer will be checked against human genome"))
            print(Colors.info("  This ensures high-quality specificity"))
            print(Colors.warning("  Note: BLAST queries take 30-60 seconds each due to NCBI rate limits\n"))
        
        results = []
        failed = []
        
        for idx, line in enumerate(lines, 1):
            print(Colors.bold(f"\n[{idx}/{len(lines)}] ") + f"{line}")
            
            # Parse variant
            variant = self.parser.parse(line)
            if not variant:
                print(Colors.error("  ✗ Parse error"))
                failed.append({'input': line, 'reason': 'Parse failed'})
                continue
            
            if not variant.is_valid:
                print(Colors.error(f"  ✗ Invalid: {', '.join(variant.warnings)}"))
                failed.append({'input': line, 'reason': ', '.join(variant.warnings)})
                continue
            
            print(Colors.info(f"  Format: {variant.format_type.value}"))
            print(Colors.info(f"  Location: {variant.chrom}:{variant.pos}"))
            if variant.gene:
                print(Colors.info(f"  Gene: {variant.gene}"))
            if variant.ref and variant.alt:
                print(Colors.info(f"  Change: {variant.ref} > {variant.alt}"))
            
            # Design primers
            primers = self.designer.design(variant)
            if not primers:
                print(Colors.error("  ✗ Primer design failed"))
                failed.append({'input': line, 'reason': 'Primer design failed'})
                continue
            
            print(Colors.success("  ✓ Primers designed"))
            print(f"    Forward: {Colors.primer_fwd(primers.forward.sequence)} " +
                  f"(Tm={primers.forward.tm:.1f}°C, GC={primers.forward.gc:.0f}%)")
            print(f"    Reverse: {Colors.primer_rev(primers.reverse.sequence)} " +
                  f"(Tm={primers.reverse.tm:.1f}°C, GC={primers.reverse.gc:.0f}%)")
            print(f"    Product: {primers.product_size} bp")
            
            # In silico PCR validation
            if primers.amplicon_validated:
                print(Colors.success("  ✓ In Silico PCR: VALID"))
            else:
                print(Colors.warning("  ⚠ In Silico PCR: Check needed"))
            
            # BLAST validation
            specificity = None
            if run_blast and self.config.blast_enabled:
                specificity = self.blast_validator.validate_primers(
                    primers, variant.chrom, variant.pos
                )
                primers.specificity = specificity
                
                if specificity.status == ValidationStatus.PASS:
                    print(Colors.success(f"  ✓ BLAST Specificity: PASS (Score: {specificity.score:.1f}%)"))
                elif specificity.status == ValidationStatus.WARN:
                    print(Colors.warning(f"  ⚠ BLAST Specificity: WARNING (Score: {specificity.score:.1f}%)"))
                elif specificity.status == ValidationStatus.FAIL:
                    print(Colors.error(f"  ✗ BLAST Specificity: FAIL (Score: {specificity.score:.1f}%)"))
                else:
                    print(Colors.info(f"  ○ BLAST Specificity: {specificity.status.value}"))
            
            # Show visualization
            if show_visualization and len(primers.product_sequence) <= 500:
                viz = ProductVisualizer.visualize(
                    primers.product_sequence,
                    primers.forward.sequence,
                    primers.reverse.sequence,
                    primers.variant_position_in_product,
                    variant.ref,
                    variant.alt
                )
                print(viz)
            
            # Show detailed BLAST results visualization
            if show_visualization and specificity and specificity.potential_amplicons:
                blast_viz = ValidationReport.visualize_blast_results(specificity, variant, primers)
                print(blast_viz)
                
                # Show detailed primer binding analysis
                binding_report = PrimerBindingAnalyzer.generate_detailed_report(
                    specificity.potential_amplicons,
                    primers.forward.sequence,
                    primers.reverse.sequence,
                    primers.product_size,
                    max_display=3  # Show top 3 off-targets in detail
                )
                print(binding_report)
            
            # Show validation report
            report = ValidationReport.generate_report(variant, primers, specificity)
            print(report)
            
            results.append({
                'variant': variant,
                'primers': primers,
                'specificity': specificity
            })
        
        # Save results
        self._save_results(results, failed)
        
        return {
            'total': len(lines),
            'successful': len(results),
            'failed': len(failed),
            'results': results
        }
    
    def _save_results(self, results: List, failed: List):
        """Save results to files"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Main results TSV
        output_file = os.path.join(self.config.output_dir, f"primers_{timestamp}.tsv")
        with open(output_file, 'w') as f:
            headers = [
                "variant_id", "gene", "chrom", "pos", "ref", "alt", "hgvs_p",
                "forward_seq", "forward_tm", "forward_gc", "forward_len",
                "reverse_seq", "reverse_tm", "reverse_gc", "reverse_len",
                "product_size", "mutation_pos_in_product",
                "insilico_pcr", "blast_status", "specificity_score",
                "on_target_hits", "off_target_hits"
            ]
            f.write("\t".join(headers) + "\n")
            
            for r in results:
                v = r['variant']
                p = r['primers']
                s = r.get('specificity')
                
                row = [
                    v.id,
                    v.gene or "",
                    v.chrom,
                    str(v.pos),
                    v.ref,
                    v.alt,
                    v.hgvs_p or "",
                    p.forward.sequence,
                    f"{p.forward.tm:.1f}",
                    f"{p.forward.gc:.1f}",
                    str(p.forward.length),
                    p.reverse.sequence,
                    f"{p.reverse.tm:.1f}",
                    f"{p.reverse.gc:.1f}",
                    str(p.reverse.length),
                    str(p.product_size),
                    str(p.variant_position_in_product),
                    "PASS" if p.amplicon_validated else "CHECK",
                    s.status.value if s else "N/A",
                    f"{s.score:.1f}" if s else "N/A",
                    str(s.on_target_hits) if s else "N/A",
                    str(s.off_target_hits) if s else "N/A"
                ]
                f.write("\t".join(row) + "\n")
        
        print(Colors.success(f"\n✓ Results saved: {output_file}"))
        
        # Detailed reports
        report_file = os.path.join(self.config.output_dir, f"reports_{timestamp}.txt")
        with open(report_file, 'w') as f:
            f.write("PRIMER VALIDATION REPORTS\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 70 + "\n\n")
            
            for r in results:
                v = r['variant']
                p = r['primers']
                s = r.get('specificity')
                
                report = ValidationReport.generate_text_report(v, p, s)
                f.write(report + "\n\n")
                
                # Add product details
                product_detail = ProductVisualizer.text_report(
                    p.product_sequence,
                    p.forward.sequence,
                    p.reverse.sequence,
                    p.variant_position_in_product,
                    v.ref, v.alt
                )
                f.write(product_detail + "\n\n")
        
        print(Colors.success(f"✓ Validation reports saved: {report_file}"))
        
        # Product sequences FASTA
        fasta_file = os.path.join(self.config.output_dir, f"products_{timestamp}.fasta")
        with open(fasta_file, 'w') as f:
            for r in results:
                v = r['variant']
                p = r['primers']
                f.write(f">{v.id}|{v.gene or 'NA'}|{p.product_size}bp\n")
                # Wrap sequence at 60 characters
                seq = p.product_sequence
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")
        
        print(Colors.success(f"✓ Product sequences saved: {fasta_file}"))
        
        # Failed variants
        if failed:
            fail_file = os.path.join(self.config.output_dir, f"failed_{timestamp}.tsv")
            with open(fail_file, 'w') as f:
                f.write("input\treason\n")
                for item in failed:
                    f.write(f"{item['input']}\t{item['reason']}\n")
            print(Colors.warning(f"⚠ Failed variants: {fail_file}"))


def create_test_file(filename: str = "test_20_variants.txt") -> str:
    """Create test file with 20 known variants"""
    with open(filename, 'w') as f:
        f.write("# 20 Known Test Variants\n")
        f.write("# Sources: ClinVar, COSMIC, dbSNP\n\n")
        for var in KNOWN_TEST_VARIANTS:
            f.write(var + "\n")
    return filename


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="VarConfirm v6.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Supported Variant Formats:
  rs121913529                      rsID (dbSNP)
  BRAF:V600E                       Protein notation (50+ common variants)
  NM_004333.6(BRAF):c.1799T>A      HGVS coding notation
  chr7:140753336:A:T               Genomic coordinates (VCF style)
  7-140753336-A-T                  Genomic coordinates (dash style)

Examples:
  python varconfirm.py --test20              # Test with 20 known variants
  python varconfirm.py variants.txt          # Process your variants
  python varconfirm.py variants.txt --no-blast  # Skip BLAST validation
  python varconfirm.py --list-variants       # List known variants

Features:
  - NCBI BLAST specificity validation
  - In silico PCR product verification
  - Specificity scoring (0-100%)
  - Colored terminal output
  - Validation reports
        """
    )
    
    parser.add_argument("input", nargs='?', help="Input variant file (one variant per line)")
    parser.add_argument("--output-dir", default="./results", help="Output directory")
    parser.add_argument("--test20", action="store_true", help="Run test with 20 known cancer variants")
    parser.add_argument("--no-blast", action="store_true", help="Skip BLAST validation (faster)")
    parser.add_argument("--no-viz", action="store_true", help="Disable product visualization")
    parser.add_argument("--list-variants", action="store_true", help="List all known variants in database")
    parser.add_argument("--ncbi-email", default="user@example.com", help="Email for NCBI (recommended)")
    parser.add_argument("--ncbi-api-key", default="", help="NCBI API key (optional, increases rate limit)")
    
    args = parser.parse_args()
    
    # List variants
    if args.list_variants:
        print(Colors.header("\nKnown Protein Variants in Database:"))
        print("-" * 60)
        genes = {}
        for var in sorted(COMMON_VARIANTS.keys()):
            info = COMMON_VARIANTS[var]
            gene = info['gene']
            if gene not in genes:
                genes[gene] = []
            genes[gene].append(var)
        
        for gene in sorted(genes.keys()):
            print(f"\n{Colors.bold(gene)}:")
            for var in genes[gene]:
                info = COMMON_VARIANTS[var]
                print(f"  {var:20} {info['chrom']}:{info['pos']} ({info.get('hgvs_p', '')})")
        
        print(f"\nTotal: {len(COMMON_VARIANTS)} variants")
        return 0
    
    # Configuration
    config = PipelineConfig(
        output_dir=args.output_dir,
        ncbi_email=args.ncbi_email,
        ncbi_api_key=args.ncbi_api_key,
        blast_enabled=not args.no_blast
    )
    
    # Test mode
    if args.test20:
        print(Colors.header("\n" + "=" * 70))
        print(Colors.header("  TESTING WITH 20 KNOWN CANCER VARIANTS"))
        print(Colors.header("=" * 70))
        test_file = create_test_file()
        args.input = test_file
        print(f"Test file created: {test_file}\n")
    
    if not args.input:
        parser.print_help()
        return 1
    
    if not os.path.exists(args.input):
        print(Colors.error(f"File not found: {args.input}"))
        return 1
    
    # Run pipeline
    pipeline = PrimerPipeline(config)
    report = pipeline.run(
        args.input,
        run_blast=not args.no_blast,
        show_visualization=not args.no_viz
    )
    
    # Summary
    print(Colors.header(f"\n{'='*70}"))
    print(Colors.header("  PIPELINE SUMMARY"))
    print(Colors.header(f"{'='*70}"))
    print(f"  Total variants:    {report['total']}")
    print(Colors.success(f"  Successful:        {report['successful']}"))
    if report['failed'] > 0:
        print(Colors.error(f"  Failed:            {report['failed']}"))
    
    # Count by validation status
    if report['results']:
        pass_count = sum(1 for r in report['results'] 
                        if r.get('specificity') and r['specificity'].status == ValidationStatus.PASS)
        warn_count = sum(1 for r in report['results']
                        if r.get('specificity') and r['specificity'].status == ValidationStatus.WARN)
        fail_count = sum(1 for r in report['results']
                        if r.get('specificity') and r['specificity'].status == ValidationStatus.FAIL)
        
        if pass_count or warn_count or fail_count:
            print(f"\n  Validation Results:")
            print(Colors.success(f"    PASS:    {pass_count}"))
            print(Colors.warning(f"    WARNING: {warn_count}"))
            print(Colors.error(f"    FAIL:    {fail_count}"))
    
    print(Colors.header(f"{'='*70}\n"))
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
