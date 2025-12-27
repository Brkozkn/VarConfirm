# VarConfirm

Primer design tool for Sanger sequencing validation of NGS-detected variants.

## Overview

NGS → Variant Detection → **VarConfirm** → Sanger Validation

## Features

- Multiple input formats (rsID, HGVS, VCF, genomic coordinates)
- Automatic primer design with optimal Tm/GC
- BLAST-based specificity validation
- Off-target amplicon analysis
- Gel electrophoresis simulation

## Installation

```bash
git clone https://github.com/Brkozkn/VarConfirm.git
cd VarConfirm

conda create -n varconfirm python=3.10 -y
conda activate varconfirm
pip install -r requirements.txt
```

### Local BLAST (Optional, faster)

```bash
sudo apt install ncbi-blast+
mkdir -p ~/blast_db && cd ~/blast_db
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
makeblastdb -in GCF_000001405.40_GRCh38.p14_genomic.fna -dbtype nucl -parse_seqids -out human_genome
```

## Usage

```bash
# Design primers for variants
python varconfirm.py variants.txt

# Skip BLAST validation
python varconfirm.py variants.txt --no-blast

# Test with examples
python varconfirm.py --test20
```

## Input Formats

| Format | Example |
|--------|---------|
| Gene:Protein | `BRAF:V600E` |
| rsID | `rs113488022` |
| HGVS | `NM_004333.4(BRAF):c.1799T>A` |
| Genomic | `chr7:140753336:A>T` |

## Output

Results saved to `./results/`:
- `primers_*.tsv` - Primer sequences
- `reports_*.txt` - Validation reports  
- `products_*.fasta` - Amplicon sequences

## License

MIT License

## Author

Burak Özkan - Ulm University
