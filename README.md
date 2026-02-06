# Ancestra  
**A lineage-explicit simulator for benchmarking B-cell receptor repertoire and lineage inference methods**

Ancestra is a biologically grounded simulator for generating affinity-matured B-cell receptor (BCR) repertoires with explicit clonal lineage structure. It is designed for **benchmarking repertoire analysis, clonal clustering, phylogenetic reconstruction, and lineage inference methods** under controlled and reproducible conditions.

Unlike purely sequence-level simulators, Ancestra jointly models:

- V(D)J recombination  
- Junctional diversity  
- Context-dependent somatic hypermutation  
- Affinity-based selection against antigen epitopes  
- Explicit genealogical lineage trees with abundances  

---

## Key Features

- **Lineage-explicit simulation**  
  Full parent–child relationships are retained, enabling direct evaluation of lineage and phylogeny inference accuracy.

- **VDJ recombination with biological constraints**  
  Random V/D/J selection, trimming, N-nucleotide insertions, and frame selection consistent with productive rearrangements.

- **Context-dependent somatic hypermutation (SHM)**  
  Implements hotspot (WRC, GYW, DGYW, WRCH) and coldspot motifs using IUPAC ambiguity codes, with tunable mutation rates and transition bias.

- **Affinity-based selection**  
  BCRs are selected based on normalized local amino-acid alignment scores between CDR3 regions and supplied antigen epitopes.

- **Explicit modeling of deleterious mutations**  
  Stop codons and low-affinity variants are penalized probabilistically rather than deterministically.

- **Tree-aware outputs**  
  Exports FASTA repertoires, Newick lineage trees, per-sequence metadata, and optional lineage visualizations.

- **Reproducible simulation metadata**  
  Each run saves parameters, statistics, and outcomes for downstream analysis.

---

## Installation

### Requirements

- Python ≥ 3.8
- Core dependencies:
  - numpy
  - biopython
  - matplotlib
  - ete3
  - mplcursors

Install dependencies via pip:

```bash
pip install numpy biopython matplotlib ete3 mplcursors

### Input Files

Ancestra expects an input directory containing:

input_dir/
├── V.fasta        # V gene segments
├── D.fasta        # D gene segments
├── J.fasta        # J gene segments
└── epitope.txt    # One amino-acid epitope per line

