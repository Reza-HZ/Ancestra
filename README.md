# Ancestra

**A lineage-explicit simulator for benchmarking B-cell receptor repertoire and lineage inference methods**

Ancestra is a biologically grounded simulator for generating affinity-matured B-cell receptor (BCR) repertoires with explicit clonal lineage structure. It is designed for **benchmarking repertoire analysis, clonal clustering, phylogenetic reconstruction, and lineage inference methods** under controlled and reproducible conditions.

Unlike purely sequence-level simulators, Ancestra jointly models:

* V(D)J recombination
* Junctional diversity
* Context-dependent somatic hypermutation
* Affinity-based selection against antigen epitopes
* Explicit genealogical lineage trees with abundances

---

## Key Features

* **Lineage-explicit simulation**
  Full parent–child relationships are retained, enabling direct evaluation of lineage and phylogeny inference accuracy.

* **VDJ recombination with biological constraints**
  Random V/D/J selection, trimming, N-nucleotide insertions, and frame selection consistent with productive rearrangements.

* **Context-dependent somatic hypermutation (SHM)**
  Implements hotspot (WRC, GYW, DGYW, WRCH) and coldspot motifs, with tunable mutation rates and transition bias.

* **Flexible affinity-based selection**
  BCRs are selected based on their affinity to the antigen epitopes. Affinity can be computed using:

  *  CDR3-only region (default; biologically focused and computationally efficient)

  *  Full BCR amino acid sequence (optional)

* **Explicit modeling of deleterious mutations**
  Stop codons and low-affinity variants are penalized probabilistically rather than deterministically.

* **Tree-aware outputs**
  Exports FASTA repertoires, Newick lineage trees, per-sequence metadata, and optional lineage visualizations.

* **Reproducible simulation metadata**
  Each run saves parameters, statistics, and outcomes for downstream analysis.

---

## Installation

### Requirements

* Python ≥ 3.8
* Core dependencies:

  * `numpy`
  * `biopython`
  * `matplotlib`
  * `ete3`
  * `mplcursors`

Install dependencies via pip:

```bash
pip install numpy biopython matplotlib ete3 mplcursors
```

---

## Input Files

Ancestra expects an input directory containing:

```
input_dir/
├── V.fasta        # V gene segments
├── D.fasta        # D gene segments
├── J.fasta        # J gene segments
└── epitope.txt    # One amino-acid epitope per line
```

---

## Usage

### Basic run (default parameters)

```bash
python Ancestra_V1.py
```

### Full parameterized example

```bash
python Ancestra_V1.py \
  --input-dir ./IGH_genes \
  --clones 5 \
  --max-gen 20 \
  --min-seq 50 \
  --t-min 0.3 \
  --t-max 0.8 \
  --mu 0.001 \
  --p-sub 0.3 \
  --p-stop 0.005 \
  --p-min 0.001 \
  --p-trans 0.7 \
  --use-full-BCR True \
  --plot-tree
```

### For full usage:

```bash
python Ancestra_V1.py --help
```


---

## Command-Line Arguments

### Simulation control

| Argument    | Description                                      | Default |
| ----------- | ------------------------------------------------ | ------- |
| `--clones`  | Number of successful clonal lineages to generate | `1`     |
| `--max-gen` | Maximum lineage depth                            | `15`    |
| `--min-seq` | Minimum unique sequences required for acceptance | `10`    |

### Affinity selection

| Argument  | Description                                    | Default |
| --------- | ---------------------------------------------- | ------- |
| `--t-min` | Minimum affinity threshold (early generations) | `0.3`   |
| `--t-max` | Maximum affinity threshold (late generations)  | `0.8`   |
| `--use-full-BCR` | If True, affinity is computed using the full amino acid sequence. If False, only the CDR3 region is used.  | `False`   |

### Mutation model

| Argument    | Description                           | Default |
| ----------- | ------------------------------------- | ------- |
| `--mu`      | Baseline per-nucleotide mutation rate | `0.001` |
| `--p-trans` | Transition vs. transversion bias      | `0.7`   |

### Survival probabilities

| Argument   | Description                                    | Default |
| ---------- | ---------------------------------------------- | ------- |
| `--p-sub`  | Survival of low-affinity, in-frame sequences   | `0.3`   |
| `--p-stop` | Survival of high-affinity sequences with stops | `0.005` |
| `--p-min`  | Survival of low-affinity sequences with stops  | `0.001` |

### Output

| Argument       | Description                         |
| -------------- | ----------------------------------- |
| `--plot-tree`  | Generate lineage tree visualization |
| `--output-dir` | Output directory                    |
| `--verbose`    | Enable debug logging                |

  By default, `--plot-tree` is set to `false`. When enabled (`--plot-tree`), the simulated lineage tree is rendered and saved as `lineage_tree.png` in the corresponding output directory.
  
---

## Output Structure

```
output_dir/
├── successful_clones/
│   └── run_1_gen15_minSeq10/
│       ├── repertoire.fasta
│       ├── repertoire.nk
│       ├── repertoire_info.csv
│       ├── lineage_tree.png
│       └── run_metadata.json
└── rejected_clones/
    └── rejected_attempt_*
```

### Output files

* **FASTA**: Nucleotide sequences annotated with abundances
* **Newick**: True lineage tree with mutation distances
* **CSV**: Per-sequence metadata (generation, affinity, mutations, parent, abundance)
* **PNG** (optional): Interactive lineage visualization
* **JSON**: Full run configuration and summary statistics


---


## Contact

For questions, bug reports, or contributions, please contact the repository maintainer, **Reza Hassanzadeh**, at **rhz.sbu@gmail.com**.
