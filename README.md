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
  Implements hotspot (WRC, GYW, DGYW, WRCH) and coldspot motifs using IUPAC ambiguity codes, with tunable mutation rates and transition bias.

* **Affinity-based selection**
  BCRs are selected based on normalized local amino-acid alignment scores between CDR3 regions and supplied antigen epitopes.

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

Epitope sequences are used to compute affinity via local amino-acid alignment (BLOSUM62).

---

## Usage

### Basic run (default parameters)

```bash
python bcr_simulator.py
```

### Full parameterized example

```bash
python bcr_simulator.py \
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
  --plot-tree
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

## Intended Use Cases

* Benchmarking clonal clustering algorithms
* Evaluating BCR phylogenetic reconstruction methods
* Stress-testing lineage inference under controlled selection pressure
* Method development for repertoire-scale immunogenomics

Ancestra is explicitly designed for **ground-truth benchmarking**, not for fitting experimental repertoires directly.

---

## Limitations

* Single-chain (heavy chain) simulation
* Simplified affinity model based on local alignment
* No explicit germinal center spatial structure
* No light-chain pairing

These choices are intentional to maintain interpretability and computational efficiency.

---

## Citation

If you use Ancestra in academic work, please cite appropriately. A manuscript or preprint reference can be added here once available.

---

## License

[Specify license here — e.g., MIT, BSD-3, GPL-3]

---

## Contact

For questions, issues, or contributions, please open a GitHub issue or contact the repository maintainer.

