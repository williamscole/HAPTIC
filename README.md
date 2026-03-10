# HAPTiC: HAPlotype Tiling and Clustering

HAPTiC (HAPlotype Tiling and Clustering) is an algorithm for inter-chromosomal haplotype phasing that leverages IBD (Identity By Descent) segments to reconstruct ancestral haplotypes.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
  - [Option 1: pip / venv](#option-1-pip--venv)
  - [Option 2: conda](#option-2-conda)
- [Clustering](#clustering)
  - [Usage](#usage)
  - [Required Arguments](#required-arguments)
  - [Optional Arguments](#optional-arguments)
    - [Input Processing](#input-processing)
    - [IBD Filtering Parameters](#ibd-filtering-parameters)
    - [Chromosome Processing](#chromosome-processing)
    - [ROH and IBD Density Processing](#roh-and-ibd-density-processing)
    - [Performance Options](#performance-options)
    - [Additional Features](#additional-features)
  - [Output](#output)
  - [Examples](#examples)
  - [Notes](#notes)
- [Working with the Results Object](#working-with-the-results-object)
  - [Loading the Results](#loading-the-results)
  - [Iterating Over Individuals](#iterating-over-individuals)
  - [Getting the Optimal Clustering](#getting-the-optimal-clustering)
  - [Getting Tiles for a Specific Threshold](#getting-tiles-for-a-specific-threshold)
  - [Clustering Summary](#clustering-summary)
  - [Accessing IBD Segments](#accessing-ibd-segments)
  - [Full Example](#full-example)
- [Phasing](#phasing)
  - [Usage](#usage-1)
  - [Required Arguments](#required-arguments-1)
- [Example Run](#example-run)
  - [Clustering Step](#clustering-step)
  - [Phasing Step](#phasing-step)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Overview

HAPTiC uses IBD segments shared between individuals to cluster and reconstruct ancestral haplotypes. It can process both autosomal and X chromosome data, handle parent-child relationships, and includes features for ROH (runs of homozygosity) detection and filtering of high IBD density regions.

## Installation

```bash
git clone https://github.com/williamscole/HAPTIC.git
cd HAPTIC
```

### Option 1: pip / venv

```bash
python -m venv haptic_env
source haptic_env/bin/activate        # On Windows: haptic_env\Scripts\activate
pip install -r requirements.txt
```

### Option 2: conda

```bash
conda env create -f environment.yaml
conda activate haptic
```

> **Note:** The conda environment installs `scikit-allel` via pip, as it is more reliably available through PyPI than conda channels. All other dependencies are installed from `conda-forge`.

HAPTiC works in two steps. The ```Clustering``` step takes IBD segments as input and, for each focal individual, clusters their IBD segments into two sets, each set corresponding to a parent of the focal individual. The ```Phasing``` step takes the output of the ```Clustering``` step and a pre-phased VCF to correct for long-range phase errors and phase inter-chromosomally.


# **Clustering**

## Usage

Basic usage:
```bash
python relative_clustering.py -focal_file individuals.txt -relative_ibd_file ibd_segments.feather -out results
```

## Required Arguments

- `-focal_file`: Path to file containing focal individuals. Format:
  - First column: Individual IDs
  - Second and third columns (optional): Parent IDs
  
- `-relative_ibd_file`: Path to file containing IBD segments between relatives.
  - A pandas feather file in the [phasedibd](https://github.com/23andMe/phasedibd) format.

  To convert a phasedibd pandas dataframe to a feather file, use:
```python
df.to_feather("path/to/ibd_segments.feather")
```

## Optional Arguments

### Input Processing

- `-batches`: Specify batch processing format: `batch,current_batch,range`
  - Example: `batch,1,1-26` processes focal individuals from batch 1, loading batches 1-26
  - Example: `batch,1,1-10,12-20,26` processes specific batch ranges
  - Default: `,,,` (no batching)

- `-parent_ibd_file`: Path to file containing IBD segments shared with parents
  - Default: Empty string

- `-single_ibd_file`: Flag to indicate all chromosomes are in one IBD file. If flag is not used, please supply the file for chromosome 1, and HAPTiC will look for the other chromosomes' files. 
  - Default: False (expects separate files per chromosome)

### IBD Filtering Parameters

- `-min_seg_length`: Minimum centiMorgan (cM) length for IBD segment inclusion
  - Default: 5.0 cM

- `-min_k`: Minimum kinship in cM with focal individual for inclusion
  - Default: 10 cM

- `-max_k`: Maximum kinship in cM with focal individual for inclusion
  - Default: 2000 cM

- `-descendant_k`: Maximum kinship threshold for descendant detection
  - Default: 1100.0 cM
  - Set to 0 to skip descendant detection

### Chromosome Processing

- `-noX`: Flag to ignore X chromosome tiles
  - Default: False

- `-x_name`: Specify name for X chromosome in input files
  - Default: "23"

### ROH and IBD Density Processing

- `-keepROH`: Flag to retain ROH and high IBD density regions
  - Default: False

- `-stdevs`: Standard deviations threshold for high-density region filtering
  - Default: 3
  - Removes IBD in regions with density > mean + (stdevs × standard_deviation)

- `-relative_roh`: Flag to compute ROH from relatives
  - Default: False

### Performance Options

- `-threads`: Number of threads for parallel processing
  - Default: 5

- `-multiprocess`: Enable multiprocessing
  - Default: False

- `-coverage_gain`: Early stopping threshold based on coverage gain
  - Default: -inf
  - Stops clustering if (current_coverage - previous_coverage)/current_coverage < coverage_gain

### Additional Features

- `-sex_file`: Path to file containing sex information. A two-column file with `[Individual ID]` and `[Sex Code]` where `1` = male, `2` = female.
  - Default: None

- `-dev_mode`: Enable developer mode
  - Default: False

- `-ibd_only`: Only output IBD segments without clustering
  - Default: False

- `-ibd_pkl`: Path to pre-pickled IBD file
  - Default: Empty string

- `-hap_annotation`: Only perform haplotype annotation. This uses focal-parent IBD segments and outputs the pre-phased regions (haplotype index and coordinates) that were inherited from each parent.
  - Default: False

- `-continue_run`: Continue previous run, append to existing results
  - Default: False

- `-t2`: Custom thresholds for tier 2 clustering
  - Default: [10, 25, 50, 75, 100]

## Output

The algorithm produces several output files with prefix specified by `-out` (default: "clustering"):

- `{out}_results.pkl`: Pickled results file containing clustering information
- Additional output files depending on selected options

## Examples

1. Basic run with minimal parameters:
```bash
python relative_clustering.py -focal_file individuals.txt -relative_ibd_file ibd_segments.feather
```

2. Run with custom IBD filtering:
```bash
python relative_clustering.py -focal_file individuals.txt -relative_ibd_file ibd_segments.feather -min_seg_length 7.5 -min_k 15 -max_k 1500
```

3. Batch processing with parent information:
```bash
python relative_clustering.py -focal_file family_data.txt -relative_ibd_file relatives_ibd.feather -parent_ibd_file parent_ibd.feather -batches batch,1,1-10
```

## Notes

- For large datasets, consider using batch processing and multiprocessing options
- Adjust `-stdevs` parameter based on your population's IBD density distribution
- Use `-dev_mode` for additional debugging information
- The `-coverage_gain` parameter can significantly speed up processing for large datasets data

---

# **Working with the Results Object**

The clustering step produces a `{out}_results.pkl` file. This is a Python pickle file containing a `StoreResults` object that holds the clustering output for every focal individual. This section explains how to load, explore, and extract useful data from this object.

## Loading the Results

```python
import pickle as pkl

with open("clustering_results.pkl", "rb") as f:
    results = pkl.load(f)
```

The `results` object has one main attribute and one method for accessing the data:

- `results.focal_ids` — a list of all focal individual IDs that were successfully clustered.
- `results.iterator()` — returns an iterator of `(focal_id, focal_cluster_object)` pairs (see below).

Each focal individual's clustering is stored as an attribute on the results object. For example, if `"sample_001"` was clustered, its data lives at `results.sample_001`. In practice, you'll almost always want to use the iterator rather than accessing these directly.

## Iterating Over Individuals

The most common pattern for working with the results is to loop over all focal individuals:

```python
for focal_id, focal_obj in results.iterator():
    print(focal_id)
    # focal_obj is a FocalCluster object (see below)
```

Each `focal_obj` is a `FocalCluster` object with the following key attributes and methods:

| Attribute / Method | Type | Description |
|---|---|---|
| `focal_obj.focal` | `str` | The focal individual's ID |
| `focal_obj.male` | `bool` | Whether the focal individual is male |
| `focal_obj.scheme` | `bool` | `True` if clustering was run successfully |
| `focal_obj.segments` | `DataFrame` | The processed IBD segments used for clustering |
| `focal_obj.tile_df` | `DataFrame` | The full tile dataframe with all clustering columns |
| `focal_obj.clusterings` | `dict` | Maps `(t1, t2)` threshold pairs to clustering results |
| `focal_obj.get_optimal(tile_df=False)` | method | Returns the best clustering (see below) |
| `focal_obj.return_tiles(t1, t2)` | method | Returns the tile dataframe for a specific threshold pair |
| `focal_obj.cluster_summary()` | method | Returns a summary dataframe of all clusterings |

## Getting the Optimal Clustering

HAPTiC tries multiple threshold combinations (`t1`, `t2`) and selects the best one automatically. To get the optimal result:

```python
for focal_id, focal_obj in results.iterator():
    # Get the optimal clustering metadata
    optimal = focal_obj.get_optimal()
    print(f"{focal_id}: eigenvalue = {optimal.eigenval1}, coverage = {optimal.tot_cov} cM")
```

The object returned by `get_optimal()` (without `tile_df=True`) is a `SpectralCluster` object with these attributes:

| Attribute | Type | Description |
|---|---|---|
| `eigenval1` | `float` | Smallest relevant eigenvalue of the signed Laplacian. Lower is better; a value near 0 indicates clean separation into two parental clusters. |
| `tot_cov` | `float` | Total genomic coverage (in cM) of the tiles used in this clustering. |
| `optimal` | `bool` | `True` if this was selected as the best clustering. |
| `worked` | `bool` | `True` if the spectral clustering succeeded. |
| `t1` | `int` | The `t1` threshold (minimum segment length in cM for edge inclusion). |
| `t2` | `int` | The `t2` threshold (minimum kinship in cM for edge inclusion). |
| `x1` | `float` | Total IBD (cM) on the X chromosome assigned to cluster 1. |
| `x2` | `float` | Total IBD (cM) on the X chromosome assigned to cluster 2. |

To get the optimal tile dataframe (which is what the phasing step uses):

```python
tile_df = focal_obj.get_optimal(tile_df=True)
```

This returns a DataFrame with the following columns:

| Column | Description |
|---|---|
| `chromosome` | Chromosome number (1–23, where 23 = X) |
| `coord` | Tuple of `(start_cM, end_cM)` for the tile |
| `cluster1` | List of relative segment IDs assigned to parental cluster 1 |
| `cluster2` | List of relative segment IDs assigned to parental cluster 2 |
| `segs1` | List of segment indices corresponding to `cluster1` |
| `segs2` | List of segment indices corresponding to `cluster2` |
| `haplotype1` | Haplotype index (0 or 1) associated with cluster 1 |
| `haplotype2` | Haplotype index (0 or 1) associated with cluster 2 |
| `sex` | `True` if the clusters are ordered as [maternal, paternal] based on X-chromosome IBD (males only); `False` otherwise |

Each relative segment ID in `cluster1`/`cluster2` is a tuple of `(relative_id, segment_number, kinship_cM, merged_segment_length)`.

### Threshold Columns in the Raw Tile DataFrame

The raw tile DataFrame (accessed via `focal_obj.tile_df`) also contains additional columns corresponding to each `(t1, t2)` threshold pair that was tested during clustering. These columns are named as tuples, e.g. `(10, 50)`, `(5, 100)`, etc. The values in these columns indicate how each tile relates to the global parental assignment for that clustering:

- **`False`** — the tile's original haplotype labeling is already consistent with the global parental assignment.
- **`True`** — the tile's `cluster1`/`cluster2` (and their associated haplotype indices) need to be swapped to align with the global parental assignment.
- **`NaN`** — the tile was not connected to the main graph under this threshold pair (the thresholds were too conservative for any qualifying edges to reach this tile).

When `return_tiles(t1, t2)` is called, it drops all tiles that are `NaN` for that threshold, then uses the `True`/`False` values to swap `cluster1`↔`cluster2`, `segs1`↔`segs2`, and `haplotype1`↔`haplotype2` where needed, so that `cluster1` consistently represents the same parent across all returned tiles.

## Getting Tiles for a Specific Threshold

If you want to inspect the clustering at a particular `(t1, t2)` threshold rather than the optimal:

```python
tile_df = focal_obj.return_tiles(t1=10, t2=50)
```

The returned DataFrame has the same columns as described above. If the requested threshold was not run, this returns an empty DataFrame.

## Clustering Summary

To get a quick overview of all threshold combinations that were tried for a focal individual:

```python
summary = focal_obj.cluster_summary()
print(summary)
```

This returns a DataFrame with one row per threshold combination:

| Column | Description |
|---|---|
| `focal` | Focal individual ID |
| `t1` | Segment length threshold |
| `t2` | Kinship threshold |
| `tot_cov` | Total tile coverage in cM |
| `eigenval1` | Smallest eigenvalue (lower = cleaner clustering) |
| `optimal` | `True` for the selected best clustering |

## Accessing IBD Segments

The processed IBD segments used for clustering are stored on each focal object:

```python
segments = focal_obj.segments
```

This is a DataFrame with columns `id2`, `id1_haplotype`, `chromosome`, `start_cm`, and `end_cm`. The `id2` column contains tuple IDs of the form `(relative_id, segment_number, kinship_cM, merged_segment_length)`.

## Full Example

```python
import pickle as pkl
import pandas as pd

# Load results
with open("clustering_results.pkl", "rb") as f:
    results = pkl.load(f)

# Summarize all focal individuals
rows = []
for focal_id, focal_obj in results.iterator():
    opt = focal_obj.get_optimal()
    rows.append({
        "focal": focal_id,
        "eigenvalue": opt.eigenval1,
        "coverage_cM": opt.tot_cov,
        "t1": opt.t1,
        "t2": opt.t2,
        "is_male": focal_obj.male
    })

summary = pd.DataFrame(rows)
print(summary)

# Get the tile dataframe for the first individual
focal_id, focal_obj = next(results.iterator())
tile_df = focal_obj.get_optimal(tile_df=True)

# Print tiles for chromosome 1
chr1_tiles = tile_df[tile_df.chromosome == 1]
for _, tile in chr1_tiles.iterrows():
    start, end = tile.coord
    n_cluster1 = len(tile.cluster1)
    n_cluster2 = len(tile.cluster2)
    print(f"  {start:.1f}-{end:.1f} cM: {n_cluster1} relatives in cluster 1, {n_cluster2} in cluster 2")
```

---

# **Phasing**

## Usage

Basic usage:
```bash
python phase_vcf.py -phase -vcf input.vcf -map chr1.map -chr 1 -results clustering_results.pkl -output corrected_input.vcf -write_phase
```

## Required Arguments

- `-phase`: this flag instructs the program to phase the input VCF.
- `-vcf`: a ```.vcf``` file that is pre-phased with Beagle, SHAPEIT, or Eagle. This must be the same file that was used to call the IBD input in the ```Clustering``` step.
- `-map`: a PLINK-formatted .map file whose sites correspond to the VCF sites. HAPTiC expects the .map file for chromosome 1 and will find the .map file for the other chromosomes.
- `-chr`: the chromosome number of the VCF.
- `-results`: the ```[output]_results.pkl``` file from the ```Clustering``` step.
- `-output`: the full path and file name of the outputted VCF.
- `-write_phase`: instructs the program to write out the phase.

# Example Run

We have provided a simple dataset for users to test HAPTIC.

## Clustering Step

```bash
python3 relative_clustering.py -focal_file example/focal_ids.txt -relative_ibd_file example/ibd_segments.feather -single_ibd_file -min_seg_length 5 -min_k 0 -keepROH -out example/clustering
```

## Phasing Step

```bash
python3 phase_vcf.py -phase -vcf example/chr1.vcf -map example/chr1.map -chr 1 -results example/clustering_results.pkl -write_phase
```

See the README in ```examples/```

# Citation

If you use HAPTiC in your research, please cite our preprint:

```
Williams, C. M., O'Connell, J., Jewett, E., Freyman, W. A., 23andMe Research Team, Gignoux, C. R., Ramachandran, S., & Williams, A. L. (2024). Phasing millions of samples achieves near perfect accuracy, enabling parent-of-origin classification of variants. bioRxiv, 2024.05.06.592816. Cold Spring Harbor Laboratory.
```

## License

HAPTiC is released under the GNU General Public License v3.0 (GPL-3.0).
This means you are free to:

Use the software for any purpose
- Change the software to suit your needs 
- Share the software with others 
- Share the changes you make 

Under the following terms:

- You must include the original copyright and license notices in any copy of the software/source 
- If you modify the code, you must release these modifications under the GPL-3.0 
- You must state significant changes made to the software 
- You must make your source code available when you distribute the software 

For the full license text, see the LICENSE file in the repository.

## Contact

cole [underscore] williams [at] brown [dot] edu