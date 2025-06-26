# HAPTiC: HAPlotype Tiling and Clustering

HAPTiC (HAPlotype Tiling and Clustering) is an algorithm for inter-chromosomal haplotype phasing that leverages IBD (Identity By Descent) segments to reconstruct ancestral haplotypes.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
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
conda create -n haptic python=3.9 "numpy>=1.21.0" "pandas>=1.3.0" "networkx>=2.6.0" matplotlib seaborn  pyarrow
conda activate haptic
```

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
python relative_clustering.py -focal_file example/focal_ids.txt -relative_ibd_file example/ibd_segments.feather -single_ibd_file -min_seg_length 5 -min_k 0 -keepROH -out example/clustering
```

## Phasing Step

```bash
python phase_vcf.py -phase -vcf example/chr1.vcf -map example/chr1.map -chr 1 -results example/clustering_results.pkl -write_phase
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