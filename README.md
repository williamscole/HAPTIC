# HAPTIC: HAPlotype TIling and Clustering

HAPTIC (HAPlotype TIling and Clustering) is an algorithm for inter-chromosomal haplotype phasing that leverages IBD (Identity By Descent) segments to reconstruct ancestral haplotypes.

## Overview

HAPTIC uses IBD segments shared between individuals to cluster and reconstruct ancestral haplotypes. It can process both autosomal and X chromosome data, handle parent-child relationships, and includes features for ROH (Runs of Homozygosity) detection and filtering of high IBD density regions.

## Installation

```bash
git clone https://github.com/williamscole/haptic.git
cd haptic
pip install -r requirements.txt
```

## Usage

Basic usage:
```bash
python haptic.py -focal_file individuals.txt -relative_ibd_file ibd_segments.txt -out results
```

## Required Arguments

- `-focal_file`: Path to file containing focal individuals. Format:
  - First column: Individual IDs
  - Second and third columns (optional): Parent IDs
  
- `-relative_ibd_file`: Path to file containing IBD segments between relatives

## Optional Arguments

### Input Processing

- `-batches`: Specify batch processing format: `batch,current_batch,range`
  - Example: `batch,1,1-26` processes focal individuals from batch 1, loading batches 1-26
  - Example: `batch,1,1-10,12-20,26` processes specific batch ranges
  - Default: `,,,` (no batching)

- `-parent_ibd_file`: Path to file containing IBD segments shared with parents
  - Default: Empty string

- `-single_ibd_file`: Flag to indicate all chromosomes are in one IBD file
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

- `-sex_file`: Path to file containing sex information
  - Default: None

- `-dev_mode`: Enable developer mode
  - Default: False

- `-ibd_only`: Only output IBD segments without clustering
  - Default: False

- `-ibd_pkl`: Path to pre-pickled IBD file
  - Default: Empty string

- `-hap_annotation`: Only perform haplotype annotation
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
python haptic.py -focal_file individuals.txt -relative_ibd_file ibd_segments.txt
```

2. Run with custom IBD filtering:
```bash
python haptic.py -focal_file individuals.txt -relative_ibd_file ibd_segments.txt -min_seg_length 7.5 -min_k 15 -max_k 1500
```

3. Batch processing with parent information:
```bash
python haptic.py -focal_file family_data.txt -relative_ibd_file relatives_ibd.txt -parent_ibd_file parent_ibd.txt -batches batch,1,1-10
```

## Notes

- For large datasets, consider using batch processing and multiprocessing options
- Adjust `-stdevs` parameter based on your population's IBD density distribution
- Use `-dev_mode` for additional debugging information
- The `-coverage_gain` parameter can significantly speed up processing for large datasets data

## Citation

If you use HAPTIC in your research, please cite our preprint:

```
Williams, C. M., O'Connell, J., Freyman, W. A., 23andMe Research Team, Gignoux, C. R., Ramachandran, S., & Williams, A. L. (2024). Phasing millions of samples achieves near perfect accuracy, enabling parent-of-origin classification of variants. bioRxiv, 2024.05.06.592816. Cold Spring Harbor Laboratory.
```

## License

[License information to be added]

## Contact

cole [underscore] williams [at] brown [dot] edu
