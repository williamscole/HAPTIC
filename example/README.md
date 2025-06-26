## Example: IBD Tiling with Phase Switch Error

In this example, we have a focal individual called **"focal"** who shares IBD segments with 6 relatives:

- **A, B, C** are paternal relatives sharing IBD with haplotype 0
- **D, E, F** are maternal relatives sharing IBD with haplotype 1

### Phase Switch Error
All relatives correctly phase until a **switch error occurs between 59-60 cM**, causing haplotypes to flip.

### Chromosome Coverage
**Chromosome length:** 100 cM

### Resulting Tiles
| Tile | Region | Haplotype 0 | Haplotype 1 |
|------|--------|-------------|-------------|
| 1    | 0-19 cM | {A} | {D} |
| 2    | 20-59 cM | {A} | {E} |
| 3    | 60-100 cM | {E, F} | {B, C} |

**Note:** In Tile 3, the haplotype assignments have switched due to the phase error at 60 cM.

### Expected Output VCF
HAPTIC will detect the phase switch and correct it in the output VCF:

- **Input VCF:** Shows `0|1` from 0-59 cM, then switches to `1|0` from 60-100 cM
- **Output VCF:** Will be corrected to show consistent `0|1` phasing across the entire chromosome

The algorithm identifies that relatives switched haplotypes at 60 cM and applies the appropriate phase correction to maintain consistent parent-of-origin assignment.