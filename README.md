# LexiCut: Verification Code for Max-Cut with Geometric Edge Weights

This repository contains the computational verification code accompanying the paper:

**"Hierarchical threshold structure in Max-Cut with geometric edge weights"**  
by Nevena Marić  


## Problem Description

We study Max-Cut on the complete graph $K_n$ where edges are weighted geometrically in lexicographic order: the $i$-th edge has weight $r^{N-i}$ where $N = \binom{n}{2}$ and $r \in (1,2)$.

The code verifies the conjecture that for $n \geq 7$, the optimal cut is always a $k$-isolated cut $C_k = \{1,\ldots,k\} \mid \{k+1,\ldots,n\}$.

## Files

- `lexicut_verify.R` — Main verification code (R, no external packages required)

## Usage

```r
source("lexicut_verify.R")

# Verify conjecture for a single n
result <- verify_conjecture(25)

# Batch verification for multiple n values
results <- verify_batch(21:30)

# Compute threshold values r_k(n)
thresholds <- compute_all_thresholds(20)

# Focused analysis near threshold boundaries
violations <- verify_near_thresholds(15)
```

## Main Functions

| Function | Description |
|----------|-------------|
| `verify_conjecture(n)` | Main verification for single $n$ |
| `verify_batch(n_values)` | Batch verification for multiple $n$ |
| `compute_all_thresholds(n)` | Compute threshold values $r_k(n)$ |
| `verify_near_thresholds(n)` | Focused analysis near thresholds |

## Verification Strategy

1. **Threshold computation**: Find $r_k(n)$ where $k$-isolated and $(k+1)$-isolated cuts have equal weight
2. **Near-isolated verification**: Exhaustively check cuts of the form $S^*_k = \{1,\ldots,k,n\}$
3. **Random sampling**: Test random cuts against the predicted optimal isolated cut

## Results

- Exhaustive enumeration for $n \leq 21$: No violations found for $n \geq 7$
- Near-isolated verification for $n \leq 100$: Conjecture supported
- Counterexamples for $n \in \{4,5,6\}$: Characterized completely in the paper

## Requirements

- R (tested on version 4.0+)
- No external packages required

## License

MIT License

## Contact

Nevena Marić  
School of Computing, Union University, Belgrade, Serbia  
nmaric@raf.rs
