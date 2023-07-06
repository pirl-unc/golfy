# golfy

Heuristic solver for peptide pool assingments

## Installation

```sh
pip install golfy
```

Also, [scikit-learn](https://scikit-learn.org/stable/index.html) is an requirement for the deconvolution module:

```sh
pip install scikit-learn
```

## Usage

### Peptide to pool assignment

```python

from golfy import init, is_valid, optimize

# create a random initial assignment of peptides to pools
s = init(num_peptides = 100, peptides_per_pool = 5, num_replicates = 3)

# the random assignment probably isn't yet a valid solution
assert not is_valid(s)

# iteratively swap peptides which violate constraints until
# a valid configuration is achieved
optimize(s)

assert is_valid(s)
```

### Deconvolution of hit peptides from ELISpot counts

```python
from golfy.deconvolution import create_linear_system, solve_linear_system

# s is a golfy.Solution object containing the mapping of peptides to pools
# counts is a dictionary from (replicate, pool) index pairs to ELISpot counts or activity values
linear_system = create_linear_system(s, counts)

# result type has an array of individual peptide activity estimates (result.activity_per_peptide)
# and a set of high confidence hit peptides (result.high_confidence_peptides)
result = solve_linear_system(linear_system)
print(result.high_confidence_hits)
```
