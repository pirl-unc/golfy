# golfy

Heuristic solver for peptide pool assignments. Golfy attempts to map peptides to pools with a fixed specified "coverage" (number of pools that each peptide occurs in) while avoiding having any pair of peptides occur in a pool twice.

## Installation

```sh
pip install golfy
```

Also, [scikit-learn](https://scikit-learn.org/stable/index.html) is an requirement for the deconvolution module:

```sh
pip install scikit-learn
```

## Usage

### High level API

Assignments of peptides to pools are called `golfy.Design` objects, which can be constructed and optimized using several different strategies. You can
ignore most of the implementation details by just calling `find_best_design`, which tries multiple different initialization strategies to create multiple designs, optimizes each one, and returns the design which fewest constraint violations and fewest number of total pools.

```python
from golfy import find_best_design

s = find_best_design(
    num_peptides=100,
    max_peptides_per_pool=5,
    num_replicates=3,
    invalid_neighbors=[(0,1), (1,2)],
    preferred_neighbors=[(0,3),(1,5)],
    allow_extra_pools=False,
    verbose=False)
```

A key parameter to `find_best_design` is `allow_extra_pools`, which determines whether Golfy is allowed to expand the number of total pools beyond the minimum by the `ceil(num_peptides * num_replicates / max_peptides_per_pool)`. If Golfy cannot add extra pools then
it may not be able to find a valid solution for every combination of parameters (but will still give you the design with the least constraint violations that it could construct and optimize).

### Initialization and optimization of designs

If you want more control over the way that desgins are initialized and optimized you can use the `golfy.init` and `golfy.optimize` functions directly.

```python

from golfy import init, is_valid, optimize

# create a random initial assignment of peptides to pools
s = init(num_peptides=100, peptides_per_pool=5, num_replicates=3, strategy='random', allow_extra_pools=False)

# the random assignment probably isn't yet a valid design
assert not is_valid(s)

# iteratively swap peptides which violate constraints until
# a valid configuration is achieved
optimize(s, allow_extra_pools=False)

assert is_valid(s)
```

### Deconvolution of hit peptides from ELISpot counts

```python
from golfy.deconvolution import create_linear_system, solve_linear_system

# s is a golfy.Design object containing the mapping of peptides to pools
# counts is a dictionary from (replicate, pool) index pairs to ELISpot counts or activity values
linear_system = create_linear_system(s, counts)

# result type has an array of individual peptide activity estimates (result.activity_per_peptide)
# and a set of high confidence hit peptides (result.high_confidence_peptides)
result = solve_linear_system(linear_system)
print(result.high_confidence_hits)
```
