# golfy

Heuristic solver for peptide pool assingments

## Usage

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
