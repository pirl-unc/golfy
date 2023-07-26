from dataclasses import dataclass

import numpy as np

from .deconvolution import create_linear_system, solve_linear_system
from .design import Design
from .simulation import simulate_number_hits_per_pool
from .validity import count_violations


@dataclass
class EvaluationResult:
    precision: float
    recall: float
    f1: float
    num_pools: int
    num_violations: int

    def sort_key(self):
        # maximize f1, minimize violations,  maximize (precision, recall), minimize pools, maximize replicates
        return (
            -round(self.f1, 2),
            self.num_violations,
            -round(self.precision, 2),
            -round(self.recall, 2),
            self.num_pools,
        )


def evaluate_design(
    s: Design,
    num_simulation_iters: int = 10,
    min_hit_fraction: float = 0.01,
    max_hit_fraction: float = 0.05,
) -> EvaluationResult:
    """
    Returns average precision, recall, F1 scores across multiple simulations of the given design for
    num_hits in range of 1% to 5% of the total number of peptides
    """
    ps = []
    rs = []
    f1s = []
    for _ in range(num_simulation_iters):
        for hit_fraction in np.arange(min_hit_fraction, max_hit_fraction + 0.001, 0.01):
            num_hits = int(np.ceil(hit_fraction * s.num_peptides))

            spot_counts, hit_peptides = simulate_number_hits_per_pool(
                s, num_hits=num_hits
            )
            linear_system = create_linear_system(s, spot_counts)
            result = solve_linear_system(
                linear_system, leave_on_out=False, min_peptide_activity=0.2
            )

            predicted_hits = result.high_confidence_hits

            tp = len(predicted_hits.intersection(hit_peptides))
            fp = len({p for p in predicted_hits if p not in hit_peptides})
            fn = len({p for p in hit_peptides if p not in predicted_hits})

            precision = tp / (tp + fp) if tp + fp > 0 else 0
            recall = tp / (tp + fn) if tp + fn > 0 else 0
            f1 = (
                2 * precision * recall / (precision + recall)
                if precision + recall > 0
                else 0
            )

            ps.append(precision)
            rs.append(recall)
            f1s.append(f1)
    return EvaluationResult(
        precision=np.mean(ps),
        recall=np.mean(rs),
        f1=np.mean(f1s),
        num_pools=s.num_pools(),
        num_violations=count_violations(
            s, error_on_duplicate=False, error_on_extra=False, error_on_missing=False
        ),
    )
