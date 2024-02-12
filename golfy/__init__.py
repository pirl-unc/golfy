# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations

from .deconvolution import (
    create_linear_system,
    solve_linear_system,
    em_deconvolve,
    deconvolve,
    DeconvolutionResult,
)
from .design import Design
from .initialization import init
from .main import find_best_design, best_design_for_pool_budget
from .optimization import optimize
from .simulation import (
    simulate_elispot_counts,
    simulate_any_hits_per_pool,
    simulate_number_hits_per_pool,
)
from .types import SpotCounts
from .validity import is_valid, count_violations, violations_per_replicate

__version__ = "2.5.2"

__all__ = [
    "__version__",
    "find_best_design",
    "best_design_for_pool_budget",
    "Design",
    "init",
    "optimize",
    "count_violations",
    "is_valid",
    "violations_per_replicate",
    "simulate_elispot_counts",
    "simulate_any_hits_per_pool",
    "simulate_number_hits_per_pool",
    "create_linear_system",
    "solve_linear_system",
    "SpotCounts",
    "deconvolve",
    "em_deconvolve",
    "DeconvolutionResult",
]
