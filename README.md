# ATLAS
code for the paper " Nonlinear model reduction for slow-fast stochastic systems near  manifolds" by Felix X.-F. Ye, Sichen Yang, Mauro Maggioni. 

The paper contains results for three example problems based on Pinched sphere(nick name: peanut), Oscillating half-moons and Butane model. Code for each example can be found in the respective subfolder in the examples folder.

In each example, there are two steps. The first is the ATLAS learning stage and the second is the mean resident time simulating stage.


## Python implementation

The repository now contains a pure Python port of the MATLAB ATLAS reference implementation under the [`pyatlas`](pyatlas/) package. The translation mirrors the original workflow while relying on NumPy and SciPy-ready data structures and simulators:

* `pyatlas.core` defines chart data models, the slow-manifold learning routine, and projection/interpolation helpers.
* `pyatlas.simulator` exposes Euler–Maruyama integrators that accumulate covariance and mean statistics in streaming fashion, making the large Monte-Carlo batches from the MATLAB code practical on CPUs.
* `pyatlas.models` provides ready-to-use configurations for the halfmoon, peanut, and butane examples, including drift/diffusion callables and the training-time heuristics from the original scripts.
* `pyatlas.workflows` orchestrates initial chart learning and weighted projections, providing a Pythonic entry-point for end-to-end experiments.

A minimal example for running the initial learning phase on the halfmoon model is shown below:

```python
import numpy as np
from pyatlas import halfmoon_configuration, initialise_atlas

config = halfmoon_configuration()
state = initialise_atlas(config, rng=np.random.default_rng(42))
print(f"Learned {len(state.charts)} charts")
```

All components are NumPy-based and can be combined with SciPy tooling for further analysis or downstream tasks.

