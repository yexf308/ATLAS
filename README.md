# ATLAS
code for the paper " Nonlinear model reduction for slow-fast stochastic systems near  manifolds" by Felix X.-F. Ye, Sichen Yang, Mauro Maggioni. 

The paper contains results for three example problems based on Pinched sphere(nick name: peanut), Oscillating half-moons and Butane model. Code for each example can be found in the respective subfolder in the examples folder.

In each example, there are two steps. The first is the ATLAS learning stage and the second is the mean resident time simulating stage.

## Parameter structure

All three drivers (`Butane.m`, `Halfmoon.m`, `Peanut.m`) now call their model-specific `set_parameter` function, which returns a structured configuration object named `params`. The struct exposes the settings that were previously injected into the base workspace so that downstream scripts can read and update values explicitly. The most relevant fields are:

* `params.RHS` – numerical parameters and function handles required by the right-hand-side integrators (e.g., time step `dt`, simulator horizon `T_max`, drift/diffusion handles, dimensionality `D`/`d`, thresholds, and any model-specific constants).
* `params.simulator` – simulator entry points, with `one`, `serial`, and `parallel` variants for single-trajectory, serial batch, and parallel batch integrations, respectively.
* `params.chart_sim_parameter` – the configuration used when invoking `ATLAS_simulator2`, including neighbourhood search settings and the evolving nearest landmark information.
* `params.relearn` – the relearning option flag (`option`), the mutable RHS configuration used during relearning (`RHS`), and the template relearning settings (`settings`).
* `params.learning` – high-level learning controls such as the number of initial landmarks (`K_int`), the long-trajectory horizon `T_one`, initial state `X_int`, and the default projection option (`option`).
* `params.simulation` – coarse controls for ATLAS sampling, such as the default number of time steps (`Nstep`) and sampling gap (`gap`).
* `params.paths` – filesystem paths for generated artefacts (chart data, transition matrices, analysis outputs, etc.), with `datapath` as the base directory.
* `params.MFPT` – metadata for mean first-passage time experiments (currently the number of initial conditions `N_IC`).

Scripts in `ATLAS/` and the model folders read from and write to these fields to keep shared state explicit and discoverable.
