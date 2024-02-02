# Risk-Aware-Mean-Field-Coupled-LQ-Simulations-and-Plots
Files for the Python simulation and LaTeX plots for various quantities related to the optimal control of risk-aware mean-field coupled linear quadratic subsystems.

## Python Simulation.

### Running the simulations.

To run the simulations, simply run `python simulation.py` in the terminal (or `python3 simulation.py`, `python3.9 simulation.py`, etc, to specify the specific version of Python, or if the `python` command in your system does not refer to Python 3). The simulation requires the `numpy` and `csv` packages.

The system parameters all have the same names in the code as in the paper, and so they can be located and changed in lines 78 to 98. To change the specific $\lambda$ for which the state energy and control effort data is generated vs. time, change the range of $\lambda$ in the for loop at line 115. To change the $\lambda$ over which the data for the time average of the state energy and control effort quantities is generated, change the parameters set on lines 132 to 135.

The simulation will take many hours to complete; in particular the second set of simulations to find the time averages of the various quantities of controller vs. $\lambda$ from line 152 onward makes up the bulk of the simulation time. In order to run the simulation in a more reasonable time, I recommend adding `N = 100` at line 164, and setting `λ_num` on line 157 to a smaller value. This will make the simulation against the various $\lambda$ run faster by orders of magnitude, and the results of the second simulation are not significantly impacted by more precision in the $\lambda$ or a larger number of simulations $N$.

### Simulation output.

The simulation data is saved to various CSV files in the `Simulation Data` subdirectory.

Within the `Simulation Data` subdirectory, there will be one or more subdirectories named `lambda_[lambda]`, with the number `[lambda]` having its decimal point replaced by an underscore. This subdirectory contains data for the $N$ trials simulated with risk-parameter `[lambda]` vs. time $t$ stored in CSV files named `time_varying_[type]_[quantity].csv`. `[quantity]` is either “state energy” or “control effort”, and `[type]` is either `max` or `avg`, indicating whether the data is for the maximum or average of the quantity over the $N$ trials. For example, the file `Simulation Data/lambda_0_01/time_varying_avg_control_effort.csv` consists of the data for the average control effort among the subsystems at each time step $t$ when the system is simulated with risk-parameter $\lambda = 0.01$. These files are formatted as:

| t  | 0 | 0.5 | 1 | 1.5 | ⋯   | 99.5 | 100 | mean | stddev |
|:--:|---|-----|---|-----|-----|------|-----|------|--------|
| 0  |   |     |   |     |     |      |     |      |        |
| 1  |   |     |   |     |     |      |     |      |        |
| 2  |   |     |   |     |     |      |     |      |        |
| ⋮  |   |     |   |     |     |      |     |      |        |
| 50 |   |     |   |     |     |      |     |      |        |

where the first row consists of the headers, and the first column consists of the time stamps. The cell in row with time stamp $t$ and header $k$, for $k \in \{0, 0.5, ..., 100\}$ contains the $\frac{k}{100}$ th percentile data for the quantity `[quantity]` of type `[type]` at time $t$ among the $N$ simulations. The cells in the columns “mean” and “stddev" consist of the mean and standard deviation of the data among the $N$ simulations at time $t$.

The CSV files named `lambda_varying_[type]_[quantity].csv` contain the data of the time-average of the quantity `[quantity]` of type `[type]` over the appropriate time horizon for the $N$ trials vs. $\lambda$. For example, the file `lambda_varying_max_state_energy.csv` consists of the data for the time-average of the maximum state energy among the subsystems over a range of $\lambda$.

These files are formatted as:

| lambda  | 0 | 0.5 | 1 | 1.5 | ⋯   | 99.5 | 100 | mean | stddev |
|:-------:|---|-----|---|-----|-----|------|-----|------|--------|
| 1e-5    |   |     |   |     |     |      |     |      |        |
| ⋮       |   |     |   |     |     |      |     |      |        |
| 1e3     |   |     |   |     |     |      |     |      |        |

where the first row consists of the headers, and the first column consists of the $\lambda$ values. The cell in row with risk-parameter $\lambda$ and header $k$, for $k \in \{0, 0.5, ..., 100\}$ contains the $\frac{k}{100}$ th percentile data for the time-average of quantity `[quantity]` of type `[type]` when among the $N$ trials simulated with risk-parameter $\lambda$. The cells in the columns labelled “mean” and “stddev" consist of the mean and standard deviation of the data among the $N$ simulations at time $t$.

If you wish to generate and save different quantities related to the simulation, write the data to a 2D `numpy` array. Then call the `write_data_to_csv` function in the Python file with your new data and the appropriate independent variable to save the quantity in a similar form to the files above.

## Plotting.

The plots are generated using $\LaTeX$ and the TikZ and PGFPlots packages. The required compiler(s) and the packages can all be installed with [TeX Live](https://www.tug.org/texlive/). To generate the plots, compile the file `generate_plots.tex` with the command line argument `-shell-escape`. Note that this compilation may fail with an error `TeX capacity exceeded, sorry [main memory size=<some number>].` You can increase the memory available to LaTeX, but I would instead recommend simply using LuaLatex instead, which does not have these memory issues. However, the compilation will still be extremely slow. In order to reduce the compilation time, replace all the lines in `generate_plots.tex` of the form `\input{TikZ Figures\...\density.tikz} \\` with `\\` to disable the density plots, since it is these specific plots that take up most of the time.