# Risk-Aware-Mean-Field-Coupled-LQ-Simulations-and-Plots
Files for the Python simulation and LaTeX plots for various quantities related to the optimal control of risk-aware mean-field coupled linear quadratic subsystems.

## Python Simulation.

### Running the simulations.

To run the simulations, simply run `python simulation.py` in the terminal (or `python3 simulation.py`, `python3.9 simulation.py`, etc, to specify the specific version of Python, or if the `python` command in your system does not refer to Python 3). The simulation requires the `numpy` and `csv` packages.

The system parameters all have the same names in the code as in the paper, and so they can be located and changed in lines 78 to 98. To change the specific $\lambda$ for which the state energy and control effort data is generated vs. time, change the range of $\lambda$ defined on line 129. To change the $\lambda$ over which the data for the time average of the state energy and control effort quantities is generated, change the parameters set on lines 164 to 167.

The simulation will take many hours to complete; in particular the second set of simulations to find the time averages of the various quantities of controller vs. $\lambda$ from line 161 onward makes up the bulk of the simulation time. In order to run the simulation in a more reasonable time, I recommend adding `N = 100` at line 168, and setting `λ_num` on line 166 to a smaller value. This will make the simulation against the various $\lambda$ run faster by orders of magnitude, and the results of the second simulation are not significantly impacted by more precision in the $\lambda$ or a larger number of simulations $N$.

### Simulation output.

The simulation data is saved to various CSV files in a `data` subdirectory. The `data` subdirectory is subdivided further into two subdirectories, `state_energy` and `control_effort`, containing the data for the state energy and control effort, respectively. Both of these are also subdivided into two subdirectories, `vs_time` and `vs_lambda`, containing the information about the quantity (state energy or control effort) as it varies with $t$ or $\lambda$.

The `vs_lambda` subdirectory will have two files, `avg.csv` and `max.csv`, containing data of the time average of the average or maximum of the quantity (state energy or control effort) across the subsystems. These CSV files are formatted as

| lambda  | 0 | 0.5 | 1 | 1.5 | ⋯   | 99.5 | 100 | mean | stddev |
|:-------:|---|-----|---|-----|-----|------|-----|------|--------|
| 1e-5    |   |     |   |     |     |      |     |      |        |
| ⋮       |   |     |   |     |     |      |     |      |        |
| 1e3     |   |     |   |     |     |      |     |      |        |

where the first row consists of the headers, and the first column consists of the $\lambda$ values. The cell in row with risk-parameter $\lambda$ and header $k$, for $k \in \{0, 0.5, ..., 100\}$ contains the $\frac{k}{100}$ th quantile data for the time-average of the average or maximum of the quantity among the $N$ trials simulated with risk-parameter $\lambda$. The cells in the columns labelled “mean” and “stddev" consist of the mean and standard deviation of the data among the $N$ simulations at time $t$.

The `vs_time` subdirectory will be subdivided into two subdirectories, `avg` and `max`, containing data of the time average of the average or maximum of the quantity (state energy or control effort) across the subsystems. These subdirectories will have files named `[lambda].csv`, where `lambda` is some decimal number (with the '.' replaced by '_') indicating the risk-parameter used to generate the data. These CSV files are formatted as:

| t  | 0 | 0.5 | 1 | 1.5 | ⋯   | 99.5 | 100 | mean | stddev |
|:--:|---|-----|---|-----|-----|------|-----|------|--------|
| 0  |   |     |   |     |     |      |     |      |        |
| 1  |   |     |   |     |     |      |     |      |        |
| 2  |   |     |   |     |     |      |     |      |        |
| ⋮  |   |     |   |     |     |      |     |      |        |
| 50 |   |     |   |     |     |      |     |      |        |

where the first row consists of the headers, and the first column consists of the time stamps. The cell in row with time stamp $t$ and header $k$, for $k \in \{0, 0.5, ..., 100\}$ contains the $\frac{k}{100}$ th percentile data for the average or maximum of quantity at time $t$ among the $N$ simulations with risk-parameter given by the file name. The cells in the columns “mean” and “stddev" consist of the mean and standard deviation of the data among the $N$ simulations at time $t$.

## Plotting.

The plots are generated using $\LaTeX$ and the TikZ and PGFPlots packages. The required compiler(s) and the packages can all be installed with [TeX Live](https://www.tug.org/texlive/). To generate the plots, compile the file `generate_plots.tex` with the command line argument `-shell-escape` (do not forget to include this argument, it will **not** work otherwise). The compilation may fail if you do not create a `figures` subdirectory within your build destination directory (instead of making a new subdirectory, you may also simply comment out line 20 of `generate_plots.tex`, however this might cause your repo to get polluted with many LaTeX build artifacts). The compilation may also fail with an error `TeX capacity exceeded, sorry [main memory size=<some number>].` You can increase the memory available to LaTeX, but I would instead recommend simply using LuaLaTeX instead, which do not have these memory issues. However, the compilation will still be extremely slow. In order to reduce the compilation time, comment out each line that is of the form `\input{templates/density.tikz}` to disable generating the density plots, since it is these specific plots that take up most of the time.

Remember to generate the data using the simulation before compiling the file to generate the plots! The simulation file will store the data into the correct directories.