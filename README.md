# FastICA - PULP
Implementation of FastICA algorithm fot PULP platform

## Usage
First of all, use `docker-compose run --rm pulp` to run the container and open a shell.
Then, move to the directory containing the project (i.e., `~/fast-ica`) and compile it using `make clean all`; you can set the following variables:

- `COMP`: number of independent components;
- `OBS`: number of observations;
- `USE_CLUSTER`: whether to offload the task to the cluster or use only the Fabric Controller;
- `CORES`: number of cores to use (from 1 to 8);
- `STRATEGY`: FastICA strategy (0 for Parallel and 1 for Deflation);
- `G_FUNC`: G function to use in FastICA (0 for LogCosh, 1 for Exp and 2 for Cube);
- `MAX_ITER`: maximum number of iterations;
- `SAMPLES`: number of samples of the signal to generate;
- `WINDOW_SIZE`: temporal resolution of the signal to generate;
- `ADD_NOISE`: whether to add Gaussian noise to the generated signal or not;
- `VERB`: whether to print the original signal, the observations and the recovered signal to standard output.

Finally, the program can be run with `make run`.

The results can be visualized using the `pyutils/plot_signals.py REPORT.LOG N_COMP N_REC` program.
