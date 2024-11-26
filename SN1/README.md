# SN1

The included code runs a simplified SN code. 
Cases given have assumed isotropic sources (in the absolute simple case, this is just 0).
However, the code has been written to show where a non-isotropic source can be inputted.
Details for the run can be found in the code's output.

## Running 
The compile the code, use the following command: `g++ -I./src SN1.cpp ./src/*.cpp -o SN1`. 
To run the code from the SN1 directory, simply run `./SN1` in a terminal.
This all assumes that `g++` compiler is up to date.
For a streamlined case, use the `run` bash script included to compile and run the code and plot the output.

## Plotting
A `make_plots.py` has been written to generate a plot of the 100-cell computation, taking the `.csv` file as input.
