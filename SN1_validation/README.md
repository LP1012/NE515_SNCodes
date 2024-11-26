# SN1 Verification
The included code runs a simplified SN code with source term determined by the Method of Manufactured Solutions (MMS). 
Details can be found in the code's output.

## 
To test the libraries, run the following command: `g++ -I./src test_c.cpp ./src/*.cpp -o test_c && ./test_c'.

## Running 
The compile the code, use the following command: `g++ -I./src SN1_mms.cpp ./src/*.cpp -o SN1_mms`. 
To run the code from the SN1 directory, simply run `./SN1` in a terminal.
This all assumes that `g++` compiler is up to date.
For a streamlined run, use the `run` bash executable provided.
It will compile and run the `SN1_mms script`, and it will plot the results.

## Plotting
A `make_plots.py` has been written to generate the plots of each run, taking the `.csv` files as input.
In addition, a scaling study has been performed, and mentioned script will export the results of this.
