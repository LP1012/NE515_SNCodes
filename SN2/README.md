# SN2 -- Final Project
The included code was developed to solve the problems outlined in the `NE515-520_F24_Project-Final.pdf` file, which is the final problem set for the course.

## Testing
To test the libraries, run the following command: `g++ -I./src test_c.cpp ./src/*.cpp -o test_c && ./test_c'.

## Cases
Two codes have been implemented.
The first, `SN2_freeSurfaces.cpp`, assumes zero flux on the boundaries and includes a normalized volumetric source (details can be found in the report).
The second, `SN2_reed.cpp`, solves the Reed problem with normalized sources.
In this, an albedo (reflected) case has been implemented along the left boundary, while the right boundary was designated as a vacuum condition.

## Running 
Code can be compiled and run in a similar manner as before.
Two bash scripts: `compile_programs` and `run` have also been created to streamline the running of the code.
Additionally, for these problems, the code has been generalized to all for an input text file (see examples included). 
For this, the slab can be made arbitrarily large and of as many sections as desirable.
However, there is an upper limit due to stack size limitations, so the user is cautioned to be wary of this when the total number of cells grows large.

## Plotting
A `make_plots.py` has been written to generate the plots of each run, taking the `.csv` files as input.
