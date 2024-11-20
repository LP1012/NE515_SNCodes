# SN1

The included code runs a simplified SN code with source term determined by the Method of Manufactured Solutions (MMS). 
Details can be found in the code's output.

## 
To test the libraries, run the following command: `g++ -I./src test_c.cpp ./src/*.cpp -o test_c && ./test_c'.

## Running 
The compile the code, use the following command: `g++ -I./src SN1_mms.cpp ./src/*.cpp -o SN1_mms`. 
To run the code from the SN1 directory, simply run `./SN1` in a terminal.
This all assumes that `g++` compiler is up to date.
