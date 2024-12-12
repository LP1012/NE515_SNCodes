# NE515_SNCodes
Contains code developed for NE515: Radiation Interaction and Transport offered at the University of New Mexico in Fall 2024.

## Codes
The codes written were prepared to submit for the class mentioned.
They are not intended to serve as production or research codes -- simply an academic exploration of a particular numerical method to solve the 1-D linear Boltzmann transport equation.
A known limitation of the code is that will fail when the number of cells exceeds a certain value. 
This is because of C++ running out of memory in the array containing the angular flux values.
This can likely be improved by transitioning from arrays to vectors.

### SN1
The `SN1.cpp` code was the first to be created, and it used for very simple cases of slab geometry.

### SN1 Verification
The `SN1_mms.cpp` code located `SN1_vericication` directory is an implementation of the Method of Manufactured Solutions to study the accuracy and convergence of the developed numerical scheme.

### SN2 
The `SN2.cpp` code greatly improves on the `SN1` by allowing for a 1-D slab with arbitrary layers and variable physical values. 
This allowed for implementation of the so-called Reed problem, which is a benchmark case for 1-D slab geometry that greatly tests the solver's capabilities.

## Report
A report has been supplied along with the given code.
Please consult this document for a general feel of the methods used in the supplied programs.