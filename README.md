

# Code for bounding steady-state observables

This is the code for a paper.

## Compiling

You need MOSEK, SCS, Eigen, Gurobi to run this code, so it's a bit of a pain to compile. Unless you're one of the referees of the paper I would strongly advise against bothering. But, if you have the above dependencies and have set the various MOSEKHOME etc. as defined in the makefile, you can compile the code by running
```
make
```

## Usage

The code can be ran using:
```bash
./run --help
```
to display the various possible Lindbladians, objectives and constraints.

