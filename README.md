

# Code for bounding steady-state observables

This is the code for a paper.

## Compiling

You need MOSEK, SCS, Eigen, Gurobi to run this code, so it's a bit of a pain to compile. If you're not one of the referees of the paper I wouldn't strongly advise against bothering. But, if you have the above dependencies, you can compile the code by running
```
make
```

## Usage

The code can be ran by running:
```bash
./run --help
```
to display the various possible Lindbladians, objectives and constraints.

