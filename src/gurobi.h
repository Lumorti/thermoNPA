#pragma once
#include "poly.h"

// Convert to Gurobi form and solve
std::pair<double,double> solveGurobi(Poly obj, std::vector<Poly> constraintsZero, int verbosity=1);
