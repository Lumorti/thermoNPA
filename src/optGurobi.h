#pragma once
#include "poly.h"

// Convert to Gurobi form and solve
std::pair<double,double> boundGurobi(Poly obj={}, std::vector<Poly> constraintsZero={}, int verbosity=1);
