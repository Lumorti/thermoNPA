#pragma once
#include <vector>
#include "poly.h"

// Convert to SCS form and solve
double solveSCS(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero, std::vector<Poly> constraintsPositive, int verbosity, std::pair<int,int> varBounds = {-1, 1}, std::map<Mon, std::complex<double>>* solution = nullptr);

// Convert to SCS form and bound
std::pair<double,double> boundSCS(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero={}, std::vector<Poly> constraintsPositive={}, int verbosity=1, std::pair<int,int> varBounds = {-1, 1}, std::map<Mon, std::complex<double>>* solution = nullptr);


