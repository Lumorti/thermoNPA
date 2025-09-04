// Standard includes
#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <unordered_set>
#include <set>
#include <chrono>
#include <map>
#include <fstream>
#include <random>

// Import OpenMP
#include <omp.h>

// Import Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

// PolyNC
#include "../../PolyNC/src/polync.h"
 
// Pauli definitions
Eigen::SparseMatrix<std::complex<double>> pauliI(2,2);
Eigen::SparseMatrix<std::complex<double>> pauliX(2,2);
Eigen::SparseMatrix<std::complex<double>> pauliY(2,2);
Eigen::SparseMatrix<std::complex<double>> pauliZ(2,2);
std::map<int, std::string> letterMap = {{0, "I"}, {1, "X"}, {2, "Y"}, {3, "Z"}};

// https://stackoverflow.com/questions/2647858/multiplying-complex-with-constant-in-c
template <typename T>
struct identity_t { typedef T type; };
#define COMPLEX_OPS(OP)                                                 \
  template <typename _Tp>                                               \
  std::complex<_Tp>                                                     \
  operator OP(std::complex<_Tp> lhs, const typename identity_t<_Tp>::type & rhs) \
  {                                                                     \
    return lhs OP rhs;                                                  \
  }                                                                     \
  template <typename _Tp>                                               \
  std::complex<_Tp>                                                     \
  operator OP(const typename identity_t<_Tp>::type & lhs, const std::complex<_Tp> & rhs) \
  {                                                                     \
    return lhs OP rhs;                                                  \
  }
COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)
#undef COMPLEX_OPS

// https://stackoverflow.com/questions/900326/how-do-i-elegantly-format-string-in-c-so-that-it-is-rounded-to-6-decimal-place
std::string strRound(double value, int decimals) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(decimals) << value;
    std::string s = ss.str();
    if (decimals > 0 && s[s.find_last_not_of('0')] == '.') {
        s.erase(s.size() - decimals + 1);
    }
    return s;
}

// Useful constants
const std::complex<double> imag(0, 1);

// Trace of an eigen sparse matrix
std::complex<double> trace(const Eigen::SparseMatrix<std::complex<double>>& mat) {
    std::complex<double> trace = 0;
    for (int k=0; k<mat.cols(); ++k) {
        trace += mat.coeff(k, k);
    }
    return trace;
}

// Generate a Pauli string matrix on demand
Eigen::SparseMatrix<std::complex<double>> generatePauliMatrix(std::vector<int> pauliInds) {

    // The matrix starts as 1x1 identity
    Eigen::SparseMatrix<std::complex<double>> mat(1, 1);
    mat.insert(0, 0) = 1;

    // For each ind, tensor product with that Pauli
    for (size_t i=0; i<pauliInds.size(); ++i) {
        Eigen::SparseMatrix<std::complex<double>> pauliMat;
        if (pauliInds[i] == 0) {
            pauliMat = pauliI;
        } else if (pauliInds[i] == 1) {
            pauliMat = pauliX;
        } else if (pauliInds[i] == 2) {
            pauliMat = pauliY;
        } else if (pauliInds[i] == 3) {
            pauliMat = pauliZ;
        } else {
            throw std::invalid_argument("Invalid Pauli index: " + std::to_string(pauliInds[i]));
        }
        mat = Eigen::kroneckerProduct(mat, pauliMat).eval();
    }
    mat.makeCompressed();

    // Return the resulting matrix
    return mat;

}

// Partial trace of a density matrix
// Here we have a matrix of size 2^N x 2^N, and we want to trace out the siteInd-th qubit
// the site index is zero-based
std::vector<std::vector<Poly>> partialTrace(std::vector<std::vector<Poly>> mat, int siteInd) {

    // Get the size of the matrix
    int size = mat.size();
    int newSize = size / 2;
    siteInd = std::log2(size) - siteInd - 1;

    // Initialize the new matrix
    std::vector<std::vector<Poly>> newMat(newSize, std::vector<Poly>(newSize, Poly()));

    // Perform partial trace
    for (int i = 0; i < newSize; ++i) {                                                         
        for (int j = 0; j < newSize; ++j) {                                                     
            int i1 = i / (1 << siteInd);
            int i2 = i % (1 << siteInd);
            int j1 = j / (1 << siteInd);
            int j2 = j % (1 << siteInd);
            for (int k = 0; k < 2; ++k) {
                int row = (i1 << (siteInd + 1)) | (k << siteInd) | i2;
                int col = (j1 << (siteInd + 1)) | (k << siteInd) | j2;
                newMat[i][j] += mat[row][col];
            }
        }                                                                                       
    }

    // Return the new matrix
    return newMat;

}

// Generic entry function
int main(int argc, char* argv[]) {

    // Random seed unless specified
    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(NULL));

    // Define the scenario
    std::string note = "";
    bool allTimeInMillis = false;
    bool usingHeatCurrent = false;
    bool benchmark = false;
    bool checkObj = false;
    bool symSample = false;
    bool useKnown = false;
    bool allSymmetries = false;
    std::vector<Poly> zeroConsForSampling;
    std::string stateFile = "state.dat";
    std::string autoType = "first";
    std::string autoMomentType = "first";
    int level = 0;
    int precision = 10;
    int lindbladLevel = 0;
    bool precompute = false;
    bool precomputed = false;
    int numQubits = 1;
    int gridWidth = 1;
    int gridHeight = 1;
    int imagType = 0;
    long int numSamplesPer = 0;
    long int numShots = 0;
    double tol = 1e-7;
    bool usePurity = false;
    bool useEnergy = false;
    bool groundStateProblem = false;
    bool outputToFile = false;
    std::string sampleChoice = "all";
    int maxSampleDegree = 2;
    int maxPaulis = 1000000;
    std::vector<int> includeOnly;
    bool excludeObjective = false;
    bool excludeX = false;
    bool excludeY = false;
    bool excludeZ = false;
    bool outputLindbladAsLaTeX = false;
    double percentile = 95;
    Poly objective("<Z1>");
    std::map<Mon, double> samples;
    Poly pseudoObjective;
    int verbosity = 1;
    std::string solver = "auto";
    std::string modelName = "default";
    std::complex<double> knownIdeal = 0.0;
    bool idealIsKnown = false;
    bool findMinimal = false;
    bool tryRemove = false;
    bool traceCons = false;
    int reconLevel = 0;
    int autoMomentAmount = 0;
    int findMinimalAmount = 0;
    int reductiveCons = 0;
    std::string seed = "";
    Poly lindbladian("<X1A0X1>+<Y1A0Y1>-<A0>");
    std::vector<std::vector<Poly>> hamiltonianInter;
    Poly lindbladianHot;
    Poly lindbladianCold;
    std::vector<std::string> extraMonomials;
    std::vector<std::string> extraMonomialsLim;
    std::vector<Poly> constraintsZero;
    std::vector<Poly> constraintsPositive;
    std::vector<std::pair<std::vector<int>,std::vector<int>>> symmetries;

    // Pauli definitions
    pauliI.insert(0,0) = 1;
    pauliI.insert(1,1) = 1;
    pauliI.makeCompressed();
    pauliX.insert(0,1) = 1;
    pauliX.insert(1,0) = 1;
    pauliX.makeCompressed();
    pauliY.insert(0,1) = std::complex<double>(0,-1);
    pauliY.insert(1,0) = std::complex<double>(0,1);
    pauliY.makeCompressed();
    pauliZ.insert(0,0) = 1;
    pauliZ.insert(1,1) = -1;
    pauliZ.makeCompressed();

    // Set the precision
    std::cout << std::setprecision(precision);

    // The max physical cores
    int numCores = omp_get_max_threads();

    // Process command-line args
    for (int i=1; i<argc; i++) {
        std::string argAsString = std::string(argv[i]);

        // Set the level of the moment matrix
        if (argAsString == "-m") {
            level = std::stoi(argv[i+1]);
            i++;

        // Manually set the objective
        } else if (argAsString == "-O" || argAsString == "--obj") {
            objective = Poly(std::string(argv[i+1]));
            i++;

        // Manually set the Lindbladian
        } else if (argAsString == "-L") {
            lindbladian = Poly(std::string(argv[i+1]));
            i++;

        // Pauli X objective
        } else if (argAsString == "--objX" || argAsString == "--objx") {
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<X" + std::to_string(i) + ">");
            }

        // Pauli Y objective
        } else if (argAsString == "--objY" || argAsString == "--objy") {
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Y" + std::to_string(i) + ">");
            }

        // Pauli Z objective
        } else if (argAsString == "--objZ" || argAsString == "--objx") {
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }

        // Combination of all Pauli objectives
        } else if (argAsString == "--objXYZ" || argAsString == "--objxyz") {
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<X" + std::to_string(i) + ">");
                objective += Poly(1.0/numQubits, "<Y" + std::to_string(i) + ">");
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }

        // Purity objective
        } else if (argAsString == "--objPurity" || argAsString == "--objpurity") {
            usePurity = true;
            solver = "mosek";

        // Energy objective
        } else if (argAsString == "--objEnergy" || argAsString == "--objenergy") {
            useEnergy = true;

        // Objective is made of random Paulis of some order
        } else if (argAsString == "--objRand" || argAsString == "--objrand") {
            int order = std::stoi(argv[i+1]);
            i++;
            objective = Poly();
            if (order >= 1) {
                for (int i=0; i<numQubits; i++) {
                    objective += Poly(1.0/numQubits, "<X" + std::to_string(i+1) + ">");
                    objective += Poly(1.0/numQubits, "<Y" + std::to_string(i+1) + ">");
                    objective += Poly(1.0/numQubits, "<Z" + std::to_string(i+1) + ">");
                }
            }
            if (order >= 2) {
                for (int i=0; i<numQubits; i++) {
                    for (int j=i+1; j<numQubits; j++) {
                        objective += Poly(1.0/numQubits, "<X" + std::to_string(i+1) + "X" + std::to_string(j+1) + ">");
                        objective += Poly(1.0/numQubits, "<Y" + std::to_string(i+1) + "Y" + std::to_string(j+1) + ">");
                        objective += Poly(1.0/numQubits, "<Z" + std::to_string(i+1) + "Z" + std::to_string(j+1) + ">");
                    }
                }
            }
            objective.randomize();
            objective.normalize(true);

        // Setting the tolerance
        } else if (argAsString == "-t") {
            tol = std::stod(argv[i+1]);
            i++;

        // Pauli Lindbladian
        } else if (argAsString == "--pauli") {
            if (i+3 >= argc) {
                std::cout << "Not enough arguments for --pauli" << std::endl;
                return 1;
            }
            double coeff1 = std::stod(argv[i+1]);
            double coeff2 = std::stod(argv[i+2]);
            double coeff3 = std::stod(argv[i+3]);
            if (coeff1 < 0 || coeff2 < 0 || coeff3 < 0) {
                std::cout << "Coefficients must be non-negative" << std::endl;
                return 1;
            }
            std::string pauliString = "";
            pauliString += "+" + std::to_string(coeff1) + "<X1A0X1>";
            pauliString += "+" + std::to_string(coeff2) + "<Y1A0Y1>";
            pauliString += "+" + std::to_string(coeff3) + "<Z1A0Z1>";
            pauliString += "-" + std::to_string(coeff1+coeff2+coeff3) + "<A0>";
            lindbladian = Poly(pauliString);
            i += 3;

        // Lindbladian "breaking" the second law
        } else if (argAsString == "--second") {

            // Params
            double omega_h = 10;
            double omega_c = 5;
            double T_h = 12;
            double T_c = 10;
            double kappa = 1e-4;
            numQubits = std::stoi(argv[i+1]);
            double epsilon = std::stod(argv[i+2]);
            i += 2;

            // Calculated quantities
            double beta_h = 1.0 / T_h;
            double beta_c = 1.0 / T_c;
            double gamma_h = kappa * std::pow(omega_h, 3) * (1 - std::exp(-beta_h * omega_h));
            double gamma_c = kappa * std::pow(omega_c, 3) * (1 - std::exp(-beta_c * omega_c));

            // Creation/annihilation operators
            std::vector<Poly> as(numQubits, Poly());
            std::vector<Poly> adags(numQubits, Poly());
            for (int i=0; i<numQubits; i++) {
                std::string XString = "<X" + std::to_string(i+1) + ">";
                std::string YString = "<Y" + std::to_string(i+1) + ">";
                for (int j=0; j<numQubits; j++) {
                    if (j != i) {
                        XString += "Z" + std::to_string(j+1);
                        YString += "Z" + std::to_string(j+1);
                    }
                }
                as[i] = Poly(0.5, XString) - Poly(0.5*imag, YString);
                adags[i] = Poly(0.5, XString) + Poly(0.5*imag, YString);
            }

            // Hamiltonian
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            for (int i=0; i<numQubits; i++) {
                hamiltonianInter[i][i] = epsilon * as[i] * adags[i];
                hamiltonianInter[i][i].reduce();
            }
            for (int i=0; i<numQubits-1; i++) {
                hamiltonianInter[i][i+1] = epsilon * (adags[i] * as[i+1] + as[i] * adags[i+1]);
                hamiltonianInter[i][i+1].reduce();
                hamiltonianInter[i+1][i] = hamiltonianInter[i][i+1];
            }
            Poly H = Poly();
            for (int i=0; i<numQubits; i++) {
                for (int j=0; j<numQubits; j++) {
                    H += hamiltonianInter[i][j];
                }
            }
            H.reduce();

            // The full Lindbladian
            Poly rho("<R1>");
            lindbladianHot = Poly();
            lindbladianHot += gamma_h * (adags[0] * rho * as[0] - 0.5*as[0]*adags[0]*rho - 0.5*rho*as[0]*adags[0]);
            lindbladianHot += gamma_h * std::exp(-beta_h*omega_h) * (as[0] * rho * adags[0] - 0.5*adags[0]*as[0]*rho - 0.5*rho*adags[0]*as[0]);
            lindbladianCold = Poly();
            lindbladianCold += gamma_c * (adags[numQubits-1] * rho * as[numQubits-1] - 0.5*as[numQubits-1]*adags[numQubits-1]*rho - 0.5*rho*as[numQubits-1]*adags[numQubits-1]);
            lindbladianCold += gamma_c * std::exp(-beta_c*omega_c) * (as[numQubits-1] * rho * adags[numQubits-1] - 0.5*adags[numQubits-1]*as[numQubits-1]*rho - 0.5*rho*adags[numQubits-1]*as[numQubits-1]);
            lindbladian = -imag*H.commutator(rho) + lindbladianHot + lindbladianCold;

            // Simplify the Lindbladians
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladianHot = Poly("<A0>") * lindbladianHot;
            lindbladianCold = Poly("<A0>") * lindbladianCold;
            lindbladianHot.cycleToAndRemove('R', 1);
            lindbladianCold.cycleToAndRemove('R', 1);
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.reduce();
            lindbladianHot.reduce();
            lindbladianCold.reduce();

            // Simple objective
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }

        // Two-body Lindbladian
        } else if (argAsString == "--two" || argAsString == "--twov") {
            modelName = argAsString;

            // Defining quantities
            double gamma_c = 1.1e-2;
            double gamma_h = 1e-3;
            double g = 1.6e-3;
            double T_h = 1.0;
            double T_c = 0.1;
            double delta = 0.005;
            double epsilon_h = 1.0;

            // If the arg is --twov, we have values
            if (argAsString == "--twov") {
                T_h = std::stod(argv[i+1]);
                delta = std::stod(argv[i+2]);
                i += 2;
            }

            // Calculated quantities
            double epsilon_c = epsilon_h + delta;
            double n_c = 1.0 / (std::exp(epsilon_c / T_c) - 1.0);
            double n_h = 1.0 / (std::exp(epsilon_h / T_h) - 1.0);
            double gamma_h_plus = gamma_h * n_h;
            double gamma_h_minus = gamma_h * (n_h + 1.0);
            double gamma_c_plus = gamma_c * n_c;
            double gamma_c_minus = gamma_c * (n_c + 1.0);
            double Gamma_h = gamma_h_plus + gamma_h_minus;
            double Gamma_c = gamma_c_plus + gamma_c_minus;
            double Gamma = Gamma_h + Gamma_c;
            double chi = (4.0*g*g+Gamma_h*Gamma_c)*Gamma*Gamma + 4.0*delta*delta*Gamma_h*Gamma_c;

            // The rho from the paper
            Eigen::MatrixXcd rho = Eigen::MatrixXcd::Zero(4,4);
            rho(0,0) = (4.0*g*g*(gamma_h_plus + gamma_c_plus)*(gamma_h_plus + gamma_c_plus)  + gamma_h_plus*gamma_c_plus*(Gamma*Gamma + 4.0*delta*delta)) / chi;
            rho(1,1) = (4.0*g*g*(gamma_h_minus + gamma_c_minus)*(gamma_h_plus + gamma_c_plus)  + gamma_h_plus*gamma_c_minus*(Gamma*Gamma + 4.0*delta*delta)) / chi;
            rho(2,2) = (4.0*g*g*(gamma_h_minus + gamma_c_minus)*(gamma_h_plus + gamma_c_plus)  + gamma_h_minus*gamma_c_plus*(Gamma*Gamma + 4.0*delta*delta)) / chi;
            rho(3,3) = (4.0*g*g*(gamma_h_minus + gamma_c_minus)*(gamma_h_minus + gamma_c_minus)  + gamma_h_minus*gamma_c_minus*(Gamma*Gamma + 4.0*delta*delta)) / chi;
            rho(1,2) = (2.0*g*(gamma_h_plus*gamma_c_minus - gamma_h_minus*gamma_c_plus)*(imag*Gamma-2.0*delta)) / chi;
            rho(2,1) = std::conj(rho(1,2));

            // Constants
            Eigen::MatrixXcd sigma_z = Eigen::MatrixXcd::Zero(2,2);
            sigma_z(0,0) = 1.0;
            sigma_z(1,1) = -1.0;
            Eigen::MatrixXcd sigma_y = Eigen::MatrixXcd::Zero(2,2);
            sigma_y(0,1) = -imag;
            sigma_y(1,0) = imag;
            Eigen::MatrixXcd sigma_x = Eigen::MatrixXcd::Zero(2,2);
            sigma_x(0,1) = 1;
            sigma_x(1,0) = 1;
            Eigen::MatrixXcd sigma_plus = Eigen::MatrixXcd::Zero(2,2);
            sigma_plus(0,1) = 1.0;
            Eigen::MatrixXcd sigma_minus = Eigen::MatrixXcd::Zero(2,2);
            sigma_minus(1,0) = 1.0;
            Eigen::MatrixXcd eye = Eigen::MatrixXcd::Identity(2,2);
            Eigen::MatrixXcd eye4 = Eigen::MatrixXcd::Identity(4,4);
            Eigen::MatrixXcd sigma_h_plus = kroneckerProduct(sigma_plus, eye);
            Eigen::MatrixXcd sigma_h_minus = kroneckerProduct(sigma_minus, eye);
            Eigen::MatrixXcd sigma_c_plus = kroneckerProduct(eye, sigma_plus);
            Eigen::MatrixXcd sigma_c_minus = kroneckerProduct(eye, sigma_minus);
            Eigen::MatrixXcd sigma_z_h = kroneckerProduct(sigma_z, eye);
            Eigen::MatrixXcd sigma_z_c = kroneckerProduct(eye, sigma_z);
            Eigen::MatrixXcd sigma_h_z = sigma_z_h;
            Eigen::MatrixXcd sigma_c_z = sigma_z_c;
            Eigen::MatrixXcd sigma_h_y = kroneckerProduct(sigma_y, eye);
            Eigen::MatrixXcd sigma_c_y = kroneckerProduct(eye, sigma_y);
            Eigen::MatrixXcd sigma_h_x = kroneckerProduct(sigma_x, eye);
            Eigen::MatrixXcd sigma_c_x = kroneckerProduct(eye, sigma_x);

            // Calculate the ideal J_ss
            Eigen::MatrixXcd H_s = epsilon_h*sigma_h_plus*sigma_h_minus 
                                   + epsilon_c*sigma_c_plus*sigma_c_minus;  
            Eigen::MatrixXcd term2 = gamma_h_plus*(sigma_h_plus*rho*sigma_h_minus 
                                                   - 0.5*sigma_h_minus*sigma_h_plus*rho 
                                                   - 0.5*rho*sigma_h_minus*sigma_h_plus) 
                                   + gamma_h_minus*(sigma_h_minus*rho*sigma_h_plus 
                                                   - 0.5*sigma_h_plus*sigma_h_minus*rho 
                                                   - 0.5*rho*sigma_h_plus*sigma_h_minus);
            double exp_sigma_z_h = tr(sigma_z_h*rho);
            double exp_sigma_z_c = tr(sigma_z_c*rho);
            double Q_ss_h = 0.5*epsilon_h*(gamma_h_plus-gamma_h_minus) - 0.5*epsilon_h*(gamma_h_plus+gamma_h_minus)*exp_sigma_z_h;
            double Q_ss_c = 0.5*epsilon_c*(gamma_c_plus-gamma_c_minus) - 0.5*epsilon_c*(gamma_c_plus+gamma_c_minus)*exp_sigma_z_c;
            double J_ss_from_rho = Q_ss_h - Q_ss_c;

            // Other system settings
            numQubits = 2;
            knownIdeal = J_ss_from_rho;
            idealIsKnown = true;

            // Construct the objective as a polynomial with plus/minus mats
            objective = Poly();
            objective += Poly(epsilon_h*gamma_h_plus,       "<M1P1M1P1>");
            objective += Poly(-0.5*epsilon_h*gamma_h_plus,  "<P1M1M1P1>");
            objective += Poly(-0.5*epsilon_h*gamma_h_plus,  "<M1P1P1M1>");
            objective += Poly(epsilon_h*gamma_h_minus,      "<P1P1M1M1>");
            objective += Poly(-0.5*epsilon_h*gamma_h_minus, "<P1M1P1M1>");
            objective += Poly(-0.5*epsilon_h*gamma_h_minus, "<P1M1P1M1>");
            objective += Poly(epsilon_c*gamma_h_plus,       "<M1P2M2P1>");
            objective += Poly(-0.5*epsilon_c*gamma_h_plus,  "<P2M2M1P1>");
            objective += Poly(-0.5*epsilon_c*gamma_h_plus,  "<M1P1P2M2>");
            objective += Poly(epsilon_c*gamma_h_minus,      "<P1P2M2M1>");
            objective += Poly(-0.5*epsilon_c*gamma_h_minus, "<P2M2P1M1>");
            objective += Poly(-0.5*epsilon_c*gamma_h_minus, "<P1M1P2M2>");
            objective += Poly(-epsilon_h*gamma_c_plus,       "<M2P1M1P2>");
            objective += Poly(0.5*epsilon_h*gamma_c_plus,  "<P1M1M2P2>");
            objective += Poly(0.5*epsilon_h*gamma_c_plus,  "<M2P2P1M1>");
            objective += Poly(-epsilon_h*gamma_c_minus,      "<P2P1M1M2>");
            objective += Poly(0.5*epsilon_h*gamma_c_minus, "<P1M1P2M2>");
            objective += Poly(0.5*epsilon_h*gamma_c_minus, "<P2M2P1M1>");
            objective += Poly(-epsilon_c*gamma_c_plus,       "<M2P2M2P2>");
            objective += Poly(0.5*epsilon_c*gamma_c_plus,  "<P2M2M2P2>");
            objective += Poly(0.5*epsilon_c*gamma_c_plus,  "<M2P2P2M2>");
            objective += Poly(-epsilon_c*gamma_c_minus,      "<P2P2M2M2>");
            objective += Poly(0.5*epsilon_c*gamma_c_minus, "<P2M2P2M2>");
            objective += Poly(0.5*epsilon_c*gamma_c_minus, "<P2M2P2M2>");
            objective.convertToPaulis();
            
            // Construct the Limbadlian as a polynomial from plus/minus
            lindbladian = Poly();
            lindbladian += Poly(-epsilon_h*imag, "<A0P1M1>");
            lindbladian += Poly(-epsilon_c*imag, "<A0P2M2>");
            lindbladian += Poly(-g*imag, "<A0P1M2>");
            lindbladian += Poly(-g*imag, "<A0M1P2>");
            lindbladian += Poly(epsilon_h*imag, "<P1M1A0>");
            lindbladian += Poly(epsilon_c*imag, "<P2M2A0>");
            lindbladian += Poly(g*imag, "<P1M2A0>");
            lindbladian += Poly(g*imag, "<M1P2A0>");
            lindbladian += Poly(gamma_h_plus, "<M1A0P1>");
            lindbladian += Poly(-0.5*gamma_h_plus, "<A0M1P1>");
            lindbladian += Poly(-0.5*gamma_h_plus, "<M1P1A0>");
            lindbladian += Poly(gamma_h_minus, "<P1A0M1>");
            lindbladian += Poly(-0.5*gamma_h_minus, "<A0P1M1>");
            lindbladian += Poly(-0.5*gamma_h_minus, "<P1M1A0>");
            lindbladian += Poly(gamma_c_plus, "<M2A0P2>");
            lindbladian += Poly(-0.5*gamma_c_plus, "<A0M2P2>");
            lindbladian += Poly(-0.5*gamma_c_plus, "<M2P2A0>");
            lindbladian += Poly(gamma_c_minus, "<P2A0M2>");
            lindbladian += Poly(-0.5*gamma_c_minus, "<A0P2M2>");
            lindbladian += Poly(-0.5*gamma_c_minus, "<P2M2A0>");
            lindbladian.convertToPaulis();

            // Note the local and non-local Hamiltonians in case we need them later
            hamiltonianInter = std::vector<std::vector<Poly>>(2, std::vector<Poly>(2, Poly()));
            hamiltonianInter[0][0] = Poly(epsilon_h, "<P1M1>");
            hamiltonianInter[0][1] = Poly(g, "<P1M2>") + Poly(g, "<M1P2>");
            hamiltonianInter[1][0] = hamiltonianInter[0][1];
            hamiltonianInter[1][1] = Poly(epsilon_c, "<P2M2>");
            for (int i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    hamiltonianInter[i][j].convertToPaulis();
                }
            }
            lindbladianHot = Poly();
            lindbladianHot += Poly(gamma_h_plus, "<M1A0P1>");
            lindbladianHot += Poly(-0.5*gamma_h_plus, "<A0M1P1>");
            lindbladianHot += Poly(-0.5*gamma_h_plus, "<M1P1A0>");
            lindbladianHot += Poly(gamma_h_minus, "<P1A0M1>");
            lindbladianHot += Poly(-0.5*gamma_h_minus, "<A0P1M1>");
            lindbladianHot += Poly(-0.5*gamma_h_minus, "<P1M1A0>");
            lindbladianHot.convertToPaulis();
            lindbladianCold = Poly();
            lindbladianCold += Poly(gamma_c_plus, "<M2A0P2>");
            lindbladianCold += Poly(-0.5*gamma_c_plus, "<A0M2P2>");
            lindbladianCold += Poly(-0.5*gamma_c_plus, "<M2P2A0>");
            lindbladianCold += Poly(gamma_c_minus, "<P2A0M2>");
            lindbladianCold += Poly(-0.5*gamma_c_minus, "<A0P2M2>");
            lindbladianCold += Poly(-0.5*gamma_c_minus, "<P2M2A0>");
            lindbladianCold.convertToPaulis();

        // Majumdar-Ghosh model
        } else if (argAsString == "--mg") {
            modelName = argAsString;

            // Defining quantities
            double J = 1;

            // We should be given the number of qubits
            numQubits = std::stoi(argv[i+1]);
            i++;
        
            // H = J \sum S_j S_{j+1} + (J/2) \sum S_j S_{j+2}
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            for (int i=0; i<numQubits; i++) {
                int ind1 = i;
                int ind2 = (i+1) % numQubits;
                hamiltonianInter[ind1][ind2] += Poly(J/4.0, "<X" + std::to_string(ind1+1) + "X" + std::to_string(ind2+1) + ">")
                                             + Poly(J/4.0, "<Y" + std::to_string(ind1+1) + "Y" + std::to_string(ind2+1) + ">")
                                             + Poly(J/4.0, "<Z" + std::to_string(ind1+1) + "Z" + std::to_string(ind2+1) + ">");
                hamiltonianInter[ind2][ind1] = hamiltonianInter[ind1][ind2];
            }
            for (int i=0; i<numQubits; i++) {
                int ind1 = i;
                int ind2 = (i+2) % numQubits;
                hamiltonianInter[ind1][ind2] += Poly(J/8.0, "<X" + std::to_string(ind1+1) + "X" + std::to_string(ind2+1) + ">")
                                             + Poly(J/8.0, "<Y" + std::to_string(ind1+1) + "Y" + std::to_string(ind2+1) + ">")
                                             + Poly(J/8.0, "<Z" + std::to_string(ind1+1) + "Z" + std::to_string(ind2+1) + ">");
                hamiltonianInter[ind2][ind1] = hamiltonianInter[ind1][ind2];
            }

            // Full Hamiltonian
            Poly H = Poly();
            for (int i=0; i<numQubits; i++) {
                for (int j=i; j<numQubits; j++) {
                    H += hamiltonianInter[i][j];
                }
            }
            H.reduce();

            // This is an energy only model
            objective = H;
            useEnergy = true;
            groundStateProblem = true;

            // We know the true ground state energy
            idealIsKnown = true;
            knownIdeal = -(3.0/8.0) * J * numQubits;

            // No baths
            lindbladianHot = Poly();
            lindbladianCold = Poly();

            // Lindbladian = -i[H, rho] + \sum (sigma_minus*rho*sigma_plus - 0.5*sigma_plus*sigma_minus*rho - 0.5*rho*sigma_plus*sigma_minus)
            Poly rho("<R1>");
            lindbladian = -imag*H.commutator(rho);
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian += lindbladianHot;
            lindbladian += lindbladianCold;
            lindbladian.convertToPaulis();
            lindbladian.reduce();

        // 2D Lindbladian test
        } else if (argAsString == "--2d" || argAsString == "--2dtfiperiodic" || argAsString == "--2dtfi" || argAsString == "--2dtwo") {
            modelName = argAsString;

            // Defining quantities
            double gamma_c = 1.1e-2;
            double gamma_h = 1e-3;
            double g = 1.6e-3;
            double T_h = 1.0;
            double T_c = 0.1;
            double delta = 0.005;
            double epsilon_h = 1.0;

            // Regardless we should be given the number of qubits
            gridWidth = std::stoi(argv[i+1]);
            gridHeight = std::stoi(argv[i+2]);
            numQubits = gridWidth * gridHeight;
            i += 2;

            // Calculated quantities
            double epsilon_c = epsilon_h + delta;
            double epsilon_other = (epsilon_h + epsilon_c) / 2.0;
            if (argAsString == "--2dtfi" || argAsString == "--2dtfiperiodic") {
                epsilon_other = epsilon_h;
                epsilon_c = epsilon_h;
            }
            std::vector<double> epsilons(numQubits, epsilon_other);
            double n_c = 1.0 / (std::exp(epsilon_c / T_c) - 1.0);
            double n_h = 1.0 / (std::exp(epsilon_h / T_h) - 1.0);
            double gamma_h_plus = gamma_h * n_h;
            double gamma_h_minus = gamma_h * (n_h + 1.0);
            double gamma_c_plus = gamma_c * n_c;
            double gamma_c_minus = gamma_c * (n_c + 1.0);

            // For the tfi, higher g TODO
            if (argAsString == "--2dtfi" || argAsString == "--2dtfiperiodic") {
                g = 0.5;
            }

            // Connectivity: 2D grid - nearest neighbours
            std::vector<double> gamma_plus(numQubits, 0);
            std::vector<double> gamma_minus(numQubits, 0);
            std::vector<std::vector<double>> gs(numQubits, std::vector<double>(numQubits, 0));
            for (int i=0; i<numQubits; i++) {

                // Get the x and y location
                int xLoc = i % gridWidth;
                int yLoc = i / gridWidth;

                // The leftmost qubits are connected to the hot bath
                if (xLoc == 0) {
                    gamma_plus[i] = gamma_h_plus;
                    gamma_minus[i] = gamma_h_minus;
                    epsilons[i] = epsilon_h;

                // The rightmost qubits are connected to the cold bath
                } else if (xLoc == gridWidth-1) {
                    gamma_plus[i] = gamma_c_plus;
                    gamma_minus[i] = gamma_c_minus;
                    epsilons[i] = epsilon_c;
                }

                // The qubit to the right
                if (xLoc < gridWidth-1) {
                    int otherInd = xLoc+1 + yLoc*gridWidth;
                    gs[i][otherInd] = g;
                    gs[otherInd][i] = g;
                }

                // The qubit to the left
                if (xLoc > 0) {
                    int otherInd = xLoc-1 + yLoc*gridWidth;
                    gs[i][otherInd] = g;
                    gs[otherInd][i] = g;
                }

                // The qubit below
                if (yLoc < gridHeight-1) {
                    int otherInd = xLoc + (yLoc+1)*gridWidth;
                    gs[i][otherInd] = g;
                    gs[otherInd][i] = g;
                }

                // The qubit above
                if (yLoc > 0) {
                    int otherInd = xLoc + (yLoc-1)*gridWidth;
                    gs[i][otherInd] = g;
                    gs[otherInd][i] = g;
                }

                // If y-periodic
                if (argAsString == "--2dtfiperiodic") {
                    if (yLoc == 0) {
                        int otherInd = xLoc + (gridHeight-1)*gridWidth;
                        gs[i][otherInd] = g;
                        gs[otherInd][i] = g;
                    } else if (yLoc == gridHeight-1) {
                        int otherInd = xLoc;
                        gs[i][otherInd] = g;
                        gs[otherInd][i] = g;
                    }
                }
                
            }

            // The Hamiltonian
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            if (argAsString == "--2dtwo") {
                for (int i=0; i<numQubits; i++) {
                    hamiltonianInter[i][i] = Poly(epsilons[i], "<P" + std::to_string(i+1) + "M" + std::to_string(i+1) + ">");
                }
                for (int i=0; i<numQubits; i++) {
                    for (int j=i+1; j<numQubits; j++) {
                        if (gs[i][j] != 0) {
                            hamiltonianInter[i][j] = Poly(gs[i][j], "<P" + std::to_string(i+1) + "M" + std::to_string(j+1) + ">") 
                                                   + Poly(gs[i][j], "<M" + std::to_string(i+1) + "P" + std::to_string(j+1) + ">");
                            hamiltonianInter[j][i] = hamiltonianInter[i][j];
                        }
                    }
                }
            } else if (argAsString == "--2dtfi" || argAsString == "--2dtfiperiodic") {
                for (int i=1; i<=numQubits; i++) {
                    hamiltonianInter[i-1][i-1] = Poly(epsilons[i-1], "<Z" + std::to_string(i) + ">");
                }
                for (int i=1; i<=numQubits; i++) {
                    for (int j=i+1; j<=numQubits; j++) {
                        hamiltonianInter[i-1][j-1] = Poly(gs[i-1][j-1], "<X" + std::to_string(i) + "X" + std::to_string(j) + ">");
                        hamiltonianInter[j-1][i-1] = hamiltonianInter[i-1][j-1];
                    }
                }
            } else {
                for (int i=1; i<=numQubits; i++) {
                    hamiltonianInter[i-1][i-1] = Poly(epsilons[i-1], "<Z" + std::to_string(i) + ">");
                }
                for (int i=1; i<=numQubits; i++) {
                    for (int j=i+1; j<=numQubits; j++) {
                        hamiltonianInter[i-1][j-1] = Poly(gs[i-1][j-1], "<X" + std::to_string(i) + "X" + std::to_string(j) + ">");
                        hamiltonianInter[j-1][i-1] = hamiltonianInter[i-1][j-1];
                    }
                }
            }
            for (int i=0; i<numQubits; i++) {
                for (int j=0; j<numQubits; j++) {
                    hamiltonianInter[i][j].convertToPaulis();
                }
            }
            Poly H = Poly();
            for (int i=0; i<numQubits; i++) {
                for (int j=i; j<numQubits; j++) {
                    H += hamiltonianInter[i][j];
                }
            }

            // The hot and cold parts of the Lindbladian
            lindbladianHot = Poly();
            lindbladianCold = Poly();
            for (int i=1; i<=numQubits; i++) {
                int xLoc = (i-1) % gridWidth;
                if (xLoc == 0) {
                    lindbladianHot += Poly(gamma_h_plus, "<M" + std::to_string(i) + "A0P" + std::to_string(i) + ">");
                    lindbladianHot += Poly(-0.5*gamma_h_plus, "<A0M" + std::to_string(i) + "P" + std::to_string(i) + ">");
                    lindbladianHot += Poly(-0.5*gamma_h_plus, "<M" + std::to_string(i) + "P" + std::to_string(i) + "A0>");
                    lindbladianHot += Poly(gamma_h_minus, "<P" + std::to_string(i) + "A0M" + std::to_string(i) + ">");
                    lindbladianHot += Poly(-0.5*gamma_h_minus, "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                    lindbladianHot += Poly(-0.5*gamma_h_minus, "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
                }
                if (xLoc == gridWidth-1) {
                    lindbladianCold += Poly(gamma_c_plus, "<M" + std::to_string(i) + "A0P" + std::to_string(i) + ">");
                    lindbladianCold += Poly(-0.5*gamma_c_plus, "<A0M" + std::to_string(i) + "P" + std::to_string(i) + ">");
                    lindbladianCold += Poly(-0.5*gamma_c_plus, "<M" + std::to_string(i) + "P" + std::to_string(i) + "A0>");
                    lindbladianCold += Poly(gamma_c_minus, "<P" + std::to_string(i) + "A0M" + std::to_string(i) + ">");
                    lindbladianCold += Poly(-0.5*gamma_c_minus, "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                    lindbladianCold += Poly(-0.5*gamma_c_minus, "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
                }
            }
            lindbladianHot.convertToPaulis();
            lindbladianCold.convertToPaulis();

            // Construct the objective as a polynomial
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }

            // Lindbladian = -i[H, rho] + \sum (sigma_minus*rho*sigma_plus - 0.5*sigma_plus*sigma_minus*rho - 0.5*rho*sigma_plus*sigma_minus)
            Poly rho("<R1>");
            lindbladian = -imag*H.commutator(rho);
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian += lindbladianHot;
            lindbladian += lindbladianCold;
            lindbladian.convertToPaulis();
            lindbladian.reduce();

            // Add symmetries
            if (allSymmetries) {

                // For the periodic case
                if (argAsString == "--2dtfiperiodic") {

                    // All rows should be equal
                    int rowInd = 0;
                    std::vector<int> rowInds;
                    for (int j=0; j<gridWidth; j++) {
                        rowInds.push_back(rowInd*gridWidth + j + 1);
                    }
                    for (int i=0; i<gridHeight; i++) {
                        if (i == rowInd) {
                            continue;
                        }
                        std::vector<int> rowIndsOther;
                        for (int j=0; j<gridWidth; j++) {
                            rowIndsOther.push_back(i*gridWidth + j + 1);
                        }
                        symmetries.push_back({rowInds, rowIndsOther});
                        if (verbosity >= 1) {
                            std::cout << "Adding symmetry: " << rowInds << " <-> " << rowIndsOther << std::endl;
                        }
                    }

                // For all other cases
                } else {

                    // Top and bottom rows should be equal
                    int rowInd1 = 0;
                    int rowInd2 = gridHeight - 1;
                    while (rowInd1 < rowInd2) {
                        std::vector<int> topInds;
                        for (int i=0; i<gridWidth; i++) {
                            topInds.push_back(rowInd1*gridWidth + i + 1);
                        }
                        std::vector<int> bottomInds;
                        for (int i=0; i<gridWidth; i++) {
                            bottomInds.push_back(rowInd2*gridWidth + i + 1);
                        }
                        symmetries.push_back({topInds, bottomInds});
                        if (verbosity >= 1) {
                            std::cout << "Adding symmetry: " << topInds << " <-> " << bottomInds << std::endl;
                        }
                        rowInd1++;
                        rowInd2--;
                    }

                }

            }

        // Simple phase transition example
        } else if (argAsString == "--phase1") {
            modelName = argAsString;

            double gamma = std::stod(argv[i+1]);
            i++;
            numQubits = 1;

            // Construct the objective as a polynomial
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }
            
            // The full Lindbladian
            // (1-gamma)*(sigma_minus*rho*sigma_plus - 0.5*sigma_plus*sigma_minus*rho - 0.5*rho*sigma_plus*sigma_minus)
            // + gamma*(sigma_z*rho*sigma_z - rho)
            Poly rho("<R1>");
            lindbladian = (1-gamma)*(Poly("<M1>")*rho*Poly("<P1>") - 0.5*Poly("<P1>")*Poly("<M1>")*rho - 0.5*rho*Poly("<P1>")*Poly("<M1>"))
                        + gamma*(Poly("<Z1>")*rho*Poly("<Z1>") - rho);
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.convertToPaulis();
            lindbladian.reduce();

        // TODO
        // Jin "Cluster Mean-Field Approach to the Steady-State Phase Diagram of Dissipative Spin Systems"
        // Lee "Unconventional magnetism via optical pumping of interacting spin systems"
        } else if (argAsString == "--phase2") {
            modelName = argAsString;

            // Parameters
            int latticeWidth = std::stoi(argv[i+1]);
            int latticeHeight = std::stoi(argv[i+2]);
            double Jx = std::stod(argv[i+3]);
            double Jy = std::stod(argv[i+4]);
            double gamma = 1.0;
            double Jz = 1.0;
            numQubits = latticeWidth * latticeHeight;
            i += 4;

            // Objective is random for now
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }

            // H = sum Jx X_i X_j + Jy Y_i Y_j + Jz Z_i Z_j
            // 2D grid
            Poly H;
            for (int x=0; x<latticeWidth; x++) {
                for (int y=0; y<latticeHeight; y++) {

                    // The index of this spin
                    int i = x + y*latticeWidth + 1;
                    int j = 0;

                    // The index of the spin x++
                    if (x < latticeWidth-1) {
                        j = (x+1) + y*latticeWidth + 1;
                        H += Poly(Jx, "<X" + std::to_string(i) + "X" + std::to_string(j) + ">");
                        H += Poly(Jy, "<Y" + std::to_string(i) + "Y" + std::to_string(j) + ">");
                        H += Poly(Jz, "<Z" + std::to_string(i) + "Z" + std::to_string(j) + ">");
                    }

                    // The index of the spin y++
                    if (y < latticeHeight-1) {
                        j = x + (y+1)*latticeWidth + 1;
                        H += Poly(Jx, "<X" + std::to_string(i) + "X" + std::to_string(j) + ">");
                        H += Poly(Jy, "<Y" + std::to_string(i) + "Y" + std::to_string(j) + ">");
                        H += Poly(Jz, "<Z" + std::to_string(i) + "Z" + std::to_string(j) + ">");
                    }

                    // The index of the spin x--
                    if (x > 0) {
                        j = (x-1) + y*latticeWidth + 1;
                        H += Poly(Jx, "<X" + std::to_string(i) + "X" + std::to_string(j) + ">");
                        H += Poly(Jy, "<Y" + std::to_string(i) + "Y" + std::to_string(j) + ">");
                        H += Poly(Jz, "<Z" + std::to_string(i) + "Z" + std::to_string(j) + ">");
                    }

                    // The index of the spin y--
                    if (y > 0) {
                        j = x + (y-1)*latticeWidth + 1;
                        H += Poly(Jx, "<X" + std::to_string(i) + "X" + std::to_string(j) + ">");
                        H += Poly(Jy, "<Y" + std::to_string(i) + "Y" + std::to_string(j) + ">");
                        H += Poly(Jz, "<Z" + std::to_string(i) + "Z" + std::to_string(j) + ">");
                    }

                }
            }

            // Lindbladian = -i[H, rho] + gamma*(sum sigma_minus*rho*sigma_plus - 0.5*sigma_plus*sigma_minus*rho - 0.5*rho*sigma_plus*sigma_minus)
            Poly rho("<R1>");
            lindbladian = -imag*H.commutator(rho);
            for (int i=1; i<=numQubits; i++) {
                lindbladian += gamma*(Poly("<M" + std::to_string(i) + ">")*rho*Poly("<P" + std::to_string(i) + ">") - 0.5*Poly("<P" + std::to_string(i) + ">")*Poly("<M" + std::to_string(i) + ">")*rho - 0.5*rho*Poly("<P" + std::to_string(i) + ">")*Poly("<M" + std::to_string(i) + ">"));
            }
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.convertToPaulis();
            lindbladian.reduce();

        // TODO https://journals.aps.org/pra/pdf/10.1103/PhysRevA.98.042118?casa_token=CFcuEJLHI8IAAAAA%3AQA26OCFU-TCiLQGncKG3TRRu9FW-Bf5mfgdspWRoCLiLa-1UJN9qLSBlSy9qqyfGQPAvs4e54FXAwcwp
        } else if (argAsString == "--phase3") {
            modelName = argAsString;

        // For an arbitrary connectivity, given as a file
        // This file should have the number of qubits on the first line
        // Then each line should have two integers detailing connectivity (e.g. "3 6")
        // For baths, use -1 for the cold and -2 for the hot
        } else if (argAsString == "--file") {
            modelName = argAsString;

            // Get the filename
            std::string filename = argv[i+1];
            i++;
            if (filename == "") {
                std::cerr << "No filename given for --file" << std::endl;
                return 1;
            }

            // Load the file
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Could not open file " << filename << std::endl;
                return 1;
            }

            // Load the connectivity
            std::set<std::pair<int, int>> connectivity;
            std::set<int> connectedToHot;
            std::set<int> connectedToCold;
            int q1, q2;
            while (file >> q1 >> q2) {
                std::pair<int, int> conn;
                if (q1 < q2) {
                    conn = std::make_pair(q1, q2);
                } else {
                    conn = std::make_pair(q2, q1);
                }
                connectivity.insert(conn);
                if (q1 == -1) {
                    connectedToCold.insert(q2);
                } else if (q2 == -1) {
                    connectedToCold.insert(q1);
                } else if (q1 == -2) {
                    connectedToHot.insert(q2);
                } else if (q2 == -2) {
                    connectedToHot.insert(q1);
                }
            }

            // Find the minimum and maximum qubit indices
            int minQubit = 1000000;
            int maxQubit = -1;
            for (auto conn : connectivity) {
                minQubit = std::min(minQubit, conn.first);
                minQubit = std::min(minQubit, conn.second);
                maxQubit = std::max(maxQubit, conn.first);
                maxQubit = std::max(maxQubit, conn.second);
            }

            // Rescale the qubits to start at 0
            std::set<std::pair<int, int>> newConnectivity;
            for (auto &conn : connectivity) {
                newConnectivity.insert(std::make_pair(conn.first - minQubit, conn.second - minQubit));
            }
            connectivity = newConnectivity;
            std::set<int> newConnectedToHot;
            for (auto &q : connectedToHot) {
                newConnectedToHot.insert(q - minQubit);
            }
            connectedToHot = newConnectedToHot;
            std::set<int> newConnectedToCold;
            for (auto &q : connectedToCold) {
                newConnectedToCold.insert(q - minQubit);
            }
            connectedToCold = newConnectedToCold;
            numQubits = maxQubit - minQubit + 1;

            // Defining quantities
            double gamma_c = 1e-3;
            double gamma_h = 1e-3;
            double g = 0.1;
            double T_h = 1.0;
            double T_c = 0.1;
            double delta = 0.005;
            double epsilon_h = 1.0;

            // Construct the objective as a polynomial
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }
            
            // The Hamiltonian
            // H = \sum_nn <X{i}X{j}> + g * \sum_i <Zi>
            Poly H;
            for (int i=1; i<=numQubits; i++) {
                H += Poly(g, "<Z" + std::to_string(i) + ">");
            }
            for (auto conn : connectivity) {
                H += Poly(g, "<X" + std::to_string(conn.first) + "X" + std::to_string(conn.second) + ">");
            }

            // Output the Hamiltonian
            if (verbosity >= 2) {
                std::cout << "Hamiltonian: " << H << std::endl;
            }

            // The jump operators
            std::vector<Poly> Gamma_k(numQubits);
            for (int i=0; i<numQubits; i++) {
                Gamma_k[i] = Poly("<X" + std::to_string(i+1) + ">") - imag*Poly("<Y" + std::to_string(i+1) + ">");
                Gamma_k[i] /= 2.0;
                if (connectedToHot.find(i) != connectedToHot.end()) {
                    Gamma_k[i] *= std::sqrt(gamma_h);
                } else if (connectedToCold.find(i) != connectedToCold.end()) {
                    Gamma_k[i] *= std::sqrt(gamma_c);
                } else {
                    Gamma_k[i] *= 0.0;
                }
            }

            // The full Lindbladian
            // -i[H, rho] + \sum_k 0.5 * (2*gamma_h * Gamma_k rho Gamma_k^dagger - Gamma_k^dagger Gamma_k rho - rho Gamma_k^dagger Gamma_k)
            Poly rho("<R1>");
            lindbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                lindbladian += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
            }
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.reduce();

            // Note the local and non-local Hamiltonians in case we need them later
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            for (int i=1; i<=numQubits; i++) {
                hamiltonianInter[i-1][i-1] = Poly(g, "<Z" + std::to_string(i) + ">");
            }
            for (auto conn : connectivity) {
                hamiltonianInter[conn.first][conn.second] = Poly(1, "<X" + std::to_string(conn.first+1) + "X" + std::to_string(conn.second+1) + ">");
                hamiltonianInter[conn.second][conn.first] = Poly(1, "<X" + std::to_string(conn.first+1) + "X" + std::to_string(conn.second+1) + ">");
            }
            lindbladianHot = Poly();
            lindbladianCold = Poly();
            for (int i=0; i<numQubits; i++) {
                if (connectedToHot.find(i) != connectedToHot.end()) {
                    lindbladianHot += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
                if (connectedToCold.find(i) != connectedToCold.end()) {
                    lindbladianCold += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            lindbladianHot = Poly("<A0>") * lindbladianHot;
            lindbladianCold = Poly("<A0>") * lindbladianCold;
            lindbladianHot.cycleToAndRemove('R', 1);
            lindbladianCold.cycleToAndRemove('R', 1);
            lindbladianHot.reduce();
            lindbladianCold.reduce();

        // The Lindbladian from David
        } else if (argAsString == "--david" || argAsString == "--davidr") {
            modelName = argAsString;

            // Parameters
            double g = std::stod(argv[i+1]);
            numQubits = std::stoi(argv[i+2]);
            i+=2;
            double gamma_h = 1.0;
            double gamma_c = 0.5;

            // Construct the objective as a polynomial
            objective = Poly("<Z1>");
            
            // The Hamiltonian
            // H = \sum_i <X{i}X{i+1}> + g * \sum_i <Zi>
            Poly H;
            for (int i=1; i<=numQubits; i++) {
                H += Poly(g, "<Z" + std::to_string(i) + ">");
            }
            for (int i=1; i<numQubits; i++) {
                H += Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(i+1) + ">");
            }
            if (argAsString == "--davidr") {
                H.randomize();
            }

            // The jump operators
            std::vector<Poly> Gamma_k(numQubits);
            for (int i=1; i<=numQubits; i++) {
                Gamma_k[i-1] = Poly("<X" + std::to_string(i) + ">") - imag*Poly("<Y" + std::to_string(i) + ">");
                Gamma_k[i-1] /= 2.0;
            }
            Gamma_k[0] *= std::sqrt(gamma_h);
            Gamma_k[numQubits-1] *= std::sqrt(gamma_c);

            // The full Lindbladian
            // -i[H, rho] + \sum_k 0.5 * (2*gamma_h * Gamma_k rho Gamma_k^dagger - Gamma_k^dagger Gamma_k rho - rho Gamma_k^dagger Gamma_k)
            Poly rho("<R1>");
            lindbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                if (i == 0 || i == numQubits-1) {
                    lindbladian += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.reduce();

            // Note the local and non-local Hamiltonians in case we need them later
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            for (int i=0; i<numQubits; i++) {
                hamiltonianInter[i][i] = Poly(g, "<Z" + std::to_string(i+1) + ">");
            }
            for (int i=1; i<numQubits; i++) {
                hamiltonianInter[i-1][i] = Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(i+1) + ">");
                hamiltonianInter[i][i-1] = Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(i+1) + ">");
            }
            lindbladianHot = Poly();
            lindbladianCold = Poly();
            for (int i=0; i<numQubits; i++) {
                if (i == 0) {
                    lindbladianHot += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
                if (i == numQubits-1) {
                    lindbladianCold += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            lindbladianHot = Poly("<A0>") * lindbladianHot;
            lindbladianCold = Poly("<A0>") * lindbladianCold;
            lindbladianHot.cycleToAndRemove('R', 1);
            lindbladianCold.cycleToAndRemove('R', 1);
            lindbladianHot.reduce();
            lindbladianCold.reduce();

        // The Lindbladian from David but 2D
        } else if (argAsString == "--david2d" || argAsString == "--david2dr") {
            modelName = argAsString;

            // Parameters
            double g = std::stod(argv[i+1]);
            gridWidth = std::stoi(argv[i+2]);
            gridHeight = std::stoi(argv[i+3]);
            i+=3;
            numQubits = gridWidth * gridHeight;
            double gamma_h = 1.0;
            double gamma_c = 0.5;

            // Construct the objective as a polynomial
            objective = Poly("<Z1>");
            
            // The Hamiltonian
            // H = \sum_nn <X{i}X{j}> + g * \sum_i <Zi>
            Poly H;
            for (int i=1; i<=numQubits; i++) {
                H += Poly(g, "<Z" + std::to_string(i) + ">");
            }
            for (int i=1; i<=numQubits; i++) {

                // Get the x and y location
                int xLoc = (i-1) % gridWidth;
                int yLoc = (i-1) / gridWidth;

                // The qubit to the right
                if (xLoc < gridWidth-1) {
                    int otherInd = xLoc+1 + yLoc*gridWidth + 1;
                    H += Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(otherInd) + ">");
                }

                // The qubit below
                if (yLoc < gridHeight-1) {
                    int otherInd = xLoc + (yLoc+1)*gridWidth + 1;
                    H += Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(otherInd) + ">");
                }

            }

            // Randomize the Hamiltonian if asked to
            if (argAsString == "--david2dr") {
                H.randomize();
            }

            // Output the Hamiltonian
            if (verbosity >= 2) {
                std::cout << "Hamiltonian: " << H << std::endl;
            }

            // The jump operators
            std::vector<Poly> Gamma_k(numQubits);
            for (int i=1; i<=numQubits; i++) {
                Gamma_k[i-1] = Poly("<X" + std::to_string(i) + ">") - imag*Poly("<Y" + std::to_string(i) + ">");
                Gamma_k[i-1] /= 2.0;
                int xLoc = (i-1) % gridWidth;
                if (xLoc != 0 && xLoc != gridWidth-1) {
                    Gamma_k[i-1] *= 0;
                } else if (xLoc == 0) {
                    Gamma_k[i-1] *= std::sqrt(gamma_h);
                } else if (xLoc == gridWidth-1) {
                    Gamma_k[i-1] *= std::sqrt(gamma_c);
                }
            }

            // The full Lindbladian
            // -i[H, rho] + \sum_k 0.5 * (2*gamma_h * Gamma_k rho Gamma_k^dagger - Gamma_k^dagger Gamma_k rho - rho Gamma_k^dagger Gamma_k)
            Poly rho("<R1>");
            lindbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                int xLoc = i % gridWidth;
                int yLoc = i / gridWidth;
                if (xLoc == 0 || xLoc == gridWidth-1) {
                    if (verbosity >= 2) {
                        std::cout << "Connecting " << (i+1) << " to the bath" << std::endl;
                    }
                    lindbladian += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.reduce();

            // Note the local and non-local Hamiltonians in case we need them later
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            for (int i=0; i<numQubits; i++) {
                hamiltonianInter[i][i] = Poly(g, "<Z" + std::to_string(i+1) + ">");
            }
            for (int i=1; i<=numQubits; i++) {
                int xLoc = (i-1) % gridWidth;
                int yLoc = (i-1) / gridWidth;
                if (xLoc < gridWidth-1) {
                    int otherInd = xLoc+1 + yLoc*gridWidth + 1;
                    hamiltonianInter[i-1][otherInd-1] = Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(otherInd) + ">");
                    hamiltonianInter[otherInd-1][i-1] = Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(otherInd) + ">");
                }
                if (yLoc < gridHeight-1) {
                    int otherInd = xLoc + (yLoc+1)*gridWidth + 1;
                    hamiltonianInter[i-1][otherInd-1] = Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(otherInd) + ">");
                    hamiltonianInter[otherInd-1][i-1] = Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(otherInd) + ">");
                }
            }
            lindbladianHot = Poly();
            lindbladianCold = Poly();
            for (int i=0; i<numQubits; i++) {
                int xLoc = i % gridWidth;
                if (xLoc == 0) {
                    lindbladianHot += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
                if (xLoc == gridWidth-1) {
                    lindbladianCold += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            lindbladianHot = Poly("<A0>") * lindbladianHot;
            lindbladianCold = Poly("<A0>") * lindbladianCold;
            lindbladianHot.cycleToAndRemove('R', 1);
            lindbladianCold.cycleToAndRemove('R', 1);
            lindbladianHot.reduce();
            lindbladianCold.reduce();

        // A 1.5D chain
        } else if (argAsString == "--1.5d" || argAsString == "--1.5dr") {
            modelName = argAsString;

            // Parameters
            double g = std::stod(argv[i+1]);
            numQubits = std::stoi(argv[i+2]);
            i+=2;
            double gamma_h = 1.0;
            double gamma_c = 0.5;

            // Make sure it's odd
            if (numQubits % 2 == 0) {
                std::cerr << "Error - Number of qubits must be odd" << std::endl;
                return 1;
            }

            // Construct the objective as a polynomial
            objective = Poly("<Z1>");

            // The Hamiltonian
            // H = \sum_i <X{i}X{i+1}> + g * \sum_i <Zi>
            Poly H;
            for (int i=1; i<=numQubits; i++) {
                H += Poly(g, "<Z" + std::to_string(i) + ">");
            }
            for (int i=1; i<numQubits-1; i+=2) {
                H += Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(i+2) + ">");
            }
            for (int i=2; i<numQubits; i+=2) {
                H += Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(i-1) + ">");
                H += Poly(1, "<X" + std::to_string(i) + "X" + std::to_string(i+1) + ">");
            }
            if (argAsString == "--1.5dr") {
                H.randomize(-1, 1);
            }

            // Output the Hamiltonian
            if (verbosity >= 2) {
                std::cout << "Hamiltonian: " << H << std::endl;
            }

            // The jump operators
            std::vector<Poly> Gamma_k(numQubits);
            for (int i=1; i<=numQubits; i++) {
                Gamma_k[i-1] = Poly("<X" + std::to_string(i) + ">") - imag*Poly("<Y" + std::to_string(i) + ">");
                Gamma_k[i-1] /= 2.0;
            }
            Gamma_k[0] *= std::sqrt(gamma_h);
            Gamma_k[numQubits-1] *= std::sqrt(gamma_c);

            // The full Lindbladian
            // -i[H, rho] + \sum_k 0.5 * (2*gamma_h * Gamma_k rho Gamma_k^dagger - Gamma_k^dagger Gamma_k rho - rho Gamma_k^dagger Gamma_k)
            Poly rho("<R1>");
            lindbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                if (i == 0 || i == numQubits-1) {
                    lindbladian += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.reduce();

        // The Lindbladian from the tensor paper
        } else if (argAsString == "--tensor") {
            modelName = argAsString;

            // Defining quantities
            double J = 0.5;
            double gamma_minus = 1;
            double h = 0.5;

            // Regardless we should be given the number of qubits
            numQubits = std::stoi(argv[i+1]);
            i++;

            // The Hamiltonian
            // H = J * \sum_i <Z{i}Z{i+1}> + h * \sum_i <Xi>
            Poly H;
            for (int i=1; i<=numQubits; i++) {
                H += Poly(h, "<X" + std::to_string(i) + ">");
            }
            for (int i=1; i<numQubits; i++) {
                H += Poly(J, "<Z" + std::to_string(i) + "Z" + std::to_string(i+1) + ">");
            }
            H += Poly(J, "<Z" + std::to_string(numQubits) + "Z1>");

            // The jump operators
            // Gamma_k = sqrt(gamma_minus)/2 * (<X{i}> - i <Y{i}>)
            std::vector<Poly> Gamma_k(numQubits);
            for (int i=1; i<=numQubits; i++) {
                Gamma_k[i-1] = Poly("<X" + std::to_string(i) + ">") - imag*Poly("<Y" + std::to_string(i) + ">");
                Gamma_k[i-1] *= std::sqrt(gamma_minus)/2.0;
            }

            // The full Lindbladian
            // -i[H, rho] + \sum_k Gamma_k rho Gamma_k^dagger - 0.5 {Gamma_k^dagger Gamma_k, rho}
            Poly rho("<R1>");
            lindbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                lindbladian += Gamma_k[i] * rho * Gamma_k[i].dagger();
                lindbladian -= 0.5 * (Gamma_k[i].dagger() * Gamma_k[i]).anticommutator(rho);
            }
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.reduce();

            // Objective is just the average magnetization
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }
            
            // Print everything for sanity
            if (verbosity >= 2) {
                std::cout << "Hamiltonian: " << H << std::endl;
                std::cout << "Jump operators: " << std::endl;
                for (int i=0; i<numQubits; i++) {
                    std::cout << "  " << Gamma_k[i] << std::endl;
                }
                std::cout << "Lindbladian: " << lindbladian << std::endl;
                std::cout << "Objective: " << objective << std::endl;
            }

            // Note the local and non-local Hamiltonians in case we need them later
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            for (int i=0; i<numQubits; i++) {
                hamiltonianInter[i][i] = Poly(h, "<X" + std::to_string(i+1) + ">");
            }
            for (int i=0; i<numQubits-1; i++) {
                hamiltonianInter[i][i+1] = Poly(J, "<Z" + std::to_string(i+1) + "Z" + std::to_string(i+2) + ">");
                hamiltonianInter[i+1][i] = hamiltonianInter[i][i+1];
            }
            hamiltonianInter[numQubits-1][0] = Poly(J, "<Z" + std::to_string(numQubits) + "Z1>");
            hamiltonianInter[0][numQubits-1] = hamiltonianInter[numQubits-1][0];
            lindbladianHot = Poly();
            for (int i=0; i<numQubits; i++) {
                lindbladianHot += Gamma_k[i] * rho * Gamma_k[i].dagger();
                lindbladianHot -= 0.5 * (Gamma_k[i].dagger() * Gamma_k[i]).anticommutator(rho);
            }
            lindbladianHot = Poly("<A0>") * lindbladian;
            lindbladianHot.cycleToAndRemove('R', 1);
            lindbladianHot.reduce();

            // Note the symmetries
            if (allSymmetries) {
                for (int i=1; i<=numQubits; i++) {
                    for (int j=i+1; j<=numQubits; j++) {
                        symmetries.push_back({{i}, {j}});
                    }
                }
            }

        // Lindbladian from quasiperiodic systems paper
        } else if (argAsString == "--quasi") {
            modelName = argAsString;

            // System params
            numQubits = std::stoi(argv[i+1]);
            i++;
            double beta = (std::sqrt(5)-1)  / 2.0;
            double U = 0.1;
            double phi = 0;
            double lambda = 1.0;
            double mu = 0.01;
            std::vector<double> h(numQubits, 0); 
            for (int i=0; i<numQubits; i++) {
                h[i] = lambda*std::cos(2.0*M_PI*beta*i + phi);
            }

            // H = XX + YY + ZZ + h_i Z_i
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            for (int i=1; i<=numQubits; i++) {
                hamiltonianInter[i-1][i-1] = Poly(1.0, "<Z" + std::to_string(i) + ">");
            }
            for (int i=1; i<numQubits; i++) {
                hamiltonianInter[i-1][i] += Poly(1.0, "<X" + std::to_string(i) + "X" + std::to_string(i+1) + ">");
                hamiltonianInter[i-1][i] += Poly(1.0, "<Y" + std::to_string(i) + "Y" + std::to_string(i+1) + ">");
                hamiltonianInter[i-1][i] += Poly(1.0, "<Z" + std::to_string(i) + "Z" + std::to_string(i+1) + ">");
                hamiltonianInter[i][i-1] = hamiltonianInter[i-1][i];
            }

            // Objective is 2XY - 2YX
            objective = Poly(2.0, "<X1Y2>") - Poly(2.0, "<Y1X2>");

            // Jump operators
            std::vector<Poly> L_k(numQubits);
            L_k[0] = Poly(std::sqrt(1+mu), "<P1>");
            L_k[1] = Poly(std::sqrt(1-mu), "<M1>");
            L_k[numQubits-2] = Poly(std::sqrt(1-mu), "<P" + std::to_string(numQubits) + ">");
            L_k[numQubits-1] = Poly(std::sqrt(1+mu), "<M" + std::to_string(numQubits) + ">");

            // The full Lindbladian
            Poly rho("<R1>");
            lindbladian = Poly();
            for (int i=0; i<numQubits; i++) {
                if (L_k[i].size() > 0) {
                    lindbladian += (L_k[i] * rho).commutator(L_k[i].dagger());
                    lindbladian += (L_k[i]).commutator(rho * L_k[i].dagger());
                }
            }
            lindbladian.convertToPaulis();
            lindbladian = Poly("<A0>") * lindbladian;
            lindbladian.cycleToAndRemove('R', 1);
            lindbladian.reduce();

        // Many-body Lindbladian
        } else if (argAsString == "--many" || argAsString == "--manyv" || argAsString == "--manyr") {
            modelName = argAsString;

            // Defining quantities
            double gamma_c = 1.1e-2;
            double gamma_h = 1e-3;
            double g = 1.6e-3;
            double T_h = 1.0;
            double T_c = 0.1;
            double delta = 0.005;
            double epsilon_h = 1.0;

            // Regardless we should be given the number of qubits
            numQubits = std::stoi(argv[i+1]);
            i++;

            // If the arg has a "v", we have values
            if (argAsString == "--manyv") {
                T_h = std::stod(argv[i+1]);
                delta = std::stod(argv[i+2]);
                i += 2;
            }

            // Calculated quantities
            double epsilon_c = epsilon_h + delta;
            double epsilon_other = (epsilon_h + epsilon_c) / 2.0;
            std::vector<double> epsilons(numQubits, epsilon_other);
            epsilons[0] = epsilon_h;
            epsilons[numQubits-1] = epsilon_c;
            double n_c = 1.0 / (std::exp(epsilon_c / T_c) - 1.0);
            double n_h = 1.0 / (std::exp(epsilon_h / T_h) - 1.0);
            double gamma_h_plus = gamma_h * n_h;
            double gamma_h_minus = gamma_h * (n_h + 1.0);
            double gamma_c_plus = gamma_c * n_c;
            double gamma_c_minus = gamma_c * (n_c + 1.0);
            std::vector<double> gamma_plus(numQubits, 0);
            gamma_plus[0] = gamma_h_plus;
            gamma_plus[numQubits-1] = gamma_c_plus;
            std::vector<double> gamma_minus(numQubits, 0);
            gamma_minus[0] = gamma_h_minus;
            gamma_minus[numQubits-1] = gamma_c_minus;
            std::vector<double> gs(numQubits, g);
            
            // If told to randomize
            if (argAsString == "--manyr") {
                for (int j=0; j<numQubits; j++) {
                    epsilons[j] = rand(epsilon_h, epsilon_c);
                }
                for (int j=0; j<numQubits; j++) {
                    gs[j] = rand(1e-5, 1e-2);
                }
            }

            // Construct the objective as a polynomial with plus/minus mats
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }
            
            // Construct the Limbadlian as a polynomial from plus/minus
            lindbladian = Poly();
            for (int i=1; i<=numQubits; i++) {
                lindbladian += Poly(-imag*epsilons[i-1], "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                lindbladian += Poly(imag*epsilons[i-1], "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
            }
            for (int i=1; i<=numQubits-1; i++) {
                lindbladian += Poly(-imag*gs[i], "<A0P" + std::to_string(i) + "M" + std::to_string(i+1) + ">");
                lindbladian += Poly(-imag*gs[i], "<A0M" + std::to_string(i) + "P" + std::to_string(i+1) + ">");
                lindbladian += Poly(imag*gs[i], "<P" + std::to_string(i) + "M" + std::to_string(i+1) + "A0>");
                lindbladian += Poly(imag*gs[i], "<M" + std::to_string(i) + "P" + std::to_string(i+1) + "A0>");
            }
            for (int i : {1, numQubits}) {
                lindbladian += Poly(gamma_plus[i-1], "<M" + std::to_string(i) + "A0P" + std::to_string(i) + ">");
                lindbladian += Poly(-0.5*gamma_plus[i-1], "<A0M" + std::to_string(i) + "P" + std::to_string(i) + ">");
                lindbladian += Poly(-0.5*gamma_plus[i-1], "<M" + std::to_string(i) + "P" + std::to_string(i) + "A0>");
                lindbladian += Poly(gamma_minus[i-1], "<P" + std::to_string(i) + "A0M" + std::to_string(i) + ">");
                lindbladian += Poly(-0.5*gamma_minus[i-1], "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                lindbladian += Poly(-0.5*gamma_minus[i-1], "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
            }
            lindbladian.convertToPaulis();

            // Note the local and non-local Hamiltonians in case we need them later
            hamiltonianInter = std::vector<std::vector<Poly>>(numQubits, std::vector<Poly>(numQubits, Poly()));
            for (int j=0; j<numQubits; j++) {
                hamiltonianInter[j][j] = Poly(epsilons[j], "<P" + std::to_string(j+1) + "M" + std::to_string(j+1) + ">");
            }
            for (int j=0; j<numQubits-1; j++) {
                hamiltonianInter[j][j+1] = Poly(gs[j], "<P" + std::to_string(j+1) + "M" + std::to_string(j+2) + ">") + Poly(gs[j], "<M" + std::to_string(j+1) + "P" + std::to_string(j+2) + ">");
                hamiltonianInter[j+1][j] = hamiltonianInter[j][j+1];
            }
            lindbladianHot = Poly();
            for (int i : {1}) {
                lindbladianHot += Poly(gamma_plus[i-1], "<M" + std::to_string(i) + "A0P" + std::to_string(i) + ">");
                lindbladianHot += Poly(-0.5*gamma_plus[i-1], "<A0M" + std::to_string(i) + "P" + std::to_string(i) + ">");
                lindbladianHot += Poly(-0.5*gamma_plus[i-1], "<M" + std::to_string(i) + "P" + std::to_string(i) + "A0>");
                lindbladianHot += Poly(gamma_minus[i-1], "<P" + std::to_string(i) + "A0M" + std::to_string(i) + ">");
                lindbladianHot += Poly(-0.5*gamma_minus[i-1], "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                lindbladianHot += Poly(-0.5*gamma_minus[i-1], "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
            }
            lindbladianHot.convertToPaulis();
            lindbladianCold = Poly();
            for (int i : {numQubits}) {
                lindbladianCold += Poly(gamma_plus[i-1], "<M" + std::to_string(i) + "A0P" + std::to_string(i) + ">");
                lindbladianCold += Poly(-0.5*gamma_plus[i-1], "<A0M" + std::to_string(i) + "P" + std::to_string(i) + ">");
                lindbladianCold += Poly(-0.5*gamma_plus[i-1], "<M" + std::to_string(i) + "P" + std::to_string(i) + "A0>");
                lindbladianCold += Poly(gamma_minus[i-1], "<P" + std::to_string(i) + "A0M" + std::to_string(i) + ">");
                lindbladianCold += Poly(-0.5*gamma_minus[i-1], "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                lindbladianCold += Poly(-0.5*gamma_minus[i-1], "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
            }
            lindbladianCold.convertToPaulis();

        // Set the seed
        } else if (argAsString == "-S") {
            seed = std::string(argv[i+1]);
            srand(std::hash<std::string>{}(seed));
            std::seed_seq seedSeq (seed.begin(), seed.end());
            gen.seed(seedSeq);
            i++;

        // If setting verbosity
        } else if (argAsString == "-v") {
            verbosity = std::stoi(argv[i+1]);
            i++;

        // If setting the level of the Lindbladian
        } else if (argAsString == "-l") {
            lindbladLevel = std::stoi(argv[i+1]);
            i++;

        // If adding an extra monomial to the top row
        } else if (argAsString == "-e") {
            extraMonomials.push_back(std::string(argv[i+1]));
            i++;

        // If adding an extra monomial to the list of Lindbladian replacements
        } else if (argAsString == "-E") {
            extraMonomialsLim.push_back(std::string(argv[i+1]));
            i++;

        // If trying to find the minimal set of linear constraints
        } else if (argAsString == "-M") {
            findMinimal = true;
            findMinimalAmount = std::stoi(argv[i+1]);
            i++;

        // Which solver to use
        } else if (argAsString == "-s") {
            if (argv[i+1][0] == 'G') {
                solver = "gurobi";
            } else if (argv[i+1][0] == 'E') {
                solver = "eigen";
            } else if (argv[i+1][0] == 'S') {
                solver = "scs";
            } else if (argv[i+1][0] == 'M') {
                solver = "mosek";
            } else if (argv[i+1][0] == 'N') {
                solver = "none";
            }
            i++;

        // If told to try removing constraints
        } else if (argAsString == "-R") {
            tryRemove = true;

        // Add some number of extra constraints
        } else if (argAsString == "-c") {
            reductiveCons = std::stoi(argv[i+1]);
            i++;

        // Setting the number of cores
        } else if (argAsString == "-C") {
            numCores = std::stoi(argv[i+1]);
            i++;

        // If adding a symmetry
        } else if (argAsString == "-y") {
            std::vector<int> group1;
            std::vector<int> group2;
            std::string groupAsString1 = std::string(argv[i+1]);
            std::string groupAsString2 = std::string(argv[i+2]);
            i+=2;
            std::string toConvert = "";
            for (size_t j=0; j<groupAsString1.size(); j++) {
                if (groupAsString1[j] != ',') {
                    toConvert += groupAsString1[j];
                }
                if (groupAsString1[j] == ',' || j == groupAsString1.size()-1) {
                    group1.push_back(std::stoi(toConvert));
                    toConvert = "";
                }
            }
            for (size_t j=0; j<groupAsString2.size(); j++) {
                if (groupAsString2[j] != ',') {
                    toConvert += groupAsString2[j];
                }
                if (groupAsString2[j] == ',' || j == groupAsString2.size()-1) {
                    group2.push_back(std::stoi(toConvert));
                    toConvert = "";
                }
            }
            symmetries.push_back({group1, group2});

        // Heat current objective
        // --objHC H C   is    hot -> cold
        // --objHC 2 C   is    spin 2 -> cold
        // --objHC 2 5   is    spin 2 -> spin 5
        // --conHC P 2 C   is    spin 2 -> cold should be positive
        // --conHC N 2 C   is    spin 2 -> cold should be negative
        } else if (argAsString == "--objHC" || argAsString == "--conHC") {
            modelName = argAsString;
            usingHeatCurrent = true;

            // Determine which spins we're dealing with
            std::string type = "obj";
            std::vector<int> spinsToUse;
            while (i+1 < argc && argv[i+1][0] != '-') {
                if (argv[i+1][0] == 'P') {
                    type = "positive";
                } else if (argv[i+1][0] == 'N') {
                    type = "negative";
                } else if (argv[i+1][0] == 'H') {
                    spinsToUse.push_back(-1);
                    if (lindbladianHot.size() == 0) {
                        std::cerr << "Error - list of spins connected to the hot bath not defined" << std::endl;
                        return 1;
                    }
                } else if (argv[i+1][0] == 'C') {
                    spinsToUse.push_back(-2);
                    if (lindbladianCold.size() == 0) {
                        std::cerr << "Error - list of spins connected to the hot bath not defined" << std::endl;
                        return 1;
                    }
                } else {
                    spinsToUse.push_back(std::stoi(argv[i+1]));
                }
                i++;
            }

            // Make sure we know the Hamiltonian
            if (hamiltonianInter.size() == 0) {
                std::cerr << "Error - hamiltonian not explicitly defined" << std::endl;
                return 1;
            }

            // Make sure at least one spin was specified
            if (spinsToUse.size() == 0) {
                std::cerr << "Error - no spins specified" << std::endl;
                return 1;
            }

            // The full Hamiltonian
            Poly hamiltonianFull;
            for (int j=0; j<numQubits; j++) {
                for (int k=j; k<numQubits; k++) {
                    hamiltonianFull += hamiltonianInter[j][k];
                }
            }

            // If calculating the spin current between two spins
            Poly tempObj;
            if (spinsToUse.size() == 2 && spinsToUse[0] != -1 && spinsToUse[0] != -2) {
                int spin1 = spinsToUse[0]-1;
                int spin2 = spinsToUse[1]-1;
                //tempObj = imag * hamiltonianInter[spin1][spin1].commutator(hamiltonianInter[spin1][spin2]);
                tempObj = imag * hamiltonianFull.commutator(hamiltonianInter[spin1][spin2]);
                
            // Determine the current to/from the hot bath
            } else if (spinsToUse[0] == -1) {
                std::pair<char,int> oldMon('A', 0);
                Poly newPoly = hamiltonianFull;
                tempObj = lindbladianHot.replaced(oldMon, newPoly);

            // Determine the current to/from the cold bath
            } else if (spinsToUse[0] == -2) {
                std::pair<char,int> oldMon('A', 0);
                Poly newPoly = hamiltonianFull;
                tempObj = lindbladianCold.replaced(oldMon, newPoly);

            // Otherwise we use J_i = i[H, H_i]
            } else {
                int spin = spinsToUse[0]-1;
                tempObj = imag * hamiltonianFull.commutator(hamiltonianInter[spin][spin]);

            }

            // Make sure it's as Paulis
            if (verbosity >= 2) {
                std::cout << "overall obj = " << tempObj << std::endl;
            }
            tempObj.convertToPaulis();
            if (verbosity >= 2) {
                std::cout << "as Paulis = " << tempObj << std::endl;
            }

            // If it's an objective
            if (type == "obj") {
                objective = tempObj;
            } else if (type == "positive") {
                constraintsPositive.push_back(tempObj);
            } else if (type == "negative") {
                constraintsPositive.push_back(-tempObj);
            } else {
                std::cerr << "Error - unknown type of heat current constraint" << std::endl;
                return 1;
            }

        // Reconstruction constraints
        } else if (argAsString == "-r") {
            reconLevel = std::stoi(argv[i+1]);
            i+=1;

        // Tracing out constraints
        } else if (argAsString == "-T") {
            traceCons = true;

        // If taking a pseudo-objective
        } else if (argAsString == "-P") {
            pseudoObjective = objective;

        // If taking samples
        } else if (argAsString == "--samples") {
            numSamplesPer = std::stoi(argv[i+1]);
            i++;

        // If taking samples
        } else if (argAsString == "--shots") {
            numShots = std::stol(argv[i+1]);
            numSamplesPer = 1;
            i++;

        // If setting the percentile
        } else if (argAsString == "-p") {
            percentile = std::stod(argv[i+1]);
            i++;

        // If excluding the objective from the sampling
        } else if (argAsString == "--noobj") {
            excludeObjective = true;

        // If excluding any X terms from the sampling
        } else if (argAsString == "--nox") {
            excludeX = true;

        // If excluding any Y terms from the sampling
        } else if (argAsString == "--noy") {
            excludeY = true;

        // If excluding any Z terms from the sampling
        } else if (argAsString == "--noz") {
            excludeZ = true;

        // If sampling from all Pauli strings evenly
        } else if (argAsString == "--all") {
            sampleChoice = "all";
            maxSampleDegree = std::stoi(argv[i+1]);
            i++;

        // If sampling from only the objective
        } else if (argAsString == "--onlyobj") {
            sampleChoice = "onlyobj";

        // If sampling from a subset of Pauli strings
        } else if (argAsString == "--first") {
            sampleChoice = "auto";
            autoType = "first";
            maxPaulis = std::stoi(argv[i+1]);
            i++;

        // If sampling from a subset of Pauli strings
        } else if (argAsString == "--auto") {
            sampleChoice = "auto";
            autoType = "common";
            maxPaulis = std::stoi(argv[i+1]);
            i++;

        // If sampling from the energy
        } else if (argAsString == "--energy") {
            sampleChoice = "energy";

        // If outputting the Lindbladian in LaTeX form
        } else if (argAsString == "-X") {
            outputLindbladAsLaTeX = true;

        // If outputting the final moments to file
        } else if (argAsString == "-f") {
            outputToFile = true;

        // If checking the final results by reconstructing the full state
        } else if (argAsString == "-o") {
            checkObj = true;
            findMinimal = true;
            findMinimalAmount = 0;

        // If precomputing the solution
        } else if (argAsString == "--precompute") {
            if (i+1 >= argc || argv[i+1][0] == '-') {
                stateFile = "state.dat";
            } else {
                stateFile = std::string(argv[i+1]);
                i++;
            }
            if (stateFile.substr(stateFile.size()-4) != ".dat") {
                stateFile += ".dat";
            }
            precompute = true;
            checkObj = true;
            //lindbladLevel = numQubits;
            findMinimal = true;
            findMinimalAmount = 0;

        // If using a precomputed solution
        } else if (argAsString == "--precomputed") {
            if (i+1 >= argc || argv[i+1][0] == '-') {
                stateFile = "state.dat";
            } else {
                stateFile = std::string(argv[i+1]);
                i++;
            }
            if (stateFile.substr(stateFile.size()-4) != ".dat") {
                stateFile += ".dat";
            }
            precomputed = true;

        // If setting the output precision
        } else if (argAsString == "--prec") {
            precision = std::stoi(argv[i+1]);
            std::cout << std::setprecision(precision);
            i++;

        // If told to only include these qubits
        } else if (argAsString == "-I") {
            std::vector<int> qubitsToUse;
            std::string qubitsAsString = std::string(argv[i+1]);
            i++;
            std::string toConvert = "";
            for (size_t j=0; j<qubitsAsString.size(); j++) {
                if (qubitsAsString[j] != ',') {
                    toConvert += qubitsAsString[j];
                }
                if (qubitsAsString[j] == ',' || j == qubitsAsString.size()-1) {
                    qubitsToUse.push_back(std::stoi(toConvert));
                    toConvert = "";
                }
            }
            for (int q : qubitsToUse) {
                if (q < 1 || q > numQubits) {
                    std::cerr << "Error - invalid qubit index " << q << " specified" << std::endl;
                    return 1;
                }
            }
            includeOnly = qubitsToUse;

        // If told to output all times in milliseconds
        } else if (argAsString == "--millis") {
            allTimeInMillis = true;

        // If a note is given
        } else if (argAsString == "-N") {
            note = std::string(argv[i+1]);
            i++;

        // If considering as a H min problem
        } else if (argAsString == "-H") {
            groundStateProblem = true;

        // If sampling from the symmetries
        } else if (argAsString == "--sym") {
            symSample = true;

        // If sampling using the known ground-state
        } else if (argAsString == "--known") {
            useKnown = true;

        // Output the help
        } else if (argAsString == "-h" || argAsString == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -S <str>            Seed for the random number generator" << std::endl;
            std::cout << "  -C <int>            Number of cores to use" << std::endl;
            std::cout << "  -c <int>            Add some number of extra constraints to reduce the num of vars" << std::endl;
            std::cout << "  -t <dbl>            Set the tolerance (used in various places)" << std::endl;
            std::cout << "  -P                  Take the current objective as a pseudo-objective" << std::endl;
            std::cout << "  -v <int>            Verbosity level" << std::endl;
            std::cout << "  -B                  Benchmarking mode" << std::endl;
            std::cout << "  -X                  Output the Lindbladian in LaTeX form" << std::endl;
            std::cout << "  -f                  Output the final moments to file (monoms.csv)" << std::endl;
            std::cout << "  -H                  Consider it as a Hamiltonian minimization problem" << std::endl;
            std::cout << "  -o                  Check the final results by reconstructing the full state" << std::endl;
            std::cout << "  --prec <int>        Set the output precision (default 10)" << std::endl;
            std::cout << "  -N <str>            Add a note to the run (use with benchmarking mode)" << std::endl;
            std::cout << "  --millis            Output all times in milliseconds (i.e. don't auto convert)" << std::endl;
            std::cout << "  --precompute  <str> Solve the system exactly, saving the state to a file" << std::endl;
            std::cout << "  --precomputed <str> Load the known optimum state from a file" << std::endl;
            std::cout << "Sampling options:" << std::endl;
            std::cout << "  --samples <int>     Solve exactly, take this many samples per Pauli string" << std::endl;
            std::cout << "  --shots   <int>     Same as above, but limit the total number of measurements" << std::endl;
            std::cout << "  -p <int>            Percentile for the error (e.g. 95)" << std::endl;
            std::cout << "  --all <int>         Sample from all Pauli strings evenly up to this degree" << std::endl;
            std::cout << "  --onlyobj           Sample only from the objective" << std::endl;
            std::cout << "  --energy            Sample only the energy" << std::endl;
            std::cout << "  --auto <int>        Sample from the most common Pauli strings" << std::endl;
            std::cout << "  --first <int>        Sample from first Pauli strings" << std::endl;
            std::cout << "  --noobj             Exclude the objective from the sampling" << std::endl;
            std::cout << "  --nox               Exclude any X terms from the sampling" << std::endl;
            std::cout << "  --noy               Exclude any Y terms from the sampling" << std::endl;
            std::cout << "  --noz               Exclude any Z terms from the sampling" << std::endl;
            std::cout << "  --sym               Sample from the symmetries" << std::endl;
            std::cout << "  --known             Sample using the known ground-state" << std::endl;
            std::cout << "Constraint options:" << std::endl;
            std::cout << "  -m <int>            Level of the moment matrix" << std::endl;
            std::cout << "  -l <int>            Level of the moments to put in the Lindbladian" << std::endl;
            std::cout << "  -e <mon>            Add an extra monomial to the top row of the moment matrix" << std::endl;
            std::cout << "  -E <mon>            Add an extra monomial to the list of Lindbladian replacements" << std::endl;
            std::cout << "  -M <int>            Try to generate the minimal set of linear constraints" << std::endl;
            std::cout << "  -A <int>            Generate a moment matrix from the most common moments" << std::endl;
            std::cout << "  -F <int>            Generate a moment matrix from the first occuring moments" << std::endl;
            std::cout << "  -r <int>            Insist that the reconstructed density matrix be positive" << std::endl;
            std::cout << "  -R                  Try removing random constraints" << std::endl;
            std::cout << "  -y <ints>           Add a symetry between two groups e.g. -y 1,2 3,4" << std::endl;
            std::cout << "  -Y                  Use all known symmetries (use before the Lindbladian flag)" << std::endl;
            std::cout << "  -I <ints>           Only put operators with these qubits into the Lindbladian" << std::endl;
            std::cout << "  -T                  Add tracing-out constraints between density mats" << std::endl;
            std::cout << "  -i <int>            1 => ignore imag, 2 => alternate imag handling" << std::endl;
            std::cout << "  --conHC P/N     same as objHC, but constraint to be positive/negative" << std::endl;
            std::cout << "Solver options:" << std::endl;
            std::cout << "  -s G                Use Gurobi as the solver" << std::endl;
            std::cout << "  -s E                Use Eigen as the solver" << std::endl;
            std::cout << "  -s M                Use MOSEK as the solver" << std::endl;
            std::cout << "  -s S                Use SCS as the solver" << std::endl;
            std::cout << "  -s N                Don't solve after generating" << std::endl;
            std::cout << "Objective options:" << std::endl;
            std::cout << "  -O --obj <str>      Manually set the objective" << std::endl;
            std::cout << "  --objX              Use avg sigma_X as the objective" << std::endl;
            std::cout << "  --objY              Use avg sigma_Y as the objective" << std::endl;
            std::cout << "  --objZ              Use avg sigma_Z as the objective" << std::endl;
            std::cout << "  --objXYZ            Sum of all three above objectives" << std::endl;
            std::cout << "  --objHC H/C         Heat current for either the hot or cold bath" << std::endl;
            std::cout << "  --objRand <int>     Random objective up to a certain order" << std::endl;
            std::cout << "  --objPurity         Try to minimize the purity of the state" << std::endl;
            std::cout << "  --objEnergy         Try to minimize the energy of the state" << std::endl;
            std::cout << "Limbliadian choices:" << std::endl;
            std::cout << "  -L <str>            Manually set the Lindbladian" << std::endl;
            std::cout << "  --file <str>        Read the Lindbladian from a file" << std::endl;
            std::cout << "  --2d <int> <int>    (n.b. width then height)" << std::endl;
            std::cout << "  --2dtwo <int> <int>" << std::endl;
            std::cout << "  --2dtfi <int> <int>" << std::endl;
            std::cout << "  --2dtfiperiodic <int> <int>" << std::endl;
            std::cout << "  --mg <int>" << std::endl;
            std::cout << "  --pauli <dbl> <dbl> <dbl>" << std::endl;
            std::cout << "  --second <int> <dbl>" << std::endl;
            std::cout << "  --two" << std::endl;
            std::cout << "  --twov <dbl> <dbl>" << std::endl;
            std::cout << "  --many  <int>" << std::endl;
            std::cout << "  --manyr <int>" << std::endl;
            std::cout << "  --manyv <int> <dbl> <dbl>" << std::endl;
            std::cout << "  --tensor <int>" << std::endl;
            std::cout << "  --david <dbl> <int>" << std::endl;
            std::cout << "  --davidr <dbl> <int>" << std::endl;
            std::cout << "  --david2d <dbl> <int> <int>" << std::endl;
            std::cout << "  --david2dr <dbl> <int> <int>" << std::endl;
            std::cout << "  --quasi" << std::endl;
            std::cout << "  --phase1" << std::endl;
            std::cout << "  --phase2" << std::endl;
            std::cout << "  --phase3" << std::endl;
            return 0;

        // Benchmarking mode
        } else if (argAsString == "-B") {
            verbosity = 0;
            benchmark = true;

        // Ignore the imaginary components of the moment matrix
        } else if (argAsString == "-i") {
            imagType = std::stoi(argv[i+1]);
            i++;

        // If using all symmetries
        } else if (argAsString == "-Y") {
            allSymmetries = true;

        // If auto generating the moment matrix
        } else if (argAsString == "-A") {
            autoMomentAmount = std::stoi(argv[i+1]);
            autoMomentType = "common";
            i++;

        // If auto generating the moment matrix
        } else if (argAsString == "-F") {
            autoMomentAmount = std::stoi(argv[i+1]);
            autoMomentType = "first";
            i++;

        // Otherwise we don't know what this is
        } else if (argAsString != "./run") {
            std::cout << "Unknown argument: " << argAsString << std::endl;
            return 1;

        }
    }

    // If it's a purity problem, use a generic objective for now
    if (usePurity) {
        objective = Poly();
        for (int i=1; i<=numQubits; i++) {
            objective += Poly(1.0/numQubits, "<X" + std::to_string(i) + ">");
            objective += Poly(1.0/numQubits, "<Y" + std::to_string(i) + ">");
            objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
        }
    }

    // If it's a Hamiltonian problem, don't do any Lindblad stuff
    if (groundStateProblem) {
        useEnergy = true;

        // Get the Hamiltonian
        Poly H = Poly();
        for (int i=0; i<numQubits; i++) {
            for (int j=i; j<numQubits; j++) {
                H += hamiltonianInter[i][j];
            }
        }
        H.reduce();

        // Lindbladian = -i[H, rho]
        Poly rho("<R1>");
        lindbladian = -imag*H.commutator(rho);
        lindbladian = Poly("<A0>") * lindbladian;
        lindbladian.cycleToAndRemove('R', 1);
        lindbladian.convertToPaulis();
        lindbladian.reduce();

    }

    // Start a timer
    std::chrono::steady_clock::time_point timeFinishedArgs = std::chrono::steady_clock::now();

    // If using energy as the objective
    if (useEnergy) {

        // The objective is the Hamiltonian
        Poly H = Poly();
        for (int i=0; i<numQubits; i++) {
            for (int j=i; j<numQubits; j++) {
                H += hamiltonianInter[i][j];
            }
        }
        H.reduce();
        objective = H;

    }

    // If we haven't given a pseudo-objective, set it to the objective
    if (pseudoObjective.size() == 0) {
        pseudoObjective = objective;
    }

    // If we're just outputting the Lindbladian in LaTeX form, do that and exit
    if (outputLindbladAsLaTeX) {
        std::cout << "Lindbladian in LaTeX form: " << std::endl;
        std::stringstream ss;
        ss << lindbladian;
        std::string lindbladianLaTeX = ss.str();

        // Replace all A0 with G
        size_t pos = 0;
        while ((pos = lindbladianLaTeX.find("A0", pos)) != std::string::npos) {
            lindbladianLaTeX.replace(pos, 2, "G");
            pos += 1;
        }

        // Replace all < with \braket{
        pos = 0;
        while ((pos = lindbladianLaTeX.find("<", pos)) != std::string::npos) {
            lindbladianLaTeX.replace(pos, 1, "\\braket{");
            pos += 7;
        }

        // Replace all > with }
        pos = 0;
        while ((pos = lindbladianLaTeX.find(">", pos)) != std::string::npos) {
            lindbladianLaTeX.replace(pos, 1, "}");
            pos += 1;
        }

        // Put a newline every 4 }
        pos = 0;
        int numBrackets = 0;
        while ((pos = lindbladianLaTeX.find("}", pos)) != std::string::npos) {
            numBrackets++;
            if (numBrackets % 4 == 0) {
                lindbladianLaTeX.replace(pos, 1, "}\\\\\n&");
                pos += 2;
            }
            pos += 1;
        }

        // Replace all X10 with X_{10} etc.
        for (int i=numQubits; i>=1; i--) {
            std::string toReplace = "X" + std::to_string(i);
            size_t pos = 0;
            while ((pos = lindbladianLaTeX.find(toReplace, pos)) != std::string::npos) {
                lindbladianLaTeX.replace(pos, toReplace.size(), "X_{" + std::to_string(i) + "}");
                pos += 3;
            }
        }

        // Replace all Y10 with Y_{10} etc.
        for (int i=numQubits; i>=1; i--) {
            std::string toReplace = "Y" + std::to_string(i);
            size_t pos = 0;
            while ((pos = lindbladianLaTeX.find(toReplace, pos)) != std::string::npos) {
                lindbladianLaTeX.replace(pos, toReplace.size(), "Y_{" + std::to_string(i) + "}");
                pos += 3;
            }
        }

        // Replace all Z10 with Z_{10} etc.
        for (int i=numQubits; i>=1; i--) {
            std::string toReplace = "Z" + std::to_string(i);
            size_t pos = 0;
            while ((pos = lindbladianLaTeX.find(toReplace, pos)) != std::string::npos) {
                lindbladianLaTeX.replace(pos, toReplace.size(), "Z_{" + std::to_string(i) + "}");
                pos += 3;
            }
        }

        // The starting align
        lindbladianLaTeX = "&" + lindbladianLaTeX;

        // Output and stop
        std::cout << lindbladianLaTeX << std::endl;
        return 0;

    }

    // Determine which monomials to put in the Lindbladian
    std::set<Mon> monomsUsed;
    std::vector<Mon> monomsUsedVec;
    std::vector<std::vector<std::vector<Poly>>> momentMatrices = {};
    std::vector<Mon> variables = {};
    for (int i=0; i<numQubits; i++) {
        variables.push_back(Mon("<X" + std::to_string(i+1) + ">"));
        variables.push_back(Mon("<Y" + std::to_string(i+1) + ">"));
        variables.push_back(Mon("<Z" + std::to_string(i+1) + ">"));
    }
    addSingleMonomials(variables, objective);
    std::vector<Poly> variablesToPut = generateMonomials(variables, lindbladLevel, verbosity);
    for (size_t i=0; i<extraMonomialsLim.size(); i++) {
        variablesToPut.push_back(Poly(extraMonomialsLim[i]));
    }
    if (includeOnly.size() > 0) {
        std::vector<Poly> variablesToPutFiltered;
        for (size_t i=0; i<variablesToPut.size(); i++) {
            Mon mon = variablesToPut[i].getKey();
            bool allGood = true;
            for (int j=0; j<mon.size(); j++) {
                if (std::find(includeOnly.begin(), includeOnly.end(), mon[j].second) == includeOnly.end()) {
                    allGood = false;
                    break;
                }
            }
            if (allGood) {
                variablesToPutFiltered.push_back(variablesToPut[i]);
            }
        }
        variablesToPut = variablesToPutFiltered;
    }
    if (verbosity >= 2) {
        std::cout << std::endl;
        std::cout << "Variables to put in Lindbladian: " << std::endl;
        for (size_t i=0; i<variablesToPut.size(); i++) {
            std::cout << variablesToPut[i] << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Adjoint Lindbladian: " << lindbladian << std::endl;
        if (usingHeatCurrent) {
            std::cout << "Lindbladian hot: " << lindbladianHot << std::endl;
            std::cout << "Lindbladian cold: " << lindbladianCold << std::endl;
            Poly hamiltonianFull;
            for (int j=0; j<numQubits; j++) {
                for (int k=j; k<numQubits; k++) {
                    hamiltonianFull += hamiltonianInter[j][k];
                }
            }
            std::cout << "Hamiltonian full: " << hamiltonianFull << std::endl;
        }
    }

    // Put them in the Lindbladian
    for (size_t i=0; i<variablesToPut.size(); i++) {
        std::pair<std::complex<double>, Mon> reducedMon = variablesToPut[i].getKey().reduce();
        if (!monomsUsed.count(reducedMon.second)) {
            std::pair<char,int> oldMon('A', 0);
            Mon newPoly(reducedMon.second);
            Poly newConstraint = lindbladian.replaced(oldMon, newPoly);
            if (verbosity >= 3) {
                std::cout << std::endl;
                std::cout << "Variable to put: " << newPoly << std::endl;
                std::cout << "New constraint: " << newConstraint << std::endl;
            }
            if (!newConstraint.isZero()) {
                constraintsZero.push_back(newConstraint);
            }
            monomsUsed.insert(variablesToPut[i].getKey());
            monomsUsedVec.push_back(variablesToPut[i].getKey());
        }
    }

    // If we don't have a Lindbladian, add things to monoms used
    if (lindbladian.size() == 0) {
        for (size_t i=0; i<variablesToPut.size(); i++) {
            monomsUsed.insert(variablesToPut[i].getKey());
            monomsUsedVec.push_back(variablesToPut[i].getKey());
        }
        std::cout << "Added " << variablesToPut.size() << " variables to the monomials used" << std::endl;
    }

    // Reduce the constraints as much as possible
    for (size_t i=0; i<constraintsZero.size(); i++) {
        if (verbosity >= 3) {
            std::cout << std::endl;
            std::cout << "Reducing constraint: " << constraintsZero[i] << std::endl;
        }
        constraintsZero[i].reduce();
        if (verbosity >= 3) {
            std::cout << "Reduced constraint: " << constraintsZero[i] << std::endl;
        }
    }

    // If asked to find the minimal set of linear constraints 
    if (findMinimal) {

        // If given zero, set to what seems to be the minimum needed
        if (findMinimalAmount == 0) {
            findMinimalAmount = std::pow(4, numQubits) - 1;
            if (verbosity >= 1) {
                std::cout << "Auto setting constraint limit to: " << findMinimalAmount << std::endl;
            }
        }

        // If a max num of constraints is given
        if (findMinimalAmount > 0) {

            // Add constraints based on the monomials we already have
            std::set<Mon> monomsInConstraints;
            std::vector<Mon> queue;
            queue.push_back(Mon());
            monomsUsed.insert(Mon());
            monomsInConstraints.insert(Mon());
            for (auto& term : pseudoObjective) {
                if (!monomsUsed.count(term.first)) {
                    monomsUsed.insert(term.first);
                    queue.push_back(term.first);
                }
            }
            for (size_t i=0; i<variablesToPut.size(); i++) {
                Mon monToAdd = variablesToPut[i].getKey();
                if (!monomsUsed.count(monToAdd)) {
                    monomsUsed.insert(monToAdd);
                    queue.push_back(monToAdd);
                }
            }
            for (size_t i=0; i<constraintsZero.size(); i++) {
                for (auto& term : constraintsZero[i]) {
                    if (!monomsUsed.count(term.first)) {
                        monomsUsed.insert(term.first);
                        monomsInConstraints.insert(term.first);
                        queue.push_back(term.first);
                    }
                }
            }
            for (size_t i=0; i<momentMatrices.size(); i++) {
                for (size_t j=0; j<momentMatrices[i].size(); j++) {
                    for (size_t k=0; k<momentMatrices[i][j].size(); k++) {
                        Mon monToAdd = momentMatrices[i][j][k].getKey();
                        if (!monomsUsed.count(monToAdd)) {
                            monomsUsed.insert(monToAdd);
                            monomsInConstraints.insert(monToAdd);
                            queue.push_back(monToAdd);
                        }
                    }
                }
            }

            // Verbose output
            if (verbosity >= 3) {
                std::cout << "Initial queue:";
                for (size_t i=0; i<queue.size(); i++) {
                    std::cout << queue[i] << ", ";
                }
                std::cout << std::endl;
                std::cout << "Monomials used: ";
                for (auto& mon : monomsUsed) {
                    std::cout << mon << ", ";
                }
                std::cout << std::endl;
                std::cout << "Monomials in constraints: ";
                for (auto& mon : monomsInConstraints) {
                    std::cout << mon << ", ";
                }
                std::cout << std::endl;
            }

            // Keep putting things back into the Lindbladian
            int nextQueueLoc = 0;
            while (int(constraintsZero.size()) < findMinimalAmount || precompute) {
                double ratio = double(constraintsZero.size()) / (double(monomsInConstraints.size())-1);
                if (verbosity >= 1) {
                    std::cout << monomsInConstraints.size() << " vars, " << constraintsZero.size() << " cons, aim: " << findMinimalAmount << " (" << ratio << ")                         \r" << std::flush;
                }

                // Stop if we're fully constrained
                if (!precompute && ratio >= 1.0 && constraintsZero.size() > 2) {
                    break;
                } else if (precompute && monomsInConstraints.size() >= 1+findMinimalAmount && 
                                         constraintsZero.size() >= findMinimalAmount) {
                    break;
                }

                // Find a monomial that hasn't been used
                Mon monToAdd;
                bool found = false;
                if (nextQueueLoc < int(queue.size())) {
                    monToAdd = queue[nextQueueLoc];
                    found = true;
                    nextQueueLoc++;

                // Otherwise we're looping, need to break out of the cycle
                } else {

                    // Try removing a term from the start of a monomial
                    for (auto& mon : monomsUsed) {
                        if (mon.size() >= 2) {
                            Mon reducedMon = mon;
                            reducedMon.monomial.erase(reducedMon.monomial.begin());
                            if (!monomsUsed.count(reducedMon)) {
                                monToAdd = reducedMon;
                                found = true;
                                break;
                            }
                        }
                    }

                    // Try combining two terms
                    if (!found) {
                        for (auto it1 = monomsUsed.begin(); it1 != monomsUsed.end(); ++it1) {
                            for (auto it2 = std::next(it1); it2 != monomsUsed.end(); ++it2) {
                                Mon combinedMon = (*it1) * (*it2);
                                std::pair<std::complex<double>, Mon> reducedMon = combinedMon.reduce();
                                if (!monomsUsed.count(reducedMon.second)) {
                                    monToAdd = reducedMon.second;
                                    found = true;
                                    break;
                                }
                            }
                            if (found) {
                                break;
                            }
                        }
                    }

                    // If we found something new, output
                    if (verbosity >= 2 && found) {
                        std::cout << std::endl;
                        std::cout << "Starting new cycle with monomial: " << monToAdd << std::endl;
                    }

                }

                // If we can't find anything else, break
                if (!found) {
                    if (verbosity >= 2) {
                        std::cout << std::endl;
                        std::cout << "Couldn't find any more monomials to add" << std::endl;
                    }
                    break;
                }

                // Put the monomial in the Lindbladian
                Mon toPut(monToAdd);
                if (verbosity >= 3) {
                    std::cout << std::endl;
                    std::cout << "Putting in monomial: " << toPut << std::endl;
                }
                std::pair<char,int> oldMon('A', 0);
                Poly newConstraint = lindbladian.replaced(oldMon, toPut);

                // Add the constraint
                if (!newConstraint.isZero()) {
                    constraintsZero.push_back(newConstraint); 
                }
                monomsUsed.insert(monToAdd);
                monomsUsedVec.push_back(monToAdd);
                std::set<Mon> newTerms;
                for (auto& term : newConstraint) {
                    newTerms.insert(term.first);
                }
                for (auto& term : newTerms) {
                    monomsInConstraints.insert(term);
                    if (!monomsUsed.count(term)) {
                        bool allGood = true;
                        if (includeOnly.size() > 0) {
                            for (int j=0; j<term.size(); j++) {
                                if (std::find(includeOnly.begin(), includeOnly.end(), term[j].second) == includeOnly.end()) {
                                    allGood = false;
                                    break;
                                }
                            }
                        }
                        if (allGood) {
                            queue.push_back(term);
                            monomsUsed.insert(term);
                        }
                    }
                }

            }

            // Final output
            double ratio = double(constraintsZero.size()) / (monomsInConstraints.size()-1);
            if (verbosity >= 1) {
                std::cout << constraintsZero.size() << " / " << findMinimalAmount << " (" << ratio << ")               " << std::endl;
            }

        // If a value not given, binary search
        } else {

            // First try adding constraints until it's exact
            std::set<Mon> monomsInConstraints;
            std::vector<Mon> queue;
            queue.push_back(Mon());
            monomsUsed.insert(Mon());
            monomsInConstraints.insert(Mon());
            for (auto& term : pseudoObjective) {
                if (!monomsUsed.count(term.first)) {
                    monomsUsed.insert(term.first);
                    queue.push_back(term.first);
                }
            }
            for (size_t i=0; i<variablesToPut.size(); i++) {
                Mon monToAdd = variablesToPut[i].getKey();
                if (!monomsUsed.count(monToAdd)) {
                    monomsUsed.insert(monToAdd);
                    queue.push_back(monToAdd);
                }
            }
            for (size_t i=0; i<constraintsZero.size(); i++) {
                for (auto& term : constraintsZero[i]) {
                    if (!monomsUsed.count(term.first)) {
                        monomsUsed.insert(term.first);
                        monomsInConstraints.insert(term.first);
                        queue.push_back(term.first);
                    }
                }
            }
            for (size_t i=0; i<momentMatrices.size(); i++) {
                for (size_t j=0; j<momentMatrices[i].size(); j++) {
                    for (size_t k=0; k<momentMatrices[i][j].size(); k++) {
                        Mon monToAdd = momentMatrices[i][j][k].getKey();
                        if (!monomsUsed.count(monToAdd)) {
                            monomsUsed.insert(monToAdd);
                            monomsInConstraints.insert(monToAdd);
                            queue.push_back(monToAdd);
                        }
                    }
                }
            }

            // Keep branching
            constraintsZero.clear();
            int amountToTest = 50;
            int nextQueueLoc = 0;
            for (int l=0; l<20; l++) {

                // Add constraints based on the monomials we already have
                while (int(constraintsZero.size()) < amountToTest && monomsUsed.size() < monomsInConstraints.size()) {

                    // Find a monomial that hasn't been used
                    Mon monToAdd;
                    if (nextQueueLoc < int(queue.size())) {
                        monToAdd = queue[nextQueueLoc];
                        nextQueueLoc++;

                    // Otherwise we're looping, need to break out of the cycle
                    } else {
                        for (auto& mon : monomsInConstraints) {
                            if (mon.size() >= 2) {
                                Mon reducedMon = mon;
                                reducedMon.monomial.erase(reducedMon.monomial.begin());
                                if (!monomsUsed.count(reducedMon)) {
                                    monToAdd = reducedMon;
                                    break;
                                }
                            }
                        }
                    }

                    // Put the monomial in the Lindbladian
                    Mon toPut(monToAdd);
                    std::pair<char,int> oldMon('A', 0);
                    Poly newConstraint = lindbladian.replaced(oldMon, toPut);

                    // Add the constraint
                    if (!newConstraint.isZero()) {
                        constraintsZero.push_back(newConstraint); 
                    }
                    monomsUsed.insert(monToAdd);
                    monomsUsedVec.push_back(monToAdd);
                    for (auto& term : newConstraint) {
                        if (!monomsInConstraints.count(term.first)) {
                            monomsInConstraints.insert(term.first);
                            queue.push_back(term.first);
                        }
                    }

                }

                // Run the SDP
                std::pair<double,double> boundsTemp = boundMOSEK(objective, momentMatrices, constraintsZero, {}, verbosity);
                double lowerBoundTemp = boundsTemp.first;
                double upperBoundTemp = boundsTemp.second;
                double diff = upperBoundTemp - lowerBoundTemp;
                if (verbosity >= 1) {
                    std::cout << "cons: " << amountToTest << ", diff: " << diff << std::endl;
                }
                if (diff < 1e-5) {
                    break;
                } else {
                    amountToTest *= 2;
                }

            }

            // Now try removing constraints until it's not exact
            int minNum = 0;
            int maxNum = int(constraintsZero.size());
            std::vector<Poly> constraintsZeroCopy = constraintsZero;
            for (int l=0; l<100; l++) {

                // Subset of constraints
                int toTest = (maxNum + minNum) / 2;
                constraintsZero.clear();
                for (int i=0; i<toTest; i++) {
                    constraintsZero.push_back(constraintsZeroCopy[i]);
                }

                // Run the SDP
                std::pair<double,double> boundsTemp = boundMOSEK(objective, momentMatrices, constraintsZero, {}, verbosity);
                double lowerBoundTemp = boundsTemp.first;
                double upperBoundTemp = boundsTemp.second;
                double diff = upperBoundTemp - lowerBoundTemp;
                if (verbosity >= 1) {
                    std::cout << "cons: " << toTest << ", diff: " << diff << std::endl;
                }
                if (diff < 1e-5) {
                    maxNum = toTest;
                } else {
                    minNum = toTest;
                }
                if (maxNum - minNum <= 1) {
                    break;
                }

            }

            // Output the minimal number of constraints
            if (verbosity >= 1) {
                std::cout << "Minimal num constraints: " << maxNum << std::endl;
            }
            constraintsZero = constraintsZeroCopy;

        }

    }

    // Generate the moment matrix from the monomsUsed
    if (autoMomentAmount > 0) {

        // Begin forming the top row of the moment matrix
        std::vector<Poly> topRow = {Poly(1)};
        std::set<Mon> monomsInMat;

        // Add the most common monoms to the top row
        if (autoMomentType == "common") { 

            // Determine the most-used monomial
            std::map<Mon, int> monCounts;
            for (auto& term : pseudoObjective) {
                monCounts[term.first]++;
            }
            for (size_t i=0; i<constraintsZero.size(); i++) {
                for (auto& term : constraintsZero[i]) {
                    monCounts[term.first]++;
                }
            }
            std::vector<Mon> sortedMonoms;
            for (auto& monCount : monCounts) {
                sortedMonoms.push_back(monCount.first);
            }
            std::sort(sortedMonoms.begin(), sortedMonoms.end(), 
                      [](const Mon& a, const Mon& b) { return a.size() < b.size(); });

            // Add them to the top row
            int added = 1;
            for (auto& mon : sortedMonoms) {
                if (mon.size() == 0) {
                    continue;
                }
                topRow.push_back(Poly(mon));
                monomsInMat.insert(mon);
                added++;
                if (added >= autoMomentAmount) {
                    break;
                }
            }

        // Otherwise just add the used monoms to the top row
        } else {
            std::set<Mon> toPut = monomsUsed;
            for (auto& term : pseudoObjective) {
                toPut.insert(term.first);
            }
            int added = 1;
            for (auto& mon : toPut) {
                if (mon.size() == 0) {
                    continue;
                }
                topRow.push_back(Poly(mon));
                monomsInMat.insert(mon);
                added++;
                if (added >= autoMomentAmount) {
                    break;
                }
            }
        }

        // If the top row isn't big enough, add products of things in it
        while (int(topRow.size()) < autoMomentAmount) {

            // Find the next monomial to add
            Mon nextMon;
            for (size_t i=0; i<topRow.size(); i++) {
                for (size_t j=i; j<topRow.size(); j++) {
                    Mon product = topRow[i].getKey() * topRow[j].getKey();
                    std::pair<std::complex<double>, Mon> reducedMon = product.reduce();
                    if (!monomsInMat.count(reducedMon.second) && reducedMon.second.size() > 0) {
                        nextMon = reducedMon.second;
                        break;
                    }
                }
                if (nextMon.size() > 0) {
                    break;
                }
            }

            // If we found a new monomial, add it
            if (nextMon.size() > 0) {
                topRow.push_back(Poly(nextMon));
                monomsInMat.insert(nextMon);
            } else {
                break;
            }

        }

        // Generate the moment matrix from the top row
        momentMatrices.push_back(generateFromTopRow(topRow, verbosity));

    }

    // If adding some reductive constraints
    int newReductConsAdded = 0;
    while (newReductConsAdded < reductiveCons) {
        if (verbosity >= 1) {
            std::cout << newReductConsAdded << " / " << reductiveCons << "                      \r" << std::flush;
        }

        // Add constraints based on the monomials we already have
        std::set<Mon> monomsInConstraints;
        std::vector<Mon> monomsInConstraintsVec;
        for (auto& term : pseudoObjective) {
            if (!monomsInConstraints.count(term.first)) {
                monomsInConstraints.insert(term.first);
                monomsInConstraintsVec.push_back(term.first);
            }
        }
        for (size_t i=0; i<constraintsZero.size(); i++) {
            for (auto& term : constraintsZero[i]) {
                if (!monomsInConstraints.count(term.first)) {
                    monomsInConstraints.insert(term.first);
                    monomsInConstraintsVec.push_back(term.first);
                }
            }
        }
        for (size_t i=0; i<variablesToPut.size(); i++) {
            Mon monToAdd = variablesToPut[i].getKey();
            if (!monomsInConstraints.count(monToAdd)) {
                monomsInConstraints.insert(monToAdd);
                monomsInConstraintsVec.push_back(monToAdd);
            }
        }

        // Find a monomial that hasn't been used
        Mon monToAdd;
        if (monomsUsed.size() < monomsInConstraints.size()) {
            for (auto& mon : monomsInConstraintsVec) {
                if (!monomsUsed.count(mon)) {
                    monToAdd = mon;
                    break;
                }
            }
        }

        // If we can't find anything else, try a new cycle
        if (monToAdd.size() == 0) {
            for (auto& mon : monomsInConstraintsVec) {
                if (mon.size() >= 2) {
                    Mon reducedMon = mon;
                    reducedMon.monomial.erase(reducedMon.monomial.begin());
                    if (!monomsUsed.count(reducedMon)) {
                        monToAdd = reducedMon;
                        break;
                    }
                }
            }
        }

        // If still nothing, break
        if (monToAdd.size() == 0 && monomsUsed.count(monToAdd)) {
            break;
        }

        // Put the monomial in the Lindbladian
        Mon toPut(monToAdd);
        std::pair<char,int> oldMon('A', 0);
        Poly newConstraint = lindbladian.replaced(oldMon, toPut);

        // Check if the new constraint adds more variables
        int numNewVars = 0;
        for (auto& term : newConstraint) {
            if (!monomsInConstraints.count(term.first)) {
                numNewVars++;
            }
        }

        // If it doesn't, add it
        double currentRatio = double(constraintsZero.size()) / double(monomsInConstraints.size());
        double newRatio = double(constraintsZero.size() + 1) / double(monomsInConstraints.size() + numNewVars);
        if (!newConstraint.isZero() && newRatio > currentRatio) {
            constraintsZero.push_back(newConstraint);
            for (auto& term : newConstraint) {
                if (!monomsInConstraints.count(term.first)) {
                    monomsInConstraints.insert(term.first);
                    monomsInConstraintsVec.push_back(term.first);
                }
            }
            newReductConsAdded++;
        }

        // Regardless, don't try to add it again
        monomsUsed.insert(monToAdd);
        monomsUsedVec.push_back(monToAdd);

    }

    // Create the moment matrices
    if (autoMomentAmount == 0) {
        momentMatrices = generateAllMomentMatrices(pseudoObjective, constraintsZero, level, verbosity);
    }

    // If told to add extra to the top row
    if (extraMonomials.size() > 0) {
        std::vector<Poly> topRow = momentMatrices[0][0];
        for (size_t i=0; i<extraMonomials.size(); i++) {
            Poly extraMonomial(Mon(extraMonomials[i]).reversed());
            topRow.push_back(extraMonomial);
            if (verbosity >= 1) {
                std::cout << "Added " << extraMonomial << " to the top row" << std::endl;
            }
        }
        momentMatrices[0] = generateFromTopRow(topRow, verbosity);
    }

    // See how big the moment matrix is
    if (verbosity >= 1) {
        int largestMomentMatrix = 0;
        for (size_t i=0; i<momentMatrices.size(); i++) {
            if (int(momentMatrices[i].size()) > largestMomentMatrix) {
                largestMomentMatrix = momentMatrices[i].size();
            }
        }
    }

    // Add symmetries
    std::set<Mon> vars;
    for (size_t i=0; i<momentMatrices.size(); i++) {
        addVariables(vars, momentMatrices[i]);
    }
    addVariables(vars, constraintsZero);
    addVariables(vars, pseudoObjective);
    addVariables(vars, constraintsPositive);
    addVariables(vars, objective);
    std::map<Mon,Mon> symmetriesMap;
    for (auto sym : symmetries) {
        
        // If verbose, output the symmetry
        if (verbosity >= 2) {
            std::cout << "Adding symmetry between " << sym.first << " and " << sym.second << std::endl;
        }

        // Handy dandy maps
        std::map<int,int> from2to1;
        std::map<int,int> from1to2;
        for (size_t i=0; i<sym.first.size(); i++) {
            from2to1[sym.second[i]] = sym.first[i];
            from1to2[sym.first[i]] = sym.second[i];
        }

        // For each monomial used
        for (auto& mon : vars) {

            // Check if the monomial only contains things in the group
            bool allInGroup = true;
            for (auto& term : mon) {
                if (!from1to2.count(term.second) && !from2to1.count(term.second)) {
                    allInGroup = false;
                    break;
                }
            }

            // Try swapping and see if it's lexically smaller
            if (allInGroup) {
                Mon newMon;
                for (auto& term : mon) {
                    if (from1to2.count(term.second)) {
                        newMon.monomial.push_back({term.first, from1to2[term.second]});
                    } else {
                        newMon.monomial.push_back({term.first, from2to1[term.second]});
                    }
                    
                }
                if (vars.count(newMon)) {
                    symmetriesMap[mon] = newMon;
                    if (verbosity >= 3) {
                        std::cout << mon << " -> " << newMon << std::endl;
                    }
                }
            }

        }

    }

    // If using the symmetries to sample from
    int numSyms = 0;
    if (symSample) {

        // Note the symmetries for later sampling
        for (auto sym : symmetriesMap) {
            Poly newCon = Poly(sym.first)-Poly(sym.second);
            if (newCon.size() > 0 && !newCon.isZero()) {
                zeroConsForSampling.push_back(newCon);
            }
        }

    } else {

        // Add symmetry constraints
        for (auto sym : symmetriesMap) {
            Poly newCon = Poly(sym.first)-Poly(sym.second);
            if (newCon.size() > 0 && !newCon.isZero()) {
                constraintsZero.push_back(newCon);
                numSyms++;
            }
        }

    }

    // Add constraints that the reconstructed density matrix is positive
    if (reconLevel != 0) {

        // Keep track of what all the various density mats are
        std::map<std::vector<int>, int> sitesToInd;

        // Single-body terms
        if (reconLevel >= 1) {

            // For each qubit
            for (int i=0; i<numQubits; i++) {
                int matSize = std::pow(2, 1);
                std::vector<std::vector<Poly>> rho1(matSize, std::vector<Poly>(matSize));
                double scaling = std::pow(2, -1);

                // For each Pauli matrix
                for (int l=0; l<4; l++) {
                    Eigen::SparseMatrix<std::complex<double>> mat = generatePauliMatrix({l});

                    // For each element of said matrix
                    for (int k=0; k<mat.outerSize(); ++k) {
                        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(mat,k); it; ++it) {
                            std::complex<double> val = it.value();
                            int j = it.row();
                            std::string monString = "<";
                            if (l > 0) {
                                monString += letterMap[l] + std::to_string(i+1);
                            }
                            monString += ">";
                            if (monString.size() > 2) {
                                rho1[j][k] += scaling * val * Mon(monString);
                            } else {
                                rho1[j][k] += scaling * val;
                            }
                        }
                    }

                }

                // Add the matrix to the list of mats that should be positive
                momentMatrices.push_back(rho1);
                sitesToInd[{i}] = momentMatrices.size()-1;

            }
        }

        // Two-body terms
        if (reconLevel >= 2) {

            // For each selection of qubits
            for (int i=0; i<numQubits; i++) {
                for (int i2=i+1; i2<numQubits; i2++) {
                    int matSize = std::pow(2, 2);
                    std::vector<std::vector<Poly>> rho2(matSize, std::vector<Poly>(matSize));
                    double scaling = std::pow(2, -2);

                    // For each Pauli matrix
                    for (int l=0; l<4; l++) {
                        for (int l2=0; l2<4; l2++) {
                        Eigen::SparseMatrix<std::complex<double>> mat = generatePauliMatrix({l,l2});

                            // For each element of said matrix
                            for (int k=0; k<mat.outerSize(); ++k) {
                                for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(mat,k); it; ++it) {
                                    std::complex<double> val = it.value();
                                    int j = it.row();
                                    std::string monString = "<";
                                    if (l > 0) {
                                        monString += letterMap[l] + std::to_string(i+1);
                                    }
                                    if (l2 > 0) {
                                        monString += letterMap[l2] + std::to_string(i2+1);
                                    }
                                    monString += ">";
                                    if (monString.size() > 2) {
                                        rho2[j][k] += scaling * val * Mon(monString);
                                    } else {
                                        rho2[j][k] += scaling * val;
                                    }
                                }
                            }

                        }
                    }

                    // Add the matrix to the list of mats that should be positive
                    momentMatrices.push_back(rho2);
                    sitesToInd[{i,i2}] = momentMatrices.size()-1;

                }
            }

        }

        // Three-body terms
        if (reconLevel >= 3) {

            // For each selection of qubits
            for (int i=0; i<numQubits; i++) {
                for (int i2=i+1; i2<numQubits; i2++) {
                    for (int i3=i2+1; i3<numQubits; i3++) {
                        int matSize = std::pow(2, 3);
                        std::vector<std::vector<Poly>> rho3(matSize, std::vector<Poly>(matSize));
                        double scaling = std::pow(2, -3);

                        // For each Pauli matrix
                        for (int l=0; l<4; l++) {
                            for (int l2=0; l2<4; l2++) {
                                for (int l3=0; l3<4; l3++) {

                                    // For each element of said matrix
                                    Eigen::SparseMatrix<std::complex<double>> mat = generatePauliMatrix({l,l2,l3});
                                    for (int k=0; k<mat.outerSize(); ++k) {
                                        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(mat,k); it; ++it) {
                                            std::complex<double> val = it.value();
                                            int j = it.row();
                                            std::string monString = "<";
                                            if (l > 0) {
                                                monString += letterMap[l] + std::to_string(i+1);
                                            }
                                            if (l2 > 0) {
                                                monString += letterMap[l2] + std::to_string(i2+1);
                                            }
                                            if (l3 > 0) {
                                                monString += letterMap[l3] + std::to_string(i3+1);
                                            }
                                            monString += ">";
                                            if (monString.size() > 2) {
                                                rho3[j][k] += scaling * val * Mon(monString);
                                            } else {
                                                rho3[j][k] += scaling * val;
                                            }
                                        }
                                    }

                                }
                            }
                        }

                        // Add the matrix to the list of mats that should be positive
                        momentMatrices.push_back(rho3);
                        sitesToInd[{i,i2,i3}] = momentMatrices.size()-1;

                    }
                }   
            }

        }

        // Four body terms
        if (reconLevel >= 4) {

            // For each selection of qubits
            for (int i=0; i<numQubits; i++) {
                for (int i2=i+1; i2<numQubits; i2++) {
                    for (int i3=i2+1; i3<numQubits; i3++) {
                        for (int i4=i3+1; i4<numQubits; i4++) {
                            int matSize = std::pow(2, 4);
                            std::vector<std::vector<Poly>> rho4(matSize, std::vector<Poly>(matSize));
                            double scaling = std::pow(2, -4);

                            // For each Pauli matrix
                            for (int l=0; l<4; l++) {
                                for (int l2=0; l2<4; l2++) {
                                    for (int l3=0; l3<4; l3++) {
                                        for (int l4=0; l4<4; l4++) {

                                            // For each element of said matrix
                                            Eigen::SparseMatrix<std::complex<double>> mat = generatePauliMatrix({l,l2,l3,l4});
                                            for (int k=0; k<mat.outerSize(); ++k) {
                                                for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(mat,k); it; ++it) {
                                                    std::complex<double> val = it.value();
                                                    int j = it.row();
                                                    std::string monString = "<";
                                                    if (l > 0) {
                                                        monString += letterMap[l] + std::to_string(i+1);
                                                    }
                                                    if (l2 > 0) {
                                                        monString += letterMap[l2] + std::to_string(i2+1);
                                                    }
                                                    if (l3 > 0) {
                                                        monString += letterMap[l3] + std::to_string(i3+1);
                                                    }
                                                    if (l4 > 0) {
                                                        monString += letterMap[l4] + std::to_string(i4+1);
                                                    }
                                                    monString += ">";
                                                    if (monString.size() > 2) {
                                                        rho4[j][k] += scaling * val * Mon(monString);
                                                    } else {
                                                        rho4[j][k] += scaling * val;
                                                    }
                                                }
                                            }

                                        }
                                    }
                                }
                            }

                            // Add the matrix to the list of mats that should be positive
                            momentMatrices.push_back(rho4);
                            sitesToInd[{i,i2,i3,i4}] = momentMatrices.size()-1;

                        }
                    }
                }
            }
        }

        // Five body terms
        if (reconLevel >= 5) {

            // For each selection of qubits
            for (int i=0; i<numQubits; i++) {
                for (int i2=i+1; i2<numQubits; i2++) {
                    for (int i3=i2+1; i3<numQubits; i3++) {
                        for (int i4=i3+1; i4<numQubits; i4++) {
                            for (int i5=i4+1; i5<numQubits; i5++) {
                                int matSize = std::pow(2, 5);
                                std::vector<std::vector<Poly>> rho5(matSize, std::vector<Poly>(matSize));
                                double scaling = std::pow(2, -5);

                                // For each Pauli matrix
                                for (int l=0; l<4; l++) {
                                    for (int l2=0; l2<4; l2++) {
                                        for (int l3=0; l3<4; l3++) {
                                            for (int l4=0; l4<4; l4++) {
                                                for (int l5=0; l5<4; l5++) {

                                                    // For each element of said matrix
                                                    Eigen::SparseMatrix<std::complex<double>> mat = generatePauliMatrix({l,l2,l3,l4,l5});
                                                    for (int k=0; k<mat.outerSize(); ++k) {
                                                        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(mat,k); it; ++it) {
                                                            std::complex<double> val = it.value();
                                                            int j = it.row();
                                                            std::string monString = "<";
                                                            if (l > 0) {
                                                                monString += letterMap[l] + std::to_string(i+1);
                                                            }
                                                            if (l2 > 0) {
                                                                monString += letterMap[l2] + std::to_string(i2+1);
                                                            }
                                                            if (l3 > 0) {
                                                                monString += letterMap[l3] + std::to_string(i3+1);
                                                            }
                                                            if (l4 > 0) {
                                                                monString += letterMap[l4] + std::to_string(i4+1);
                                                            }
                                                            if (l5 > 0) {
                                                                monString += letterMap[l5] + std::to_string(i5+1);
                                                            }
                                                            monString += ">";
                                                            if (monString.size() > 2) {
                                                                rho5[j][k] += scaling * val * Mon(monString);
                                                            } else {
                                                                rho5[j][k] += scaling * val;
                                                            }
                                                        }
                                                    }
                                                    
                                                }
                                            }
                                        }
                                    }
                                }

                                // Add the matrix to the list of mats that should be positive
                                momentMatrices.push_back(rho5);
                                sitesToInd[{i,i2,i3,i4,i5}] = momentMatrices.size()-1;

                            }
                        }
                    }
                }
            }
        }

        // Add trace out constraints
        if (traceCons) {

            // Output all indices 
            if (verbosity >= 1) {
                for (auto& site : sitesToInd) {
                    std::cout << site.first << " -> " << site.second << std::endl;
                }
            }

            // Determine which qubits to trace out
            std::vector<std::pair<std::vector<int>, int>> thingsToDo;
            if (reconLevel >= 3) {
                for (int i=0; i<numQubits; i++) {
                    for (int i2=i+1; i2<numQubits; i2++) {
                        for (int i3=i2+1; i3<numQubits; i3++) {
                            std::vector<int> sites = {i,i2,i3};
                            for (int k=0; k<sites.size(); k++) {
                                thingsToDo.push_back({sites, k});
                            }
                        }
                    }
                }
            }
            if (reconLevel >= 4) {
                for (int i=0; i<numQubits; i++) {
                    for (int i2=i+1; i2<numQubits; i2++) {
                        for (int i3=i2+1; i3<numQubits; i3++) {
                            for (int i4=i3+1; i4<numQubits; i4++) {
                                std::vector<int> sites = {i,i2,i3,i4};
                                for (int k=0; k<sites.size(); k++) {
                                    thingsToDo.push_back({sites, k});
                                }
                            }
                        }
                    }
                }
            }
            if (reconLevel >= 5) {
                for (int i=0; i<numQubits; i++) {
                    for (int i2=i+1; i2<numQubits; i2++) {
                        for (int i3=i2+1; i3<numQubits; i3++) {
                            for (int i4=i3+1; i4<numQubits; i4++) {
                                for (int i5=i4+1; i5<numQubits; i5++) {
                                    std::vector<int> sites = {i,i2,i3,i4,i5};
                                    for (int k=0; k<sites.size(); k++) {
                                        thingsToDo.push_back({sites, k});
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Generate the constraints
            for (auto& thing : thingsToDo) {
                if (verbosity >= 1) {
                    std::cout << "Tracing out " << thing.first << " at " << thing.second << std::endl;
                }
                std::vector<int> fullSites = thing.first;
                int toTrace = thing.second;
                std::vector<int> reducedSites;
                for (int i=0; i<fullSites.size(); i++) {
                    if (i != toTrace) {
                        reducedSites.push_back(fullSites[i]);
                    }
                }
                std::vector<std::vector<Poly>> mat12 = momentMatrices[sitesToInd[reducedSites]];
                std::vector<std::vector<Poly>> mat12T = partialTrace(momentMatrices[sitesToInd[fullSites]], toTrace);
                for (int i=0; i<mat12.size(); i++) {
                    for (int j=0; j<mat12.size(); j++) {
                        constraintsZero.push_back(mat12[i][j] - mat12T[i][j]);
                    }
                }
            }

        }

    }

    // Timing
    std::chrono::steady_clock::time_point timeFinishedGeneration = std::chrono::steady_clock::now();

    // If removing constraints
    if (tryRemove) {

        // Initial run
        std::pair<double,double> boundsStart = boundMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity);
        double currentBoundDiff = boundsStart.second - boundsStart.first;
        if (verbosity >= 1) {
            std::cout << "Original diff: " << currentBoundDiff << std::endl;
        }
        std::vector<Poly> constraintsZeroCopy = constraintsZero;
        std::vector<std::vector<std::vector<Poly>>> momentMatricesCopy = momentMatrices;

        // Per pass
        for (int l=0; l<100; l++) {
            bool removedSomething = false;
            
            // Check each linear constraint
            for (size_t i=0; i<constraintsZero.size(); i++) {

                // Remove the constraint
                constraintsZero = constraintsZeroCopy;
                constraintsZero.erase(constraintsZero.begin() + i);

                // Get the bounds
                std::pair<double,double> boundsTemp;
                boundsTemp = boundMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity);
                double lowerBoundTemp = boundsTemp.first;
                double upperBoundTemp = boundsTemp.second;
                double diff = upperBoundTemp - lowerBoundTemp;
                if (verbosity >= 1) {
                    std::cout << "Removed " << i << ", diff: " << diff << std::endl;
                }

                // If it's better, remove it
                if (diff - currentBoundDiff <= tol) {
                    if (verbosity >= 1) {
                        std::cout << "Removing constraint " << i << std::endl;
                    }
                    constraintsZeroCopy = constraintsZero;
                    removedSomething = true;
                    break;
                }

            }

            // If nothing was removed, break
            if (!removedSomething) {
                break;
            }

        }

        // Restore the constraints
        if (verbosity >= 1) {
            std::cout << "Final constraints: " << constraintsZeroCopy.size() << std::endl;
        }
        constraintsZero = constraintsZeroCopy;
        momentMatrices = momentMatricesCopy;

    }

    // If we take fake samples
    if (numSamplesPer != 0 && !precompute) {

        // If simulating samples from symmetries
        if (symSample) {

            // Get the list of Pauli strings that we would actually measure
            std::set<Mon> sampleOperators;
            for (auto& con : zeroConsForSampling) {
                for (auto& term : con) {
                    sampleOperators.insert(term.first);
                }
            }

            // If we're taking samples to a total number of shots
            if (numShots == -1) {
                numSamplesPer = -1;
            } else if (numShots != 0) {
                numSamplesPer = numShots / sampleOperators.size();
            } else {
                numShots = numSamplesPer * sampleOperators.size();
            }

            // Verbose output
            int totalSamples = sampleOperators.size() * numSamplesPer;
            totalSamples = std::abs(totalSamples);
            if (verbosity >= 1) {
                std::cout << "Taking " << numSamplesPer << " samples of " << sampleOperators.size() << " operators" << std::endl;
                std::cout << "Total measurements: " << totalSamples << std::endl;
            }
            if (verbosity >= 3) {
                std::cout << "Sample operators: " << std::endl;
                for (auto mon : sampleOperators) {
                    std::cout << mon << std::endl;
                }
            }

            // Errors
            double K = sampleOperators.size();
            double delta = 1.0 - percentile / 100.0;
            double epsilon = std::sqrt(2 * std::log((2.0 * K) / delta) / numSamplesPer);
            if (verbosity >= 2) {
                std::cout << "K = " << K << std::endl;
                std::cout << "delta = " << delta << std::endl;
                std::cout << "epsilon = " << epsilon << std::endl;
            }

            // Add the constraint that each symmetry should be between -epsilon and epsilon
            if (numSamplesPer == -1) {
                for (auto& con : zeroConsForSampling) {
                    Poly newCon1 = con - Poly(epsilon);
                    constraintsZero.push_back(newCon1);
                }
            } else {
                for (auto& con : zeroConsForSampling) {
                    Poly newCon1 = con + Poly(epsilon);
                    Poly newCon2 = Poly(epsilon) - con;
                    constraintsPositive.push_back(newCon1);
                    constraintsPositive.push_back(newCon2);
                }
            }

        // Otherwise need to take samples
        } else {

            // Pauli matrices
            int matSize = 1 << numQubits;
            Eigen::SparseMatrix<std::complex<double>> X(2, 2);
            Eigen::SparseMatrix<std::complex<double>> Y(2, 2);
            Eigen::SparseMatrix<std::complex<double>> Z(2, 2);
            Eigen::SparseMatrix<std::complex<double>> iden1(1, 1);
            Eigen::SparseMatrix<std::complex<double>> iden2(2, 2);
            X.insert(0, 1) = 1;
            X.insert(1, 0) = 1;
            Y.insert(0, 1) = std::complex<double>(0, -1);
            Y.insert(1, 0) = std::complex<double>(0, 1);
            Z.insert(0, 0) = 1;
            Z.insert(1, 1) = -1;
            iden1.insert(0, 0) = 1;
            iden2.insert(0, 0) = 1;
            iden2.insert(1, 1) = 1;
            X.makeCompressed();
            Y.makeCompressed();
            Z.makeCompressed();
            iden1.makeCompressed();
            iden2.makeCompressed();

            // Declare here, init later only if needed
            Eigen::SparseMatrix<std::complex<double>> groundTruth;
            std::vector<Eigen::SparseMatrix<std::complex<double>>> Xs(numQubits);
            std::vector<Eigen::SparseMatrix<std::complex<double>>> Ys(numQubits);
            std::vector<Eigen::SparseMatrix<std::complex<double>>> Zs(numQubits);

            // Only create these if we need to
            if (!useKnown) {

                // The state we want to find
                groundTruth = Eigen::SparseMatrix<std::complex<double>>(matSize, matSize);

                // Pauli operators on the state
                for (int i = 0; i < numQubits; i++) {
                    Xs[i] = iden1;
                    Ys[i] = iden1;
                    Zs[i] = iden1;
                    for (int j = 0; j < numQubits; j++) {
                        if (j == i) {
                            Xs[i] = kroneckerProduct(Xs[i], X).eval();
                            Ys[i] = kroneckerProduct(Ys[i], Y).eval();
                            Zs[i] = kroneckerProduct(Zs[i], Z).eval();
                        } else {
                            Xs[i] = kroneckerProduct(Xs[i], iden2).eval();
                            Ys[i] = kroneckerProduct(Ys[i], iden2).eval();
                            Zs[i] = kroneckerProduct(Zs[i], iden2).eval();
                        }

                    }
                    Xs[i].makeCompressed();
                    Ys[i].makeCompressed();
                    Zs[i].makeCompressed();
                }
            }

            // If precomputed, load the ground truth
            if (precomputed) {
                if (verbosity >= 1) {
                    std::cout << "Loading precomputed ground truth" << std::endl;
                }
                std::string filename = stateFile;
                std::ifstream inFile(filename);
                if (!inFile.is_open()) {
                    throw std::runtime_error("Could not open ground truth file: " + filename);
                }
                std::string line;
                while (std::getline(inFile, line)) {
                    std::istringstream iss(line);
                    int row, col;
                    double real, imag;
                    if (!(iss >> row >> col >> real >> imag)) {break;}
                    groundTruth.insert(row, col) = std::complex<double>(real, imag);
                }
                if (usePurity) {
                    double purity = 0;
                    for (int k=0; k<matSize; k++) {
                        for (int l=0; l<matSize; l++) {
                            purity += std::norm(groundTruth.coeff(k,l));
                        }
                    }
                    knownIdeal = purity;
                    idealIsKnown = true;
                    if (verbosity >= 1) {
                        std::cout << "Loaded state has purity " << purity << std::endl;
                    }
                }
                if (verbosity >= 3) {
                    std::cout << "Ground truth matrix:" << std::endl;
                    std::cout << Eigen::MatrixXcd(groundTruth) << std::endl;
                }
                inFile.close();

            // Otherwise we need to solve for the ground truth
            } else if (!useKnown) {

                // If it's an energy problem
                if (groundStateProblem) {

                    // Construct the Hamiltonian
                    Poly H = Poly();
                    for (int i=0; i<numQubits; i++) {
                        for (int j=i; j<numQubits; j++) {
                            H += hamiltonianInter[i][j];
                        }
                    }
                    H.reduce();

                    // Form the explicit Hamiltonian
                    Eigen::SparseMatrix<std::complex<double>> HMat(matSize, matSize);
                    for (auto& term : H) {
                        Mon monom = term.first;
                        std::complex<double> coeff = term.second;
                        Eigen::SparseMatrix<std::complex<double>> tempOp(matSize, matSize);
                        for (int j = 0; j < matSize; j++) {
                            tempOp.insert(j, j) = 1;
                        }
                        for (int j = monom.size() - 1; j >= 0; j--) {
                            char letter = monom[j].first;
                            int index = monom[j].second - 1;
                            if (letter == 'X') {
                                tempOp = (Xs[index] * tempOp).eval();
                            } else if (letter == 'Y') {
                                tempOp = (Ys[index] * tempOp).eval();
                            } else if (letter == 'Z') {
                                tempOp = (Zs[index] * tempOp).eval();
                            }
                        }
                        tempOp.makeCompressed();
                        HMat += coeff * tempOp;
                    }

                    // If verbose, output the Hamiltonian
                    if (verbosity >= 3) {
                        std::cout << "Hamiltonian:" << H << std::endl;
                        std::cout << "Hamiltonian matrix:" << std::endl;
                        Eigen::MatrixXcd HMatDense = HMat;
                        std::cout << HMatDense << std::endl;
                    }

                    // Solve the eigenvalue problem
                    Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<std::complex<double>>> es(HMat);
                    Eigen::VectorXcd eigenValues = es.eigenvalues();
                    Eigen::MatrixXcd eigenVectors = es.eigenvectors();
                    std::complex<double> smallestEigenvalue = 100000;
                    int bestInd = -1;
                    for (int i = 0; i < eigenValues.size(); i++) {
                        if (eigenValues[i].real() < smallestEigenvalue.real()) {
                            smallestEigenvalue = eigenValues[i];
                            bestInd = i;
                        }
                    }
                    Eigen::VectorXcd groundTruthVec = eigenVectors.col(bestInd);
                    groundTruth = Eigen::SparseMatrix<std::complex<double>>(matSize, matSize);
                    for (int i = 0; i < matSize; i++) {
                        for (int j = 0; j < matSize; j++) {
                            groundTruth.insert(i, j) = groundTruthVec[i] * std::conj(groundTruthVec[j]);
                        }
                    }
                    groundTruth.makeCompressed();
                    if (verbosity >= 1) {
                        std::cout << "Smallest eigenvalue: " << smallestEigenvalue << std::endl;
                    }
                    if (verbosity >= 3) {
                        std::cout << "Ground truth matrix:" << std::endl;
                        std::cout << Eigen::MatrixXcd(groundTruth) << std::endl;
                    }

                } else {

                    // Save the old problem
                    std::vector<Poly> prevConsZero = constraintsZero;

                    // First we need to solve exactly the problem
                    if (verbosity >= 1) {
                        std::cout << "Genenerating variable set for exact solution" << std::endl;
                    }
                    std::set<Mon> monomsUsed;
                    std::vector<Mon> monomsUsedVec;
                    std::vector<std::vector<std::vector<Poly>>> momentMatrices = {};
                    std::vector<Mon> variables = {};
                    for (int i=0; i<numQubits; i++) {
                        variables.push_back(Mon("<X" + std::to_string(i+1) + ">"));
                        variables.push_back(Mon("<Y" + std::to_string(i+1) + ">"));
                        variables.push_back(Mon("<Z" + std::to_string(i+1) + ">"));
                    }
                    addSingleMonomials(variables, objective);
                    int fullLevel = numQubits;
                    std::vector<Poly> variablesToPut = generateMonomials(variables, fullLevel, verbosity);
                    for (size_t i=0; i<extraMonomialsLim.size(); i++) {
                        variablesToPut.push_back(Poly(extraMonomialsLim[i]));
                    }
                    if (verbosity >= 1) {
                        std::cout << "Generating equation set for exact solution with " << variablesToPut.size() << " variables" << std::endl;
                    }
                    for (size_t i=0; i<variablesToPut.size(); i++) {
                        std::pair<std::complex<double>, Mon> reducedMon = variablesToPut[i].getKey().reduce();
                        if (!monomsUsed.count(reducedMon.second)) {
                            std::pair<char,int> oldMon('A', 0);
                            Mon newPoly(reducedMon.second);
                            Poly newConstraint = lindbladian.replaced(oldMon, newPoly);
                            if (!newConstraint.isZero()) {
                                constraintsZero.push_back(newConstraint);
                            }
                            monomsUsed.insert(variablesToPut[i].getKey());
                            monomsUsedVec.push_back(variablesToPut[i].getKey());
                        }
                    }

                    // Solve
                    std::map<Mon, std::complex<double>> results;
                    if (verbosity >= 1) {
                        std::cout << "Solving exact problem with " << constraintsZero.size() << " constraints" << std::endl;
                    }
                    double res = solveEigen(objective, constraintsZero, 0, numCores, &results);
                    if (verbosity >= 2) {
                        std::cout << "Exact solver result: " << res << std::endl;
                    }
                    if (verbosity >= 3) {
                        for (const auto& [monom, value] : results) {
                            std::cout << monom << " = " << value << std::endl;
                        }
                    }
                    
                    // Restore the old problem
                    constraintsZero = prevConsZero;

                    // Form the ground truth from the results
                    if (verbosity >= 1) {
                        std::cout << "Forming ground truth from results" << std::endl;
                    }
                    for (const auto& [monom, value] : results) {

                        // If it's non-zero
                        if (std::abs(value) > 1e-10) {

                            // Loop over the monomial to form the operator
                            Eigen::SparseMatrix<std::complex<double>> tempOp(matSize, matSize);
                            for (int j = 0; j < matSize; j++) {
                                tempOp.insert(j, j) = 1;
                            }
                            for (int j = monom.size() - 1; j >= 0; j--) {
                                char letter = monom[j].first;
                                int index = monom[j].second - 1;
                                if (letter == 'X') {
                                    tempOp = (Xs[index] * tempOp).eval();
                                } else if (letter == 'Y') {
                                    tempOp = (Ys[index] * tempOp).eval();
                                } else if (letter == 'Z') {
                                    tempOp = (Zs[index] * tempOp).eval();
                                }
                            }
                            tempOp.makeCompressed();

                            // Add the operator to the ground truth
                            groundTruth += value * tempOp;

                        }

                    }

                    // Add the identity and normalize
                    Eigen::SparseMatrix<std::complex<double>> iden(matSize, matSize);
                    for (int i = 0; i < matSize; i++) {
                        iden.insert(i, i) = 1;
                    }
                    groundTruth += iden;
                    groundTruth /= (1 << (numQubits));
                    groundTruth.makeCompressed();

                }

            }

            // The samples we should take
            if (verbosity >= 1) {
                std::cout << "Generating sample operators" << std::endl;
            }
            std::set<Mon> sampleOperators;
            if (sampleChoice == "all") {
                std::vector<Mon> variables = {};
                for (int i=0; i<numQubits; i++) {
                    variables.push_back(Mon("<X" + std::to_string(i+1) + ">"));
                    variables.push_back(Mon("<Y" + std::to_string(i+1) + ">"));
                    variables.push_back(Mon("<Z" + std::to_string(i+1) + ">"));
                }
                std::vector<Poly> variablesToPut = generateMonomials(variables, maxSampleDegree, verbosity);
                for (size_t i=0; i<variablesToPut.size(); i++) {
                    sampleOperators.insert(variablesToPut[i].getKey());
                }

            // If we should sample only the energy
            } else if (sampleChoice == "energy") {

                // Construct the Hamiltonian
                Poly H = Poly();
                for (int i=0; i<numQubits; i++) {
                    for (int j=i; j<numQubits; j++) {
                        H += hamiltonianInter[i][j];
                    }
                }
                H.reduce();

                // Add all of these terms to the sample operators
                for (auto term : H) {
                    Mon mon = term.first;
                    sampleOperators.insert(mon);
                }

            // If we should try to determine the best set
            } else if (sampleChoice == "auto") {

                // Count how many times each Pauli string appears in the constraints
                std::map<Mon, int> monCount;
                std::vector<Mon> monList;
                for (auto term : objective) {
                    if (monCount.count(term.first)) {
                        monCount[term.first]++;
                    } else {
                        monCount[term.first] = 1;
                        monList.push_back(term.first);
                    }
                }
                for (auto con : constraintsZero) {
                    for (auto term : con) {
                        if (monCount.count(term.first)) {
                            monCount[term.first]++;
                        } else {
                            monCount[term.first] = 1;
                            monList.push_back(term.first);
                        }
                    }
                }
                for (auto con : constraintsPositive) {
                    for (auto term : con) {
                        if (monCount.count(term.first)) {
                            monCount[term.first]++;
                        } else {
                            monCount[term.first] = 1;
                            monList.push_back(term.first);
                        }
                    }
                }
                for (auto mat : momentMatrices) {
                    for (int i=0; i<mat.size(); i++) {
                        for (int j=0; j<mat.size(); j++) {
                            for (auto term : mat[i][j]) {
                                if (monCount.count(term.first)) {
                                    monCount[term.first]++;
                                } else {
                                    monCount[term.first] = 1;
                                    monList.push_back(term.first);
                                }
                            }
                        }
                    }
                }
                if (verbosity >= 3) {
                    std::cout << "Monomial counts size: " << monCount.size() << std::endl;
                }

                // Delete the trivial monomial
                monCount.erase(Mon());

                // Sort the above map and take the top N monomials
                if (autoType == "common") { 
                    std::vector<std::pair<Mon, int>> monCountVec(monCount.begin(), monCount.end());
                    std::sort(monCountVec.begin(), monCountVec.end(), [](const std::pair<Mon, int>& a, const std::pair<Mon, int>& b) {
                        return a.second > b.second;
                    });
                    for (int i = 0; i < std::min(maxPaulis, (int)monCountVec.size()); i++) {
                        sampleOperators.insert(monCountVec[i].first);
                    }

                // Otherwise just take the first N monomials
                } else {
                    for (int i = 0; i < std::min(maxPaulis, (int)monList.size()); i++) {
                        sampleOperators.insert(monList[i]);
                    }

                }

                // While we haven't reached the full size, keeping adding products TODO
                while (int(sampleOperators.size()) < maxPaulis) {

                    // Find the next monomial to add
                    Mon nextMon;
                    //for (size_t i=0; i<topRow.size(); i++) {
                        //for (size_t j=i; j<topRow.size(); j++) {
                    for (auto it1 = sampleOperators.begin(); it1 != sampleOperators.end(); ++it1) {
                        for (auto it2 = it1; it2 != sampleOperators.end(); ++it2) {
                            Mon product = (*it1) * (*it2);
                            std::pair<std::complex<double>, Mon> reducedMon = product.reduce();
                            if (!sampleOperators.count(reducedMon.second) && reducedMon.second.size() > 0) {
                                nextMon = reducedMon.second;
                                break;
                            }
                        }
                        if (nextMon.size() > 0) {
                            break;
                        }
                    }

                    // If we found a new monomial, add it
                    if (nextMon.size() > 0) {
                        sampleOperators.insert(nextMon);
                    } else {
                        break;
                    }

                }

            // Just samples from the objective
            } else if (sampleChoice == "onlyobj") {
                for (auto term : objective) {
                    Mon mon = term.first;
                    sampleOperators.insert(mon);
                }

            }

            // Simplfy all, also remove any T monomials
            std::set<Mon> sampleOperatorsNew;
            for (auto mon : sampleOperators) {
                std::pair<std::complex<double>, Mon> reducedMon = mon.reduce();
                if (!reducedMon.second.contains('T')) {
                    sampleOperatorsNew.insert(reducedMon.second);
                }
            }
            sampleOperators = sampleOperatorsNew;

            // If exluding the objective
            if (excludeObjective) {
                for (auto term : objective) {
                    Mon mon = term.first;
                    sampleOperators.erase(mon);
                }
            }

            // If exluding X terms
            if (excludeX) {
                std::set<Mon> newOps;
                for (auto mon : sampleOperators) {
                    if (!mon.contains('X')) {
                        newOps.insert(mon);
                    }
                }
                sampleOperators = newOps;
            }

            // If exluding Y terms
            if (excludeY) {
                std::set<Mon> newOps;
                for (auto mon : sampleOperators) {
                    if (!mon.contains('Y')) {
                        newOps.insert(mon);
                    }
                }
                sampleOperators = newOps;
            }

            // If exluding Z terms
            if (excludeZ) {
                std::set<Mon> newOps;
                for (auto mon : sampleOperators) {
                    if (!mon.contains('Z')) {
                        newOps.insert(mon);
                    }
                }
                sampleOperators = newOps;
            }

            // Remove the trivial monomial
            if (sampleOperators.count(Mon())) {
                sampleOperators.erase(Mon());
            }

            // Make sure sample operators doesn't go above the number of total shots
            while (numShots != 0 && sampleOperators.size() > size_t(numShots)) {
                sampleOperators.erase(sampleOperators.begin());
            }

            // If we're taking samples to a total number of shots
            if (numShots == -1) {
                numSamplesPer = -1;
            } else if (numShots != 0) {
                numSamplesPer = numShots / sampleOperators.size();
            } else {
                numShots = numSamplesPer * sampleOperators.size();
            }

            // Verbose output
            int totalSamples = sampleOperators.size() * numSamplesPer;
            totalSamples = std::abs(totalSamples);
            if (verbosity >= 1) {
                std::cout << "Taking " << numSamplesPer << " samples of " << sampleOperators.size() << " operators" << std::endl;
                std::cout << "Total measurements: " << totalSamples << std::endl;
            }
            if (verbosity >= 3) {
                std::cout << "Sample operators: " << std::endl;
                for (auto mon : sampleOperators) {
                    std::cout << mon << std::endl;
                }
            }

            // Errors
            double K = sampleOperators.size();
            double delta = 1.0 - percentile / 100.0;
            double epsilon = std::sqrt(2 * std::log((2.0 * K) / delta) / numSamplesPer);
            if (verbosity >= 2) {
                std::cout << "K = " << K << std::endl;
                std::cout << "delta = " << delta << std::endl;
                std::cout << "epsilon = " << epsilon << std::endl;
            }

            // Take samples
            samples = {};
            for (auto mon : sampleOperators) {

                // The true value and the probability of getting a 1
                double trueExpectation = 0.0;
                double prob1 = 0.0;

                // If we know the true solution
                if (modelName == "--mg" && useKnown) {
                    bool allNonZero = true;
                    if (mon.size() % 2 != 0) {
                        allNonZero = false;
                    } else {
                        for (int i = 0; i < mon.size(); i+=2) {
                            if (mon[i].first != mon[i+1].first || (mon[i].second % 2) == 0 || std::abs(mon[i].second - mon[i+1].second) != 1) {
                                allNonZero = false;
                                break;
                            }
                            
                        }
                    }
                    if (allNonZero) {
                        trueExpectation = std::pow(-1, mon.size() / 2);
                    }

                // If we don't
                } else {

                    // Construct the operator
                    Eigen::SparseMatrix<std::complex<double>> op = Eigen::SparseMatrix<std::complex<double>>(matSize, matSize);
                    for (int i = 0; i < matSize; i++) {
                        op.insert(i, i) = 1;
                    }
                    for (int i = mon.size()-1; i >= 0; i--) {
                        char pauli = mon[i].first;
                        int ind = mon[i].second-1;
                        if (pauli == 'X') {
                            op = op * Xs[ind];
                        } else if (pauli == 'Y') {
                            op = op * Ys[ind];
                        } else if (pauli == 'Z') {
                            op = op * Zs[ind];
                        }
                    }
                    op.makeCompressed();
                    
                    // Get the true expectation value
                    trueExpectation = Eigen::MatrixXcd(op * groundTruth).trace().real();

                }

                // The probability of measuring 1 from the positive part of the operator
                prob1 = (trueExpectation + 1) / 2.0;

                // If -1 given as the number of samples, use the exact
                if (numSamplesPer == -1) {
                    constraintsZero.push_back(Poly(1, mon) - Poly(trueExpectation));
                    samples[mon] = trueExpectation;
                    if (verbosity >= 3) {
                        std::cout << "True expectation value of " << mon << " = " << trueExpectation << std::endl;
                    }
                    continue;
                }

                // Determine the average from this many samples
                double expFromProb = 2 * prob1 - 1;
                std::binomial_distribution<> binom(numSamplesPer, prob1);
                int success_count = binom(gen);
                double avg = static_cast<double>(success_count) / numSamplesPer;

                // Scale and save this value
                avg = 2 * avg - 1;
                samples[mon] = avg;

                // Bound each quantity
                double lower = avg - epsilon;
                double upper = avg + epsilon;
                if (std::abs(lower) < 1) {
                    constraintsPositive.push_back(Poly(1, mon) - Poly(lower));
                }
                if (std::abs(upper) < 1) {
                    constraintsPositive.push_back(Poly(-1, mon) - Poly(-upper));
                }

                // Verbose output
                if (verbosity >= 3) {
                    std::cout << "True expectation value of " << mon << " = " << trueExpectation << std::endl;
                    std::cout << "Probability of 1 = " << prob1 << std::endl;
                    std::cout << "Expectation from probability = " << expFromProb << std::endl;
                    std::cout << "Expectation value of " << mon << " from samples = " << avg << std::endl;
                    std::cout << "Bounding " << mon << " to [" << lower << ", " << upper << "]" << std::endl;
                }

            }

            // Trivially bound the objective using the above samples
            double minVal = 0;
            double maxVal = 0;
            for (auto term : objective) {
                Mon mon = term.first;
                double coeff = std::real(term.second);
                if (mon.size() == 0) {
                    minVal += coeff;
                    maxVal += coeff;
                } else if (samples.find(mon) != samples.end()) {
                    double sample = samples[mon];
                    if (coeff > 0) {
                        maxVal += coeff * std::min(sample + epsilon, 1.0);
                        minVal += coeff * std::max(sample - epsilon, -1.0);
                    } else {
                        maxVal += coeff * std::max(sample - epsilon, -1.0);
                        minVal += coeff * std::min(sample + epsilon, 1.0);
                    }
                } else {
                    minVal -= std::abs(coeff);
                    maxVal += std::abs(coeff);
                }
            }
            if (verbosity >= 1) {
                std::cout << "Objective bounds from measurement data: [" << minVal << ", " << maxVal << "]" << std::endl;
            }

        }

    }

    // Count the size of the problem
    if (verbosity >= 2) {
        std::cout << std::endl;
    }
    int maxMatSize = 0;
    for (size_t i=0; i<momentMatrices.size(); i++) {
        if (momentMatrices[i].size() > maxMatSize) {
            maxMatSize = momentMatrices[i].size();
        }
    }
    int numCons = constraintsZero.size();
    int numMats = momentMatrices.size();
    std::set<Mon> variableSet;
    variableSet.insert(Mon());
    for (size_t i=0; i<momentMatrices.size(); i++) {
        addVariables(variableSet, momentMatrices[i]);
    }
    for (size_t i=0; i<constraintsZero.size(); i++) {
        addVariables(variableSet, constraintsZero[i]);
    }
    for (size_t i=0; i<constraintsPositive.size(); i++) {
        addVariables(variableSet, constraintsPositive[i]);
    }

    // If using purity as the objective
    std::vector<Mon> quadCone;
    if (usePurity) {

        // Add all the Pauli strings in variableSet
        quadCone.push_back(Mon("<T1>"));
        quadCone.push_back(Mon());
        for (auto mon : variableSet) {
            if (mon != Mon() && !mon.contains('T')) {
                quadCone.push_back(mon);
            }
        }

        // Objective is the purity
        objective = Poly(1, Mon("<T1>"));

    }

    // Also add the variables from the objective
    addVariables(variableSet, objective);
    int numVars = variableSet.size()-1;

    // Output the problem
    if (verbosity >= 2) {
        std::cout << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "            PROBLEM TO SOLVE          " << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << std::endl;
        if (objective.size() > 0) {
            std::cout << "Objective: " << std::endl;
            std::cout << objective << std::endl << std::endl;
            std::cout << "Objective has size " << objective.size() << std::endl;
        }
        if (momentMatrices.size() > 0) {
            for (size_t i=0; i<momentMatrices.size(); i++) {
                std::cout << "Moment matrix " << i  << " (" << momentMatrices[i].size() << "x" << momentMatrices[i].size() << "): " << std::endl;
                if (momentMatrices[i].size() < 20 || verbosity >= 3) {
                    std::cout << momentMatrices[i] << std::endl;
                } else {
                    std::cout << " - Has size " << momentMatrices[i].size() << ", set verbosity 3 to show" << std::endl << std::endl;
                }
            }
        }
        if (constraintsZero.size() > 0) {
            std::cout << "Zero constraints (" << constraintsZero.size() << "): " << std::endl;
            std::cout << constraintsZero << std::endl;
        }
        if (constraintsPositive.size() > 0) {
            std::cout << "Positive constraints (" << constraintsPositive.size() << "): " << std::endl;
            std::cout << constraintsPositive << std::endl;
        }
        if (quadCone.size() > 0) {
            std::cout << "Quadratic cone constraints (" << quadCone.size() << "): " << std::endl;
            for (const auto& mon : quadCone) {
                std::cout << mon << std::endl;
            }
        }
        std::cout << "----------------------------------------" << std::endl;
    }

    // Auto detect the best solver
    if (solver == "auto") {
        if (numCons >= numVars) {
            solver = "eigen";
        } else {
            solver = "mosek";
        }
    }

    // Output what we're doing
    if (verbosity >= 1) {
        std::string optType = "LP";
        if (maxMatSize > 1) {
            optType = "SDP";
        } else if (quadCone.size() > 0) {
            optType = "QP";
        }
        if (solver == "scs") {
            std::cout << "Solving " << optType << " using SCS with " << numCons << " constraints, " << numMats << " moment mats with a max size of " << maxMatSize << " and " << numVars << " variables..." << std::endl;
        } else if (maxMatSize > 1 && solver == "mosek") {
            std::cout << "Solving " << optType << " using MOSEK with " << numCons << " constraints, " << numMats << " moment mats with a max size of " << maxMatSize << " and " << numVars << " variables..." << std::endl;
        } else if (solver == "eigen") {
            std::cout << "Solving " << optType << " using Eigen with " << numCons << " constraints and " << numVars << " variables..." << std::endl;
        } else if (solver == "mosek") {
            std::cout << "Solving " << optType << " using MOSEK with " << numCons << " constraints and " << numVars << " variables (ratio: " << double(numCons)/numVars << ")..." << std::endl;
        } else if (solver == "gurobi") {
            std::cout << "Solving " << optType << " using Gurobi with " << numCons << " constraints and " << numVars << " variables (ratio: " << double(numCons)/numVars << ")..." << std::endl;
        }
    }

    // Compute trivial bounds
    double trivialMax = 0;
    double trivialMin = 0;
    for (auto& term : objective) {
        double coeff = std::real(term.second);
        if (term.first.size() == 0) {
            trivialMax += coeff;
            trivialMin += coeff;
        } else {
            if (coeff > 0) {
                trivialMax += coeff;
                trivialMin -= coeff;
            } else {
                trivialMax -= coeff;
                trivialMin += coeff;
            }
        }
    }
    if (usePurity) {
        trivialMin = 0.0;
        trivialMax = 1.0;
    }
    if (verbosity >= 1) {
        std::cout << "Trivial bounds: " << trivialMin << " <= obj <= " << trivialMax << std::endl;
    }

    // Solve the relaxation
    std::pair<double,double> bounds = {-1000000, 1000000};
    std::map<Mon, std::complex<double>> results;
    if (!(precompute && groundStateProblem)) {
        if (solver == "scs") {
            bounds = boundSCS(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, {-1,1}, &results);
        } else if (solver == "mosek" || maxMatSize > 1) {
            bounds = boundMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, {-1,1}, imagType, &results, quadCone);
        } else if (solver == "gurobi") {
            bounds = boundGurobi(objective, constraintsZero, verbosity, &results);
        } else if (solver == "eigen") {
            double res = solveEigen(objective, constraintsZero, verbosity, numCores, &results);
            bounds = {res, res};
        }
    }
    double lowerBound = bounds.first;
    double upperBound = bounds.second;
    if (usePurity) {
        lowerBound = lowerBound / std::pow(2, numQubits);
        upperBound = 0.0;
        double absSum = 0.0;
        for (const auto& val : results) {
            if (!val.first.contains('T')) {
                upperBound += std::norm(val.second);
                absSum += std::abs(val.second);
            }
        }
        if (verbosity >= 1) {
            std::cout << "sum of pauli squares = " << upperBound << std::endl;
        }
        upperBound /= std::pow(2, numQubits);
        if (verbosity >= 1) {
            std::cout << "T = " << results[Mon("<T1>")] << std::endl;
        }
        std::swap(lowerBound, upperBound);
    }
    double diff = std::abs(upperBound - lowerBound);
    double error = (diff / std::abs(trivialMax - trivialMin)) * 100;

    // Output the results
    if (idealIsKnown && verbosity >= 1) {
        std::cout << "Known True Optimum: " << knownIdeal << std::endl;
        std::cout << "Relative Error: " << diff / std::abs(knownIdeal) * 100 << "%" << std::endl;
    }
    if (verbosity >= 1) {
        std::cout << "Bounds: " << lowerBound << "  <  " << upperBound << std::endl;
        std::cout << "Difference: " << upperBound - lowerBound << std::endl;
        std::cout << "As Percentage: " << strRound(error, 2) << "%" << std::endl;
    }

    // Timing
    std::chrono::steady_clock::time_point timeFinishedSolving = std::chrono::steady_clock::now();

    // If told to output final moments to file
    if (outputToFile && results.size() > 0) {
        std::string filename = "monoms.csv";
        std::ofstream outFile(filename);
        if (outFile.is_open()) {
            for (const auto& pair : results) {
                outFile << pair.first.size() << "," << pair.first << "," << pair.second.real() << "," << pair.second.imag() << "\n";
            }
            outFile.close();
        } else {
            std::cerr << "Error opening file: " << filename << std::endl;
        }
    }

    // Output the timings
    if (verbosity >= 1) { 
        std::cout << "Time to generate: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedGeneration - timeFinishedArgs).count() / 1000.0 << "s" << std::endl;
        std::cout << "Time to solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedSolving - timeFinishedGeneration).count() / 1000.0 << "s" << std::endl;
    }

    // If told to check the final results by reconstructing the full matrix
    // (also used when precompute is called)
    if (checkObj) {

        // Reconstruct the full matrix
        int fullMatSize = 1 << numQubits;
        Eigen::MatrixXcd reconMatrix = Eigen::MatrixXcd::Zero(fullMatSize, fullMatSize);

        // If precomputing a Hamiltonian we should diagonalize it here
        if (precompute && groundStateProblem) {

            // Construct the Hamiltonian
            Poly H = Poly();
            for (int i=0; i<numQubits; i++) {
                for (int j=i; j<numQubits; j++) {
                    H += hamiltonianInter[i][j];
                }
            }
            H.reduce();

            // Form the explicit Hamiltonian
            Eigen::SparseMatrix<std::complex<double>> HMat(fullMatSize, fullMatSize);
            for (auto& term : H) {
                Mon monom = term.first;
                std::complex<double> coeff = term.second;

                // Turn X1Y2I3 -> {1, 2, 0}
                std::vector<int> inds(numQubits, 0);
                for (int i = 0; i < monom.size(); i++) {
                    char letter = monom[i].first;
                    int index = monom[i].second - 1;
                    if (letter == 'X') {
                        inds[index] = 1;
                    } else if (letter == 'Y') {
                        inds[index] = 2;
                    } else if (letter == 'Z') {
                        inds[index] = 3;
                    }
                }

                // Verbose output
                if (verbosity >= 3) {
                    std::cout << monom << " = " << coeff  << " (";
                    for (int i = 0; i < numQubits; i++) {
                        std::cout << inds[i];
                        if (i < numQubits - 1) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << ")" << std::endl;
                }

                // The corresponding matrix
                Eigen::SparseMatrix<std::complex<double>> mat = generatePauliMatrix(inds);

                // Add it, scaled by the value
                HMat += coeff * mat;

            }
            HMat.makeCompressed();

            // If verbose, output the Hamiltonian
            if (verbosity >= 3) {
                std::cout << "Hamiltonian:" << H << std::endl;
                std::cout << "Hamiltonian matrix:" << std::endl;
                Eigen::MatrixXcd HMatDense = HMat;
                std::cout << HMatDense << std::endl;
            }

            // Solve the eigenvalue problem
            if (verbosity >= 1) {
                std::cout << "Diagonalizing the Hamiltonian (" << fullMatSize << "x" << fullMatSize << ")..." << std::endl;
            }
            omp_set_num_threads(numCores);
            Eigen::setNbThreads(numCores);
            Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<std::complex<double>>> es(HMat);
            Eigen::VectorXcd eigenValues = es.eigenvalues();
            Eigen::MatrixXcd eigenVectors = es.eigenvectors();

            // Get the smallest eigenvalue and its corresponding eigenvector
            std::complex<double> smallestEigenvalue = 100000;
            int bestInd = -1;
            for (int i = 0; i < eigenValues.size(); i++) {
                if (eigenValues[i].real() < smallestEigenvalue.real()) {
                    smallestEigenvalue = eigenValues[i];
                    bestInd = i;
                }
            }
            lowerBound = smallestEigenvalue.real();
            upperBound = lowerBound;
            bounds = {lowerBound, upperBound};
            diff = 0;
            error = 0;
            if (verbosity >= 1) {
                std::cout << "Minimum energy of H: " << lowerBound << std::endl;
                if (modelName == "--mg") {
                    knownIdeal = -(3.0/8.0)*numQubits;
                    idealIsKnown = true;
                    std::cout << "Known minimum: " << knownIdeal << std::endl;
                }
            }

            // Form the full matrix
            if (verbosity >= 1) {
                std::cout << "Reconstructing the full matrix..." << std::endl;
            }
            Eigen::VectorXcd groundTruthVec = eigenVectors.col(bestInd);
            reconMatrix = groundTruthVec * groundTruthVec.adjoint();

        } else {

            // Make sure we have enough vars to reconstruct the matrix
            int numNeeded = std::pow(4, numQubits) - 1;
            if (results.size() < numNeeded) {
                std::cout << "Warning: not enough variables to reconstruct the full matrix (" << results.size() << " < " << numNeeded << ")" << std::endl;
            }
            if (results.find(Mon()) == results.end()) {
                results[Mon()] = std::complex<double>(1.0, 0.0);
            }

            // For each Pauli string in the results
            for (auto& [monom, value] : results) {

                // If the monomial is not trivial
                if (!monom.contains('T')) {

                    // Turn X1Y2I3 -> {1, 2, 0}
                    std::vector<int> inds(numQubits, 0);
                    for (int i = 0; i < monom.size(); i++) {
                        char letter = monom[i].first;
                        int index = monom[i].second - 1;
                        if (letter == 'X') {
                            inds[index] = 1;
                        } else if (letter == 'Y') {
                            inds[index] = 2;
                        } else if (letter == 'Z') {
                            inds[index] = 3;
                        }
                    }

                    // Verbose output
                    if (verbosity >= 3) {
                        std::cout << monom << " = " << value  << " (";
                        for (int i = 0; i < numQubits; i++) {
                            std::cout << inds[i];
                            if (i < numQubits - 1) {
                                std::cout << ", ";
                            }
                        }
                        std::cout << ")" << std::endl;
                    }

                    // The corresponding matrix
                    Eigen::SparseMatrix<std::complex<double>> mat = generatePauliMatrix(inds);

                    // Add it, scaled by the value
                    reconMatrix += value * mat;

                }

            }

            // Normalize
            reconMatrix /= fullMatSize;

        }

        // Output the full matrix
        if (verbosity >= 3) {
            std::cout << "Full matrix:" << std::endl;
            std::cout << Eigen::MatrixXcd(reconMatrix) << std::endl;
        }

        // Ensure that the result is positive and has trace 1
        if (verbosity >= 1) {
            std::cout << "Checking the eigenspectra of the full density matrix..." << std::endl;
        }
        omp_set_num_threads(numCores);
        Eigen::setNbThreads(numCores);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(reconMatrix);
        Eigen::VectorXcd eigenValues = es.eigenvalues();
        std::complex<double> smallestEigenvalue = eigenValues(0).real();
        std::complex<double> matTrace = 0;
        for (int i = 0; i < fullMatSize; i++) {
            matTrace += reconMatrix(i, i);
        }
        if (verbosity >= 1) {
            std::cout << "Smallest eigenvalue: " << eigenValues(0).real() << std::endl;
            std::cout << "Trace: " << matTrace.real() << std::endl;
        }

        // Check the purity
        Eigen::MatrixXcd reconMatrixSquared = reconMatrix * reconMatrix.adjoint();
        double purity = 0;
        for (int i = 0; i < fullMatSize; i++) {
            purity += reconMatrixSquared(i, i).real();
        }
        if (verbosity >= 1) {
            std::cout << "Purity: " << purity << std::endl;
            std::cout << "Minimum possible purity: " << 1.0 / fullMatSize << std::endl;
        }

        // If precomputing, save it to a file
        if (precompute) {
            if (verbosity >= 1) {
                std::cout << "Saving the optimum state to file..." << std::endl;
            }
            std::string filename = stateFile;
            std::ofstream outFile(filename);
            outFile.precision(15);
            if (outFile.is_open()) {
                for (int i = 0; i < fullMatSize; i++) {
                    for (int j = 0; j < fullMatSize; j++) {
                        outFile << i << " " << j << " " << reconMatrix(i, j).real() << " " << reconMatrix(i, j).imag() << "\n";
                    }
                }
                outFile.close();
                if (verbosity >= 1) {
                    std::cout << "Optimum state saved to: " << filename << std::endl;
                }
            } else {
                std::cerr << "Error opening file for writing: " << filename << std::endl;
            }
            return 0;
        }

    }

    // Benchmarking output
    if (benchmark) {
        std::string outputMat = "None";
        std::string outputLin = "None ";
        std::string outputRecon = "None";
        std::string outputSym = "None";
        std::string outputDiff = "None";
        std::string outputTime = "None";
        std::string outputShots = "None";
        std::string outputNote = "";
        if (level > 0) {
            std::string matWidth = std::to_string(momentMatrices[0].size());
            outputMat = "level " + std::to_string(level) + " (" + matWidth + "x" + matWidth + ")";
        }
        if (autoMomentAmount) {
            int width = autoMomentAmount+1;
            outputMat = "auto (" + std::to_string(width) + "x" + std::to_string(width) + ")";
        }
        if (lindbladLevel > 0) {
            outputLin = "level " + std::to_string(lindbladLevel) + " (" + std::to_string(constraintsZero.size()-numSyms) + ")";
        } else if (findMinimal && findMinimalAmount >= 0) {
            int width = constraintsZero.size() - numSyms;
            outputLin = "auto (" + std::to_string(width) + ")";
        }
        if (numSamplesPer != 0) {
            outputShots = "";
            outputShots += sampleChoice;
            if (sampleChoice == "all") {
                outputShots += " order-" + std::to_string(maxSampleDegree);
            }
            outputShots += " (" + std::to_string(numShots) + ")";
        }
        if (reconLevel > 0) {
            std::string numMats = std::to_string(momentMatrices.size());
            std::string matSize = std::to_string(momentMatrices[momentMatrices.size()-1].size());
            outputRecon = "all " + std::to_string(reconLevel) + "-site (" + numMats + "x" + matSize + "x" + matSize + ")";
        }
        if (symmetries.size() > 0) {
            outputSym = "yes (" + std::to_string(symmetries.size()) + ")";
        }
        outputDiff = "[" + strRound(lowerBound, precision) + ", " + strRound(upperBound, precision) + "] (" + strRound(error, 2) + "\\%)";
        int timeInMillis = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedSolving - timeFinishedArgs).count();
        if (allTimeInMillis) {
            outputTime = std::to_string(timeInMillis) + "ms";
        } else {
            int minutes = timeInMillis / 60000;
            int seconds = (timeInMillis % 60000) / 1000;
            if (minutes > 0) {
                outputTime = std::to_string(minutes) + "m" + std::to_string(seconds) + "s";
            } else {
                outputTime = strRound(timeInMillis / 1000.0, 2) + "s";
            }
        }
        if (note.size() > 0) {
            outputNote = " & " + note;
        }
        if (numSamplesPer != 0) {
            std::cout << outputMat << " & " << outputLin << " & " << outputRecon << " & " << outputSym << " & " << outputShots << " & " << outputDiff << " & " << outputTime << outputNote << " \\\\ \\hline" << std::endl;
        } else {
            std::cout << outputMat << " & " << outputLin << " & " << outputRecon << " & " << outputSym << " & " << outputDiff << " & " << outputTime << outputNote << " \\\\ \\hline" << std::endl;
        }
    }

    // Exit without errors
    return 0;

}
