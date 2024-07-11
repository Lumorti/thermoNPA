// Standard includes
#include <iostream>
#include <vector>
#include <complex>
#include <set>
#include <chrono>
#include <map>

// Import OpenMP
#include <omp.h>

// Import Eigen
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

// Local files
#include "poly.h"
#include "printing.h"
#include "utils.h"
#include "optMOSEK.h"
#include "optGurobi.h"
#include "optSCS.h"
 
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

// Useful constants
const std::complex<double> imag(0, 1);

// Generic entry function
int main(int argc, char* argv[]) {

    // Random seed unless specified
    srand(time(NULL));

    // Define the scenario
    int level = 0;
    int limbladLevel = 0;
    int numQubits = 1;
    int gridWidth = 1;
    int gridHeight = 1;
    double tol = 1e-7;
    Poly objective("<Z1>");
    int verbosity = 1;
    std::string solver = "auto";
    std::complex<double> knownIdeal = 0.0;
    bool idealIsKnown = false;
    bool findMinimal = false;
    bool tryRemove = false;
    int autoMomentAmount = 0;
    int findMinimalAmount = 0;
    int reductiveCons = 0;
    std::string seed = "";
    Poly limbladian("<X1A0X1>+<Y1A0Y1>-<A0>");
    std::vector<std::string> extraMonomials;
    std::vector<std::string> extraMonomialsLim;
    std::vector<Poly> constraintsZero;
    std::vector<std::pair<std::vector<int>,std::vector<int>>> symmetries;

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

        // Manually set the Limbladian
        } else if (argAsString == "-L") {
            limbladian = Poly(std::string(argv[i+1]));
            i++;

        // Pauli X objective
        } else if (argAsString == "--objX") {
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<X" + std::to_string(i) + ">");
            }

        // Pauli Y objective
        } else if (argAsString == "--objY") {
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Y" + std::to_string(i) + ">");
            }

        // Pauli Z objective
        } else if (argAsString == "--objZ") {
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }

        // Setting the tolerance
        } else if (argAsString == "-t") {
            tol = std::stod(argv[i+1]);
            i++;

        // Pauli Limbladian
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
            limbladian = Poly(pauliString);
            i += 3;

        // Two-body Limbladian
        } else if (argAsString == "--two" || argAsString == "--twov") {

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
            limbladian = Poly();
            limbladian += Poly(-epsilon_h*imag, "<A0P1M1>");
            limbladian += Poly(-epsilon_c*imag, "<A0P2M2>");
            limbladian += Poly(-g*imag, "<A0P1M2>");
            limbladian += Poly(-g*imag, "<A0M1P2>");
            limbladian += Poly(epsilon_h*imag, "<P1M1A0>");
            limbladian += Poly(epsilon_c*imag, "<P2M2A0>");
            limbladian += Poly(g*imag, "<P1M2A0>");
            limbladian += Poly(g*imag, "<M1P2A0>");
            limbladian += Poly(gamma_h_plus, "<M1A0P1>");
            limbladian += Poly(-0.5*gamma_h_plus, "<A0M1P1>");
            limbladian += Poly(-0.5*gamma_h_plus, "<M1P1A0>");
            limbladian += Poly(gamma_h_minus, "<P1A0M1>");
            limbladian += Poly(-0.5*gamma_h_minus, "<A0P1M1>");
            limbladian += Poly(-0.5*gamma_h_minus, "<P1M1A0>");
            limbladian += Poly(gamma_c_plus, "<M2A0P2>");
            limbladian += Poly(-0.5*gamma_c_plus, "<A0M2P2>");
            limbladian += Poly(-0.5*gamma_c_plus, "<M2P2A0>");
            limbladian += Poly(gamma_c_minus, "<P2A0M2>");
            limbladian += Poly(-0.5*gamma_c_minus, "<A0P2M2>");
            limbladian += Poly(-0.5*gamma_c_minus, "<P2M2A0>");
            limbladian.convertToPaulis();

        // 2D Limbladian test
        } else if (argAsString == "--2d") {

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
            std::vector<double> epsilons(numQubits, epsilon_other);
            double n_c = 1.0 / (std::exp(epsilon_c / T_c) - 1.0);
            double n_h = 1.0 / (std::exp(epsilon_h / T_h) - 1.0);
            double gamma_h_plus = gamma_h * n_h;
            double gamma_h_minus = gamma_h * (n_h + 1.0);
            double gamma_c_plus = gamma_c * n_c;
            double gamma_c_minus = gamma_c * (n_c + 1.0);

            // Connectivity: 2D grid - nearest neighbours
            std::vector<double> gamma_plus(numQubits, 0);
            std::vector<double> gamma_minus(numQubits, 0);
            std::vector<std::vector<double>> gs(numQubits, std::vector<double>(numQubits, 0));
            for (int i=0; i<numQubits; i++) {

                // Get the x and y location
                int xLoc = i % gridWidth;
                int yLoc = i / gridWidth;

                // The leftmost qubits are connected to the hot bath
                //if (xLoc == 0) {
                if (xLoc == 0 && yLoc == std::floor(gridHeight/2)) {
                    gamma_plus[i] = gamma_h_plus;
                    gamma_minus[i] = gamma_h_minus;
                    epsilons[i] = epsilon_h;

                // The rightmost qubits are connected to the cold bath
                //} else if (xLoc == gridWidth-1) {
                } else if (xLoc == gridWidth-1 && yLoc == std::floor(gridHeight/2)) {
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
                
            }
            
            // Construct the objective as a polynomial
            objective = Poly();
            for (int i=1; i<=numQubits; i++) {
                objective += Poly(1.0/numQubits, "<Z" + std::to_string(i) + ">");
            }
            
            // Construct the Limbadlian as a polynomial from plus/minus
            limbladian = Poly();
            for (int i=1; i<=numQubits; i++) {
                limbladian += Poly(-imag*epsilons[i-1], "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                limbladian += Poly(imag*epsilons[i-1], "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
            }
            for (int i=1; i<=numQubits; i++) {
                for (int j=i+1; j<=numQubits; j++) {
                    double coeff = gs[i-1][j-1];
                    limbladian += Poly(-imag*coeff, "<A0P" + std::to_string(i) + "M" + std::to_string(j) + ">");
                    limbladian += Poly(-imag*coeff, "<A0M" + std::to_string(i) + "P" + std::to_string(j) + ">");
                    limbladian += Poly(imag*coeff, "<P" + std::to_string(i) + "M" + std::to_string(j) + "A0>");
                    limbladian += Poly(imag*coeff, "<M" + std::to_string(i) + "P" + std::to_string(j) + "A0>");
                }
            }
            for (int i=1; i<=numQubits; i++) {
                double coeffPlus = gamma_plus[i-1];
                double coeffMinus = gamma_minus[i-1];
                limbladian += Poly(coeffPlus, "<M" + std::to_string(i) + "A0P" + std::to_string(i) + ">");
                limbladian += Poly(-0.5*coeffPlus, "<A0M" + std::to_string(i) + "P" + std::to_string(i) + ">");
                limbladian += Poly(-0.5*coeffPlus, "<M" + std::to_string(i) + "P" + std::to_string(i) + "A0>");
                limbladian += Poly(coeffMinus, "<P" + std::to_string(i) + "A0M" + std::to_string(i) + ">");
                limbladian += Poly(-0.5*coeffMinus, "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                limbladian += Poly(-0.5*coeffMinus, "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
            }
            limbladian.clean();
            limbladian.convertToPaulis();

        // The Limbladian from David
        } else if (argAsString == "--david" || argAsString == "--davidr") {

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

            // The full Limbladian
            // -i[H, rho] + \sum_k 0.5 * (2*gamma_h * Gamma_k rho Gamma_k^dagger - Gamma_k^dagger Gamma_k rho - rho Gamma_k^dagger Gamma_k)
            Poly rho("<R1>");
            limbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                if (i == 0 || i == numQubits-1) {
                    limbladian += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            limbladian = Poly("<A0>") * limbladian;
            limbladian.cycleToAndRemove('R', 1);
            limbladian.reduce();

        // The Limbladian from David but 2D
        } else if (argAsString == "--david2d" || argAsString == "--david2dr") {

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
                int thisInd = xLoc + yLoc*gridWidth + 1;

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
                if (xLoc == 0 && xLoc == gridWidth-1) {
                    Gamma_k[i-1] *= (std::sqrt(gamma_h) + std::sqrt(gamma_c));
                } else if (xLoc == 0) {
                    Gamma_k[i-1] *= std::sqrt(gamma_h);
                } else if (xLoc == gridWidth-1) {
                    Gamma_k[i-1] *= std::sqrt(gamma_c);
                }
            }

            // The full Limbladian
            // -i[H, rho] + \sum_k 0.5 * (2*gamma_h * Gamma_k rho Gamma_k^dagger - Gamma_k^dagger Gamma_k rho - rho Gamma_k^dagger Gamma_k)
            Poly rho("<R1>");
            limbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                int xLoc = i % gridWidth;
                int yLoc = i / gridWidth;
                if (xLoc == 0 || xLoc == gridWidth-1) {
                    if (verbosity >= 2) {
                        std::cout << "Connecting " << (i+1) << " to the bath" << std::endl;
                    }
                    limbladian += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            limbladian = Poly("<A0>") * limbladian;
            limbladian.cycleToAndRemove('R', 1);
            limbladian.reduce();

        // A 1.5D chain TODO
        } else if (argAsString == "--1.5d" || argAsString == "--1.5dr") {

            // Parameters
            double g = std::stod(argv[i+1]);
            numQubits = std::stoi(argv[i+2]);
            i+=2;
            double gamma_h = 1.0;
            double gamma_c = 0.5;

            // Make sure it's odd
            if (numQubits % 2 == 0) {
                std::cout << "Error - Number of qubits must be odd" << std::endl;
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

            // The full Limbladian
            // -i[H, rho] + \sum_k 0.5 * (2*gamma_h * Gamma_k rho Gamma_k^dagger - Gamma_k^dagger Gamma_k rho - rho Gamma_k^dagger Gamma_k)
            Poly rho("<R1>");
            limbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                if (i == 0 || i == numQubits-1) {
                    limbladian += 0.5 * (2 * Gamma_k[i] * rho * Gamma_k[i].dagger() - Gamma_k[i].dagger() * Gamma_k[i] * rho - rho * Gamma_k[i].dagger() * Gamma_k[i]);
                }
            }
            limbladian = Poly("<A0>") * limbladian;
            limbladian.cycleToAndRemove('R', 1);
            limbladian.reduce();

        // The Limbladian from the tensor paper
        } else if (argAsString == "--tensor") {

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

            // The full Limbladian
            // -i[H, rho] + \sum_k Gamma_k rho Gamma_k^dagger - 0.5 {Gamma_k^dagger Gamma_k, rho}
            Poly rho("<R1>");
            limbladian = -imag*H.commutator(rho);
            for (int i=0; i<numQubits; i++) {
                limbladian += Gamma_k[i] * rho * Gamma_k[i].dagger();
                limbladian -= 0.5 * (Gamma_k[i].dagger() * Gamma_k[i]).anticommutator(rho);
            }
            limbladian = Poly("<A0>") * limbladian;
            limbladian.cycleToAndRemove('R', 1);
            limbladian.reduce();

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
                std::cout << "Limbladian: " << limbladian << std::endl;
                std::cout << "Objective: " << objective << std::endl;
            }

        // Many-body Limbladian
        } else if (argAsString == "--many" || argAsString == "--manyv" || argAsString == "--manyr") {

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
            limbladian = Poly();
            for (int i=1; i<=numQubits; i++) {
                limbladian += Poly(-imag*epsilons[i-1], "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                limbladian += Poly(imag*epsilons[i-1], "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
            }
            for (int i=1; i<=numQubits-1; i++) {
                limbladian += Poly(-imag*gs[i], "<A0P" + std::to_string(i) + "M" + std::to_string(i+1) + ">");
                limbladian += Poly(-imag*gs[i], "<A0M" + std::to_string(i) + "P" + std::to_string(i+1) + ">");
                limbladian += Poly(imag*gs[i], "<P" + std::to_string(i) + "M" + std::to_string(i+1) + "A0>");
                limbladian += Poly(imag*gs[i], "<M" + std::to_string(i) + "P" + std::to_string(i+1) + "A0>");
            }
            for (int i : {1, numQubits}) {
                limbladian += Poly(gamma_plus[i-1], "<M" + std::to_string(i) + "A0P" + std::to_string(i) + ">");
                limbladian += Poly(-0.5*gamma_plus[i-1], "<A0M" + std::to_string(i) + "P" + std::to_string(i) + ">");
                limbladian += Poly(-0.5*gamma_plus[i-1], "<M" + std::to_string(i) + "P" + std::to_string(i) + "A0>");
                limbladian += Poly(gamma_minus[i-1], "<P" + std::to_string(i) + "A0M" + std::to_string(i) + ">");
                limbladian += Poly(-0.5*gamma_minus[i-1], "<A0P" + std::to_string(i) + "M" + std::to_string(i) + ">");
                limbladian += Poly(-0.5*gamma_minus[i-1], "<P" + std::to_string(i) + "M" + std::to_string(i) + "A0>");
            }
            limbladian.convertToPaulis();

        // Set the seed
        } else if (argAsString == "-S") {
            seed = std::string(argv[i+1]);
            srand(std::hash<std::string>{}(seed));
            i++;

        // If setting verbosity
        } else if (argAsString == "-v") {
            verbosity = std::stoi(argv[i+1]);
            i++;

        // If setting the level of the Limbladian
        } else if (argAsString == "-l") {
            limbladLevel = std::stoi(argv[i+1]);
            i++;

        // If adding an extra monomial to the top row
        } else if (argAsString == "-e") {
            extraMonomials.push_back(std::string(argv[i+1]));
            i++;

        // If adding an extra monomial to the list of Limbladian replacements
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

        // Output the help
        } else if (argAsString == "-h" || argAsString == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -m <int>        Level of the moment matrix" << std::endl;
            std::cout << "  -l <int>        Level of the moments to put in the Limbladian" << std::endl;
            std::cout << "  -e <mon>        Add an extra monomial to the top row of the moment matrix" << std::endl;
            std::cout << "  -E <mon>        Add an extra monomial to the list of Limbladian replacements" << std::endl;
            std::cout << "  -S <str>        Seed for the random number generator" << std::endl;
            std::cout << "  -M <int>        Try to generate the minimal set of linear constraints" << std::endl;
            std::cout << "  -A <int>        Try to generate the minimal moment matrix" << std::endl;
            std::cout << "  -R              Try removing random constraints" << std::endl;
            std::cout << "  -v <int>        Verbosity level" << std::endl;
            std::cout << "  -C <int>        Number of cores to use" << std::endl;
            std::cout << "  -c <int>        Add some number of extra constraints to reduce the num of vars" << std::endl;
            std::cout << "  -t <dbl>        Set the tolerance (used in various places)" << std::endl;
            std::cout << "  -y <ints>       Add a symetry between two groups e.g. -y 1,2 3,4" << std::endl;
            std::cout << "Solver options:" << std::endl;
            std::cout << "  -s G            Use Gurobi as the solver" << std::endl;
            std::cout << "  -s E            Use Eigen as the solver" << std::endl;
            std::cout << "  -s M            Use MOSEK as the solver" << std::endl;
            std::cout << "  -s S            Use SCS as the solver" << std::endl;
            std::cout << "  -s N            Don't solve after generating" << std::endl;
            std::cout << "Objective options:" << std::endl;
            std::cout << "  -O --obj <str>  Manually set the objective" << std::endl;
            std::cout << "  --objX          Use avg sigma_X as the objective" << std::endl;
            std::cout << "  --objY          Use avg sigma_Y as the objective" << std::endl;
            std::cout << "  --objZ          Use avg sigma_Z as the objective" << std::endl;
            std::cout << "Limbliadian options:" << std::endl;
            std::cout << "  -L <str>        Manually set the Limbladian" << std::endl;
            std::cout << "  --pauli <dbl> <dbl> <dbl>" << std::endl;
            std::cout << "  --two" << std::endl;
            std::cout << "  --twov <dbl> <dbl>" << std::endl;
            std::cout << "  --many  <int>" << std::endl;
            std::cout << "  --manyr <int>" << std::endl;
            std::cout << "  --manyv <int> <dbl> <dbl>" << std::endl;
            std::cout << "  --2d <int> <int>" << std::endl;
            std::cout << "  --tensor <int>" << std::endl;
            std::cout << "  --david <dbl> <int>" << std::endl;
            std::cout << "  --davidr <dbl> <int>" << std::endl;
            std::cout << "  --david2d <dbl> <int> <int>" << std::endl;
            std::cout << "  --david2dr <dbl> <int> <int>" << std::endl;
            return 0;

        // If auto generating the moment matrix
        } else if (argAsString == "-A") {
            autoMomentAmount = std::stoi(argv[i+1]);
            i++;

        // Otherwise we don't know what this is
        } else if (argAsString != "./run") {
            std::cout << "Unknown argument: " << argAsString << std::endl;
            return 1;

        }
    }

    // Start a timer
    std::chrono::steady_clock::time_point timeFinishedArgs = std::chrono::steady_clock::now();

    // Create the Limbladian applied to many different operators
    std::vector<std::vector<std::vector<Poly>>> momentMatrices = {};
    std::vector<Mon> variables = {};
    for (int i=0; i<numQubits; i++) {
        variables.push_back(Mon("<X" + std::to_string(i+1) + ">"));
        variables.push_back(Mon("<Y" + std::to_string(i+1) + ">"));
        variables.push_back(Mon("<Z" + std::to_string(i+1) + ">"));
    }
    std::vector<Poly> variablesToPut = generateMonomials(variables, limbladLevel, verbosity);
    for (size_t i=0; i<extraMonomialsLim.size(); i++) {
        variablesToPut.push_back(Poly(extraMonomialsLim[i]));
    }
    if (verbosity >= 2) {
        std::cout << std::endl;
        std::cout << "Variables to put in Limbladian: " << std::endl;
        for (size_t i=0; i<variablesToPut.size(); i++) {
            std::cout << variablesToPut[i] << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Original Limbladian: " << limbladian << std::endl;
    }
    if (!findMinimal) {
        for (size_t i=0; i<variablesToPut.size(); i++) {
            std::pair<char,int> oldMon('A', 0);
            Mon newPoly(variablesToPut[i].getKey());
            Poly newConstraint = limbladian.replaced(oldMon, newPoly);
            if (verbosity >= 3) {
                std::cout << std::endl;
                std::cout << "Variable to put: " << newPoly << std::endl;
                std::cout << "New constraint: " << newConstraint << std::endl;
            }
            if (!newConstraint.isZero()) {
                constraintsZero.push_back(newConstraint);
            }
        }
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
    std::set<Mon> monomsUsed;
    if (findMinimal) {

        // If given zero, set to what seems to be the minimum needed
        if (findMinimalAmount == 0) {
            if (gridHeight == 1) {
                findMinimalAmount = 2*numQubits*numQubits - numQubits;
            } else {
                findMinimalAmount = std::pow(4, numQubits)/2 - 1;
            }
            std::cout << "Auto setting constraint limit to: " << findMinimalAmount << std::endl;
        }

        // If a max num of constraints is given
        if (findMinimalAmount > 0) {

            // Add constraints based on the monomials we already have
            std::set<Mon> monomsInConstraints;
            std::vector<Mon> queue;
            queue.push_back(Mon());
            monomsInConstraints.insert(Mon());
            for (auto& term : objective) {
                if (!monomsInConstraints.count(term.first)) {
                    monomsInConstraints.insert(term.first);
                    queue.push_back(term.first);
                }
            }
            for (size_t i=0; i<variablesToPut.size(); i++) {
                Mon monToAdd = variablesToPut[i].getKey();
                if (!monomsInConstraints.count(monToAdd)) {
                    monomsInConstraints.insert(monToAdd);
                    queue.push_back(monToAdd);
                }
            }
            for (size_t i=0; i<momentMatrices.size(); i++) {
                for (size_t j=0; j<momentMatrices[i].size(); j++) {
                    for (size_t k=0; k<momentMatrices[i][j].size(); k++) {
                        Mon monToAdd = momentMatrices[i][j][k].getKey();
                        if (!monomsInConstraints.count(monToAdd)) {
                            monomsInConstraints.insert(monToAdd);
                            queue.push_back(monToAdd);
                        }
                    }
                }
            }
            int nextQueueLoc = 0;
            while (int(constraintsZero.size()) < findMinimalAmount) {
                double ratio = double(constraintsZero.size()) / (monomsInConstraints.size()-1);
                std::cout << constraintsZero.size() << " / " << findMinimalAmount << " (" << ratio << ")        \r" << std::flush;

                // Stop if we're fully constrained
                if (monomsInConstraints.size()-1 == constraintsZero.size() && constraintsZero.size() > 2) {
                    break;
                }

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
                    std::cout << std::endl;
                    std::cout << "Starting new cycle with monomial: " << monToAdd << std::endl;
                }

                // If we can't find anything else, break
                if (monToAdd.size() == 0 && monomsUsed.count(monToAdd)) {
                    std::cout << std::endl;
                    std::cout << "Couldn't find any more monomials to add" << std::endl;
                    break;
                }

                // Put the monomial in the Limbladian
                Mon toPut(monToAdd);
                if (verbosity >= 3) {
                    std::cout << std::endl;
                    std::cout << "Putting in monomial: " << toPut << std::endl;
                }
                std::pair<char,int> oldMon('A', 0);
                Poly newConstraint = limbladian.replaced(oldMon, toPut);

                // Add the constraint
                if (!newConstraint.isZero()) {
                    constraintsZero.push_back(newConstraint); 
                }
                monomsUsed.insert(monToAdd);
                std::set<Mon> newTerms;
                for (auto& term : newConstraint) {
                    newTerms.insert(term.first);
                }
                for (auto& term : newTerms) {
                    if (!monomsInConstraints.count(term)) {
                        monomsInConstraints.insert(term);
                        queue.push_back(term);
                    }
                }

            }

            // Generate the moment matrix from the monomsUsed TODO
            // ./run --tensor 12 -M 1000 -A 50
            // -0.929473  <  -0.405878      39.2103%
            if (autoMomentAmount > 0) {

                // Remove the first 10
                //monomsUsed.erase(monomsUsed.begin(), monomsUsed.begin()+10);
                //std::erase(monomsUsed, monomsUsed.begin(), monomsUsed.begin()+10);
                // need to use advance
                //monomsUsed.erase(monomsUsed.begin(), std::next(monomsUsed.begin(), 10));

                // Add the first used monoms to the top row
                std::vector<Poly> topRow = {Poly(1)};
                int added = 0;
                for (auto& mon : monomsUsed) {
                    topRow.push_back(Poly(mon));
                    added++;
                    if (added >= autoMomentAmount) {
                        break;
                    }
                }
                
                // Add the first used monoms to the top row, checking they aren't already there
                //std::set<Mon> monomsInMatrix = {Mon()};
                //for (auto& mon : monomsUsed) {
                    //std::cout << "new monom: " << mon << std::endl;
                    //if (!monomsInMatrix.count(mon)) {
                        //std::cout << "adding" << std::endl;
                        //topRow.push_back(Poly(mon));
                        //std::vector<std::vector<Poly>> momentMatrixTemp = generateFromTopRow(topRow, verbosity);
                        //std::cout << momentMatrixTemp << std::endl;
                        //addVariables(monomsInMatrix, momentMatrixTemp);
                        //added++;
                        //if (added >= autoMomentAmount) {
                            //break;
                        //}
                    //} else {
                        //std::cout << "skipping" << std::endl;
                    //}
                //}

                // Add the most common moments in the top row
                //std::map<Mon,int> monomsUsedCount;
                //for (auto& mon : monomsUsed) {
                    //for (auto& poly : constraintsZero) {
                        //if (poly.contains(mon)) {
                            //monomsUsedCount[mon]++;
                        //}
                    //}
                //}
                //std::vector<Mon> monomsUsedSorted;
                //for (auto& mon : monomsUsedCount) {
                    //monomsUsedSorted.push_back(mon.first);
                //}
                //std::sort(monomsUsedSorted.begin(), monomsUsedSorted.end(), [&monomsUsedCount](Mon a, Mon b) {
                    //return monomsUsedCount[a] > monomsUsedCount[b];
                //});
                //for (auto& mon : monomsUsedSorted) {
                    //topRow.push_back(Poly(mon));
                    //added++;
                    //if (added >= autoMomentAmount) {
                        //break;
                    //}
                //}

                // Add first order moments to the top, split bigger ones
                //for (auto& mon : monomsUsed) {
                    //if (mon.size() <= 2) {
                        //Poly toAdd = Poly(mon);
                        //if (std::find(topRow.begin(), topRow.end(), toAdd) == topRow.end()) {
                            //topRow.push_back(Poly(toAdd));
                            //added++;
                        //}
                        //if (added >= autoMomentAmount) {
                            //break;
                        //}
                    //} else {
                        //int monSize = mon.size();
                        //Poly toAdd1 = Poly(mon.first(monSize/2));
                        //Poly toAdd2 = Poly(mon.last(monSize-(monSize/2)));
                        //if (std::find(topRow.begin(), topRow.end(), toAdd1) == topRow.end()) {
                            //topRow.push_back(Poly(toAdd1));
                            //added++;
                        //}
                        //if (std::find(topRow.begin(), topRow.end(), toAdd2) == topRow.end()) {
                            //topRow.push_back(Poly(toAdd2));
                            //added++;
                        //}
                        //if (added >= autoMomentAmount) {
                            //break;
                        //}
                    //}
                //}

                momentMatrices = {generateFromTopRow(topRow, verbosity)};

            }

        // If a value not given, binary search
        } else {

            // First try adding constraints until it's exact
            std::set<Mon> monomsInConstraints;
            std::vector<Mon> queue;
            queue.push_back(Mon());
            monomsInConstraints.insert(Mon());
            for (auto& term : objective) {
                if (!monomsInConstraints.count(term.first)) {
                    monomsInConstraints.insert(term.first);
                    queue.push_back(term.first);
                }
            }
            for (size_t i=0; i<variablesToPut.size(); i++) {
                Mon monToAdd = variablesToPut[i].getKey();
                if (!monomsInConstraints.count(monToAdd)) {
                    monomsInConstraints.insert(monToAdd);
                    queue.push_back(monToAdd);
                }
            }
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

                    // Put the monomial in the Limbladian
                    Mon toPut(monToAdd);
                    std::pair<char,int> oldMon('A', 0);
                    Poly newConstraint = limbladian.replaced(oldMon, toPut);

                    // Add the constraint
                    if (!newConstraint.isZero()) {
                        constraintsZero.push_back(newConstraint); 
                    }
                    monomsUsed.insert(monToAdd);
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
                std::cout << "cons: " << amountToTest << ", diff: " << diff << std::endl;
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
                std::cout << "cons: " << toTest << ", diff: " << diff << std::endl;
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
            std::cout << "Minimal num constraints: " << maxNum << std::endl;
            constraintsZero = constraintsZeroCopy;

        }

        // Final output
        double ratio = double(constraintsZero.size()) / (monomsUsed.size()-1);
        std::cout << constraintsZero.size() << " / " << findMinimalAmount << " (" << ratio << ")               " << std::endl;

    }

    // If adding some reductive constraints
    int newReductConsAdded = 0;
    while (newReductConsAdded < reductiveCons) {
        std::cout << newReductConsAdded << " / " << reductiveCons << "        \r" << std::flush;

        // Add constraints based on the monomials we already have
        std::set<Mon> monomsInConstraints;
        std::vector<Mon> monomsInConstraintsVec;
        for (auto& term : objective) {
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

        // Put the monomial in the Limbladian
        Mon toPut(monToAdd);
        std::pair<char,int> oldMon('A', 0);
        Poly newConstraint = limbladian.replaced(oldMon, toPut);

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

    }

    // Create the moment matrices
    if (autoMomentAmount == 0) {
        momentMatrices = generateAllMomentMatrices(objective, constraintsZero, level, verbosity);
    }

    // If told to add extra to the top row
    if (extraMonomials.size() > 0) {
        std::vector<Poly> topRow = momentMatrices[0][0];
        for (size_t i=0; i<extraMonomials.size(); i++) {
            Poly extraMonomial(Mon(extraMonomials[i]).reversed());
            topRow.push_back(extraMonomial);
            std::cout << "Added " << extraMonomial << " to the top row" << std::endl;
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

    // Add symmetries TODO
    std::set<Mon> vars;
    for (size_t i=0; i<momentMatrices.size(); i++) {
        addVariables(vars, momentMatrices[i]);
    }
    addVariables(vars, constraintsZero);
    std::map<Mon,Mon> symmetriesMap;
    for (auto sym : symmetries) {
        
        // If verbose, output the symmetry
        if (verbosity >= 1) {
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
                if (newMon < mon) {
                    symmetriesMap[mon] = newMon;
                    if (verbosity >= 3) {
                        std::cout << mon << " -> " << newMon << std::endl;
                    }
                }
            }

        }

    }

    // Replace everything using symmetries
    for (auto sym : symmetriesMap) {
        constraintsZero.push_back(Poly(sym.first)-Poly(sym.second));
    }
    //for (size_t i=0; i<momentMatrices.size(); i++) {
        //for (size_t j=0; j<momentMatrices[i].size(); j++) {
            //for (size_t k=0; k<momentMatrices[i][j].size(); k++) {
                //momentMatrices[i][j][k] = momentMatrices[i][j][k].applyMap(symmetriesMap);
            //}
        //}
    //}
    //for (size_t i=0; i<constraintsZero.size(); i++) {
        //constraintsZero[i] = constraintsZero[i].applyMap(symmetriesMap);
    //}
    //objective = objective.applyMap(symmetriesMap);

    // Timing
    std::chrono::steady_clock::time_point timeFinishedGeneration = std::chrono::steady_clock::now();

    // If removing constraints
    if (tryRemove) {

        // Initial run
        std::pair<double,double> boundsStart = boundMOSEK(objective, momentMatrices, constraintsZero, {}, verbosity);
        double currentBoundDiff = boundsStart.second - boundsStart.first;
        std::cout << "Original diff: " << currentBoundDiff << std::endl;
        std::vector<Poly> constraintsZeroCopy = constraintsZero;
        std::vector<std::vector<std::vector<Poly>>> momentMatricesCopy = momentMatrices;

        // Per pass
        for (int l=0; l<100; l++) {
            bool removedSomething = false;
            
            // Check each linear constraint
            //for (size_t i=0; i<constraintsZero.size(); i++) {

                //// Remove the constraint
                //constraintsZero = constraintsZeroCopy;
                //constraintsZero.erase(constraintsZero.begin() + i);

                //// Get the bounds
                //std::pair<double,double> boundsTemp;
                //boundsTemp = boundMOSEK(objective, momentMatrices, constraintsZero, {}, verbosity);
                //double lowerBoundTemp = boundsTemp.first;
                //double upperBoundTemp = boundsTemp.second;
                //double diff = upperBoundTemp - lowerBoundTemp;
                //std::cout << "Removed " << i << ", diff: " << diff << std::endl;

                //// If it's better, remove it
                //if (diff - currentBoundDiff <= tol) {
                    //std::cout << "Removing constraint " << i << std::endl;
                    //constraintsZeroCopy = constraintsZero;
                    //removedSomething = true;
                    //break;
                //}

            //}

            // Check each column of the moment matrix
            for (size_t i=0; i<momentMatrices[0].size(); i++) {

                    // Remove the column
                    momentMatrices = momentMatricesCopy;
                    momentMatrices[0].erase(momentMatrices[0].begin() + i);
                    for (size_t j=0; j<momentMatrices[0].size(); j++) {
                        momentMatrices[0][j].erase(momentMatrices[0][j].begin() + i);
                    }

                    // Get the bounds
                    std::pair<double,double> boundsTemp;
                    boundsTemp = boundMOSEK(objective, momentMatrices, constraintsZero, {}, verbosity);
                    double lowerBoundTemp = boundsTemp.first;
                    double upperBoundTemp = boundsTemp.second;
                    double diff = upperBoundTemp - lowerBoundTemp;
                    std::cout << "Removed row/col " << i << ", diff: " << diff << ", change: " << diff - currentBoundDiff << std::endl;

                    // If it's better, remove it
                    if (diff - currentBoundDiff <= tol) {
                        std::cout << "Removing row/col " << i << std::endl;
                        momentMatricesCopy = momentMatrices;
                        removedSomething = true;
                        break;
                    }

            }

            // If nothing was removed, break
            if (!removedSomething) {
                break;
            }

        }

        std::cout << "Final constraints: " << constraintsZeroCopy.size() << std::endl;
        constraintsZero = constraintsZeroCopy;
        momentMatrices = momentMatricesCopy;

    }

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
        std::cout << "----------------------------------------" << std::endl;
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
    std::set<Mon> variableSet;
    variableSet.insert(Mon());
    for (size_t i=0; i<momentMatrices.size(); i++) {
        addVariables(variableSet, momentMatrices[i]);
    }
    for (size_t i=0; i<constraintsZero.size(); i++) {
        addVariables(variableSet, constraintsZero[i]);
    }
    addVariables(variableSet, objective);
    int numVars = variableSet.size()-1;

    // Auto detect the best solver
    if (solver == "auto") {
        if (maxMatSize > 1) {
            solver = "mosek";
        } else if (numCons == numVars) {
            solver = "eigen";
        } else {
            solver = "gurobi";
        }
    }

    // Output what we're doing
    if (verbosity >= 1) {
        if (solver == "scs") {
            std::cout << "Solving SDP using SCS with " << numCons << " constraints, max moment mat size of " << maxMatSize << " and " << numVars << " variables..." << std::endl;
        } else if (maxMatSize > 1 && solver == "mosek") {
            std::cout << "Solving SDP using MOSEK with " << numCons << " constraints, max moment mat size of " << maxMatSize << " and " << numVars << " variables..." << std::endl;
        } else if (solver == "eigen") {
            std::cout << "Solving linear system using Eigen with " << numCons << " constraints and " << numVars << " variables..." << std::endl;
        } else if (solver == "mosek") {
            std::cout << "Solving LP using MOSEK with " << numCons << " constraints and " << numVars << " variables (ratio: " << double(numCons)/numVars << ")..." << std::endl;
        } else if (solver == "gurobi") {
            std::cout << "Solving LP using Gurobi with " << numCons << " constraints and " << numVars << " variables (ratio: " << double(numCons)/numVars << ")..." << std::endl;
        }
    }

    // Solve
    std::pair<double,double> bounds;
    if (solver == "scs") {
        bounds = boundSCS(objective, momentMatrices, constraintsZero, {}, verbosity);
    } else if (solver == "mosek" || maxMatSize > 1) {
        bounds = boundMOSEK(objective, momentMatrices, constraintsZero, {}, verbosity);
    } else if (solver == "gurobi") {
        bounds = boundGurobi(objective, constraintsZero, verbosity);
    } else if (solver == "eigen") {
        double res = solveEigen(objective, constraintsZero, verbosity, numCores);
        bounds = {res, res};
    }
    double lowerBound = bounds.first;
    double upperBound = bounds.second;
    if (verbosity >= 1) {
        std::cout << "Bounds: " << lowerBound << "  <  " << upperBound << std::endl;
        std::cout << "Difference: " << upperBound - lowerBound << std::endl;
        if (idealIsKnown) {
            std::cout << "Known ideal: " << knownIdeal << std::endl;
            std::cout << "Error: " << 100*(upperBound-lowerBound)/knownIdeal << "%" << std::endl;
        } else if (upperBound * lowerBound > 0) {
            knownIdeal = (upperBound + lowerBound) / 2;
            std::cout << "Averaged: " << knownIdeal << " +/- " << 100*(upperBound-lowerBound)/std::abs(2*knownIdeal) << "%" << std::endl;
        }
    }

    // Timing
    std::chrono::steady_clock::time_point timeFinishedSolving = std::chrono::steady_clock::now();

    // Output the timings
    if (verbosity >= 1) { 
        std::cout << "Time to generate: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedGeneration - timeFinishedArgs).count() / 1000.0 << "s" << std::endl;
        std::cout << "Time to solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedSolving - timeFinishedGeneration).count() / 1000.0 << "s" << std::endl;
    }

    // Exit without errors
    return 0;

}
