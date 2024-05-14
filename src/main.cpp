// Standard includes
#include <iostream>
#include <vector>
#include <complex>
#include <unordered_map>
#include <unordered_set>

// Import Eigen
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

// Local files
#include "poly.h"
#include "printing.h"
#include "utils.h"
#include "solver.h"
 
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

    // Define the scenario
    int level = 1;
    int limbladLevel = 1;
    int numQubits = 1;
    Poly objective("<Z1>");
    int verbosity = 1;
    std::complex<double> knownIdeal = 0.0;
    bool idealIsKnown = false;
    bool findMinimal = false;
    int findMinimalAmount = 0;
    std::string seed = "";
    Poly limbladian("<X1A0X1>+<Y1A0Y1>-<A0>");
    std::vector<std::string> extraMonomials;
    std::vector<std::string> extraMonomialsLim;
    std::vector<Poly> constraintsZero;
    std::vector<int> reductionsToIgnore = {};

    // Process command-line args
    for (int i=1; i<argc; i++) {
        std::string argAsString = std::string(argv[i]);

        // Set the level of the moment matrix
        if (argAsString == "-m") {
            level = std::stoi(argv[i+1]);
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

            //double J_ss = (8.0*g*g*(gamma_h_plus*gamma_c_minus - gamma_h_minus*gamma_c_plus) / chi) * (epsilon_h*Gamma_c + epsilon_c*Gamma_h);
            //std::cout << "ideal J_ss: " << J_ss << std::endl;

            Eigen::MatrixXcd rho = Eigen::MatrixXcd::Zero(4,4);
            rho(0,0) = (4.0*g*g*(gamma_h_plus + gamma_c_plus)*(gamma_h_plus + gamma_c_plus)  + gamma_h_plus*gamma_c_plus*(Gamma*Gamma + 4.0*delta*delta)) / chi;
            rho(1,1) = (4.0*g*g*(gamma_h_minus + gamma_c_minus)*(gamma_h_plus + gamma_c_plus)  + gamma_h_plus*gamma_c_minus*(Gamma*Gamma + 4.0*delta*delta)) / chi;
            rho(2,2) = (4.0*g*g*(gamma_h_minus + gamma_c_minus)*(gamma_h_plus + gamma_c_plus)  + gamma_h_minus*gamma_c_plus*(Gamma*Gamma + 4.0*delta*delta)) / chi;
            rho(3,3) = (4.0*g*g*(gamma_h_minus + gamma_c_minus)*(gamma_h_minus + gamma_c_minus)  + gamma_h_minus*gamma_c_minus*(Gamma*Gamma + 4.0*delta*delta)) / chi;
            rho(1,2) = (2.0*g*(gamma_h_plus*gamma_c_minus - gamma_h_minus*gamma_c_plus)*(imag*Gamma-2.0*delta)) / chi;
            rho(2,1) = std::conj(rho(1,2));

            //std::cout << "ideal rho: " << std::endl;
            //std::cout << rho << std::endl;
            
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

            Eigen::MatrixXcd H_s = epsilon_h*sigma_h_plus*sigma_h_minus 
                                   + epsilon_c*sigma_c_plus*sigma_c_minus;  
            Eigen::MatrixXcd term2 = gamma_h_plus*(sigma_h_plus*rho*sigma_h_minus 
                                                   - 0.5*sigma_h_minus*sigma_h_plus*rho 
                                                   - 0.5*rho*sigma_h_minus*sigma_h_plus) 
                                   + gamma_h_minus*(sigma_h_minus*rho*sigma_h_plus 
                                                   - 0.5*sigma_h_plus*sigma_h_minus*rho 
                                                   - 0.5*rho*sigma_h_plus*sigma_h_minus);
            double Q_ss_h_direct = tr(H_s*term2);
            
            double exp_hmhphmhp = tr(sigma_h_minus*sigma_h_plus*sigma_h_minus*sigma_h_plus*rho);
            double exp_hphmhphm = tr(sigma_h_plus*sigma_h_minus*sigma_h_plus*sigma_h_minus*rho);
            double exp_hmhp = tr(sigma_h_minus*sigma_h_plus*rho);
            double exp_hphm = tr(sigma_h_plus*sigma_h_minus*rho);
            
            double Q_ss_h_reduced = epsilon_h*gamma_h_plus*exp_hmhp - epsilon_h*gamma_h_minus*exp_hphm;
                
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
            objective.sort();
            
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
            limbladian.sort();

        // Many-body Limbladian TODO
        } else if (argAsString == "--many" || argAsString == "--manyv") {

            // Defining quantities
            double gamma_c = 1.1e-2;
            double gamma_h = 1e-3;
            double g = 1.6e-3;
            double T_h = 1.0;
            double T_c = 0.1;
            double delta = 0.005;
            double epsilon_h = 1.0;

            // Regardless we have an
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
                limbladian += Poly(-imag*g, "<A0P" + std::to_string(i) + "M" + std::to_string(i+1) + ">");
                limbladian += Poly(-imag*g, "<A0M" + std::to_string(i) + "P" + std::to_string(i+1) + ">");
                limbladian += Poly(imag*g, "<P" + std::to_string(i) + "M" + std::to_string(i+1) + "A0>");
                limbladian += Poly(imag*g, "<M" + std::to_string(i) + "P" + std::to_string(i+1) + "A0>");
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
            limbladian.sort();

        // Set the seed
        } else if (argAsString == "-S") {
            seed = std::string(argv[i+1]);
            i++;

        // If told to ignore certain Pauli reductions
        } else if (argAsString == "-2") {
            reductionsToIgnore.push_back(2);
        } else if (argAsString == "-3") {
            reductionsToIgnore.push_back(3);

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

        // Output the help
        } else if (argAsString == "-h" || argAsString == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -O <str>        Manually set the objective" << std::endl;
            std::cout << "  -L <str>        Manually set the Limbladian" << std::endl;
            std::cout << "  --objX          Use sigma_X as the objective" << std::endl;
            std::cout << "  --objY          Use sigma_Y as the objective" << std::endl;
            std::cout << "  --objZ          Use sigma_Z as the objective" << std::endl;
            std::cout << "  -2              Don't use second order Pauli reductions" << std::endl;
            std::cout << "  -3              Don't use third order Pauli reductions" << std::endl;
            std::cout << "  --pauli <num> <num> <num>" << std::endl;
            std::cout << "                  Use the Pauli Limbladian with coeffs" << std::endl;
            std::cout << "  --two" << std::endl;
            std::cout << "  --twov <num> <num>" << std::endl;
            std::cout << "                  Use the two-body Limbladian with coeffs" << std::endl;
            std::cout << "  --many <num>" << std::endl;
            std::cout << "  --manyv <num> <num> <num>" << std::endl;
            std::cout << "                  Use the many-body Limbladian with coeffs" << std::endl;
            std::cout << "  -m <num>        Level of the moment matrix" << std::endl;
            std::cout << "  -l <num>        Level of the moments to put in the Limbladian" << std::endl;
            std::cout << "  -e <monom>      Add an extra monomial to the top row of the moment matrix" << std::endl;
            std::cout << "  -E <monom>      Add an extra monomial to the list of Limbladian replacements" << std::endl;
            std::cout << "  -S <str>        Seed for the random number generator" << std::endl;
            std::cout << "  -M              Try to find the minimal set of linear constraints" << std::endl;
            std::cout << "  -v <num>        Verbosity level" << std::endl;
            std::cout << "  -t <num>        Run a section of not-yet-finished code" << std::endl;
            return 0;

        // Otherwise we don't know what this is
        } else if (argAsString != "./run") {
            std::cout << "Unknown argument: " << argAsString << std::endl;
            return 1;

        }
    }

    // If the seed isn't set
    if (seed == "") {
        srand(time(NULL));
    } else {
        srand(std::stoi(seed));
    }

    // Create the Limbladian applied to many different operators
    std::vector<Mon> variables = {};
    for (int i=0; i<numQubits; i++) {
        variables.push_back(Mon("<X" + std::to_string(i+1) + ">"));
        variables.push_back(Mon("<Y" + std::to_string(i+1) + ">"));
        variables.push_back(Mon("<Z" + std::to_string(i+1) + ">"));
    }
    std::vector<Poly> variablesToPut = generateMonomials(variables, limbladLevel, verbosity);
    for (int i=0; i<extraMonomialsLim.size(); i++) {
        variablesToPut.push_back(Poly(extraMonomialsLim[i]));
    }
    if (verbosity >= 2) {
        std::cout << std::endl;
        std::cout << "Variables to put in Limbladian: " << std::endl;
        for (int i=0; i<variablesToPut.size(); i++) {
            std::cout << variablesToPut[i] << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Original Limbladian: " << limbladian << std::endl;
    }
    for (int i=0; i<variablesToPut.size(); i++) {
        std::pair<char,int> oldMon('A', 0);
        Poly newPoly(variablesToPut[i]);
        Poly newConstraint = limbladian.replaced(oldMon, newPoly);
        if (verbosity >= 3) {
            std::cout << std::endl;
            std::cout << "Variable to put: " << newPoly << std::endl;
            std::cout << "New constraint: " << newConstraint << std::endl;
        }
        constraintsZero.push_back(newConstraint);
    }

    // Reduce the constraints as much as possible
    for (int i=0; i<constraintsZero.size(); i++) {
        if (verbosity >= 3) {
            std::cout << std::endl;
            std::cout << "Reducing constraint: " << constraintsZero[i] << std::endl;
        }
        constraintsZero[i].reduce();
        if (verbosity >= 3) {
            std::cout << "Reduced constraint: " << constraintsZero[i] << std::endl;
        }
    }

    // Define the moment matrix
    std::vector<std::vector<std::vector<Poly>>> momentMatrices = generateAllMomentMatrices(objective, constraintsZero, level, verbosity, reductionsToIgnore);

    // If told to add extra to the top row
    if (extraMonomials.size() > 0) {
        std::vector<Poly> topRow = momentMatrices[0][0];
        for (int i=0; i<extraMonomials.size(); i++) {
            Poly extraMonomial(extraMonomials[i]);
            for (int j=0; j<extraMonomial.size(); j++) {
                extraMonomial[j].second.reverse();
            }
            topRow.push_back(extraMonomial);
            std::cout << "Added " << extraMonomial << " to the top row" << std::endl;
        }
        momentMatrices[0] = generateFromTopRow(topRow, verbosity);
    }

    // See how big the moment matrix is
    if (verbosity >= 1) {
        int largestMomentMatrix = 0;
        for (int i=0; i<momentMatrices.size(); i++) {
            if (momentMatrices[i].size() > largestMomentMatrix) {
                largestMomentMatrix = momentMatrices[i].size();
            }
        }
    }

    // If asked to try removing constraints
    if (findMinimal) {

        // Add constraints based on the monomials we already have TODO
        std::unordered_set<Mon> monomsUsed;
        std::unordered_set<Mon> monomsInConstraints;
        std::vector<Mon> monomsInCon = objective.monomials();
        for (int i=0; i<monomsInCon.size(); i++) {
            monomsInConstraints.insert(monomsInCon[i]);
        }
        for (int i=0; i<constraintsZero.size(); i++) {
            monomsInCon = constraintsZero[i].monomials();
            for (int j=0; j<monomsInCon.size(); j++) {
                monomsInConstraints.insert(monomsInCon[j]);
            }
        }
        for (int i=0; i<variablesToPut.size(); i++) {
            monomsUsed.insert(variablesToPut[i][0].second);
        }
        while (constraintsZero.size() < findMinimalAmount) {

            // Find a monomial that hasn't been used
            Mon monToAdd;
            for (auto& mon : monomsInConstraints) {
                if (monomsUsed.find(mon) == monomsUsed.end()) {
                    monToAdd = mon;
                    break;
                }
            }

            // Put the monomial in the Limbladian
            Poly toPut(monToAdd);
            std::pair<char,int> oldMon('A', 0);
            Poly newConstraint = limbladian.replaced(oldMon, toPut);

            // Add the constraint
            constraintsZero.push_back(newConstraint); 
            monomsUsed.insert(monToAdd);
            monomsInCon = newConstraint.monomials();
            for (int j=0; j<monomsInCon.size(); j++) {
                monomsInConstraints.insert(monomsInCon[j]);
            }

            // Run the SDP
            std::vector<Mon> varNames2;
            std::vector<std::complex<double>> varVals2;
            double upperBoundTemp = solveMOSEK(objective, momentMatrices, constraintsZero, verbosity, varNames2, varVals2);
            for (int i=0; i<objective.size(); i++) {
                objective[i].first *= -1;
            }
            double lowerBoundTemp = -solveMOSEK(objective, momentMatrices, constraintsZero, verbosity, varNames2, varVals2);
            std::cout << "Constraint: " << monToAdd << std::endl;
            std::cout << "Lower bound: " << lowerBoundTemp << std::endl;
            std::cout << "Upper bound: " << upperBoundTemp << std::endl;
            if (upperBoundTemp - lowerBoundTemp < 1e-4) {
                break;
            }


        }

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
            for (int i=0; i<momentMatrices.size(); i++) {
                std::cout << "Moment matrix " << i  << " (" << momentMatrices[i].size() << "x" << momentMatrices[i].size() << "): " << std::endl;
                if (momentMatrices[i].size() < 10 || verbosity >= 3) {
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

    // Convert to MOSEK form and solve
    std::vector<Mon> varNames;
    std::vector<std::complex<double>> varVals;
    double upperBound = solveMOSEK(objective, momentMatrices, constraintsZero, verbosity, varNames, varVals);
    for (int i=0; i<objective.size(); i++) {
        objective[i].first *= -1;
    }
    double lowerBound = -solveMOSEK(objective, momentMatrices, constraintsZero, verbosity, varNames, varVals);
    if (verbosity >= 2) {
        std::cout << std::endl;
    }
    if (verbosity >= 1) {
        if (idealIsKnown) {
            std::cout << "Lower bound: " << lowerBound << " (" << 100*(lowerBound-knownIdeal)/knownIdeal << "%)" << std::endl;
            std::cout << "Upper bound: " << upperBound << " (" << 100*(upperBound-knownIdeal)/knownIdeal << "%)" << std::endl;
            std::cout << "Known ideal: " << knownIdeal << std::endl;
            std::cout << "Difference: " << upperBound - lowerBound << std::endl;
            std::cout << "Error: " << 100*(upperBound-lowerBound)/knownIdeal << "%" << std::endl;
        } else {
            std::cout << "Lower bound: " << lowerBound << std::endl;
            std::cout << "Upper bound: " << upperBound << std::endl;
            std::cout << "Difference: " << upperBound - lowerBound << std::endl;

        }
    }

    // Exit without errors
    return 0;

}
