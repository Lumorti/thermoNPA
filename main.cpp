// Standard includes
#include <iostream>
#include <vector>

// Import MOSEK
#include "fusion.h"

// Import Eigen
#include <Eigen/Dense>

// Define the type for the polynomials
typedef std::vector<std::pair<char, int>> monomial;
typedef std::vector<std::pair<double, monomial>> polynomialReal;
typedef std::vector<std::pair<std::complex<double>, monomial>> polynomial;
typedef std::vector<polynomial> polynomialVector;
typedef std::vector<std::vector<polynomial>> polynomialMatrix;

// Given a polynomial, combine all of the terms of the same monomial
polynomial simplify(polynomial p) {

    // Create the return polynomial
    polynomial toReturn;

    // Iterate through the polynomial
    for (unsigned long int i=0; i<p.size(); i++) {

        // Check if this monomial is already in the return polynomial
        bool found = false;
        for (unsigned long int j=0; j<toReturn.size(); j++) {
            if (toReturn[j].second == p[i].second) {
                toReturn[j].first += p[i].first;
                found = true;
                break;
            }
        }

        // If not, add it
        if (!found) {
            toReturn.push_back(p[i]);
        }

    }

    // Return the simplified polynomial
    return toReturn;

}

// Generate a monomial from a given string (e.g. "<A1B2>")
monomial stringToMonomial(std::string asString) {

    // If the string is empty, return an empty monomial
    if (asString == "" || asString.find("<") == std::string::npos) {
        return monomial();
    }

    // Iterate through the string
    monomial toReturn;
    char currentVariable = ' ';
    std::string currentIndex;
    for (unsigned long int i=0; i<asString.size(); i++) {

        // Check if this char is an integer
        if (asString[i] >= '0' && asString[i] <= '9') {
            currentIndex += asString[i];

        // Otherwise it's a variable
        } else if (asString[i] >= 'A' && asString[i] <= 'Z') {
            currentIndex += asString[i];
            if (currentVariable != ' ') {
                toReturn.push_back(std::make_pair(currentVariable, std::stoi(currentIndex)));
            }
            currentVariable = asString[i];
            currentIndex = "";
        }

    }

    // Add the last variable
    toReturn.push_back(std::make_pair(currentVariable, std::stoi(currentIndex)));

    // Return the monomial
    return toReturn;

}

// Pretty print complex numbers
std::ostream& operator<<(std::ostream& os, const std::complex<double>& c) {
    if (c.imag() == 0) {
        os << c.real();
    } else if (c.real() == 0) {
        if (c.imag() == 1) {
            os << "i";
        } else if (c.imag() == -1) {
            os << "-i";
        } else {
            os << c.imag() << "i";
        }
    } else {
        os << c.real() << "+" << c.imag() << "i";
    }
    return os;
}
        
// When printing a monomial, print it as a string
std::ostream& operator<<(std::ostream& os, const monomial& m) {

    // Iterate through the monomial
    if (m.size() > 0) {
        os << "<";
        for (unsigned long int i=0; i<m.size(); i++) {
            os << m[i].first << m[i].second;
        }
        os << ">";
    } else {
        os << "1";
    }

    // Return the stream
    return os;

}

// When printing a polynomial, print it as a string
std::ostream& operator<<(std::ostream& os, const polynomial& p) {

    // Check if it's zero
    if (p.size() == 0 || (p.size() == 1 && p[0].first == std::complex<double>(0,0))) {
        os << "0";
        return os;
    }

    // Iterate through the polynomial
    for (unsigned long int i=0; i<p.size(); i++) {
        double realPart = p[i].first.real();
        double imagPart = p[i].first.imag();
        if (imagPart == 0) {
            if (realPart == std::complex<double>(-1, 0)) {
                os << "-" << p[i].second;
            } else if (realPart == std::complex<double>(1, 0)) {
                if (i == 0) {
                    os << p[i].second;
                } else {
                    os << "+" << p[i].second;
                }
            } else if (p[i].second.size() == 0 && realPart > 0) {
                if (i == 0) {
                    os << realPart;
                } else {
                    os << "+" << realPart;
                }
            } else if (p[i].second.size() == 0 && realPart < 0) {
                os << realPart;
            } else if (realPart > 0) {
                if (i == 0) {
                    os << realPart << p[i].second;
                } else {
                    os << "+" << realPart << p[i].second;
                }
            } else if (realPart != 0) {
                os << realPart << p[i].second;
            }
        } else {
            os << p[i].first << p[i].second;
        }
    }

    // Return the stream
    return os;

}

// When printing a polynomial matrix, print it as a string
std::ostream& operator<<(std::ostream& os, const polynomialMatrix& m) {

    // Determine the maximum width of each column
    std::vector<int> columnWidths(m[0].size(), 0);
    for (unsigned long int i=0; i<m.size(); i++) {
        for (unsigned long int j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            columnWidths[j] = std::max(columnWidths[j], sizeWhenWritten);
        }
    }

    // Iterate through the matrix
    for (unsigned long int i=0; i<m.size(); i++) {
        for (unsigned long int j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            for (unsigned long int k=0; k<columnWidths[j]-sizeWhenWritten; k++) {
                ss << " ";
            }
            ss << " ";
            os << ss.str();
        }
        os << std::endl;
    }

    // Return the stream
    return os;

}

// When printing a matrix of doubles, print it as a string
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& m) {

    // Determine the maximum width of each column
    std::vector<int> columnWidths(m[0].size(), 0);
    for (unsigned long int i=0; i<m.size(); i++) {
        for (unsigned long int j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            columnWidths[j] = std::max(columnWidths[j], sizeWhenWritten);
        }
    }

    // Iterate through the matrix
    for (unsigned long int i=0; i<m.size(); i++) {
        for (unsigned long int j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            for (unsigned long int k=0; k<columnWidths[j]-sizeWhenWritten; k++) {
                ss << " ";
            }
            ss << " ";
            os << ss.str();
        }
        os << std::endl;
    }

    // Return the stream
    return os;

}

// When printing a vector of polynomials
std::ostream& operator<<(std::ostream& os, const std::vector<polynomial>& v) {

    // Iterate through the vector
    for (unsigned long int i=0; i<v.size(); i++) {
        os << v[i] << std::endl;
    }

    // Return the stream
    return os;

}

// When printing a vector of something
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {

    // Iterate through the vector
    os << "{";
    for (unsigned long int i=0; i<v.size(); i++) {
        os << v[i];
        if (i < v.size()-1) {
            os << ", ";
        }
    }
    os << "}";

    // Return the stream
    return os;

}

// Generate a polynomial from a given string (e.g. "2<A1B2>-<A3A1>")
polynomial stringToPolynomial(std::string asString) {

    // Remove all spaces and *
    std::string newString = "";
    for (unsigned long int i=0; i<asString.size(); i++) {
        if (asString[i] != ' ' && asString[i] != '*') {
            newString += asString[i];
        }
    }
    asString = newString;

    // Iterate through the string
    polynomial toReturn;
    std::string currentCoefficient;
    std::string currentMonomial;
    for (unsigned long int i=0; i<asString.size(); i++) {

        // If finished a monomial
        if (i > 0 && (asString[i] == '>' || ((asString[i] == '+' || asString[i] == '-') && asString[i-1] != '>'))) {

            // Ensure the coefficient is convertable
            if (currentCoefficient == "" || currentCoefficient == "+") {
                currentCoefficient = "1";
            } else if (currentCoefficient == "-") {
                currentCoefficient = "-1";
            }

            // Add the term
            if (asString[i] == '>') {
                currentMonomial += asString[i];
            }
            toReturn.push_back(std::make_pair(std::stod(currentCoefficient), stringToMonomial(currentMonomial)));
            currentMonomial = "";
            currentCoefficient = "";

            // The plus and the minus are for the next term
            if (asString[i] == '+' || asString[i] == '-') {
                currentCoefficient += asString[i];
            }

        // If starting or continuing a monomial
        } else if (asString[i] == '<' || currentMonomial.size() > 0) {
            currentMonomial += asString[i];

        // Otherwise it's for the coefficient
        } else if (asString[i] != ' ') {
            currentCoefficient += asString[i];

        }

    }

    // If there's a coefficient left over, add it
    if (currentCoefficient != "") {
        toReturn.push_back(std::make_pair(std::stod(currentCoefficient), stringToMonomial(currentMonomial)));
    }

    // Return the polynomial
    return toReturn;

}

// Given a matrix and a variable list, return the matrix with the variables replaced
Eigen::MatrixXcd replaceVariables(polynomialMatrix& momentMatrix, const std::vector<monomial>& variables, const std::vector<std::complex<double>>& varVals) {

    // Replace each variable with its value
    Eigen::MatrixXcd momentMatrixEigen = Eigen::MatrixXcd::Zero(momentMatrix.size(), momentMatrix.size());
    for (int i=0; i<momentMatrix.size(); i++) {
        for (int j=0; j<momentMatrix[i].size(); j++) {
            for (int k=0; k<momentMatrix[i][j].size(); k++) {
                for (int l=0; l<variables.size(); l++) {
                    if (momentMatrix[i][j][k].second == variables[l]) {
                        if (i > j) {
                            momentMatrixEigen(i, j) += momentMatrix[i][j][k].first * varVals[l];
                        } else {
                            momentMatrixEigen(i, j) += momentMatrix[i][j][k].first * std::conj(varVals[l]);
                        }
                        break;
                    }
                }
            }
        }
    }
    return momentMatrixEigen;

}

// Get the eigenvalues and vectors of a matrix after replacement
void getEigens(polynomialMatrix& momentMatrix, const std::vector<monomial>& variables, const std::vector<std::complex<double>>& varVals, std::vector<std::vector<std::complex<double>>>& eigenvectors, std::vector<std::complex<double>>& eigenvalues) {

    // Replace each variable with its value
    Eigen::MatrixXcd momentMatrixEigen = replaceVariables(momentMatrix, variables, varVals);

    // Get the eigenvalues and vectors of this
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(momentMatrixEigen);
    Eigen::MatrixXcd eigenVectorsEigen = es.eigenvectors();
    Eigen::VectorXcd eigenValuesEigen = es.eigenvalues();

    // Copy into the output vectors
    eigenvalues.clear();
    eigenvectors.clear();
    for (int i=0; i<eigenVectorsEigen.cols(); i++) {
        std::vector<std::complex<double>> eigenVector;
        for (int j=0; j<eigenVectorsEigen.rows(); j++) {
            eigenVector.push_back(eigenVectorsEigen(j, i));
        }
        eigenvectors.push_back(eigenVector);
    }
    for (int i=0; i<eigenValuesEigen.size(); i++) {
        eigenvalues.push_back(eigenValuesEigen(i));
    }

}

// Given a monomial, reduce it to its simplest form
std::pair<std::complex<double>, monomial> reduceMonomial(const monomial& mon_, int verbosity, bool trySwap=true, bool diffLettersCommute=false, bool pauliReductions=true, std::vector<int> reductionsToIgnore={}) {

    // Sort the monomial as much as we can
    monomial mon = mon_;
    if (diffLettersCommute) {
        for (int i=0; i<mon.size(); i++) {
            for (int j=0; j<int(mon.size())-1; j++) {
                if (mon[j].first != mon[j+1].first && mon[j] > mon[j+1]) {
                    std::swap(mon[j], mon[j+1]);
                }
            }
        }
    }

    // <A1A1> = <1>
    int i = 0;
    while (i < int(mon.size())-1) {
        if (mon[i] == mon[i+1]) {
            mon.erase(mon.begin()+i+1);
            mon.erase(mon.begin()+i);
            i = -1;
        }
        i++;
    }

    // Simplify Pauli strings using the following rules (plus conjugates):
    // XY = iZ
    // ZX = iY
    // YZ = iX
    std::complex<double> imag(0, 1);
    std::complex<double> coeff(1, 0);
    if (std::find(reductionsToIgnore.begin(), reductionsToIgnore.end(), mon.size()) == reductionsToIgnore.end() && pauliReductions) {
        for (i=mon.size()-1; i>0; i--) {
            if (mon[i-1].second == mon[i].second) {
                if (mon[i-1].first == 'X' && mon[i].first == 'Y') {
                    coeff *= imag;
                    mon[i-1] = std::make_pair('Z', mon[i-1].second);
                    mon.erase(mon.begin()+i);
                } else if (mon[i-1].first == 'X' && mon[i].first == 'Z') {
                    coeff *= -imag;
                    mon[i-1] = std::make_pair('Y', mon[i-1].second);
                    mon.erase(mon.begin()+i);
                } else if (mon[i-1].first == 'Y' && mon[i].first == 'Z') {
                    coeff *= imag;
                    mon[i-1] = std::make_pair('X', mon[i-1].second);
                    mon.erase(mon.begin()+i);
                } else if (mon[i-1].first == 'Y' && mon[i].first == 'X') {
                    coeff *= -imag;
                    mon[i-1] = std::make_pair('Z', mon[i-1].second);
                    mon.erase(mon.begin()+i);
                } else if (mon[i-1].first == 'Z' && mon[i].first == 'X') {
                    coeff *= imag;
                    mon[i-1] = std::make_pair('Y', mon[i-1].second);
                    mon.erase(mon.begin()+i);
                } else if (mon[i-1].first == 'Z' && mon[i].first == 'Y') {
                    coeff *= -imag;
                    mon[i-1] = std::make_pair('X', mon[i-1].second);
                    mon.erase(mon.begin()+i);
                }
            }
        }
    }

    // Flip it to see if it's smaller
    if (trySwap) {
        monomial monFlipped = mon;
        std::reverse(monFlipped.begin(), monFlipped.end());
        if (monFlipped < mon) {
            mon = monFlipped;
        }
    }

    // Verbose output
    if (verbosity >= 3) {
        std::cout << "Reduced monomial: " << mon_ << " -> " << coeff << " " << mon << std::endl;
    }

    return {coeff, mon};

}

// Given a polynomial, reduce each monomial and combine
polynomial reducePolynomial(polynomial p, int verbosity, bool trySwap=true, bool diffLettersCommute=false, bool pauliReductions=true, std::vector<int> reductionsToIgnore={}) {

    // Apply the reduction to each monomial
    for (int j=0; j<p.size(); j++) {
        std::pair<std::complex<double>, monomial> reducedMonomial = reduceMonomial(p[j].second, verbosity, trySwap, diffLettersCommute, pauliReductions, reductionsToIgnore);
        p[j].first *= reducedMonomial.first;
        p[j].second = reducedMonomial.second;
    }

    // Combine like terms
    p = simplify(p);

    return p;

}

// Add all single order monomials from a functional to a list of variables
void addSingleMonomials(std::vector<monomial>& variables, polynomial functional) {

    // Iterate through the functional
    for (long unsigned int i=0; i<functional.size(); i++) {

        // Iterate through the monomial
        for (long unsigned int j=0; j<functional[i].second.size(); j++) {
            monomial currentMonomial = {functional[i].second[j]};

            // Check if this moment is already in the list
            bool found = false;
            for (long unsigned int k=0; k<variables.size(); k++) {
                if (variables[k] == currentMonomial) {
                    found = true;
                    break;
                }
            }

            // If not, add it in the correct place
            if (!found) {
                bool hasBeenAdded = false;
                for (long unsigned int k=0; k<variables.size(); k++) {
                    if (variables[k][0] > currentMonomial[0]) {
                        variables.insert(variables.begin()+k, currentMonomial);
                        hasBeenAdded = true;
                        break;
                    }
                }
                if (!hasBeenAdded) {
                    variables.push_back(currentMonomial);
                }
            }

        }

    }

}

// Generate a moment matrix given the top row
polynomialMatrix generateFromTopRow(std::vector<monomial> monomsInTopRow, int verbosity) {

    // Generate all combinations of the top row
    polynomialMatrix matrixToReturn = std::vector<std::vector<polynomial>>(monomsInTopRow.size(), std::vector<polynomial>(monomsInTopRow.size()));
    for (long unsigned int i=0; i<monomsInTopRow.size(); i++) {
        for (long unsigned int j=i; j<monomsInTopRow.size(); j++) {

            // Form the new moment
            monomial newMonomial;
            for (long unsigned int k=0; k<monomsInTopRow[j].size(); k++) {
                newMonomial.push_back(monomsInTopRow[j][k]);
            }
            for (int k=int(monomsInTopRow[i].size())-1; k>=0; k--) {
                newMonomial.push_back(monomsInTopRow[i][k]);
            }

            // Reduce the monomial
            std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(newMonomial, verbosity);

            // Set the matrix elements
            matrixToReturn[i][j] = polynomial(1, std::make_pair(monomCoeff.first, monomCoeff.second));
            matrixToReturn[j][i] = polynomial(1, std::make_pair(std::conj(monomCoeff.first), monomCoeff.second));

        }
    }

    // Return the matrix
    return matrixToReturn;

}

// Generate a moment matrix given the top row as polynomials
polynomialMatrix generateFromTopRow(std::vector<polynomial> monomsInTopRow, int verbosity) {

    // Generate all combinations of the top row
    polynomialMatrix matrixToReturn = std::vector<std::vector<polynomial>>(monomsInTopRow.size(), std::vector<polynomial>(monomsInTopRow.size()));
    for (long unsigned int i=0; i<monomsInTopRow.size(); i++) {
        for (long unsigned int j=i; j<monomsInTopRow.size(); j++) {

            // Form the new polynomial
            polynomial newPolynomial;
            for (long unsigned int k=0; k<monomsInTopRow[i].size(); k++) {
                for (long unsigned int l=0; l<monomsInTopRow[j].size(); l++) {

                    // Form the new monomial
                    monomial newMonomial;
                    for (long unsigned int m=0; m<monomsInTopRow[i][k].second.size(); m++) {
                        newMonomial.push_back(monomsInTopRow[i][k].second[m]);
                    }
                    for (int m=int(monomsInTopRow[j][l].second.size())-1; m>=0; m--) {
                        newMonomial.push_back(monomsInTopRow[j][l].second[m]);
                    }

                    // Reduce the monomial
                    std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(newMonomial, verbosity);
                    newMonomial = monomCoeff.second;

                    // Add to the polynomial
                    bool found = false;
                    std::complex<double> newCoefficient = monomCoeff.first * monomsInTopRow[i][k].first * std::conj(monomsInTopRow[j][l].first);
                    for (long unsigned int m=0; m<newPolynomial.size(); m++) {
                        if (newPolynomial[m].second == newMonomial) {
                            newPolynomial[m].first += newCoefficient;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        newPolynomial.push_back(std::make_pair(newCoefficient, newMonomial));
                    }

                }
            }

            // The conjugate of the polynomial
            polynomial newPolynomialConj;
            for (long unsigned int k=0; k<newPolynomial.size(); k++) {
                newPolynomialConj.push_back(std::make_pair(std::conj(newPolynomial[k].first), newPolynomial[k].second));
            }

            // Set the matrix elements
            matrixToReturn[i][j] = newPolynomial;
            matrixToReturn[j][i] = newPolynomialConj;

        }
    }

    // Return the matrix
    return matrixToReturn;

}

// Generate a list of monomials of a certain level
std::vector<polynomial> generateMonomials(std::vector<monomial> variables, int level, int verbosity) {

    // Generate all monomials up to the given level
    std::vector<polynomial> monomsInTopRow = {};
    if (level >= 1) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            monomial currentMonomial = {variables[i][0]};
            std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(currentMonomial, verbosity, false);
            polynomial currentPolynomial = {monomCoeff};
            if (std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentPolynomial) == monomsInTopRow.end()) {
                monomsInTopRow.push_back(currentPolynomial);
            }
        }
    }
    if (level >= 2) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                monomial currentMonomial = {variables[i][0], variables[j][0]};
                std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(currentMonomial, verbosity, false);
                polynomial currentPolynomial = {monomCoeff};
                if (std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentPolynomial) == monomsInTopRow.end()) {
                    monomsInTopRow.push_back(currentPolynomial);
                }
            }
        }
    }
    if (level >= 3) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0]};
                    std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(currentMonomial, verbosity, false);
                    polynomial currentPolynomial = {monomCoeff};
                    if (std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentPolynomial) == monomsInTopRow.end()) {
                        monomsInTopRow.push_back(currentPolynomial);
                    }
                }
            }
        }
    }
    if (level >= 4) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    for (long unsigned int l=0; l<variables.size(); l++) {
                        monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0], variables[l][0]};
                        std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(currentMonomial, verbosity, false);
                        polynomial currentPolynomial = {monomCoeff};
                        if (std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentPolynomial) == monomsInTopRow.end()) {
                            monomsInTopRow.push_back(currentPolynomial);
                        }
                    }
                }
            }
        }
    }
    if (level >= 5) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    for (long unsigned int l=0; l<variables.size(); l++) {
                        for (long unsigned int m=0; m<variables.size(); m++) {
                            monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0], variables[l][0], variables[m][0]};
                            std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(currentMonomial, verbosity, false);
                            polynomial currentPolynomial = {monomCoeff};
                            if (std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentPolynomial) == monomsInTopRow.end()) {
                                monomsInTopRow.push_back(currentPolynomial);
                            }
                        }
                    }
                }
            }
        }
    }
    if (level >= 6) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    for (long unsigned int l=0; l<variables.size(); l++) {
                        for (long unsigned int m=0; m<variables.size(); m++) {
                            for (long unsigned int n=0; n<variables.size(); n++) {
                                monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0], variables[l][0], variables[m][0], variables[n][0]};
                                std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(currentMonomial, verbosity, false);
                                polynomial currentPolynomial = {monomCoeff};
                                if (std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentPolynomial) == monomsInTopRow.end()) {
                                    monomsInTopRow.push_back(currentPolynomial);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if  (level >= 7) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    for (long unsigned int l=0; l<variables.size(); l++) {
                        for (long unsigned int m=0; m<variables.size(); m++) {
                            for (long unsigned int n=0; n<variables.size(); n++) {
                                for (long unsigned int o=0; o<variables.size(); o++) {
                                    monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0], variables[l][0], variables[m][0], variables[n][0], variables[o][0]};
                                    std::pair<std::complex<double>, monomial> monomCoeff = reduceMonomial(currentMonomial, verbosity, false);
                                    polynomial currentPolynomial = {monomCoeff};
                                    if (std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentPolynomial) == monomsInTopRow.end()) {
                                        monomsInTopRow.push_back(currentPolynomial);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Return the monomials
    return monomsInTopRow;

}

// Generate all moment matrices for a given level given a polynomial
std::vector<polynomialMatrix> generateAllMomentMatrices(const polynomial& functional, std::vector<polynomial> zeroCons, int level, int verbosity, std::vector<int> reductionsToIgnore={}) {

    // First get the list of all monomials used by iterating through the polynomial
    std::vector<monomial> variables;
    addSingleMonomials(variables, functional);
    for (int i=0; i<zeroCons.size(); i++) {
        addSingleMonomials(variables, zeroCons[i]);
    }

    // Generate all monomials up to the given level
    std::vector<polynomial> monomsInTopRow = generateMonomials(variables, level, verbosity);

    // Create the index array to be used for combinations
    std::vector<int> indices;
    for (int i=0; i<monomsInTopRow.size(); i++) {
        indices.push_back(i+1);
    }
    std::vector<std::vector<int>> combinations;

    // In this code we just want one big moment matrix
    combinations = {indices};

    // Each moment mat should start with 1
    monomsInTopRow.insert(monomsInTopRow.begin(), stringToPolynomial("1"));
    for (long unsigned int i=0; i<combinations.size(); i++) {
        combinations[i].insert(combinations[i].begin(), 0);
    }

    // Form the moment matrices
    std::vector<polynomialMatrix> matricesToReturn;
    for (int k=0; k<combinations.size(); k++) {

        // Get this top row combination
        std::vector<polynomial> monomsInTopRowComb;
        for (long unsigned int i=0; i<combinations[k].size(); i++) {
            monomsInTopRowComb.push_back(monomsInTopRow[combinations[k][i]]);
        }

        // Form that matrix
        polynomialMatrix newMatrix = generateFromTopRow(monomsInTopRowComb, verbosity);
        matricesToReturn.push_back(newMatrix);

    }

    // Return the moment matrix
    return matricesToReturn;

}

// Add variables from a moment matrix
void addVariables(std::vector<monomial>& variables, polynomialMatrix toAdd) {

    // Iterative through the matrix
    for (long unsigned int i=0; i<toAdd.size(); i++) {
        for (long unsigned int j=0; j<toAdd[i].size(); j++) {

            // Iterate through the polynomial
            for (long unsigned int k=0; k<toAdd[i][j].size(); k++) {
                monomial currentMonomial = toAdd[i][j][k].second;

                // Check if this monomial is already in the list
                bool found = false;
                for (long unsigned int l=0; l<variables.size(); l++) {
                    if (variables[l] == currentMonomial) {
                        found = true;
                        break;
                    }
                }

                // If not, add it
                if (!found) {
                    variables.push_back(currentMonomial);
                }

            }

        }
    }

}

// Add variables from a polynomial
void addVariables(std::vector<monomial>& variables, polynomial toAdd) {

    // Iterate through the monomial
    for (long unsigned int i=0; i<toAdd.size(); i++) {
        monomial currentMonomial = {toAdd[i].second};

        // Check if this monomial is already in the list
        bool found = false;
        for (long unsigned int j=0; j<variables.size(); j++) {
            if (variables[j] == currentMonomial) {
                found = true;
                break;
            }
        }

        // If not, add it
        if (!found) {
            variables.push_back(currentMonomial);
        }

    }

}

// Convert from a matrix location to an svec location
int matLocToVecLoc(int i, int j, int n) {
    return i*n + j - i*(i+1)/2;
}

// Convert to MOSEK form and solve
double solveMOSEK(polynomial obj, std::vector<polynomialMatrix>& psd, std::vector<polynomial> constraintsZero, int verbosity, std::vector<monomial>& variables, std::vector<std::complex<double>>& variableValues) {

    // Get the list of variables
    int oneIndex = 0;
    variables = {monomial()};
    for (int i=0; i<psd.size(); i++) {
        addVariables(variables, psd[i]);
    }
    for (int i=0; i<constraintsZero.size(); i++) {
        addVariables(variables, constraintsZero[i]);
    }
    addVariables(variables, obj);

    // Add an imaginary part for each variable
    std::vector<monomial> newVarList;
    for (int i=0; i<variables.size(); i++) {
        newVarList.push_back(variables[i]);
        newVarList.push_back(variables[i]);
    }
    variables = newVarList;
    int oneIndexImag = 1;

    // Output the variable list
    if (verbosity >= 3) {
        std::cout << "Variables:" << std::endl;
        for (int i=0; i<variables.size(); i++) {
            if (i % 2 == 0) {
                std::cout << i << " " << variables[i] << std::endl;
            } else {
                std::cout << i << " " << variables[i] << " (imag)" << std::endl;
            }
        }
        std::cout << std::endl;
    }

    // The c vector defining the objective
    std::vector<double> c(variables.size());
    for (int i=0; i<obj.size(); i++) {

        // Find the location of this variable
        monomial toFind = {obj[i].second};
        for (int j=0; j<variables.size(); j++) {
            if (variables[j] == toFind) {
                c[j] += std::real(obj[i].first);
                break;
            }
        }

    }
    auto cM = monty::new_array_ptr<double>(c);

    // The A matrix defining the equality constraints
    // f(x) = (c_r + i c_i)*(x_r + i x_i) = 0
    std::vector<int> ARows;
    std::vector<int> ACols;
    std::vector<double> AVals;
    for (int i=0; i<constraintsZero.size(); i++) {
        for (int j=0; j<constraintsZero[i].size(); j++) {

            // Find the location of this variable
            monomial toFind = {constraintsZero[i][j].second};
            int realLoc = -1;
            int imagLoc = -1;
            for (int k=0; k<variables.size(); k++) {
                if (variables[k] == toFind) {
                    realLoc = k;
                    imagLoc = k+1;
                    break;
                }
            }

            // c_r*x_r - c_i*x_i = 0
            ARows.push_back(2*i);
            ACols.push_back(realLoc);
            AVals.push_back(std::real(constraintsZero[i][j].first));
            ARows.push_back(2*i);
            ACols.push_back(imagLoc);
            AVals.push_back(-std::imag(constraintsZero[i][j].first));

            // c_r*x_i + c_i*x_r = 0
            ARows.push_back(2*i+1);
            ACols.push_back(realLoc);
            AVals.push_back(std::imag(constraintsZero[i][j].first));
            ARows.push_back(2*i+1);
            ACols.push_back(imagLoc);
            AVals.push_back(std::real(constraintsZero[i][j].first));

        }
    }
    if (verbosity >= 3) {
        std::cout << "ARows: " << ARows << std::endl;
        std::cout << "ACols: " << ACols << std::endl;
        std::cout << "AVals: " << AVals << std::endl;
    }
    auto AM = mosek::fusion::Matrix::sparse(2*constraintsZero.size(), variables.size(), monty::new_array_ptr<int>(ARows), monty::new_array_ptr<int>(ACols), monty::new_array_ptr<double>(AVals));

    // The vectors defining the PSD constraints
    std::vector<std::shared_ptr<monty::ndarray<int,1>>> indicesPSDPerMat;
    std::vector<std::shared_ptr<monty::ndarray<double,1>>> coeffsPSDPerMat;
    std::vector<std::pair<int,int>> matDims;
    for (int k=0; k<psd.size(); k++) {

        // The indices and coefficients for the svec
        int fullMatSize = 2*psd[k].size();
        int sVecSize = fullMatSize * (fullMatSize + 1) / 2;
        int imagOffset = psd[k].size();
        std::vector<int> indicesPSD(sVecSize);
        std::vector<double> coeffsPSD(sVecSize);
        for (int i=0; i<psd[k].size(); i++) {
            for (int j=i; j<psd[k][i].size(); j++) {

                // Find this in the variable list
                int realLoc = -1;
                int imagLoc = -1;
                monomial toFind = {psd[k][i][j][0].second};
                for (int k2=0; k2<variables.size(); k2++) {
                    if (variables[k2] == toFind) {
                        realLoc = k2;
                        imagLoc = k2+1;
                        break;
                    }
                }

                // Locations in the svec for the real and imaginary parts
                int realInd = matLocToVecLoc(i, j, fullMatSize);
                int realInd2 = matLocToVecLoc(i+imagOffset, j+imagOffset, fullMatSize);
                int imagInd = matLocToVecLoc(i, j+imagOffset, fullMatSize);
                int imagInd2 = matLocToVecLoc(j, i+imagOffset, fullMatSize);

                // If it's a real coeff
                if (std::imag(psd[k][i][j][0].first) == 0) {
                    double coeff = std::real(psd[k][i][j][0].first);
                    indicesPSD[realInd] = realLoc;
                    coeffsPSD[realInd] = coeff;
                    indicesPSD[realInd2] = realLoc;
                    coeffsPSD[realInd2] = coeff;
                    indicesPSD[imagInd] = imagLoc;
                    coeffsPSD[imagInd] = coeff;
                    indicesPSD[imagInd2] = imagLoc;
                    coeffsPSD[imagInd2] = -coeff;

                // If it's an imaginary coeff
                } else {
                    double coeff = std::imag(psd[k][i][j][0].first);
                    indicesPSD[realInd] = imagLoc;
                    coeffsPSD[realInd] = coeff;
                    indicesPSD[realInd2] = imagLoc;
                    coeffsPSD[realInd2] = coeff;
                    indicesPSD[imagInd] = realLoc;
                    coeffsPSD[imagInd] = -coeff;
                    indicesPSD[imagInd2] = realLoc;
                    coeffsPSD[imagInd2] = coeff;

                }

                // The off-diagonals are multiplied by sqrt(2)
                if (i != j) {
                    coeffsPSD[realInd] *= std::sqrt(2.0);
                    coeffsPSD[realInd2] *= std::sqrt(2.0);
                    coeffsPSD[imagInd] *= std::sqrt(2.0);
                    coeffsPSD[imagInd2] *= std::sqrt(2.0);
                }


            }

        }

        // Construct the dense matrix for debugging
        if (verbosity >= 3) {
            std::cout << "Indices: " << indicesPSD << std::endl;
            std::cout << "Coeffs: " << coeffsPSD << std::endl;
            Eigen::MatrixXi indMat = Eigen::MatrixXi::Zero(fullMatSize, fullMatSize);
            int nextX = 0;
            int nextY = 0;
            for (int i=0; i<sVecSize; i++) {
                indMat(nextX, nextY) = indicesPSD[i];
                nextY++;
                if (nextY == fullMatSize) {
                    nextX++;
                    nextY = nextX;
                }
            }
            std::cout << "Indices matrix: " << std::endl;
            std::cout << indMat << std::endl;
        }

        // Convert to MOSEK form
        auto indicesPSDM = monty::new_array_ptr<int>(indicesPSD);
        auto coeffsPSDM = monty::new_array_ptr<double>(coeffsPSD);

        // Add to the list
        indicesPSDPerMat.push_back(indicesPSDM);
        coeffsPSDPerMat.push_back(coeffsPSDM);
        matDims.push_back({1, sVecSize});

    }

    // Create a model
    mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});
    if (verbosity >= 3) {
        M->setLogHandler([=](const std::string & msg) {std::cout << msg << std::flush;});
    }

    // Create the main variable vector
    mosek::fusion::Variable::t xM = M->variable(variables.size());

    // The objective function
    M->objective(mosek::fusion::ObjectiveSense::Maximize, mosek::fusion::Expr::dot(cM, xM));

    // The one variable should be fixed
    M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));
    M->constraint(xM->index(oneIndexImag), mosek::fusion::Domain::equalsTo(0.0));

    // The matrix of this should be PSD
    for (int k=0; k<psd.size(); k++) {
        M->constraint(
            mosek::fusion::Expr::sum(
                mosek::fusion::Expr::reshape(
                    mosek::fusion::Expr::mulElm(
                        coeffsPSDPerMat[k], 
                        xM->pick(indicesPSDPerMat[k])
                    ),
                    matDims[k].first,
                    matDims[k].second
                ),
                0
            ), 
            mosek::fusion::Domain::inSVecPSDCone()
        );
    }

    // Linear equality constraints
    M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

    // Solve the problem
    M->solve();

    // Output the primal objective value
    double objPrimal = M->primalObjValue();

    // Get all of the variable values
    auto xMLevel = *(xM->level());
    variableValues.resize(variables.size());
    for (int i=0; i<variables.size(); i+=2) {
        variableValues[i] = std::complex<double>(xMLevel[i], xMLevel[i+1]);
    }

    // Check the eigenvalues of each moment matrix
    if (verbosity >= 2) {
        std::cout << std::endl;
        std::cout << "Objective value: " << objPrimal << std::endl;
        for (int i=0; i<psd.size(); i++) {
            std::vector<std::vector<std::complex<double>>> eigenvectors;
            std::vector<std::complex<double>> eigenvalues;
            getEigens(psd[i], variables, variableValues, eigenvectors, eigenvalues);
            double minEig = 1e10;
            for (int j=0; j<eigenvalues.size(); j++) {
                if (std::imag(eigenvalues[j]) > 1e-5) {
                    std::cout << "ERROR - Eigenvalue " << j << " has imaginary part: " << std::imag(eigenvalues[j]) << std::endl;
                }
                if (std::real(eigenvalues[j]) < minEig) {
                    minEig = std::real(eigenvalues[j]);
                }
            }
            std::cout << "Min eigenvalue: " << minEig << std::endl;
        }
    }

    // If superverbose, output all monomials
    if (verbosity >= 3) {
        std::cout << "Solution: " << std::endl;
        for (int i=0; i<variables.size(); i++) {
            std::cout << variables[i] << ": " << variableValues[i] << std::endl;
        }
        std::cout << "There are " << variables.size() << " variables." << std::endl;

        Eigen::MatrixXcd A = replaceVariables(psd[0], variables, variableValues);
        std::cout << "Moment matrix:" << std::endl;
        std::cout << psd[0] << std::endl;
        std::cout << "Moment matrix with vars replaced:" << std::endl;
        std::cout << A << std::endl;

    // If verbose, just output monomials that are in the objective
    } else if (verbosity >= 2) {
        std::cout << "Solution: " << std::endl;
        std::vector<monomial> variablesInObj = {};
        addVariables(variablesInObj, obj);
        for (int i=0; i<variablesInObj.size(); i++) {
            for (int j=0; j<variables.size(); j++) {
                if (variablesInObj[i] == variables[j]) {
                    std::cout << variablesInObj[i] << ": " << variableValues[j] << std::endl;
                    break;
                }
            }
        }
    } 

    // Return the objective
    return objPrimal;

}

// Generic entry function
int main(int argc, char* argv[]) {

    // Define the scenario
    int level = 1;
    int limbladLevel = 1;
    polynomial objective = stringToPolynomial("<Z1>");
    int verbosity = 1;
    std::string seed = "";
    polynomial limbladian = stringToPolynomial("<X1A1X1>+<Y1A1Y1>-<A1>");
    std::vector<std::string> extraMonomials;
    std::vector<std::string> extraMonomialsLim;
    std::vector<polynomial> constraintsZero;
    std::vector<int> reductionsToIgnore = {};

    // Process command-line args
    for (int i=1; i<argc; i++) {
        std::string argAsString = std::string(argv[i]);

        // Set the level
        if (argAsString == "-l") {
            level = std::stoi(argv[i+1]);
            i++;

        // Pauli X objective
        } else if (argAsString == "--objX") {
            objective = stringToPolynomial("<X1>");

        // Pauli Y objective
        } else if (argAsString == "--objY") {
            objective = stringToPolynomial("<Y1>");

        // Pauli Z objective
        } else if (argAsString == "--objZ") {
            objective = stringToPolynomial("<Z1>");

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
            pauliString += "+" + std::to_string(coeff1) + "<X1A1X1>";
            pauliString += "+" + std::to_string(coeff2) + "<Y1A1Y1>";
            pauliString += "+" + std::to_string(coeff3) + "<Z1A1Z1>";
            pauliString += "-" + std::to_string(coeff1+coeff2+coeff3) + "<A1>";
            limbladian = stringToPolynomial(pauliString);
            i += 3;

        // Two-body Limbladian TODO
        } else if (argAsString == "--two") {

            // Defining quantities
            double gamma_c = 1.1e-2;
            double gamma_h = 1e-3;
            double g = 1.6e-3;
            double T_h = 1.0;
            double T_c = 0.1;
            double delta = 0.005;
            double epsilon_h = 1.0;

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
            double chi = (4*g*g+Gamma_h*Gamma_c)*Gamma*Gamma + 4*delta*delta*Gamma_h*Gamma_c;

            double J_ss = (8*g*g*(gamma_h_plus*gamma_c_minus - gamma_h_minus*gamma_c_plus) / chi) * (epsilon_h*Gamma_c + epsilon_c*Gamma_h);

            std::cout << "ideal J_ss: " << J_ss << std::endl;

            return 0;

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
        } else if (argAsString == "-m") {
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
            std::cout << "  --two <num> <num> <num>" << std::endl;
            std::cout << "                  Use the two-body Limbladian with coeffs" << std::endl;
            std::cout << "  -l <num>        Level of the moment matrix" << std::endl;
            std::cout << "  -m <num>        Level of the moments to put in the Limbladian" << std::endl;
            std::cout << "  -e <monom>      Add an extra monomial to the top row" << std::endl;
            std::cout << "  -E <monom>      Add an extra monomial to the list of Limbladian replacements" << std::endl;
            std::cout << "  -S <str>        Seed for the random number generator" << std::endl;
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
    // TODO seg fault on level 2
    std::vector<monomial> variables = {stringToMonomial("<X1>"), stringToMonomial("<Y1>"), stringToMonomial("<Z1>")};
    std::vector<polynomial> variablesToPut = generateMonomials(variables, limbladLevel, verbosity);
    for (int i=0; i<extraMonomialsLim.size(); i++) {
        variablesToPut.push_back(stringToPolynomial(extraMonomialsLim[i]));
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
        polynomial newConstraint = limbladian;
        monomial monomToPut = variablesToPut[i][0].second;
        if (monomToPut.size() > 0) {
            for (int j=0; j<newConstraint.size(); j++) {
                for (int k=0; k<newConstraint[j].second.size(); k++) {
                    if (newConstraint[j].second[k].first == 'A') {
                        newConstraint[j].second[k] = monomToPut[0];
                        for (int l=1; l<monomToPut.size(); l++) {
                            newConstraint[j].second.insert(newConstraint[j].second.begin() + k + l, monomToPut[l]);
                        }
                    }
                }
            }
        }
        if (verbosity >= 2) {
            std::cout << std::endl;
            std::cout << "Variable to put: " << variablesToPut[i] << std::endl;
            std::cout << "New constraint: " << newConstraint << std::endl;
        }
        constraintsZero.push_back(newConstraint);
    }

    // Reduce the constraints as much as possible
    for (int i=0; i<constraintsZero.size(); i++) {
        if (verbosity >= 2) {
            std::cout << std::endl;
            std::cout << "Reducing constraint: " << constraintsZero[i] << std::endl;
        }
        constraintsZero[i] = reducePolynomial(constraintsZero[i], verbosity, true, false, true, reductionsToIgnore);
        if (verbosity >= 2) {
            std::cout << "Reduced constraint: " << constraintsZero[i] << std::endl;
        }
    }

    // Define the moment matrix
    if (verbosity >= 2) {
        std::cout << std::endl;
    }
    std::vector<polynomialMatrix> momentMatrices = generateAllMomentMatrices(objective, constraintsZero, level, verbosity, reductionsToIgnore);

    // If told to add extra to the top row
    if (extraMonomials.size() > 0) {
        std::vector<polynomial> topRow = momentMatrices[0][0];
        for (int i=0; i<extraMonomials.size(); i++) {
            polynomial extraMonomial = stringToPolynomial(extraMonomials[i]);
            for (int j=0; j<extraMonomial.size(); j++) {
                std::reverse(extraMonomial[j].second.begin(), extraMonomial[j].second.end());
            }
            topRow.push_back(extraMonomial);
            std::cout << "Added " << extraMonomial << " to the top row" << std::endl;
        }
        momentMatrices[0] = generateFromTopRow(topRow, verbosity);
    }

    // See how big the moment matrix is
    if (verbosity >= 1) {
        if (verbosity >= 2) {
            std::cout << std::endl;
        }
        std::cout << "Generated " << momentMatrices.size() << " moment matrices" << std::endl;
        int largestMomentMatrix = 0;
        for (int i=0; i<momentMatrices.size(); i++) {
            if (momentMatrices[i].size() > largestMomentMatrix) {
                largestMomentMatrix = momentMatrices[i].size();
            }
        }
        std::cout << "Largest moment matrix has size " << largestMomentMatrix << std::endl;
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
                std::cout << "Moment matrix " << i << ": " << std::endl;
                if (momentMatrices[i].size() < 10 || verbosity >= 3) {
                    std::cout << momentMatrices[i] << std::endl;
                } else {
                    std::cout << " - Has size " << momentMatrices[i].size() << ", set verbosity 3 to show" << std::endl << std::endl;
                }
            }
        }
        if (constraintsZero.size() > 0) {
            std::cout << "Zero constraints: " << std::endl;
            std::cout << constraintsZero << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
    }

    // Convert to MOSEK form and solve
    std::vector<monomial> varNames;
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
        std::cout << "Upper bound: " << upperBound << std::endl;
        std::cout << "Lower bound: " << lowerBound << std::endl;
    }

    // Exit without errors
    return 0;

}
