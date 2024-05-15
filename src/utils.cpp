#include "utils.h"
#include <iostream>

// Given a matrix and a variable list, return the matrix with the variables replaced
Eigen::MatrixXcd replaceVariables(std::vector<std::vector<Poly>>& momentMatrix, const std::vector<Mon>& variables, const std::vector<std::complex<double>>& varVals) {

    // Replace each variable with its value
    Eigen::MatrixXcd momentMatrixEigen = Eigen::MatrixXcd::Zero(momentMatrix.size(), momentMatrix.size());
    for (int i=0; i<momentMatrix.size(); i++) {
        for (int j=0; j<momentMatrix[i].size(); j++) {
            //for (int k=0; k<momentMatrix[i][j].size(); k++) {
            for (auto& term : momentMatrix[i][j]) {
                for (int l=0; l<variables.size(); l++) {
                    if (term.first == variables[l]) {
                        if (i > j) {
                            momentMatrixEigen(i, j) += term.second * varVals[l];
                        } else {
                            momentMatrixEigen(i, j) += term.second * std::conj(varVals[l]);
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
void getEigens(std::vector<std::vector<Poly>>& momentMatrix, const std::vector<Mon>& variables, const std::vector<std::complex<double>>& varVals, std::vector<std::vector<std::complex<double>>>& eigenvectors, std::vector<std::complex<double>>& eigenvalues) {

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

// Add all single order monomials from a functional to a list of variables
void addSingleMonomials(std::vector<Mon>& variables, Poly functional) {

    // Iterate through the functional
    //for (long unsigned int i=0; i<functional.size(); i++) {
    for (auto& term : functional) {

        // Iterate through the monomial
        for (long unsigned int j=0; j<term.first.size(); j++) {
            Mon currentMonomial(term.first[j]);

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
std::vector<std::vector<Poly>> generateFromTopRow(std::vector<Mon> monomsInTopRow, int verbosity) {

    // Generate all combinations of the top row
    std::vector<std::vector<Poly>> matrixToReturn(monomsInTopRow.size(), std::vector<Poly>(monomsInTopRow.size()));
    for (long unsigned int i=0; i<monomsInTopRow.size(); i++) {
        for (long unsigned int j=i; j<monomsInTopRow.size(); j++) {

            // Form the new moment
            Mon newMonomial;
            for (long unsigned int k=0; k<monomsInTopRow[j].size(); k++) {
                newMonomial.monomial.push_back(monomsInTopRow[j][k]);
            }
            for (int k=int(monomsInTopRow[i].size())-1; k>=0; k--) {
                newMonomial.monomial.push_back(monomsInTopRow[i][k]);
            }

            // Reduce the monomial
            std::pair<std::complex<double>, Mon> monomCoeff = newMonomial.reduce();

            // Set the matrix elements
            matrixToReturn[i][j] = Poly(monomCoeff.first, monomCoeff.second);
            matrixToReturn[j][i] = Poly(std::conj(monomCoeff.first), monomCoeff.second);

        }
    }

    // Return the matrix
    return matrixToReturn;

}

// Generate a moment matrix given the top row as polynomials
std::vector<std::vector<Poly>> generateFromTopRow(std::vector<Poly> monomsInTopRow, int verbosity) {

    // Generate all combinations of the top row
    std::vector<std::vector<Poly>> matrixToReturn = std::vector<std::vector<Poly>>(monomsInTopRow.size(), std::vector<Poly>(monomsInTopRow.size()));
    for (long unsigned int i=0; i<monomsInTopRow.size(); i++) {
        for (long unsigned int j=i; j<monomsInTopRow.size(); j++) {

            // Form the new polynomial
            Poly newPolynomial = monomsInTopRow[j] * monomsInTopRow[i].conj();
            //for (long unsigned int k=0; k<monomsInTopRow[i].size(); k++) {
                //for (long unsigned int l=0; l<monomsInTopRow[j].size(); l++) {

                    //// Form the new monomial
                    //Mon newMonomial= monomsInTopRow[i][k].second * monomsInTopRow[j][l].second;

                    //// Reduce the monomial
                    //std::pair<std::complex<double>, Mon> monomCoeff = newMonomial.reduce();
                    //newMonomial = monomCoeff.second;

                    //// Add to the polynomial
                    //std::complex<double> newCoefficient = monomCoeff.first * monomsInTopRow[i][k].first * std::conj(monomsInTopRow[j][l].first);
                    //newPolynomial += Poly(monomCoeff.first, newMonomial);

                //}
            //}

            // Set the matrix elements
            matrixToReturn[i][j] = newPolynomial;
            matrixToReturn[j][i] = newPolynomial.conj();

        }
    }

    // Return the matrix
    return matrixToReturn;

}

// Generate a list of monomials of a certain level
std::vector<Poly> generateMonomials(std::vector<Mon> variables, int level, int verbosity) {

    // Generate all monomials up to the given level
    std::vector<Poly> monomsInTopRow = {};
    if (level >= 1) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            Mon currentMonomial = variables[i];
            std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
            Poly currentPolynomial(monomCoeff.second);
            if (std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentPolynomial) == monomsInTopRow.end()) {
                monomsInTopRow.push_back(currentPolynomial);
            }
        }
    }
    if (level >= 2) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                Mon currentMonomial = variables[i] * variables[j];
                std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                Poly currentPolynomial(monomCoeff.second);
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
                    Mon currentMonomial = variables[i] * variables[j] * variables[k];
                    std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                    Poly currentPolynomial(monomCoeff.second);
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
                        Mon currentMonomial = variables[i] * variables[j] * variables[k] * variables[l];
                        std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                        Poly currentPolynomial(monomCoeff.second);
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
                            Mon currentMonomial = variables[i] * variables[j] * variables[k] * variables[l] * variables[m];
                            std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                            Poly currentPolynomial(monomCoeff.second);
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
                                Mon currentMonomial = variables[i] * variables[j] * variables[k] * variables[l] * variables[m] * variables[n];
                                std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                                Poly currentPolynomial(monomCoeff.second);
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
                                    Mon currentMonomial = variables[i] * variables[j] * variables[k] * variables[l] * variables[m] * variables[n] * variables[o];
                                    std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                                    Poly currentPolynomial(monomCoeff.second);
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

    // Remove the 1 if there is
    for (long unsigned int i=0; i<monomsInTopRow.size(); i++) {
        if (monomsInTopRow[i].size() == 1 && monomsInTopRow[i].getKey() == Mon()) {
            monomsInTopRow.erase(monomsInTopRow.begin()+i);
            break;
        }
    }

    // Return the monomials
    return monomsInTopRow;

}

// Generate all moment matrices for a given level given a polynomial
std::vector<std::vector<std::vector<Poly>>> generateAllMomentMatrices(const Poly& functional, std::vector<Poly> zeroCons, int level, int verbosity, std::vector<int> reductionsToIgnore={}) {

    // First get the list of all monomials used by iterating through the polynomial
    std::vector<Mon> variables;
    addSingleMonomials(variables, functional);
    for (int i=0; i<zeroCons.size(); i++) {
        addSingleMonomials(variables, zeroCons[i]);
    }

    // Generate all monomials up to the given level
    std::vector<Poly> monomsInTopRow = generateMonomials(variables, level, verbosity);

    // Create the index array to be used for combinations
    std::vector<int> indices;
    for (int i=0; i<monomsInTopRow.size(); i++) {
        indices.push_back(i+1);
    }
    std::vector<std::vector<int>> combinations;

    // In this code we just want one big moment matrix
    combinations = {indices};

    // Each moment mat should start with 1
    monomsInTopRow.insert(monomsInTopRow.begin(), Poly(1));
    for (long unsigned int i=0; i<combinations.size(); i++) {
        combinations[i].insert(combinations[i].begin(), 0);
    }

    // Form the moment matrices
    std::vector<std::vector<std::vector<Poly>>> matricesToReturn;
    for (int k=0; k<combinations.size(); k++) {

        // Get this top row combination
        std::vector<Poly> monomsInTopRowComb;
        for (long unsigned int i=0; i<combinations[k].size(); i++) {
            monomsInTopRowComb.push_back(monomsInTopRow[combinations[k][i]]);
        }

        // Form that matrix
        std::vector<std::vector<Poly>> newMatrix = generateFromTopRow(monomsInTopRowComb, verbosity);
        matricesToReturn.push_back(newMatrix);

    }

    // Return the moment matrix
    return matricesToReturn;

}

// Add variables from a moment matrix
void addVariables(std::vector<Mon>& variables, std::vector<std::vector<Poly>> toAdd) {

    // Iterative through the matrix
    for (long unsigned int i=0; i<toAdd.size(); i++) {
        for (long unsigned int j=0; j<toAdd[i].size(); j++) {

            // Iterate through the polynomial
            for (auto& term : toAdd[i][j]) {
                Mon currentMonomial(term.first);

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
void addVariables(std::vector<Mon>& variables, Poly toAdd) {

    // Iterate through the polynomial
    //for (long unsigned int i=0; i<toAdd.size(); i++) {
    for (auto& term : toAdd) {
        Mon currentMonomial(term.first);

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

// Take the trace of a matrix, assuming it's real
double tr(Eigen::MatrixXcd A) {
    std::complex<double> trace = A.trace();
    if (std::abs(std::imag(trace)) > 1e-5) {
        std::cout << "WARNING - trace of matrix has non-zero imaginary part" << std::endl;
    }
    return std::real(trace);
}

