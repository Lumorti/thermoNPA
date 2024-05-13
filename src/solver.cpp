#include "solver.h"
#include "utils.h"
#include "printing.h"
#include <iostream>

// Import MOSEK
#include "fusion.h"

// Convert to MOSEK form and solve
double solveMOSEK(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero, int verbosity, std::vector<Mon>& variables, std::vector<std::complex<double>>& variableValues) {

    // Get the list of variables
    int oneIndex = 0;
    variables = {Mon()};
    for (int i=0; i<psd.size(); i++) {
        addVariables(variables, psd[i]);
    }
    for (int i=0; i<constraintsZero.size(); i++) {
        addVariables(variables, constraintsZero[i]);
    }
    addVariables(variables, obj);

    // Add an imaginary part for each variable
    std::vector<Mon> newVarList;
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
        Mon toFind(obj[i].second);
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
            Mon toFind(constraintsZero[i][j].second);
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
            ACols.push_back(imagLoc);
            AVals.push_back(std::real(constraintsZero[i][j].first));
            ARows.push_back(2*i+1);
            ACols.push_back(realLoc);
            AVals.push_back(std::imag(constraintsZero[i][j].first));

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
                Mon toFind(psd[k][i][j][0].second);
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
            Eigen::MatrixXi indMat = Eigen::MatrixXi::Zero(fullMatSize, fullMatSize);
            Eigen::MatrixXd coeffMat = Eigen::MatrixXd::Zero(fullMatSize, fullMatSize);
            int nextX = 0;
            int nextY = 0;
            for (int i=0; i<sVecSize; i++) {
                indMat(nextX, nextY) = indicesPSD[i];
                coeffMat(nextX, nextY) = coeffsPSD[i];
                nextY++;
                if (nextY == fullMatSize) {
                    nextX++;
                    nextY = nextX;
                }
            }
            std::cout << "Indices matrix: " << std::endl;
            std::cout << indMat << std::endl;
            std::cout << "Coeffs matrix: " << std::endl;
            std::cout << coeffMat << std::endl;
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
    mosek::fusion::Variable::t xM = M->variable(variables.size(), mosek::fusion::Domain::inRange(-1.0, 1.0));

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
        for (int i=0; i<variables.size(); i+=2) {
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
        std::vector<Mon> variablesInObj = {};
        addVariables(variablesInObj, obj);
        for (int i=0; i<variablesInObj.size(); i++) {
            for (int j=0; j<variables.size(); j+=2) {
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

