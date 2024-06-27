#include "optGurobi.h"
#include "utils.h"
#include "printing.h"
#include <gurobi_c++.h>

// Convert to Gurobi form and solve
std::pair<double,double> boundGurobi(Poly obj, std::vector<Poly> constraintsZero, int verbosity) {

    // Get the list of variables
    std::set<Mon> variableSet;
    variableSet.insert(Mon());
    for (size_t i=0; i<constraintsZero.size(); i++) {
        addVariables(variableSet, constraintsZero[i]);
    }
    addVariables(variableSet, obj);
    std::vector<Mon> variables = toVector(variableSet);

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.start();
    GRBModel model = GRBModel(env);
    if (verbosity < 3) {
        model.set(GRB_IntParam_OutputFlag, 0);
    }

    // Gurobi variable array
    std::vector<GRBVar> gurobiVars(variables.size());
    for (size_t i=0; i<variables.size(); i++) {
        gurobiVars[i] = model.addVar(-1.0, 1.0, 0.0, GRB_CONTINUOUS, "x" + std::to_string(i));
    }

    // Cache the variable locations
    std::map<Mon, int> variableLocs;
    for (size_t i=0; i<variables.size(); i++) {
        variableLocs[variables[i]] = i;
    }

    // Add constraints
    for (size_t i=0; i<constraintsZero.size(); i++) {
        GRBLinExpr exprReal;
        GRBLinExpr exprImag;
        for (auto& term : constraintsZero[i]) {
            if (std::abs(std::real(term.second)) > 1e-14) {
                exprReal += std::real(term.second) * gurobiVars[variableLocs[term.first]];
            }
            if (std::abs(std::imag(term.second)) > 1e-14) {
                exprImag += std::imag(term.second) * gurobiVars[variableLocs[term.first]];
            }
        }
        model.addConstr(exprReal == 0, "c" + std::to_string(i));
        model.addConstr(exprImag == 0, "c" + std::to_string(i));
    }
    model.addConstr(gurobiVars[variableLocs[Mon()]] == 1, "c" + std::to_string(constraintsZero.size()));

    // Set objective
    GRBLinExpr expr;
    for (auto& term : obj) {
        if (std::abs(std::real(term.second)) > 1e-14) {
            expr += std::real(term.second) * gurobiVars[variableLocs[term.first]];
        }
    }

    // Set the method (-1 is automatic, 2 is barrier)
    model.set(GRB_IntParam_Method, -1);

    // Maximize objective
    model.setObjective(expr, GRB_MAXIMIZE);
    model.optimize();
    double maxVal = model.get(GRB_DoubleAttr_ObjVal);

    // Minimize objective
    model.setObjective(expr, GRB_MINIMIZE);
    model.optimize();
    double minVal = model.get(GRB_DoubleAttr_ObjVal);

    // Return the bounds
    return {minVal, maxVal};

}


