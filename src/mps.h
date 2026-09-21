#pragma once

// A minimal left-canonical MPS, used as the ground truth when sampling on
// systems far too large to hold a full 2^n x 2^n density matrix. Files in this
// format are produced by src/dmrg.py.
//
// Because the tensors are left-canonical, everything to the left of a Pauli
// string's support contracts to the identity for free, and everything to the
// right contracts to a precomputed environment. Only the span of the support
// is actually swept, so short-range operators are cheap regardless of n.

#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

class MPS {
public:

    int numSites = 0;
    int physDim = 2;
    double energy = 0.0;

    // tensors[i][s] is the chiL[i] x chiR[i] matrix for physical index s
    std::vector<std::vector<Eigen::MatrixXcd>> tensors;
    std::vector<int> chiL;
    std::vector<int> chiR;

    // rightEnv[i] is the contraction of everything strictly right of site i
    std::vector<Eigen::MatrixXcd> rightEnv;
    double normSq = 1.0;

    // Read the next whitespace-separated token, skipping '#' comment lines
    static std::string nextToken(std::ifstream& f) {
        std::string tok;
        while (f >> tok) {
            if (!tok.empty() && tok[0] == '#') {
                f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                continue;
            }
            return tok;
        }
        return "";
    }

    // Load an MPS from the plain-text format written by src/dmrg.py
    void load(const std::string& filename, int verbosity = 0) {

        std::ifstream f(filename);
        if (!f.is_open()) {
            throw std::runtime_error("Could not open MPS file: " + filename);
        }

        numSites = std::stoi(nextToken(f));
        physDim = std::stoi(nextToken(f));
        energy = std::stod(nextToken(f));

        tensors.resize(numSites);
        chiL.resize(numSites);
        chiR.resize(numSites);
        for (int i = 0; i < numSites; i++) {
            chiL[i] = std::stoi(nextToken(f));
            chiR[i] = std::stoi(nextToken(f));
            tensors[i].assign(physDim, Eigen::MatrixXcd::Zero(chiL[i], chiR[i]));

            // Values are ordered (l, s, r), matching numpy's C ordering
            for (int l = 0; l < chiL[i]; l++) {
                for (int s = 0; s < physDim; s++) {
                    for (int r = 0; r < chiR[i]; r++) {
                        double re = std::stod(nextToken(f));
                        double im = std::stod(nextToken(f));
                        tensors[i][s](l, r) = std::complex<double>(re, im);
                    }
                }
            }

        }
        f.close();

        // The bonds have to line up, otherwise the sweep below is nonsense
        for (int i = 0; i + 1 < numSites; i++) {
            if (chiR[i] != chiL[i+1]) {
                throw std::runtime_error("MPS bond mismatch between sites "
                                         + std::to_string(i) + " and " + std::to_string(i+1));
            }
        }

        buildEnvironments();

        if (verbosity >= 1) {
            int maxChi = 0;
            for (int i = 0; i < numSites; i++) {
                maxChi = std::max(maxChi, chiR[i]);
            }
            std::cout << "Loaded MPS with " << numSites << " sites, max bond dimension "
                      << maxChi << ", energy " << energy << std::endl;
            std::cout << "MPS norm^2 = " << normSq
                      << ", left-canonical deviation = " << canonicalError() << std::endl;
        }

    }

    // Largest deviation from sum_s A_s^dag A_s = I, which the sweep relies on
    double canonicalError() const {
        double worst = 0.0;
        for (int i = 0; i < numSites; i++) {
            Eigen::MatrixXcd gram = Eigen::MatrixXcd::Zero(chiR[i], chiR[i]);
            for (int s = 0; s < physDim; s++) {
                gram += tensors[i][s].adjoint() * tensors[i][s];
            }
            gram -= Eigen::MatrixXcd::Identity(chiR[i], chiR[i]);
            worst = std::max(worst, gram.cwiseAbs().maxCoeff());
        }
        return worst;
    }

    // Contract everything to the right of each site, once, up front
    void buildEnvironments() {

        rightEnv.assign(numSites, Eigen::MatrixXcd());
        rightEnv[numSites-1] = Eigen::MatrixXcd::Identity(chiR[numSites-1], chiR[numSites-1]);
        for (int i = numSites-1; i >= 1; i--) {
            Eigen::MatrixXcd next = Eigen::MatrixXcd::Zero(chiL[i], chiL[i]);
            for (int s = 0; s < physDim; s++) {
                next += tensors[i][s].conjugate() * rightEnv[i] * tensors[i][s].transpose();
            }
            rightEnv[i-1] = next;
        }

        // Folding in site 0 as well gives the norm
        Eigen::MatrixXcd full = Eigen::MatrixXcd::Zero(chiL[0], chiL[0]);
        for (int s = 0; s < physDim; s++) {
            full += tensors[0][s].conjugate() * rightEnv[0] * tensors[0][s].transpose();
        }
        normSq = std::real(full(0, 0));

    }

    // Matrix element <s|P|t> of a single-site Pauli
    static std::complex<double> pauliElem(char pauli, int s, int t) {
        switch (pauli) {
            case 'X': return (s != t) ? std::complex<double>(1, 0) : std::complex<double>(0, 0);
            case 'Y': if (s == 0 && t == 1) return std::complex<double>(0, -1);
                      if (s == 1 && t == 0) return std::complex<double>(0, 1);
                      return std::complex<double>(0, 0);
            case 'Z': if (s != t) return std::complex<double>(0, 0);
                      return (s == 0) ? std::complex<double>(1, 0) : std::complex<double>(-1, 0);
            default:  return (s == t) ? std::complex<double>(1, 0) : std::complex<double>(0, 0);
        }
    }

    // Expectation value of a Pauli string, given as (letter, one-based site) pairs
    double expectation(const std::vector<std::pair<char, int>>& mon) const {

        if (mon.size() == 0) {
            return 1.0;
        }

        // Collect the support, converting to zero-based indices
        std::map<int, char> ops;
        for (const auto& part : mon) {
            int site = part.second - 1;
            if (site < 0 || site >= numSites) {
                throw std::runtime_error("Pauli acts on site " + std::to_string(part.second)
                                         + ", outside the " + std::to_string(numSites) + "-site MPS");
            }
            ops[site] = part.first;
        }
        int lo = ops.begin()->first;
        int hi = ops.rbegin()->first;

        // Left-canonical form means everything before the support is the identity
        Eigen::MatrixXcd left = Eigen::MatrixXcd::Identity(chiL[lo], chiL[lo]);
        for (int i = lo; i <= hi; i++) {
            auto it = ops.find(i);
            char pauli = (it == ops.end()) ? 'I' : it->second;
            Eigen::MatrixXcd next = Eigen::MatrixXcd::Zero(chiR[i], chiR[i]);
            for (int s = 0; s < physDim; s++) {
                for (int t = 0; t < physDim; t++) {
                    std::complex<double> coeff = pauliElem(pauli, s, t);
                    if (coeff == std::complex<double>(0, 0)) {
                        continue;
                    }
                    next += coeff * (tensors[i][s].adjoint() * left * tensors[i][t]);
                }
            }
            left = next;
        }

        // And everything after it is the precomputed environment
        std::complex<double> val = (left.array() * rightEnv[hi].array()).sum();
        return std::real(val) / normSq;

    }

};
