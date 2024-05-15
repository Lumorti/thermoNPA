#include "poly.h"
#include <algorithm>

// Useful constants
const std::complex<double> imag(0, 1);
const double zeroTol = 1e-13;

// If initialized from nothing
Poly::Poly() {}

// If initialized from just a coeff
Poly::Poly(std::complex<double> coeff) {
    polynomial[Mon()] = coeff;
}

// If initialized from just a monomial
Poly::Poly(Mon mon) {
    polynomial[mon] = 1;
}

// If initialized from a map
Poly::Poly (std::unordered_map<Mon, std::complex<double>> poly) {
    polynomial = polynomial;
}

// If initialized from a single coeff + monomial
Poly::Poly(std::pair<std::complex<double>, Mon> pair) {
    polynomial[pair.second] = pair.first;
}

// If initialized from a coeff and a monomial
Poly::Poly(std::complex<double> coeff, Mon mon) {
    polynomial[mon] = coeff;
}

// If initialized from a single part of monomial
Poly::Poly(std::pair<char, int> part) {
    polynomial[Mon(part)] = 1;
}

// If initialized from a coeff and a string
Poly::Poly(std::complex<double> coeff, std::string asString) {
    polynomial[Mon(asString)] = coeff;
}

// If initialized from a string
Poly::Poly(std::string asString) {

    // Remove all spaces and *
    std::string newString = "";
    for (size_t i=0; i<asString.size(); i++) {
        if (asString[i] != ' ' && asString[i] != '*') {
            newString += asString[i];
        }
    }
    asString = newString;

    // Iterate through the string
    std::unordered_map<Mon, std::complex<double>> toReturn;
    std::string currentCoefficient;
    std::string currentMonomial;
    for (size_t i=0; i<asString.size(); i++) {

        // If finished a monomial
        if (i > 0 && asString[i] == '>') {

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
            if (toReturn.find(Mon(currentMonomial)) == toReturn.end()) {
                toReturn[Mon(currentMonomial)] = std::stod(currentCoefficient);
            } else {
                toReturn[Mon(currentMonomial)] += std::stod(currentCoefficient);
            }
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
        if (toReturn.find(Mon(currentMonomial)) == toReturn.end()) {
            toReturn[Mon(currentMonomial)] = std::stod(currentCoefficient);
        } else {
            toReturn[Mon(currentMonomial)] += std::stod(currentCoefficient);
        }
    }

    // Set the polynomial
    polynomial = toReturn;

}

// When assigning manually with a vector
Poly& Poly::operator=(const std::unordered_map<Mon, std::complex<double>>& other) {
    polynomial = other;
    return *this;
}

// Checking equality of two polynomials
bool Poly::operator==(const Poly& other) {
    return polynomial == other.polynomial;
}

// When summing two polynomials
Poly Poly::operator+(const Poly& other) {
    Poly toReturn;
    toReturn.polynomial = polynomial;
    //for (size_t i=0; i<other.polynomial.size(); i++) {
        //toReturn.polynomial.push_back(other.polynomial[i]);
    //}
    for (auto& term : other.polynomial) {
        if (toReturn.polynomial.find(term.first) == toReturn.polynomial.end()) {
            toReturn.polynomial[term.first] = term.second;
        } else {
            toReturn.polynomial[term.first] += term.second;
        }
    }
    return toReturn;
}

// When summing in-place
Poly& Poly::operator+=(const Poly& other) {
    //for (size_t i=0; i<other.polynomial.size(); i++) {
        //polynomial.push_back(other.polynomial[i]);
    //}
    for (auto& term : other.polynomial) {
        if (polynomial.find(term.first) == polynomial.end()) {
            polynomial[term.first] = term.second;
        } else {
            polynomial[term.first] += term.second;
        }
    }
    return *this;
}

// When subtracting two polynomials
Poly Poly::operator-(const Poly& other) {
    Poly toReturn;
    toReturn.polynomial = polynomial;
    //for (size_t i=0; i<other.polynomial.size(); i++) {
        //toReturn.polynomial.push_back(std::make_pair(-other.polynomial[i].first, other.polynomial[i].second));
    //}
    for (auto& term : other.polynomial) {
        if (toReturn.polynomial.find(term.first) == toReturn.polynomial.end()) {
            toReturn.polynomial[term.first] = -term.second;
        } else {
            toReturn.polynomial[term.first] -= term.second;
        }
    }
    return toReturn;
}

// When subtracting in-place
Poly& Poly::operator-=(const Poly& other) {
    //for (size_t i=0; i<other.polynomial.size(); i++) {
        //polynomial.push_back(std::make_pair(-other.polynomial[i].first, other.polynomial[i].second));
    //}
    for (auto& term : other.polynomial) {
        if (polynomial.find(term.first) == polynomial.end()) {
            polynomial[term.first] = -term.second;
        } else {
            polynomial[term.first] -= term.second;
        }
    }
    return *this;
}

// When multiplying two polynomials
Poly Poly::operator*(const Poly& other) const {
    Poly toReturn;
    //for (size_t i=0; i<polynomial.size(); i++) {
        //for (size_t j=0; j<other.polynomial.size(); j++) {
            //std::complex<double> newCoefficient = polynomial[i].first * other.polynomial[j].first;
            //Mon newMonomial = polynomial[i].second * other.polynomial[j].second;
            //toReturn += Poly(newCoefficient, newMonomial);
        //}
    //}
    for (const auto& term1 : polynomial) {
        for (const auto& term2 : other.polynomial) {
            std::complex<double> newCoefficient = term1.second * term2.second;
            Mon newMonomial = term1.first * term2.first;
            if (toReturn.polynomial.find(newMonomial) == toReturn.polynomial.end()) {
                toReturn.polynomial[newMonomial] = newCoefficient;
            } else {
                toReturn.polynomial[newMonomial] += newCoefficient;
            }
        }
    }
    return toReturn;
}

// When multiplying by a constant
Poly Poly::operator*(const std::complex<double>& other) {
    Poly toReturn;
    //for (size_t i=0; i<polynomial.size(); i++) {
        //toReturn.polynomial.push_back(std::make_pair(polynomial[i].first * other, polynomial[i].second));
    //}
    for (auto& term : polynomial) {
        toReturn.polynomial[term.first] = term.second * other;
    }
    return toReturn;
}

// When multiplying by a constant in-place
Poly& Poly::operator*=(const std::complex<double>& other) {
    //for (size_t i=0; i<polynomial.size(); i++) {
        //polynomial[i].first *= other;
    //}
    for (auto& term : polynomial) {
        term.second *= other;
    }
    return *this;
}

// When dividing by a constant
Poly Poly::operator/(const std::complex<double>& other) {
    Poly toReturn;
    //for (size_t i=0; i<polynomial.size(); i++) {
        //toReturn.polynomial.push_back(std::make_pair(polynomial[i].first / other, polynomial[i].second));
    //}
    for (auto& term : polynomial) {
        toReturn.polynomial[term.first] = term.second / other;
    }
    return toReturn;
}

// When dividing by a constant in-place
Poly& Poly::operator/=(const std::complex<double>& other) {
    //for (size_t i=0; i<polynomial.size(); i++) {
        //polynomial[i].first /= other;
    //}
    for (auto& term : polynomial) {
        term.second /= other;
    }
    return *this;
}

// Size of the polynomial
const size_t Poly::size() const {
    return polynomial.size();
}

// Allow bracket access
const std::complex<double> Poly::operator[](Mon mon) const {
    auto it = polynomial.find(mon);
    return it != polynomial.end() ? it->second : std::complex<double>(0,0);
}

// Allow bracket access
std::complex<double>& Poly::operator[](Mon mon) {
    return polynomial[mon];
}

// Self-negative
Poly Poly::operator-() {
    Poly toReturn = *this;
    //for (size_t i=0; i<toReturn.size(); i++) {
        //toReturn.polynomial[i].first = -toReturn.polynomial[i].first;
    //}
    for (auto& term : toReturn.polynomial) {
        term.second = -term.second;
    }
    return toReturn;
}

// Self-conjugate
Poly Poly::conj() {
    Poly toReturn = *this;
    //for (size_t i=0; i<toReturn.size(); i++) {
        //toReturn.polynomial[i].first = std::conj(toReturn.polynomial[i].first);
    //}
    for (auto& term : toReturn.polynomial) {
        term.second = std::conj(term.second);
    }
    return toReturn;
}

// Given a polynomial, reduce each monomial and combine
void Poly::reduce() {

    // Apply the reduction to each monomial
    //for (int j=0; j<size(); j++) {
        //std::pair<std::complex<double>, Mon> reducedMonomial = polynomial[j].second.reduce();
        //polynomial[j].first *= reducedMonomial.first;
        //polynomial[j].second = reducedMonomial.second;
    //}

    Poly toReturn;
    for (auto& term : polynomial) {
        std::pair<std::complex<double>, Mon> reducedMonomial = term.first.reduce();
        //term.first = reducedMonomial.second;
        //tfirsterm.second *= reducedMonomial.first;
        if (toReturn.polynomial.find(reducedMonomial.second) == toReturn.polynomial.end()) {
            toReturn.polynomial[reducedMonomial.second] = term.second * reducedMonomial.first;
        } else {
            toReturn.polynomial[reducedMonomial.second] += term.second * reducedMonomial.first;
        }
    }
    *this = toReturn;
}

// Pretty print
std::ostream& operator<<(std::ostream& os, const Poly& p) {

    // Check if it's zero
    if (p.size() == 0 || (p.size() == 1 && p[Mon()] == std::complex<double>(0,0))) {
        os << "0";
        return os;
    }

    // Iterate through the polynomial
    for (auto& term : p.polynomial) {
        double realPart = term.second.real();
        double imagPart = term.second.imag();
        if (imagPart == 0) {
            if (realPart == std::complex<double>(-1, 0)) {
                os << "-" << term.first;
            } else if (realPart == std::complex<double>(1, 0)) {
                if (term.first == Mon()) {
                    os << term.second;
                } else {
                    os << "+" << term.second;
                }
            } else if (term.first.size() == 0 && realPart > 0) {
                if (term.first == Mon()) {
                    os << realPart;
                } else {
                    os << "+" << realPart;
                }
            } else if (term.first.size() == 0 && realPart < 0) {
                os << realPart;
            } else if (realPart > 0) {
                if (term.first == Mon()) {
                    os << realPart;
                } else {
                    os << "+" << realPart;
                }
            } else if (realPart != 0) {
                os << realPart;
            }
        } else {
            os << term.second << term.first;
        }
    }

    // Return the stream
    return os;

}

// Replace a monomial by a polynomial
void Poly::replace(std::pair<char,int> mon, Poly replacement) {

    // Iterate through the polynomial
    Poly pNew;
    //for (int i=0; i<size(); i++) {
        //Poly newTerm = Poly(polynomial[i].first);

        //// For each of the terms in the monomial
        //for (int j=0; j<polynomial[i].second.size(); j++) {

            //// If we have a match
            //Poly toMultiply(polynomial[i].second[j]);
            //if (polynomial[i].second[j] == mon) {
                //toMultiply = replacement;
            //} else if (polynomial[i].second[j].first == mon.first && mon.second == 0) {
                //toMultiply = replacement;
                //for (int k=0; k<toMultiply.size(); k++) {
                    //for (int l=0; l<toMultiply[k].second.size(); l++) {
                        //toMultiply[k].second[l].second = polynomial[i].second[j].second;
                    //}
                //}
            //}

            //// Perform the multiplication and reduce
            //newTerm = newTerm * toMultiply;
            //newTerm.reduce();


        //}

        //// Add this term to the polynomial
        //pNew += newTerm;

    //}

    // For each term in the polynomial
    for (auto& term : polynomial) {
        Poly newTerm = Poly(term.second);

        // For each of the sections in this monomial
        for (int j=0; j<term.first.size(); j++) {

            // Determine the monomial to multiply
            Poly toMultiply(term.first[j]);
            if (term.first[j] == mon) {
                toMultiply = replacement;
            } else if (term.first[j].first == mon.first && mon.second == 0) {
                toMultiply = replacement;
                // TODO
            }

            // Perform the multiplication
            newTerm = newTerm * toMultiply;

        }

        // Reduce the new term and add it to the polynomial
        newTerm.reduce();
        pNew += newTerm;

    }

    // Set the new polynomial
    *this = pNew;

}

// Replace a monomial by a polynomial
Poly Poly::replaced(std::pair<char,int> mon, Poly replacement) {
    Poly toReturn = *this;
    toReturn.replace(mon, replacement);
    return toReturn;
}

// Replace a monomial by a polynomial
void Poly::replace(Mon mon, Poly replacement) {

    // Iterate through the polynomial
    Poly pNew;
    //for (int i=0; i<size(); i++) {
        //if (polynomial[i].second == mon) {
            //pNew += replacement * polynomial[i].first;
            //continue;
        //} else {
            //pNew += polynomial[i];
        //}
    //}
    for (auto& term : polynomial) {
        if (term.first == mon) {
            pNew += replacement * term.second;
        } else {
            pNew += term.second;
        }
    }

    // Set the new polynomial
    *this = pNew;

}

// Replace a monomial by a polynomial
Poly Poly::replaced(Mon mon, Poly replacement) {
    Poly toReturn = *this;
    toReturn.replace(mon, replacement);
    return toReturn;
}

// Reduce a polynomial using plus/minus operators to Paulis
void Poly::convertToPaulis() {
    Poly toReplaceP = Poly(0.5, "<X0>") + Poly(0.5*imag,  "<Y0>");
    Poly toReplaceM = Poly(0.5, "<X0>") + Poly(-0.5*imag, "<Y0>");
    replace(std::make_pair('P', 0), toReplaceP);
    replace(std::make_pair('M', 0), toReplaceM);
    reduce();
}

// Convert a polynomial to Paulis
Poly Poly::convertedToPaulis() {
    Poly toReturn = *this;
    toReturn.convertToPaulis();
    return toReturn;
}

// When printing a polynomial matrix, print it as a string
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<Poly>>& m) {

    // Determine the maximum width of each column
    std::vector<int> columnWidths(m[0].size(), 0);
    for (size_t i=0; i<m.size(); i++) {
        for (size_t j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            columnWidths[j] = std::max(columnWidths[j], sizeWhenWritten);
        }
    }

    // Iterate through the matrix
    for (size_t i=0; i<m.size(); i++) {
        for (size_t j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            for (size_t k=0; k<columnWidths[j]-sizeWhenWritten; k++) {
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
std::ostream& operator<<(std::ostream& os, const std::vector<Poly>& v) {

    // Iterate through the vector
    for (size_t i=0; i<v.size(); i++) {
        os << v[i] << std::endl;
    }

    // Return the stream
    return os;

}

// Get the list of monomials
std::vector<Mon> Poly::monomials() {

    // Get the list of monomials
    std::vector<Mon> toReturn;
    for (auto& term : polynomial) {
        toReturn.push_back(term.first);
    }
    return toReturn;

}

// Get the polynomial rearranged to equal a certain monomial
// e.g. X1 + 2*Y1 = 0 -> X1 = -2*Y1
Poly Poly::rearranged(Mon mon) {

    // Find this monomial
    Poly returnPoly;
    std::complex<double> coeff = 1;
    for (auto& term : polynomial) {
        if (term.first == mon) {
            coeff = term.second;
        } else {
            returnPoly += Poly(term.second, term.first);
        }
    }

    // Divide by the coefficient
    returnPoly /= -coeff;
    return returnPoly;

}

// Check if a polynomial has a monomial
bool Poly::hasMonomial(Mon mon) {
    return polynomial.find(mon) != polynomial.end();
}

// If asking for just the first monomial
std::complex<double>& Poly::getValue() {
    return polynomial.begin()->second;
}
Mon Poly::getKey() {
    return polynomial.begin()->first;
}


