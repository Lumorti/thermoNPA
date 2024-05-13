#include "poly.h"
#include <algorithm>

// Useful constants
const std::complex<double> imag(0, 1);
const double zeroTol = 1e-13;

// If initialized from nothing
Poly::Poly() {
    polynomial = {};
}

// If initialized from just a coeff
Poly::Poly(std::complex<double> coeff) {
    polynomial = {{coeff, Mon()}};
}

// If initialized from just a monomial
Poly::Poly(Mon mon) {
    polynomial = {{1, mon}};
}

// If initialized from a vector
Poly::Poly (std::vector<std::pair<std::complex<double>, Mon>> polynomial) {
    polynomial = polynomial;
}

// If initialized from a single coeff + monomial
Poly::Poly(std::pair<std::complex<double>, Mon> pair) {
    polynomial = {pair};
}

// If initialized from a coeff and a monomial
Poly::Poly(std::complex<double> coeff, Mon mon) {
    polynomial = {{coeff, mon}};
}

// If initialized from a single part of monomial
Poly::Poly(std::pair<char, int> part) {
    polynomial = {{1, Mon(part)}};
}

// If initialized from a coeff and a string
Poly::Poly(std::complex<double> coeff, std::string asString) {
    polynomial = {{coeff, Mon(asString)}};
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
    std::vector<std::pair<std::complex<double>, Mon>> toReturn;
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
            toReturn.push_back(std::make_pair(std::stod(currentCoefficient), Mon(currentMonomial)));
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
        toReturn.push_back(std::make_pair(std::stod(currentCoefficient), Mon(currentMonomial)));
    }

    // Set the polynomial
    polynomial = toReturn;

}

// When assigning manually with a vector
Poly& Poly::operator=(const std::vector<std::pair<std::complex<double>, Mon>>& other) {
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
    for (size_t i=0; i<other.polynomial.size(); i++) {
        toReturn.polynomial.push_back(other.polynomial[i]);
    }
    return toReturn;
}

// When summing in-place
Poly& Poly::operator+=(const Poly& other) {
    for (size_t i=0; i<other.polynomial.size(); i++) {
        polynomial.push_back(other.polynomial[i]);
    }
    return *this;
}

// When subtracting two polynomials
Poly Poly::operator-(const Poly& other) {
    Poly toReturn;
    toReturn.polynomial = polynomial;
    for (size_t i=0; i<other.polynomial.size(); i++) {
        toReturn.polynomial.push_back(std::make_pair(-other.polynomial[i].first, other.polynomial[i].second));
    }
    return toReturn;
}

// When subtracting in-place
Poly& Poly::operator-=(const Poly& other) {
    for (size_t i=0; i<other.polynomial.size(); i++) {
        polynomial.push_back(std::make_pair(-other.polynomial[i].first, other.polynomial[i].second));
    }
    return *this;
}

// When multiplying two polynomials
Poly Poly::operator*(const Poly& other) {
    Poly toReturn;
    for (size_t i=0; i<polynomial.size(); i++) {
        for (size_t j=0; j<other.polynomial.size(); j++) {
            std::complex<double> newCoefficient = polynomial[i].first * other.polynomial[j].first;
            Mon newMonomial = polynomial[i].second * other.polynomial[j].second;
            toReturn += Poly(newCoefficient, newMonomial);
        }
    }
    return toReturn;
}

// When multiplying by a constant
Poly Poly::operator*(const std::complex<double>& other) {
    Poly toReturn;
    for (size_t i=0; i<polynomial.size(); i++) {
        toReturn.polynomial.push_back(std::make_pair(polynomial[i].first * other, polynomial[i].second));
    }
    return toReturn;
}

// When multiplying by a constant in-place
Poly& Poly::operator*=(const std::complex<double>& other) {
    for (size_t i=0; i<polynomial.size(); i++) {
        polynomial[i].first *= other;
    }
    return *this;
}

// When dividing by a constant
Poly Poly::operator/(const std::complex<double>& other) {
    Poly toReturn;
    for (size_t i=0; i<polynomial.size(); i++) {
        toReturn.polynomial.push_back(std::make_pair(polynomial[i].first / other, polynomial[i].second));
    }
    return toReturn;
}

// When dividing by a constant in-place
Poly& Poly::operator/=(const std::complex<double>& other) {
    for (size_t i=0; i<polynomial.size(); i++) {
        polynomial[i].first /= other;
    }
    return *this;
}

// Sort the terms by monomial
void Poly::sort() {
    std::sort(polynomial.begin(), polynomial.end(), [](const std::pair<std::complex<double>, Mon>& a, const std::pair<std::complex<double>, Mon>& b) {
        return a.second < b.second;
    });
}

// Size of the polynomial
const size_t Poly::size() const {
    return polynomial.size();
}

// Allow bracket access
const std::pair<std::complex<double>, Mon>& Poly::operator[](size_t index) const {
    return polynomial[index];
}

// Allow bracket access
std::pair<std::complex<double>, Mon>& Poly::operator[](size_t index) {
    return polynomial[index];
}

// Self-negative
Poly Poly::operator-() {
    Poly toReturn = *this;
    for (size_t i=0; i<toReturn.size(); i++) {
        toReturn.polynomial[i].first = -toReturn.polynomial[i].first;
    }
    return toReturn;
}

// Self-conjugate
Poly Poly::conj() {
    Poly toReturn = *this;
    for (size_t i=0; i<toReturn.size(); i++) {
        toReturn.polynomial[i].first = std::conj(toReturn.polynomial[i].first);
    }
    return toReturn;
}

// Combine like terms
void Poly::simplify() {

    // Iterate through the polynomial
    for (size_t i=0; i<size(); i++) {

        // Check if this monomial is already in the return polynomial
        bool found = false;
        for (size_t j=0; j<size(); j++) {
            if (i == j) {
                continue;
            }
            if (polynomial[i].second == polynomial[j].second) {
                polynomial[j].first += polynomial[i].first;
                polynomial.erase(polynomial.begin()+i);
                i--;
                break;
            }
        }

    }

    // Remove any zero terms
    for (int i=size()-1; i>=0; i--) {
        if (std::abs(polynomial[i].first) < zeroTol) {
            polynomial.erase(polynomial.begin()+i);
        }
    }

}

// Given a polynomial, reduce each monomial and combine
void Poly::reduce() {

    // Apply the reduction to each monomial
    for (int j=0; j<size(); j++) {
        std::pair<std::complex<double>, Mon> reducedMonomial = polynomial[j].second.reduce();
        polynomial[j].first *= reducedMonomial.first;
        polynomial[j].second = reducedMonomial.second;
    }

}

// Pretty print
std::ostream& operator<<(std::ostream& os, const Poly& p) {

    // Check if it's zero
    if (p.size() == 0 || (p.size() == 1 && p.polynomial[0].first == std::complex<double>(0,0))) {
        os << "0";
        return os;
    }

    // Iterate through the polynomial
    for (size_t i=0; i<p.size(); i++) {
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

// Replace a monomial by a polynomial
void Poly::replace(std::pair<char,int> mon, Poly replacement) {

    // Iterate through the polynomial
    Poly pNew;
    for (int i=0; i<size(); i++) {
        Poly newTerm = Poly(polynomial[i].first);

        // For each of the terms in the monomial
        for (int j=0; j<polynomial[i].second.size(); j++) {

            // If we have a match
            Poly toMultiply(polynomial[i].second[j]);
            if (polynomial[i].second[j] == mon) {
                toMultiply = replacement;
            } else if (polynomial[i].second[j].first == mon.first && mon.second == 0) {
                toMultiply = replacement;
                for (int k=0; k<toMultiply.size(); k++) {
                    for (int l=0; l<toMultiply[k].second.size(); l++) {
                        toMultiply[k].second[l].second = polynomial[i].second[j].second;
                    }
                }
            }

            // Perform the multiplication and reduce
            newTerm = newTerm * toMultiply;
            newTerm.reduce();


        }

        // Add this term to the polynomial
        pNew += newTerm;
        pNew.simplify();

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
    for (int i=0; i<size(); i++) {
        if (polynomial[i].second == mon) {
            pNew += replacement * polynomial[i].first;
            continue;
        } else {
            pNew += polynomial[i];
        }
    }

    // Set the new polynomial
    pNew.simplify();
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
    simplify();
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
    for (size_t i=0; i<size(); i++) {
        toReturn.push_back(polynomial[i].second);
    }
    return toReturn;

}

// Get the polynomial rearranged to equal a certain monomial
// e.g. X1 + 2*Y1 = 0 -> X1 = -2*Y1
Poly Poly::rearranged(Mon mon) {

    // Find this monomial
    Poly returnPoly;
    std::complex<double> coeff = 1;
    for (size_t i=0; i<size(); i++) {
        if (polynomial[i].second == mon) {
            coeff = polynomial[i].first;
        } else {
            returnPoly += Poly(polynomial[i].first, polynomial[i].second);
        }
    }

    // Divide by the coefficient
    returnPoly /= -coeff;
    return returnPoly;

}

// Check if a polynomial has a monomial
bool Poly::hasMonomial(Mon mon) {
    for (size_t i=0; i<size(); i++) {
        if (polynomial[i].second == mon) {
            return true;
        }
    }
    return false;
}


