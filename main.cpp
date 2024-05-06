// Standard includes
#include <iostream>
#include <vector>

// Import MOSEK
#include "fusion.h"

// Import Eigen
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

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
const double zeroTol = 1e-10;

// Pretty print part of monomial
std::ostream& operator<<(std::ostream& os, const std::pair<char, int>& p) {
    os << p.first << p.second;
    return os;
}

// Pretty print complex numbers
std::ostream& operator<<(std::ostream& os, const std::complex<double>& c) {
    if (c.imag() == 0) {
        if (c.real() >= 0) {
            os << "+" << c.real();
        } else {
            os << c.real();
        }
    } else if (c.real() == 0) {
        if (c.imag() == 1) {
            os << "+i";
        } else if (c.imag() == -1) {
            os << "-i";
        } else {
            if (c.imag() > 0) {
                os << "+" << c.imag() << "i";
            } else {
                os << c.imag() << "i";
            }
        }
    } else {
        if (c.imag() >= 0) {
            os << "+(" << c.real() << "+" << c.imag() << "i)";
        } else {
            os << "+(" << c.real() << c.imag() << "i)";
        }
    }
    return os;
}
 
// Given a polynomial, combine all of the terms of the same monomial
//polynomial simplify(polynomial p) {

    //// Create the return polynomial
    //polynomial toReturn;

    //// Iterate through the polynomial
    //for (size_t i=0; i<p.size(); i++) {

        //// Check if this monomial is already in the return polynomial
        //bool found = false;
        //for (size_t j=0; j<toReturn.size(); j++) {
            //if (toReturn[j].second == p[i].second) {
                //toReturn[j].first += p[i].first;
                //found = true;
                //break;
            //}
        //}

        //// If not, add it
        //if (!found) {
            //toReturn.push_back(p[i]);
        //}

    //}

    //// Return the simplified polynomial
    //return toReturn;

//}

// Non-commuting monomial class
class Mon {
public:
    std::vector<std::pair<char, int>> monomial;
    int verbosity = 1;

    // If initialized from nothing
    Mon() {
        monomial = {};
    }

    // If initialized from a vector
    Mon(std::vector<std::pair<char, int>> monomial) {
        monomial = monomial;
    }

    // If initialized from a char and a number
    Mon(char variable, int index) {
        monomial = {{variable, index}};
    }

    // If initialized from a single part of monomial
    Mon(std::pair<char, int> part) {
        monomial = {part};
    }

    // If initialized from a string
    Mon(std::string asString) {

        // If the string is empty, return an empty monomial
        if (asString == "" || asString.find("<") == std::string::npos) {
            return;
        }

        // Iterate through the string
        std::vector<std::pair<char, int>> toReturn;
        char currentVariable = ' ';
        std::string currentIndex;
        for (size_t i=0; i<asString.size(); i++) {

            // Check if this char is an integer
            if (asString[i] >= '0' && asString[i] <= '9' || asString[i] == '-') {
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

        // Set the monomial
        monomial = toReturn;

    }

    // If asked to reverse the monomial in place
    void reverse() {
        std::reverse(monomial.begin(), monomial.end());
    }
    
    // If asked to reverse the monomial
    Mon reversed() {
        Mon toReturn = *this;
        std::reverse(toReturn.monomial.begin(), toReturn.monomial.end());
        return toReturn;
    }

    // If multiplied by another monomial
    Mon operator*(const Mon& other) {
        Mon toReturn;
        toReturn.monomial = monomial;
        for (size_t i=0; i<other.monomial.size(); i++) {
            toReturn.monomial.push_back(other.monomial[i]);
        }
        return toReturn;
    }

    // Bracket access
    std::pair<char, int>& operator[](size_t index) {
        return monomial[index];
    }

    // Bracket access
    const std::pair<char, int>& operator[](size_t index) const {
        return monomial[index];
    }

    // Check equality
    bool operator==(const Mon& other) const {
        return monomial == other.monomial;
    }

    // Size of the monomial
    size_t size() const {
        return monomial.size();
    }

    // Compare two sections of a monomial, but reversed
    bool compareReversed(const std::pair<char, int>& a, const std::pair<char, int>& b) const {
        if (a.second == b.second) {
            return a.first < b.first;
        } else {
            return a.second < b.second;
        }
    }

    // Comparison between this and another monomial, but reversed
    bool compareReversed(Mon& a, Mon& b) {
        for (size_t i=0; i<std::min(a.size(), b.size()); i++) {
            if (a[i].second != b[i].second) {
                return a[i].second < b[i].second;
            }
        }
        return a.size() < b.size();
    }

    // Comparison between this and another monomial
    bool operator<(const Mon& other) const {
        return monomial < other.monomial;
    }
    bool operator>(const Mon& other) const {
        return monomial > other.monomial;
    }

    // Given a monomial, reduce it to its simplest form
    std::pair<std::complex<double>, Mon> reduce(std::string swapType = "numFirst", bool diffLettersCommute=false, bool diffNumbersCommute=true, bool pauliReductions=true, std::vector<int> reductionsToIgnore={}) {

        // Sort the monomial as much as we can
        Mon mon = *this;
        if (diffLettersCommute) {
            for (int i=0; i<mon.size(); i++) {
                for (int j=0; j<int(mon.size())-1; j++) {
                    if (mon[j].first != mon[j+1].first && mon[j] > mon[j+1]) {
                        std::swap(mon[j], mon[j+1]);
                    }
                }
            }
        }
        if (diffNumbersCommute) {
            for (int i=0; i<mon.size(); i++) {
                for (int j=0; j<int(mon.size())-1; j++) {
                    if (mon[j].second != mon[j+1].second && !compareReversed(mon[j], mon[j+1])) {
                        std::swap(mon[j], mon[j+1]);
                    }
                }
            }
        }

        // Simplify Pauli strings using the following rules (plus conjugates):
        // XY = iZ
        // ZX = iY
        // YZ = iX
        std::complex<double> coeff(1, 0);
        if (pauliReductions) {
            if (std::find(reductionsToIgnore.begin(), reductionsToIgnore.end(), mon.size()) == reductionsToIgnore.end() && pauliReductions) {
                for (int i=mon.size()-1; i>0; i--) {
                    if (mon[i-1].second == mon[i].second) {
                        if (mon[i-1].first == 'X' && mon[i].first == 'Y') {
                            coeff *= imag;
                            mon[i-1] = std::make_pair('Z', mon[i-1].second);
                            mon.monomial.erase(mon.monomial.begin()+i);
                        } else if (mon[i-1].first == 'X' && mon[i].first == 'Z') {
                            coeff *= -imag;
                            mon[i-1] = std::make_pair('Y', mon[i-1].second);
                            mon.monomial.erase(mon.monomial.begin()+i);
                        } else if (mon[i-1].first == 'Y' && mon[i].first == 'Z') {
                            coeff *= imag;
                            mon[i-1] = std::make_pair('X', mon[i-1].second);
                            mon.monomial.erase(mon.monomial.begin()+i);
                        } else if (mon[i-1].first == 'Y' && mon[i].first == 'X') {
                            coeff *= -imag;
                            mon[i-1] = std::make_pair('Z', mon[i-1].second);
                            mon.monomial.erase(mon.monomial.begin()+i);
                        } else if (mon[i-1].first == 'Z' && mon[i].first == 'X') {
                            coeff *= imag;
                            mon[i-1] = std::make_pair('Y', mon[i-1].second);
                            mon.monomial.erase(mon.monomial.begin()+i);
                        } else if (mon[i-1].first == 'Z' && mon[i].first == 'Y') {
                            coeff *= -imag;
                            mon[i-1] = std::make_pair('X', mon[i-1].second);
                            mon.monomial.erase(mon.monomial.begin()+i);
                        } else if (mon[i-1].first == 'X' && mon[i].first == 'X') {
                            mon.monomial.erase(mon.monomial.begin()+i);
                            mon.monomial.erase(mon.monomial.begin()+i-1);
                            i--;
                        } else if (mon[i-1].first == 'Y' && mon[i].first == 'Y') {
                            mon.monomial.erase(mon.monomial.begin()+i);
                            mon.monomial.erase(mon.monomial.begin()+i-1);
                            i--;
                        } else if (mon[i-1].first == 'Z' && mon[i].first == 'Z') {
                            mon.monomial.erase(mon.monomial.begin()+i);
                            mon.monomial.erase(mon.monomial.begin()+i-1);
                            i--;
                        }
                    }
                }
            }
        }

        // <A1A1> = <1>
        int i = 0;
        while (i < int(mon.size())-1) {
            if (mon[i] == mon[i+1]) {
                mon.monomial.erase(mon.monomial.begin()+i+1);
                mon.monomial.erase(mon.monomial.begin()+i);
                i = -1;
            }
            i++;
        }

        // Flip it to see if it's smaller
        if (swapType == "letFirst") {
            Mon monFlipped = mon.reversed();
            if (monFlipped < mon) {
                mon = monFlipped;
            }
        } else if (swapType == "numFirst") {
            Mon monFlipped = mon.reversed();
            if (compareReversed(monFlipped, mon)) {
                mon = monFlipped;
            }
        }

        // Verbose output
        if (verbosity >= 3) {
            std::cout << "Reduced monomial: " << *this << " -> " << coeff << " " << mon << std::endl;
        }

        return {coeff, mon};

    }

    // Pretty printing
    friend std::ostream& operator<<(std::ostream& os, const Mon& m) {

        // Iterate through the monomial
        if (m.size() > 0) {
            os << "<";
            for (size_t i=0; i<m.size(); i++) {
                os << m.monomial[i].first << m.monomial[i].second;
            }
            os << ">";
        } else {
            os << "1";
        }

        // Return the stream
        return os;

    }

    // Convert a monomial to a string
    std::string toString() {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

};

// Non-commuting polynomial class
class Poly {
public:
    std::vector<std::pair<std::complex<double>, Mon>> polynomial;
    int verbosity = 1;

    // If initialized from nothing
    Poly() {
        polynomial = {};
    }

    // If initialized from just a coeff
    Poly(std::complex<double> coeff) {
        polynomial = {{coeff, Mon()}};
    }

    // If initialized from a vector
    Poly (std::vector<std::pair<std::complex<double>, Mon>> polynomial) {
        polynomial = polynomial;
    }

    // If initialized from a single coeff + monomial
    Poly(std::pair<std::complex<double>, Mon> pair) {
        polynomial = {pair};
    }

    // If initialized from a coeff and a monomial
    Poly(std::complex<double> coeff, Mon mon) {
        polynomial = {{coeff, mon}};
    }

    // If initialized from a single part of monomial
    Poly(std::pair<char, int> part) {
        polynomial = {{1, Mon(part)}};
    }

    // If initialized from a coeff and a string
    Poly(std::complex<double> coeff, std::string asString) {
        polynomial = {{coeff, Mon(asString)}};
    }

    // If initialized from a string
    Poly(std::string asString) {

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
    Poly& operator=(const std::vector<std::pair<std::complex<double>, Mon>>& other) {
        polynomial = other;
        return *this;
    }

    // Checking equality of two polynomials
    bool operator==(const Poly& other) {
        return polynomial == other.polynomial;
    }

    // When summing two polynomials
    Poly operator+(const Poly& other) {
        Poly toReturn;
        toReturn.polynomial = polynomial;
        for (size_t i=0; i<other.polynomial.size(); i++) {
            toReturn.polynomial.push_back(other.polynomial[i]);
        }
        return toReturn;
    }

    // When summing in-place
    Poly& operator+=(const Poly& other) {
        for (size_t i=0; i<other.polynomial.size(); i++) {
            polynomial.push_back(other.polynomial[i]);
        }
        return *this;
    }

    // When subtracting two polynomials
    Poly operator-(const Poly& other) {
        Poly toReturn;
        toReturn.polynomial = polynomial;
        for (size_t i=0; i<other.polynomial.size(); i++) {
            toReturn.polynomial.push_back(std::make_pair(-other.polynomial[i].first, other.polynomial[i].second));
        }
        return toReturn;
    }

    // When multiplying two polynomials
    Poly operator*(const Poly& other) {
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

    // Size of the polynomial
    const size_t size() const {
        return polynomial.size();
    }

    // Allow bracket access
    const std::pair<std::complex<double>, Mon>& operator[](size_t index) const {
        return polynomial[index];
    }

    // Allow bracket access
    std::pair<std::complex<double>, Mon>& operator[](size_t index) {
        return polynomial[index];
    }

    // Self-negative
    Poly operator-() {
        Poly toReturn = *this;
        for (size_t i=0; i<toReturn.size(); i++) {
            toReturn.polynomial[i].first = -toReturn.polynomial[i].first;
        }
        return toReturn;
    }

    // Self-conjugate
    Poly conj() {
        Poly toReturn = *this;
        for (size_t i=0; i<toReturn.size(); i++) {
            toReturn.polynomial[i].first = std::conj(toReturn.polynomial[i].first);
        }
        return toReturn;
    }

    // Combine like terms
    void simplify() {

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
    void reduce() {

        // Apply the reduction to each monomial
        for (int j=0; j<size(); j++) {
            std::pair<std::complex<double>, Mon> reducedMonomial = polynomial[j].second.reduce();
            polynomial[j].first *= reducedMonomial.first;
            polynomial[j].second = reducedMonomial.second;
        }

    }

    // Pretty print
    friend std::ostream& operator<<(std::ostream& os, const Poly& p) {

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

    // Replace a monomial by a polynomial TODO
    void replace(Mon mon, Poly replacement) {

        // Iterate through the polynomial
        Poly pNew;
        for (int i=0; i<size(); i++) {
            Poly newTerm = Poly(polynomial[i].first);

            // For each of the terms in the monomial
            for (int j=0; j<polynomial[i].second.size(); j++) {

                // If we have a match
                Poly toMultiply(polynomial[i].second[j]);
                if (polynomial[i].second[j] == mon[0]) {
                    toMultiply = replacement;
                } else if (polynomial[i].second[j].first == mon[0].first && mon[0].second == -1) {
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
    Poly replaced(Mon mon, Poly replacement) {
        Poly toReturn = *this;
        toReturn.replace(mon, replacement);
        return toReturn;
    }

    // Reduce a polynomial using plus/minus operators to Paulis TODO
    void convertToPaulis() {
        Poly toReplaceP = Poly(0.5, "<X-1>") + Poly(0.5*imag, "<Y-1>");
        Poly toReplaceM = Poly(0.5, "<X-1>") + Poly(-0.5*imag, "<Y-1>");
        replace(Mon('P', -1), toReplaceP);
        replace(Mon('M', -1), toReplaceM);
        simplify();
        reduce();
    }

    // Convert a polynomial to Paulis
    Poly convertedToPaulis() {
        Poly toReturn = *this;
        toReturn.convertToPaulis();
        return toReturn;
    }

};

// Generate a monomial from a given string (e.g. "<A1B2>")
//monomial stringToMonomial(std::string asString) {

    //// If the string is empty, return an empty monomial
    //if (asString == "" || asString.find("<") == std::string::npos) {
        //return monomial();
    //}

    //// Iterate through the string
    //monomial toReturn;
    //char currentVariable = ' ';
    //std::string currentIndex;
    //for (size_t i=0; i<asString.size(); i++) {

        //// Check if this char is an integer
        //if (asString[i] >= '0' && asString[i] <= '9') {
            //currentIndex += asString[i];

        //// Otherwise it's a variable
        //} else if (asString[i] >= 'A' && asString[i] <= 'Z') {
            //currentIndex += asString[i];
            //if (currentVariable != ' ') {
                //toReturn.push_back(std::make_pair(currentVariable, std::stoi(currentIndex)));
            //}
            //currentVariable = asString[i];
            //currentIndex = "";
        //}

    //}

    //// Add the last variable
    //toReturn.push_back(std::make_pair(currentVariable, std::stoi(currentIndex)));

    //// Return the monomial
    //return toReturn;

//}

       
// When printing a monomial, print it as a string
//std::ostream& operator<<(std::ostream& os, const monomial& m) {

    //// Iterate through the monomial
    //if (m.size() > 0) {
        //os << "<";
        //for (size_t i=0; i<m.size(); i++) {
            //os << m[i].first << m[i].second;
        //}
        //os << ">";
    //} else {
        //os << "1";
    //}

    //// Return the stream
    //return os;

//}

// When printing a polynomial, print it as a string
//std::ostream& operator<<(std::ostream& os, const polynomial& p) {

    //// Check if it's zero
    //if (p.size() == 0 || (p.size() == 1 && p[0].first == std::complex<double>(0,0))) {
        //os << "0";
        //return os;
    //}

    //// Iterate through the polynomial
    //for (size_t i=0; i<p.size(); i++) {
        //double realPart = p[i].first.real();
        //double imagPart = p[i].first.imag();
        //if (imagPart == 0) {
            //if (realPart == std::complex<double>(-1, 0)) {
                //os << "-" << p[i].second;
            //} else if (realPart == std::complex<double>(1, 0)) {
                //if (i == 0) {
                    //os << p[i].second;
                //} else {
                    //os << "+" << p[i].second;
                //}
            //} else if (p[i].second.size() == 0 && realPart > 0) {
                //if (i == 0) {
                    //os << realPart;
                //} else {
                    //os << "+" << realPart;
                //}
            //} else if (p[i].second.size() == 0 && realPart < 0) {
                //os << realPart;
            //} else if (realPart > 0) {
                //if (i == 0) {
                    //os << realPart << p[i].second;
                //} else {
                    //os << "+" << realPart << p[i].second;
                //}
            //} else if (realPart != 0) {
                //os << realPart << p[i].second;
            //}
        //} else {
            //os << p[i].first << p[i].second;
        //}
    //}

    //// Return the stream
    //return os;

//}

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

// When printing a matrix of doubles, print it as a string
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& m) {

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

// When printing a vector of something
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {

    // Iterate through the vector
    os << "{";
    for (size_t i=0; i<v.size(); i++) {
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
//polynomial stringToPolynomial(std::string asString) {

    //// Remove all spaces and *
    //std::string newString = "";
    //for (size_t i=0; i<asString.size(); i++) {
        //if (asString[i] != ' ' && asString[i] != '*') {
            //newString += asString[i];
        //}
    //}
    //asString = newString;

    //// Iterate through the string
    //polynomial toReturn;
    //std::string currentCoefficient;
    //std::string currentMonomial;
    //for (size_t i=0; i<asString.size(); i++) {

        //// If finished a monomial
        //if (i > 0 && (asString[i] == '>' || ((asString[i] == '+' || asString[i] == '-') && asString[i-1] != '>'))) {

            //// Ensure the coefficient is convertable
            //if (currentCoefficient == "" || currentCoefficient == "+") {
                //currentCoefficient = "1";
            //} else if (currentCoefficient == "-") {
                //currentCoefficient = "-1";
            //}

            //// Add the term
            //if (asString[i] == '>') {
                //currentMonomial += asString[i];
            //}
            //toReturn.push_back(std::make_pair(std::stod(currentCoefficient), stringToMonomial(currentMonomial)));
            //currentMonomial = "";
            //currentCoefficient = "";

            //// The plus and the minus are for the next term
            //if (asString[i] == '+' || asString[i] == '-') {
                //currentCoefficient += asString[i];
            //}

        //// If starting or continuing a monomial
        //} else if (asString[i] == '<' || currentMonomial.size() > 0) {
            //currentMonomial += asString[i];

        //// Otherwise it's for the coefficient
        //} else if (asString[i] != ' ') {
            //currentCoefficient += asString[i];

        //}

    //}

    //// If there's a coefficient left over, add it
    //if (currentCoefficient != "") {
        //toReturn.push_back(std::make_pair(std::stod(currentCoefficient), stringToMonomial(currentMonomial)));
    //}

    //// Return the polynomial
    //return toReturn;

//}

// Given a matrix and a variable list, return the matrix with the variables replaced
Eigen::MatrixXcd replaceVariables(std::vector<std::vector<Poly>>& momentMatrix, const std::vector<Mon>& variables, const std::vector<std::complex<double>>& varVals) {

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

// Compare two parts of monomials, such that "Y1" < "X2", "Y1" < "Y2"
//bool compareReversed(const std::pair<char, int>& a, const std::pair<char, int>& b) {
    //if (a.second == b.second) {
        //return a.first < b.first;
    //} else {
        //return a.second < b.second;
    //}
//}

// Same as above but for full monomials
//bool compareReversed(const monomial& a, const monomial& b) {
    //for (int i=0; i<std::min(a.size(), b.size()); i++) {
        //if (a[i] != b[i]) {
            //return compareReversed(a[i], b[i]);
        //}
    //}
    //return a.size() < b.size();
//}

// Given a monomial, reduce it to its simplest form
//std::pair<std::complex<double>, monomial> reduceMonomial(const monomial& mon_, int verbosity, std::string swapType = "numFirst", bool diffLettersCommute=false, bool diffNumbersCommute=true, bool pauliReductions=true, std::vector<int> reductionsToIgnore={}) {

    //// Sort the monomial as much as we can
    //monomial mon = mon_;
    //if (diffLettersCommute) {
        //for (int i=0; i<mon.size(); i++) {
            //for (int j=0; j<int(mon.size())-1; j++) {
                //if (mon[j].first != mon[j+1].first && mon[j] > mon[j+1]) {
                    //std::swap(mon[j], mon[j+1]);
                //}
            //}
        //}
    //}
    //if (diffNumbersCommute) {
        //for (int i=0; i<mon.size(); i++) {
            //for (int j=0; j<int(mon.size())-1; j++) {
                //if (mon[j].second != mon[j+1].second && !compareReversed(mon[j], mon[j+1])) {
                    //std::swap(mon[j], mon[j+1]);
                //}
            //}
        //}
    //}

    //// Simplify Pauli strings using the following rules (plus conjugates):
    //// XY = iZ
    //// ZX = iY
    //// YZ = iX
    //std::complex<double> coeff(1, 0);
    //if (pauliReductions) {
        //if (std::find(reductionsToIgnore.begin(), reductionsToIgnore.end(), mon.size()) == reductionsToIgnore.end() && pauliReductions) {
            //for (int i=mon.size()-1; i>0; i--) {
                //if (mon[i-1].second == mon[i].second) {
                    //if (mon[i-1].first == 'X' && mon[i].first == 'Y') {
                        //coeff *= imag;
                        //mon[i-1] = std::make_pair('Z', mon[i-1].second);
                        //mon.erase(mon.begin()+i);
                    //} else if (mon[i-1].first == 'X' && mon[i].first == 'Z') {
                        //coeff *= -imag;
                        //mon[i-1] = std::make_pair('Y', mon[i-1].second);
                        //mon.erase(mon.begin()+i);
                    //} else if (mon[i-1].first == 'Y' && mon[i].first == 'Z') {
                        //coeff *= imag;
                        //mon[i-1] = std::make_pair('X', mon[i-1].second);
                        //mon.erase(mon.begin()+i);
                    //} else if (mon[i-1].first == 'Y' && mon[i].first == 'X') {
                        //coeff *= -imag;
                        //mon[i-1] = std::make_pair('Z', mon[i-1].second);
                        //mon.erase(mon.begin()+i);
                    //} else if (mon[i-1].first == 'Z' && mon[i].first == 'X') {
                        //coeff *= imag;
                        //mon[i-1] = std::make_pair('Y', mon[i-1].second);
                        //mon.erase(mon.begin()+i);
                    //} else if (mon[i-1].first == 'Z' && mon[i].first == 'Y') {
                        //coeff *= -imag;
                        //mon[i-1] = std::make_pair('X', mon[i-1].second);
                        //mon.erase(mon.begin()+i);
                    //} else if (mon[i-1].first == 'X' && mon[i].first == 'X') {
                        //mon.erase(mon.begin()+i);
                        //mon.erase(mon.begin()+i-1);
                        //i--;
                    //} else if (mon[i-1].first == 'Y' && mon[i].first == 'Y') {
                        //mon.erase(mon.begin()+i);
                        //mon.erase(mon.begin()+i-1);
                        //i--;
                    //} else if (mon[i-1].first == 'Z' && mon[i].first == 'Z') {
                        //mon.erase(mon.begin()+i);
                        //mon.erase(mon.begin()+i-1);
                        //i--;
                    //}
                //}
            //}
        //}
    //}

    //// <A1A1> = <1>
    //int i = 0;
    //while (i < int(mon.size())-1) {
        //if (mon[i] == mon[i+1]) {
            //mon.erase(mon.begin()+i+1);
            //mon.erase(mon.begin()+i);
            //i = -1;
        //}
        //i++;
    //}

    //// Flip it to see if it's smaller
    //if (swapType == "letFirst") {
        //monomial monFlipped = mon;
        //std::reverse(monFlipped.begin(), monFlipped.end());
        //if (monFlipped < mon) {
            //mon = monFlipped;
        //}
    //} else if (swapType == "numFirst") {
        //monomial monFlipped = mon;
        //std::reverse(monFlipped.begin(), monFlipped.end());
        //if (compareReversed(monFlipped, mon)) {
            //mon = monFlipped;
        //}
    //}

    //// Verbose output
    //if (verbosity >= 3) {
        //std::cout << "Reduced monomial: " << mon_ << " -> " << coeff << " " << mon << std::endl;
    //}

    //return {coeff, mon};

//}

// Given a polynomial, reduce each monomial and combine
//polynomial reducePolynomial(polynomial p, int verbosity) {

    //// Apply the reduction to each monomial
    //for (int j=0; j<p.size(); j++) {
        //std::pair<std::complex<double>, monomial> reducedMonomial = reduceMonomial(p[j].second, verbosity);
        //p[j].first *= reducedMonomial.first;
        //p[j].second = reducedMonomial.second;
    //}

    //// Combine like terms
    //p = simplify(p);

    //return p;

//}

// Reduce a polynomial using plus/minus operators to Paulis TODO
//polynomial removePlusMinus(polynomial p, int verbosity) {

    //std::cout << "recieved polynomial: " << p << std::endl;

    //// Iterate through the polynomial
    //polynomial pNew;
    //for (int i=0; i<p.size(); i++) {
        
        //std::cout << "current term: " << p[i].first << " " << p[i].second << std::endl;
        //polynomial newTerm = {{p[i].first, monomial()}};

        //for (int j=0; j<p[i].second.size(); j++) {

            //polynomial toMultiply = {{1, monomial({p[i].second[j]})}};

            //// P -> (X + iY) / 2
            //if (p[i].second[j].first == 'P') {
                //monomial monX = monomial({{'X', p[i].second[j].second}});
                //monomial monY = monomial({{'Y', p[i].second[j].second}});
                //toMultiply = {{0.5, monX}, 
                              //{0.5*imag, monY}};

            //// M -> (X - iY) / 2
            //} else if (p[i].second[j].first == 'M') {
                //monomial monX = monomial({{'X', p[i].second[j].second}});
                //monomial monY = monomial({{'Y', p[i].second[j].second}});
                //toMultiply = {{0.5, monX}, 
                              //{-0.5*imag, monY}};

            //}

            //std::cout << "    toMultiply: " << toMultiply << std::endl;
            //std::cout << "    newTerm: " << newTerm << std::endl;

            //// Perform the multiplication
            //polynomial newTermCopy = {};
            //for (int k=0; k<toMultiply.size(); k++) {
                //for (int l=0; l<newTerm.size(); l++) {
                    //newTermCopy.push_back(newTerm[l]);
                    //int newInd = newTermCopy.size()-1;
                    //newTermCopy[newInd].first *= toMultiply[k].first;
                    //newTermCopy[newInd].second.insert(newTermCopy[newInd].second.end(), toMultiply[k].second.begin(), toMultiply[k].second.end());
                //}
            //}

            //std::cout << "    after multiply: " << newTermCopy << std::endl;
            //newTerm = reducePolynomial(newTermCopy, verbosity);
            //std::cout << "    after simplify: " << newTerm << std::endl;

        //}

        //std::cout << "result: " << newTerm << std::endl;

        //// Add this term to the polynomial
        //for (int j=0; j<newTerm.size(); j++) {
            //pNew.push_back(newTerm[j]);
        //}

    //}

    //// Simplify the polynomial
    //std::cout << "before simplify: " << pNew << std::endl;
    //pNew = simplify(pNew);
    //std::cout << "after simplify: " << pNew << std::endl;
    //pNew = reducePolynomial(pNew, verbosity);
    //std::cout << "after reduce: " << pNew << std::endl;

    //// Return the new polynomial
    //return pNew;

//}

// Add all single order monomials from a functional to a list of variables
void addSingleMonomials(std::vector<Mon>& variables, Poly functional) {

    // Iterate through the functional
    for (long unsigned int i=0; i<functional.size(); i++) {

        // Iterate through the monomial
        for (long unsigned int j=0; j<functional[i].second.size(); j++) {
            Mon currentMonomial(functional[i].second[j]);

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
            Poly newPolynomial;
            for (long unsigned int k=0; k<monomsInTopRow[i].size(); k++) {
                for (long unsigned int l=0; l<monomsInTopRow[j].size(); l++) {

                    // Form the new monomial
                    Mon newMonomial;
                    for (long unsigned int m=0; m<monomsInTopRow[i][k].second.size(); m++) {
                        newMonomial.monomial.push_back(monomsInTopRow[i][k].second[m]);
                    }
                    for (int m=int(monomsInTopRow[j][l].second.size())-1; m>=0; m--) {
                        newMonomial.monomial.push_back(monomsInTopRow[j][l].second[m]);
                    }

                    // Reduce the monomial
                    std::pair<std::complex<double>, Mon> monomCoeff = newMonomial.reduce();
                    newMonomial = monomCoeff.second;

                    // Add to the polynomial
                    std::complex<double> newCoefficient = monomCoeff.first * monomsInTopRow[i][k].first * std::conj(monomsInTopRow[j][l].first);
                    newPolynomial += Poly(monomCoeff.first, newMonomial);
                    //bool found = false;
                    //for (long unsigned int m=0; m<newPolynomial.size(); m++) {
                        //if (newPolynomial[m].second == newMonomial) {
                            //newPolynomial[m].first += newCoefficient;
                            //found = true;
                            //break;
                        //}
                    //}
                    //if (!found) {
                        //newPolynomial.push_back(std::make_pair(newCoefficient, newMonomial));
                    //}

                }
            }

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
            Poly currentPolynomial(monomCoeff);
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
                Poly currentPolynomial(monomCoeff);
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
                    Poly currentPolynomial(monomCoeff);
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
                        Poly currentPolynomial(monomCoeff);
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
                            Poly currentPolynomial(monomCoeff);
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
                                Poly currentPolynomial(monomCoeff);
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
                                    Poly currentPolynomial(monomCoeff);
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
            for (long unsigned int k=0; k<toAdd[i][j].size(); k++) {
                Mon currentMonomial(toAdd[i][j][k].second);

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

    // Iterate through the monomial
    for (long unsigned int i=0; i<toAdd.size(); i++) {
        Mon currentMonomial(toAdd[i].second);

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

// Take the trace of a matrix, assuming it's real
double tr(Eigen::MatrixXcd A) {
    std::complex<double> trace = A.trace();
    if (std::abs(std::imag(trace)) > 1e-5) {
        std::cout << "WARNING - trace of matrix has non-zero imaginary part" << std::endl;
    }
    return std::real(trace);
}

// Generic entry function
int main(int argc, char* argv[]) {

    // Define the scenario
    int level = 1;
    int limbladLevel = 1;
    int numQubits = 1;
    Poly objective("<Z1>");
    int verbosity = 1;
    std::complex<double> knownIdeal = 0.0;
    std::string seed = "";
    Poly limbladian("<X1A1X1>+<Y1A1Y1>-<A1>");
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
            objective = Poly("<X1>");

        // Pauli Y objective
        } else if (argAsString == "--objY") {
            objective = Poly("<Y1>");

        // Pauli Z objective
        } else if (argAsString == "--objZ") {
            objective = Poly("<Z1>");

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

            //std::cout << "J_ss from rho: " << J_ss_from_rho << std::endl;
            
            //Eigen::MatrixXcd A = sigma_z_c;
            
            //std::complex<double> exp_hphmA = (sigma_h_plus*sigma_h_minus*A*rho).trace();
            //std::complex<double> exp_cpcmA = (sigma_c_plus*sigma_c_minus*A*rho).trace();
            //std::complex<double> exp_hpcmA = (sigma_h_plus*sigma_c_minus*A*rho).trace();
            //std::complex<double> exp_hmcpA = (sigma_h_minus*sigma_c_plus*A*rho).trace();
            //std::complex<double> exp_Ahphm = (A*sigma_h_plus*sigma_h_minus*rho).trace();
            //std::complex<double> exp_Acpcm = (A*sigma_c_plus*sigma_c_minus*rho).trace();
            //std::complex<double> exp_Ahpcm = (A*sigma_h_plus*sigma_c_minus*rho).trace();
            //std::complex<double> exp_Ahmcp = (A*sigma_h_minus*sigma_c_plus*rho).trace();
            //std::complex<double> exp_hpAhm = (sigma_h_plus*A*sigma_h_minus*rho).trace();
            //std::complex<double> exp_hmAhp = (sigma_h_minus*A*sigma_h_plus*rho).trace();
            //std::complex<double> exp_cpAcm = (sigma_c_plus*A*sigma_c_minus*rho).trace();
            //std::complex<double> exp_cmAcp = (sigma_c_minus*A*sigma_c_plus*rho).trace();
            //std::complex<double> exp_hmhpA = (sigma_h_minus*sigma_h_plus*A*rho).trace();
            //std::complex<double> exp_cmcpA = (sigma_c_minus*sigma_c_plus*A*rho).trace();
            //std::complex<double> exp_Ahmhp = (A*sigma_h_minus*sigma_h_plus*rho).trace();
            //std::complex<double> exp_Acmcp = (A*sigma_c_minus*sigma_c_plus*rho).trace();
            
            //std::complex<double> L_A = - imag*epsilon_h*exp_Ahphm - imag*epsilon_c*exp_Acpcm 
                                       //- imag*g*exp_Ahpcm - imag*g*exp_Ahmcp
                                       //+ imag*epsilon_h*exp_hphmA + imag*epsilon_c*exp_cpcmA
                                       //+ imag*g*exp_hpcmA + imag*g*exp_hmcpA
                                       //+ gamma_h_plus*exp_hmAhp - 0.5*gamma_h_plus*exp_Ahmhp - 0.5*gamma_h_plus*exp_hmhpA
                                       //+ gamma_h_minus*exp_hpAhm - 0.5*gamma_h_minus*exp_Ahphm - 0.5*gamma_h_minus*exp_hphmA
                                       //+ gamma_c_plus*exp_cmAcp - 0.5*gamma_c_plus*exp_Acmcp - 0.5*gamma_c_plus*exp_cmcpA
                                       //+ gamma_c_minus*exp_cpAcm - 0.5*gamma_c_minus*exp_Acpcm - 0.5*gamma_c_minus*exp_cpcmA;
            //std::cout << "L_A: " << L_A << std::endl;

            //Eigen::MatrixXcd hphmA = sigma_h_plus*sigma_h_minus*A;
            //Eigen::MatrixXcd cpcmA = sigma_c_plus*sigma_c_minus*A;
            //Eigen::MatrixXcd hpcmA = sigma_h_plus*sigma_c_minus*A;
            //Eigen::MatrixXcd hmcpA = sigma_h_minus*sigma_c_plus*A;
            //Eigen::MatrixXcd Ahphm = A*sigma_h_plus*sigma_h_minus;
            //Eigen::MatrixXcd Acpcm = A*sigma_c_plus*sigma_c_minus;
            //Eigen::MatrixXcd Ahpcm = A*sigma_h_plus*sigma_c_minus;
            //Eigen::MatrixXcd Ahmcp = A*sigma_h_minus*sigma_c_plus;
            //Eigen::MatrixXcd hpAhm = sigma_h_plus*A*sigma_h_minus;
            //Eigen::MatrixXcd hmAhp = sigma_h_minus*A*sigma_h_plus;
            //Eigen::MatrixXcd cpAcm = sigma_c_plus*A*sigma_c_minus;
            //Eigen::MatrixXcd cmAcp = sigma_c_minus*A*sigma_c_plus;
            //Eigen::MatrixXcd hmhpA = sigma_h_minus*sigma_h_plus*A;
            //Eigen::MatrixXcd cmcpA = sigma_c_minus*sigma_c_plus*A;
            //Eigen::MatrixXcd Ahmhp = A*sigma_h_minus*sigma_h_plus;
            //Eigen::MatrixXcd Acmcp = A*sigma_c_minus*sigma_c_plus;
            
            //Eigen::MatrixXcd LMat = - imag*epsilon_h*Ahphm - imag*epsilon_c*Acpcm 
                                       //- imag*g*Ahpcm - imag*g*Ahmcp
                                       //+ imag*epsilon_h*hphmA + imag*epsilon_c*cpcmA
                                       //+ imag*g*hpcmA + imag*g*hmcpA
                                       //+ gamma_h_plus*hmAhp - 0.5*gamma_h_plus*Ahmhp - 0.5*gamma_h_plus*hmhpA
                                       //+ gamma_h_minus*hpAhm - 0.5*gamma_h_minus*Ahphm - 0.5*gamma_h_minus*hphmA
                                       //+ gamma_c_plus*cmAcp - 0.5*gamma_c_plus*Acmcp - 0.5*gamma_c_plus*cmcpA
                                       //+ gamma_c_minus*cpAcm - 0.5*gamma_c_minus*Acpcm - 0.5*gamma_c_minus*cpcmA;
            
            //std::complex<double> exp_A = (A*rho).trace();
            //std::complex<double> exp_Ahz = (A*sigma_h_z*rho).trace();
            //std::complex<double> exp_Acz = (A*sigma_c_z*rho).trace();
            //std::complex<double> exp_hzA = (sigma_h_z*A*rho).trace();
            //std::complex<double> exp_czA = (sigma_c_z*A*rho).trace();
            //std::complex<double> exp_Ahxcx = (A*sigma_h_x*sigma_c_x*rho).trace();
            //std::complex<double> exp_Ahycy = (A*sigma_h_y*sigma_c_y*rho).trace();
            //std::complex<double> exp_hxcxA = (sigma_h_x*sigma_c_x*A*rho).trace();
            //std::complex<double> exp_hycyA = (sigma_h_y*sigma_c_y*A*rho).trace();
            //std::complex<double> exp_hxAhx = (sigma_h_x*A*sigma_h_x*rho).trace();
            //std::complex<double> exp_hxAhy = (sigma_h_x*A*sigma_h_y*rho).trace();
            //std::complex<double> exp_hyAhx = (sigma_h_y*A*sigma_h_x*rho).trace();
            //std::complex<double> exp_hyAhy = (sigma_h_y*A*sigma_h_y*rho).trace();
            //std::complex<double> exp_cxAcx = (sigma_c_x*A*sigma_c_x*rho).trace();
            //std::complex<double> exp_cxAcy = (sigma_c_x*A*sigma_c_y*rho).trace();
            //std::complex<double> exp_cyAcx = (sigma_c_y*A*sigma_c_x*rho).trace();
            //std::complex<double> exp_cyAcy = (sigma_c_y*A*sigma_c_y*rho).trace();
            
            //Eigen::MatrixXcd Ahz = A*sigma_h_z;
            //Eigen::MatrixXcd Acz = A*sigma_c_z;
            //Eigen::MatrixXcd hzA = sigma_h_z*A;
            //Eigen::MatrixXcd czA = sigma_c_z*A;
            //Eigen::MatrixXcd Ahxcx = A*sigma_h_x*sigma_c_x;
            //Eigen::MatrixXcd Ahycy = A*sigma_h_y*sigma_c_y;
            //Eigen::MatrixXcd hxcxA = sigma_h_x*sigma_c_x*A;
            //Eigen::MatrixXcd hycyA = sigma_h_y*sigma_c_y*A;
            //Eigen::MatrixXcd hxAhx = sigma_h_x*A*sigma_h_x;
            //Eigen::MatrixXcd hxAhy = sigma_h_x*A*sigma_h_y;
            //Eigen::MatrixXcd hyAhx = sigma_h_y*A*sigma_h_x;
            //Eigen::MatrixXcd hyAhy = sigma_h_y*A*sigma_h_y;
            //Eigen::MatrixXcd cxAcx = sigma_c_x*A*sigma_c_x;
            //Eigen::MatrixXcd cxAcy = sigma_c_x*A*sigma_c_y;
            //Eigen::MatrixXcd cyAcx = sigma_c_y*A*sigma_c_x;
            //Eigen::MatrixXcd cyAcy = sigma_c_y*A*sigma_c_y;
            
            //std::complex<double> L_A_pauli = (-2*gamma_h_plus-2*gamma_h_minus-2*gamma_c_plus-2*gamma_c_minus)*exp_A
                                           //+ (-2*imag*epsilon_h-gamma_h_minus+gamma_h_plus)*exp_Ahz
                                           //+ (-2*imag*epsilon_c-gamma_c_minus+gamma_c_plus)*exp_Acz 
                                           //+ (2*imag*epsilon_h-gamma_h_minus+gamma_h_plus)*exp_hzA
                                           //+ (2*imag*epsilon_c-gamma_c_minus+gamma_c_plus)*exp_czA
                                           //+ (-2*imag*g)*exp_Ahxcx
                                           //+ (-2*imag*g)*exp_Ahycy
                                           //+ (2*imag*g)*exp_hxcxA
                                           //+ (2*imag*g)*exp_hycyA
                                           //+ (gamma_h_plus+gamma_h_minus)*exp_hxAhx
                                           //+ (imag*gamma_h_plus-imag*gamma_h_minus)*exp_hxAhy
                                           //+ (-imag*gamma_h_plus+imag*gamma_h_minus)*exp_hyAhx
                                           //+ (gamma_h_plus+gamma_h_minus)*exp_hyAhy
                                           //+ (gamma_c_plus+gamma_c_minus)*exp_cxAcx
                                           //+ (imag*gamma_c_plus-imag*gamma_c_minus)*exp_cxAcy
                                           //+ (-imag*gamma_c_plus+imag*gamma_c_minus)*exp_cyAcx
                                           //+ (gamma_c_plus+gamma_c_minus)*exp_cyAcy;
                                           
            //Eigen::MatrixXcd LMat_pauli = (-2*gamma_h_plus-2*gamma_h_minus-2*gamma_c_plus-2*gamma_c_minus)*A
                                           //+ (-2*imag*epsilon_h-gamma_h_minus+gamma_h_plus)*Ahz
                                           //+ (-2*imag*epsilon_c-gamma_c_minus+gamma_c_plus)*Acz 
                                           //+ (2*imag*epsilon_h-gamma_h_minus+gamma_h_plus)*hzA
                                           //+ (2*imag*epsilon_c-gamma_c_minus+gamma_c_plus)*czA
                                           //+ (-2*imag*g)*Ahxcx
                                           //+ (-2*imag*g)*Ahycy
                                           //+ (2*imag*g)*hxcxA
                                           //+ (2*imag*g)*hycyA
                                           //+ (gamma_h_plus+gamma_h_minus)*hxAhx
                                           //+ (imag*gamma_h_plus-imag*gamma_h_minus)*hxAhy
                                           //+ (-imag*gamma_h_plus+imag*gamma_h_minus)*hyAhx
                                           //+ (gamma_h_plus+gamma_h_minus)*hyAhy
                                           //+ (gamma_c_plus+gamma_c_minus)*cxAcx
                                           //+ (imag*gamma_c_plus-imag*gamma_c_minus)*cxAcy
                                           //+ (-imag*gamma_c_plus+imag*gamma_c_minus)*cyAcx
                                           //+ (gamma_c_plus+gamma_c_minus)*cyAcy;

            //std::cout << "L_A_pauli: " << L_A_pauli << std::endl;

            // Other system settings
            numQubits = 2;
            knownIdeal = J_ss_from_rho;

            // Construct the objective as a polynomial with plus/minus mats TODO
            objective = Poly();

            // Q_ss hot
            //objective.push_back({epsilon_h*gamma_h_plus,       Mon("<M1P1M1P1>")});
            //objective.push_back({-0.5*epsilon_h*gamma_h_plus,  Mon("<P1M1M1P1>")});
            //objective.push_back({-0.5*epsilon_h*gamma_h_plus,  Mon("<M1P1P1M1>")});
            //objective.push_back({epsilon_h*gamma_h_minus,      Mon("<P1P1M1M1>")});
            //objective.push_back({-0.5*epsilon_h*gamma_h_minus, Mon("<P1M1P1M1>")});
            //objective.push_back({-0.5*epsilon_h*gamma_h_minus, Mon("<P1M1P1M1>")});
            //objective.push_back({epsilon_c*gamma_h_plus,       Mon("<M1P2M2P1>")});
            //objective.push_back({-0.5*epsilon_c*gamma_h_plus,  Mon("<P2M2M1P1>")});
            //objective.push_back({-0.5*epsilon_c*gamma_h_plus,  Mon("<M1P1P2M2>")});
            //objective.push_back({epsilon_c*gamma_h_minus,      Mon("<P1P2M2M1>")});
            //objective.push_back({-0.5*epsilon_c*gamma_h_minus, Mon("<P2M2P1M1>")});
            //objective.push_back({-0.5*epsilon_c*gamma_h_minus, Mon("<P1M1P2M2>")});
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

            // Q_ss cold
            //objective.push_back({-epsilon_h*gamma_c_plus,       stringToMonomial("<M2P1M1P2>")});
            //objective.push_back({0.5*epsilon_h*gamma_c_plus,  stringToMonomial("<P1M1M2P2>")});
            //objective.push_back({0.5*epsilon_h*gamma_c_plus,  stringToMonomial("<M2P2P1M1>")});
            //objective.push_back({-epsilon_h*gamma_c_minus,      stringToMonomial("<P2P1M1M2>")});
            //objective.push_back({0.5*epsilon_h*gamma_c_minus, stringToMonomial("<P1M1P2M2>")});
            //objective.push_back({0.5*epsilon_h*gamma_c_minus, stringToMonomial("<P2M2P1M1>")});
            //objective.push_back({-epsilon_c*gamma_c_plus,       stringToMonomial("<M2P2M2P2>")});
            //objective.push_back({0.5*epsilon_c*gamma_c_plus,  stringToMonomial("<P2M2M2P2>")});
            //objective.push_back({0.5*epsilon_c*gamma_c_plus,  stringToMonomial("<M2P2P2M2>")});
            //objective.push_back({-epsilon_c*gamma_c_minus,      stringToMonomial("<P2P2M2M2>")});
            //objective.push_back({0.5*epsilon_c*gamma_c_minus, stringToMonomial("<P2M2P2M2>")});
            //objective.push_back({0.5*epsilon_c*gamma_c_minus, stringToMonomial("<P2M2P2M2>")});
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

            // Simplify the objective
            objective.convertToPaulis();
            std::cout << "From plus minus:" << std::endl;
            std::cout << objective << std::endl;
            
            // Construct the objective as a polynomial
            objective = Poly();
            //objective.push_back({0.5*epsilon_h*(gamma_h_plus-gamma_h_minus)-0.5*epsilon_c*(gamma_c_plus-gamma_c_minus), stringToMonomial("")});
            //objective.push_back({-0.5*epsilon_h*(gamma_h_plus+gamma_h_minus), stringToMonomial("<Z1>")});
            //objective.push_back({0.5*epsilon_c*(gamma_c_plus+gamma_c_minus), stringToMonomial("<Z2>")});
            objective += Poly(0.5*epsilon_h*(gamma_h_plus-gamma_h_minus)-0.5*epsilon_c*(gamma_c_plus-gamma_c_minus), "");
            objective += Poly(-0.5*epsilon_h*(gamma_h_plus+gamma_h_minus), "<Z1>");
            objective += Poly(0.5*epsilon_c*(gamma_c_plus+gamma_c_minus), "<Z2>");
            std::cout << "From paulis:" << std::endl;
            std::cout << objective << std::endl;
                        
            // Construct the Limbadlian as a polynomial
            limbladian = Poly();
            //limbladian.push_back({(-2*gamma_h_plus-2*gamma_h_minus-2*gamma_c_plus-2*gamma_c_minus), stringToMonomial("<A1>")});
            //limbladian.push_back({(-2*imag*epsilon_h-gamma_h_minus+gamma_h_plus), stringToMonomial("<A1Z1>")});
            //limbladian.push_back({(-2*imag*epsilon_c-gamma_c_minus+gamma_c_plus), stringToMonomial("<A1Z2>")});
            //limbladian.push_back({(2*imag*epsilon_h-gamma_h_minus+gamma_h_plus), stringToMonomial("<Z1A1>")});
            //limbladian.push_back({(2*imag*epsilon_c-gamma_c_minus+gamma_c_plus), stringToMonomial("<Z2A1>")});
            //limbladian.push_back({(-2*imag*g), stringToMonomial("<A1X1X2>")});
            //limbladian.push_back({(-2*imag*g), stringToMonomial("<A1Y1Y2>")});
            //limbladian.push_back({(2*imag*g), stringToMonomial("<X1X2A1>")});
            //limbladian.push_back({(2*imag*g), stringToMonomial("<Y1Y2A1>")});
            //limbladian.push_back({(gamma_h_plus+gamma_h_minus), stringToMonomial("<X1A1X1>")});
            //limbladian.push_back({(imag*gamma_h_plus-imag*gamma_h_minus), stringToMonomial("<X1A1Y1>")});
            //limbladian.push_back({(-imag*gamma_h_plus+imag*gamma_h_minus), stringToMonomial("<Y1A1X1>")});
            //limbladian.push_back({(gamma_h_plus+gamma_h_minus), stringToMonomial("<Y1A1Y1>")});
            //limbladian.push_back({(gamma_c_plus+gamma_c_minus), stringToMonomial("<X2A1X2>")});
            //limbladian.push_back({(imag*gamma_c_plus-imag*gamma_c_minus), stringToMonomial("<X2A1Y2>")});
            //limbladian.push_back({(-imag*gamma_c_plus+imag*gamma_c_minus), stringToMonomial("<Y2A1X2>")});
            //limbladian.push_back({(gamma_c_plus+gamma_c_minus), stringToMonomial("<Y2A1Y2>")});
            limbladian += Poly(-2*gamma_h_plus-2*gamma_h_minus-2*gamma_c_plus-2*gamma_c_minus, "<A1>");
            limbladian += Poly(-2*imag*epsilon_h-gamma_h_minus+gamma_h_plus, "<A1Z1>");
            limbladian += Poly(-2*imag*epsilon_c-gamma_c_minus+gamma_c_plus, "<A1Z2>");
            limbladian += Poly(2*imag*epsilon_h-gamma_h_minus+gamma_h_plus, "<Z1A1>");
            limbladian += Poly(2*imag*epsilon_c-gamma_c_minus+gamma_c_plus, "<Z2A1>");
            limbladian += Poly(-2*imag*g, "<A1X1X2>");
            limbladian += Poly(-2*imag*g, "<A1Y1Y2>");
            limbladian += Poly(2*imag*g, "<X1X2A1>");
            limbladian += Poly(2*imag*g, "<Y1Y2A1>");
            limbladian += Poly(gamma_h_plus+gamma_h_minus, "<X1A1X1>");
            limbladian += Poly(imag*gamma_h_plus-imag*gamma_h_minus, "<X1A1Y1>");
            limbladian += Poly(-imag*gamma_h_plus+imag*gamma_h_minus, "<Y1A1X1>");
            limbladian += Poly(gamma_h_plus+gamma_h_minus, "<Y1A1Y1>");
            limbladian += Poly(gamma_c_plus+gamma_c_minus, "<X2A1X2>");
            limbladian += Poly(imag*gamma_c_plus-imag*gamma_c_minus, "<X2A1Y2>");
            limbladian += Poly(-imag*gamma_c_plus+imag*gamma_c_minus, "<Y2A1X2>");
            limbladian += Poly(gamma_c_plus+gamma_c_minus, "<Y2A1Y2>");
                                           
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
            std::cout << "  -m <num>        Level of the moment matrix" << std::endl;
            std::cout << "  -l <num>        Level of the moments to put in the Limbladian" << std::endl;
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
        Mon oldMon("<A1>");
        Poly newPoly(variablesToPut[i]);
        Poly newConstraint = limbladian.replaced(oldMon, newPoly);
        //if (monomToPut.size() > 0) {
            //for (int j=0; j<newConstraint.size(); j++) {
                //for (int k=0; k<newConstraint[j].second.size(); k++) {
                    //if (newConstraint[j].second[k].first == 'A') {
                        //newConstraint[j].second[k] = monomToPut[0];
                        //for (int l=1; l<monomToPut.size(); l++) {
                            //newConstraint[j].second.insert(newConstraint[j].second.begin() + k + l, monomToPut[l]);
                        //}
                    //}
                //}
            //}
        //}
        if (verbosity >= 3) {
            std::cout << std::endl;
            std::cout << "Variable to put: " << newPoly << std::endl;
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
        constraintsZero[i].reduce();
        if (verbosity >= 2) {
            std::cout << "Reduced constraint: " << constraintsZero[i] << std::endl;
        }
    }

    // Define the moment matrix
    if (verbosity >= 2) {
        std::cout << std::endl;
    }
    std::vector<std::vector<std::vector<Poly>>> momentMatrices = generateAllMomentMatrices(objective, constraintsZero, level, verbosity, reductionsToIgnore);

    // If told to add extra to the top row
    if (extraMonomials.size() > 0) {
        std::vector<Poly> topRow = momentMatrices[0][0];
        for (int i=0; i<extraMonomials.size(); i++) {
            Poly extraMonomial(extraMonomials[i]);
            for (int j=0; j<extraMonomial.size(); j++) {
                //std::reverse(extraMonomial[j].second.begin(), extraMonomial[j].second.end());
                extraMonomial[j].second.reverse();
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
        std::cout << "Known ideal: " << knownIdeal << std::endl;
        std::cout << "Upper bound: " << upperBound << " (" << 100*(upperBound-knownIdeal)/knownIdeal << "%)" << std::endl;
        std::cout << "Lower bound: " << lowerBound << " (" << 100*(lowerBound-knownIdeal)/knownIdeal << "%)" << std::endl;
        std::cout << "Difference: " << upperBound - lowerBound << std::endl;
        std::cout << "Error: " << 100*(upperBound-lowerBound)/knownIdeal << "%" << std::endl;
    }

    // Exit without errors
    return 0;

}
