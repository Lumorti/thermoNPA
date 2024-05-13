#pragma once
#include <vector>
#include <complex>
#include <string>

// Non-commuting monomial class
class Mon {
public:

    // Vars
    std::vector<std::pair<char, int>> monomial;
    int verbosity = 1;

    // If initialized from nothing
    Mon();

    // If initialized from a vector
    Mon(std::vector<std::pair<char, int>> monomial);

    // If initialized from a char and a number
    Mon(char variable, int index);

    // If initialized from a single part of monomial
    Mon(std::pair<char, int> part);

    // If initialized from a string
    Mon(std::string asString);

    // If asked to reverse the monomial in place
    void reverse();
    
    // If asked to reverse the monomial
    Mon reversed();

    // If multiplied by another monomial
    Mon operator*(const Mon& other);

    // Bracket access
    std::pair<char, int>& operator[](size_t index);

    // Bracket access
    const std::pair<char, int>& operator[](size_t index) const;

    // Check equality
    bool operator==(const Mon& other) const;

    // Size of the monomial
    size_t size() const;

    // Compare two sections of a monomial, but reversed
    bool compareReversed(const std::pair<char, int>& a, const std::pair<char, int>& b) const;

    // Comparison between this and another monomial, but reversed
    bool compareReversed(Mon& a, Mon& b);

    // Comparison between this and another monomial
    bool operator<(const Mon& other) const;
    bool operator>(const Mon& other) const;

    // Given a monomial, reduce it to its simplest form
    std::pair<std::complex<double>, Mon> reduce(std::string swapType = "none", bool diffLettersCommute=false, bool diffNumbersCommute=true, bool pauliReductions=true, std::vector<int> reductionsToIgnore={});

    // Pretty printing
    friend std::ostream& operator<<(std::ostream& os, const Mon& m);
    friend std::ostream& operator<<(std::ostream& os, const std::pair<char, int>& p);
    friend std::ostream& operator<<(std::ostream& os, const std::vector<Mon>& m);

    // Convert a monomial to a string
    std::string toString() const;

};

template <>
struct std::hash<Mon>
{
  std::size_t operator()(const Mon& k) const
  {
    using std::size_t;
    using std::hash;
    using std::string;
    return hash<string>()(k.toString());
  }
};

