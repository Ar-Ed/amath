#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#define ROW 0
#define COL 1
#define PI 3.141592653589793
#define PI180 0.017453292519943

struct rational
{
    rational(unsigned long long numerator, unsigned long long denominator);

    double to_double() const;

    rational operator+(int constant) const;
    rational operator-(int constant) const;
    rational operator*(int constant) const;
    rational operator/(int constant) const;

    rational operator+(const rational &rat) const;
    rational operator-(const rational &rat) const;
    rational operator*(const rational &rat) const;
    rational operator/(const rational &rat) const;

private:
    unsigned long long numerator; // prime array???
    unsigned long long denominator;
    char sign;
};

struct complex
{
    complex(double real, double imag);
    complex(double angle);

    /*     complex operator^(int power) const; */
    complex operator~() const;

    complex operator+(double constant) const;
    complex operator-(double constant) const;
    complex operator*(double constant) const;
    complex operator/(double constant) const;

    complex operator+(const complex &number) const;
    complex operator-(const complex &number) const;
    complex operator*(const complex &number) const;
    complex operator/(const complex &number) const;

    double get_real() const;
    double get_imag() const;
    double get_angle() const;

private:
    std::vector<double> angle_to_complex(double angle);
    double real;
    double imag;
    double angle;
};

struct array
{
    array(std::vector<double> v, int rows, int cols);
    array(std::string file_path);

    void read_file(std::string file_path);

    void write_file(std::string file_path);
    void print(std::string delimeter = " ", int precision = 6) const;
    void print(int precision) const;

    array aggregate(int axis, double (*function)(double, double)) const;
    array apply(double (*function)(double)) const;

    array solveSquare(array matrix);
    array Utriangular() const;
    array Ltriangular() const;
    array transpose() const;
    double det() const;

    array operator>(double number) const;
    array operator<(double number) const;
    array operator<=(double number) const;
    array operator>=(double number) const;
    array operator==(double number) const;
    array operator!=(double number) const;

    array operator>(const array &matrix) const;
    array operator<(const array &matrix) const;
    array operator<=(const array &matrix) const;
    array operator>=(const array &matrix) const;
    array operator==(const array &matrix) const;
    array operator!=(const array &matrix) const;

    array operator*(double number) const;
    array operator/(double number) const;
    array operator-(double number) const;
    array operator+(double number) const;

    array operator*(const array &matrix) const;
    array operator/(const array &matrix) const;
    array operator-(const array &matrix) const;
    array operator+(const array &matrix) const;

    array copy() const;

    void reShape(int rows, int cols);

    double operator()(int row, int col) const;
    double operator()(int index) const;
    array operator()(std::vector<int> row_range, std::vector<int> col_range) const;

    double at(int row, int col) const;
    double at(int index) const;

    std::vector<double> get_vector() const;
    int get_size() const;
    int get_cols() const;
    int get_rows() const;

private:
    std::vector<double> vector;
    int cols;
    int rows;
    int size;
};

std::ostream &operator<<(std::ostream &out, const array &arr);
std::ostream &operator<<(std::ostream &out, const complex &number);
array leastSquares(array data, int order_of_polynomial);


/*
    templates, complex class
    isSquare, minors, cofactors
    row col and range indexing with step size
    row-echelon form and solving Ax = b
    upper triangular
    Least squares method : best fit line, best fit quadratic, best-fit exponential
    cout << output formatting --> cout is pprint, print is the normal one
    reading from file and file read constructor doesn't work

    complex power, and sqrt
*/