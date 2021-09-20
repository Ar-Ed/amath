#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <algorithm>
#include <random>
#include <exception>
#include <thread>

#define WITHOUT_NUMPY 1
#include "matplotlibcpp.h"


#define ROW 0
#define COL 1
#define BOTH 2
#define PI 3.141592653589793
#define PI180 0.017453292519943

struct array
{
    array(std::vector<double> v, int rows, int cols);
    array(double uniform_val, int rows, int cols);
    array(std::string file_path);

    template <typename... args>
    array(int rows, int cols, args... v) : rows(rows), cols(cols), size(cols * rows), vector{static_cast<double>(v)...} {};

    void write_file(std::string file_path);
    void print(std::string delimeter = " ", int precision = 6) const;
    void print(int precision) const;

    array aggregate(int axis, double (*function)(double, double)) const;
    array apply(double (*function)(double)) const;
    array pairWise(double (*function)(double, double), array second_array) const;

    array solveSquare(array matrix);
    array inverse();
    std::vector<array> LUFactor() const;
    array Utriangular() const;
    array Ltriangular() const;
    array transpose() const;
    array minorMatrix() const;
    array cofactorMatrix() const;

    array trace() const;
    double minor(int i, int j) const;
    double cofactor(int i, int j) const;
    array eigBiggerThan() const;
    array eigLessThan() const;
    double det() const;

    bool isSquare() const;
    bool isSkewSymmetric() const;
    bool isSymmetric() const;
    bool isBoolean() const;
    bool isDiagonal() const;

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

    array operator()(std::vector<int> row_range, std::vector<int> col_range) const;

    inline void assign(double val, int row, int col) { this->vector[this->cols * row + col] = val; };
    inline double at(int index) const { return this->vector[index]; }
    inline double at(int row, int col) const { return this->vector[cols * row + col]; }
    inline double operator()(int index) const { return this->vector[index]; }
    inline double operator()(int row, int col) const { return this->vector[cols * row + col]; }

    inline void reShape(int rows, int cols)
    {
        this->rows = rows;
        this->cols = cols;
    };

    inline array copy() const { return array(this->vector, this->rows, this->cols); };
    inline int get_size() const { return this->size; };
    inline int get_cols() const { return this->cols; };
    inline int get_rows() const { return this->rows; };
    inline std::vector<double> get_vector() const { return this->vector; };

private:
    std::vector<double> vector;
    int cols;
    int rows;
    int size;
};

array read_file(std::string file_path);
std::ostream &operator<<(std::ostream &out, const array &arr);
array leastSquares(array data, int order_of_polynomial);
array diagonal(double value, int rows, int cols);
array random(double start, double _end, int rows, int cols); // range implementation and intrandom version
void plot(std::vector<double> xlim, std::vector<double> ylim, array array, std::string file_name); // plot class 
void plotPoly(double start, double end, array array, std::string file_name); // plot class 

/*
    
    ********* linspace, logspace, arange 
    assigning values to indeces
    QR factorization
    eig vals
    templates, complex class
    partial solutions 
    append new columns and rows
    row col and range indexing with step size
    Least squares method : best fit line, best fit quadratic, best-fit exponential
    cout << output formatting --> cout is pprint, print is the normal one : I need to come up with solutions

    complex power, and sqrt
*/
