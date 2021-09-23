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
#include <mutex>

#define WITHOUT_NUMPY 1
#include "matplotlibcpp.h"

#define ROW 0
#define COL 1
#define BOTH 2
#define PI 3.141592653589793
#define PI180 0.017453292519943

struct Matrix
{
    Matrix(const std::vector<double> &v, const int rows, const int cols);
    Matrix(const double uniform_val, const int rows, const int cols);
    Matrix(const std::string &file_path);

    template <typename... args>
    Matrix(const int rows, const int cols, args... v) : rows(rows), cols(cols), size(cols * rows), vector{static_cast<double>(v)...} {};

    void write_file(const std::string &file_path) const;
    void print(const std::string delimeter = " ", const int precision = 6) const;
    void print(const int precision) const;

    Matrix aggregate(int axis, double (*function)(double, double)) const;
    Matrix apply(double (*function)(double)) const;
    Matrix pairWise(double (*function)(double, double), const Matrix &second_Matrix) const;

    Matrix solveSquare(const Matrix &matrix) const;
    Matrix inverse() const;
    std::vector<Matrix> LUFactor() const;
    Matrix Utriangular() const;
    Matrix Ltriangular() const;
    Matrix transpose() const;
    Matrix minorMatrix() const;
    Matrix cofactorMatrix() const;

    Matrix trace() const;
    double minor(int i, int j) const;
    double cofactor(int i, int j) const;
    double det() const;
    Matrix eigVals() const;

    bool isSquare() const;
    bool isSkewSymmetric() const;
    bool isSymmetric() const;
    bool isBoolean() const;
    bool isDiagonal() const;

    Matrix operator>(const double &number) const;
    Matrix operator<(const double &number) const;
    Matrix operator<=(const double &number) const;
    Matrix operator>=(const double &number) const;
    Matrix operator==(const double &number) const;
    Matrix operator!=(const double &number) const;

    Matrix operator>(const Matrix &matrix) const;
    Matrix operator<(const Matrix &matrix) const;
    Matrix operator<=(const Matrix &matrix) const;
    Matrix operator>=(const Matrix &matrix) const;
    Matrix operator==(const Matrix &matrix) const;
    Matrix operator!=(const Matrix &matrix) const;

    Matrix operator*(const double &number) const;
    Matrix operator/(const double &number) const;
    Matrix operator-(const double &number) const;
    Matrix operator+(const double &number) const;

    Matrix operator*(const Matrix &matrix) const;
    Matrix operator/(const Matrix &matrix) const;
    Matrix operator-(const Matrix &matrix) const;
    Matrix operator+(const Matrix &matrix) const;

    const int argmax() const;
    const int argmin() const;
    const double max() const;
    const double min() const;

    Matrix operator()(const std::vector<int> &row_range, const std::vector<int> &col_range) const;

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

    inline Matrix copy() const { return Matrix(this->vector, this->rows, this->cols); };
    inline int get_size() const { return this->size; };
    inline int get_cols() const { return this->cols; };
    inline int get_rows() const { return this->rows; };
    inline const std::vector<double> &get_vector() const { return this->vector; };

private:
    std::vector<double> vector;
    int cols;
    int rows;
    int size;

private:
    static void dotProduct(int i, const Matrix *_this, const Matrix *matrix, std::vector<double> *res);
};

Matrix read_file(const std::string &file_path);
std::ostream &operator<<(std::ostream &out, const Matrix &arr);
Matrix linspace(const double start, const double end, const int partition_count);
Matrix logspace(const double start, const double end, const int partition_count);
Matrix arange(const double start, const double end, const int step_size = 1);
Matrix leastSquares(const Matrix &data, const int order_of_polynomial);
Matrix diagonal(const double value, int rows, int cols);
Matrix random(double start, double _end, int rows, int cols);                                                                                                                      // range implementation and intrandom version
void plot(const Matrix &Matrix, const std::vector<double> &xlim, const std::vector<double> &ylim, const std::string &file_name = "");                                              // plot class
void plotModel(const Matrix &model, const Matrix &data, const Matrix &range, const std::vector<double> &xlim, const std::vector<double> &ylim, const std::string &file_name = ""); // plot class

/*
    Class Structure: Matrix --> matrix, Array --> Poly
    macros --> enums
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