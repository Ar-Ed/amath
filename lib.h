#pragma once
#include <iostream>
#include <vector>
#include <string>


struct array
{
    array(std::vector<double> v, int rows, int cols);

    void print(std::string delimeter = " ", int precision=6) const;
    void print(int precision) const;

    array transpose() const;
    array apply(double (*function) (double)) const;

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

    double at(int row, int col) const;

    int get_size() const;
    int get_cols() const;
    int get_rows() const;

private:
    std::vector<double> vector;
    int cols;
    int rows;
    int size;
};


/*
    row-echelon form Ax = b
    Determinant
    Dimension aggregation
    apply
    Least squares method : best fit line, best fit quadratic, best-fit exponential


*/