#pragma once
#include <iostream>
#include <vector>
#include <string>

#define ROW 0
#define COL 1

struct array
{
    array(std::vector<double> v, int rows, int cols);

    void print(std::string delimeter = " ", int precision = 6) const;
    void print(int precision) const;

    array transpose() const;
    array apply(double (*function)(double)) const;
    array Utriangular() const;
    array aggregate(int axis, double (*function)(double, double)) const;
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

/*
    row col and range indexing with step size
    row-echelon form and solvinf Ax = b
    Dimension aggregation   
    Least squares method : best fit line, best fit quadratic, best-fit exponential
*/