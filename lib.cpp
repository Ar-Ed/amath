#include "lib.h"

array::array(std::vector<double> v, int rows, int cols) : vector(v), cols(cols), rows(rows), size(cols * rows) {}

void array::print(std::string delimeter, int precision) const
{
    for (int i = 0; i < this->rows; i++)
    {   
        std::cout.precision(precision);
        for (int j = 0; j < this->cols; j++)
            std::cout << this->at(i, j) << delimeter;
        std::cout << "\n";
    }
}

void array::print(int precision) const
{
    for (int i = 0; i < this->rows; i++)
    {   
        std::cout.precision(precision);
        for (int j = 0; j < this->cols; j++)
            std::cout << this->at(i, j) << " ";
        std::cout << "\n";
    }
}

array array::transpose() const
{
    return array(this->vector, this->cols, this->rows);
}

array array::apply(double (*function)(double)) const
{
    std::vector<double> res;
    for (int i = 0; i < this->size; i++)
        res.push_back(function(this->vector[i]));
    return array(res, this->rows, this->cols);
}

array array::operator*(double number) const
{
    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] * number);
    return array(res, this->rows, this->cols);
}

array array::operator/(double number) const
{
    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] / number);
    return array(res, this->rows, this->cols);
}

array array::operator-(double number) const
{
    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] - number);
    return array(res, this->rows, this->cols);
}

array array::operator+(double number) const
{
    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] + number);
    return array(res, this->rows, this->cols);
}

array array::operator*(const array &matrix) const
{
    if (this->cols != matrix.rows)
    {
        std::cerr << "number of cols of the first one should be equal to the number of rows of the second one.";
        exit(-1);
    }
    std::vector<double> res;
    double product = 0;

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            for (int k = 0; k < this->cols; k++)
            {
                product += this->at(i, k) * matrix.at(k, j);
            }
            res.push_back(product);
            product = 0;
        }
    }

    return array(res, this->rows, matrix.cols);
}

array array::operator/(const array &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
    {
        std::cerr << "Number of cols and rows should match";
        exit(-1);
    }
    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] / matrix.vector[i]);
    return array(res, this->rows, this->cols);
}

array array::operator-(const array &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
    {
        std::cerr << "Number of cols and rows should match";
        exit(-1);
    }
    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] - matrix.vector[i]);
    return array(res, this->rows, this->cols);
}

array array::operator+(const array &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
    {
        std::cerr << "Number of cols and rows should match";
        exit(-1);
    }
    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] + matrix.vector[i]);
    return array(res, this->rows, this->cols);
}
array array::copy() const
{
    return array(this->vector, this->rows, this->cols);
}

void array::reShape(int rows, int cols)
{
    this->rows = rows;
    this->cols = cols;
}
double array::at(int row, int col) const
{
    return this->vector[cols * row + col];
}
int array::get_size() const
{
    return this->size;
}
int array::get_cols() const
{
    return this->cols;
}
int array::get_rows() const
{
    return this->rows;
}
