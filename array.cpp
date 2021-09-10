#include "amatrix.h"

array::array(std::vector<double> v, int rows, int cols) : vector(v), cols(cols), rows(rows), size(cols * rows) {}
array::array(std::string file_path)
{
    this->read_file(file_path);
}

void array::write_file(std::string file_path)
{
    std::ofstream f;
    f.open(file_path);

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
            f << (*this)(i, j) << ",";
        f << "\n";
    }

    f.close();
}

void array::read_file(std::string file_path)
{
    std::ifstream f;
    std::vector<double> vec;
    std::string temp, word;
    int _rows = 0;

    f.open(file_path);

    while (f >> temp)
    {
        std::stringstream s(temp);
        _rows++;

        while (std::getline(s, word, ','))
            vec.push_back(std::stod(word));
    }

    this->rows = _rows;
    this->vector = vec;
    this->size = this->vector.size();
    this->cols = this->size / _rows;
    f.close();
}

array leastSquares(array data, int order_of_polynomial)
{
    std::vector<double> x_vals(data.get_rows()), y_vals(data.get_rows()), a((order_of_polynomial + 1) * (order_of_polynomial + 1)), _b(order_of_polynomial + 1), sums(order_of_polynomial * 2 + 1);

    for (int i = 0; i < data.get_rows(); i++)
    {
        x_vals[i] = data.at(i, 0);
        y_vals[i] = data.at(i, 1);
    }

    for (int i = 0; i < order_of_polynomial + 1; i++)
    {
        double sum = 0;
        for (int j = 0; j < data.get_rows(); j++)
            sum += y_vals[j] * pow(x_vals[j], i);
        _b[i] = sum;
    }

    for (int i = 0; i < (order_of_polynomial * 2 + 1); i++)
    {
        double sum = 0;
        for (double &num : x_vals)
            sum += pow(num, i);
        sums[i] = sum;
    }

    for (int i = 0; i < order_of_polynomial + 1; i++)
    {
        for (int j = 0; j < order_of_polynomial + 1; j++)
            a[j * (order_of_polynomial + 1) + i] = sums[order_of_polynomial - i + j];
    }

    array b(_b, order_of_polynomial + 1, 1);
    array A(a, order_of_polynomial + 1, order_of_polynomial + 1);

    return b.solveSquare(A);
}

array array::solveSquare(array matrix)
{
    if (matrix.get_rows() != matrix.get_cols())
    {
        std::cerr << "Matrix should be square.";
        exit(-1);
    }
    else if (this->rows != matrix.get_rows())
    {
        std::cerr << "Number of rows of b should be equal to the number of columns of the matrix";
        exit(-1);
    }

    std::vector<double> vec{matrix.get_vector()};
    std::vector<double> x{this->get_vector()};

    for (int i = 0; i < matrix.get_rows(); i++)
    {
        for (int j = i + 1; j < matrix.get_rows(); j++)
        {
            double constant = vec[matrix.get_cols() * j + i] / vec[(matrix.get_cols() + 1) * i];
            x[j] -= constant * x[i];
            for (int k = 0; k < matrix.get_cols(); k++)
                vec[matrix.get_cols() * j + k] -= constant * vec[matrix.get_cols() * i + k];
        }
    }

    for (int i = matrix.get_rows() - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            double constant = vec[matrix.get_cols() * j + i] / vec[(matrix.get_cols() + 1) * i];
            x[j] -= constant * x[i];
            for (int k = 0; k < matrix.get_cols(); k++)
                vec[matrix.get_cols() * j + k] = vec[matrix.get_cols() * j + k] - constant * vec[matrix.get_cols() * i + k];
        }
    }

    for (int i = 0; i < this->rows; i++)
        x[i] /= vec[(matrix.cols + 1) * i];

    return array(x, this->get_rows(), 1);
}

void array::print(std::string delimeter, int precision) const
{
    for (int i = 0; i < this->rows; i++)
    {
        std::cout.precision(precision);
        for (int j = 0; j < this->cols; j++)
            std::cout << (*this)(i, j) << delimeter;
        std::cout << "\n";
    }
}

std::ostream &operator<<(std::ostream &out, const array &arr)
{
    for (int i = 0; i < arr.get_rows(); i++)
    {
        for (int j = 0; j < arr.get_cols(); j++)
            out << std::setw(12) << arr(i, j);
        out << "\n";
    }
    return out;
}

void array::print(int precision) const
{
    for (int i = 0; i < this->rows; i++)
    {
        std::cout.precision(precision);
        for (int j = 0; j < this->cols; j++)
            std::cout << (*this)(i, j) << " ";
        std::cout << "\n";
    }
}

double array::det() const
{
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix should be square.";
        exit(-1);
    }

    array upper = this->Utriangular();

    double determinant = 1;
    for (int i = 0; i < this->cols; i++)
        determinant *= upper(i, i);

    return determinant;
}

array array::Utriangular() const
{
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix should be square.";
        exit(-1);
    }

    std::vector<double> vec{this->vector};

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = i + 1; j < this->rows; j++)
        {
            double constant = vec[this->cols * j + i] / vec[(this->cols + 1) * i];
            for (int k = 0; k < this->cols; k++)
                vec[this->cols * j + k] -= constant * vec[this->cols * i + k];
        }
    }
    return array(vec, this->rows, this->cols);
}

array array::Ltriangular() const
{
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix should be square.";
        exit(-1);
    }

    std::vector<double> vec{this->vector};

    for (int i = this->rows - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            double constant = vec[this->cols * j + i] / vec[(this->cols + 1) * i];
            for (int k = 0; k < this->cols; k++)
                vec[this->cols * j + k] = vec[this->cols * j + k] - constant * vec[this->cols * i + k];
        }
    }
    return array(vec, this->rows, this->cols);
}

array array::aggregate(int axis, double (*function)(double, double)) const
{
    if (axis == ROW)
    {
        std::vector<double> vec(this->rows);
        for (int i = 0; i < this->rows; i++)
        {
            double res = (*this)(i, 0);
            for (int j = 1; j < this->cols; j++)
            {
                res = function(res, (*this)(i, j));
            }
            vec[i] = res;
        }
        return array(vec, this->rows, 1);
    }
    else if (axis == COL)
    {
        std::vector<double> vec(this->cols);
        for (int i = 0; i < this->cols; i++)
        {
            double res = (*this)(0, i);
            for (int j = 1; j < this->rows; j++)
            {
                res = function(res, (*this)(j, i));
            }
            vec[i] = res;
        }
        return array(vec, 1, this->cols);
    }
    else
    {
        std::cerr << "axis should be either COL or ROW";
        exit(-1);
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
                product += (*this)(i, k) * matrix(k, j);
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

array array::operator>(double number) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] > number;
    }
    return array(vec, this->rows, this->cols);
}

array array::operator<(double number) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] < number;
    }
    return array(vec, this->rows, this->cols);
}

array array::operator>=(double number) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] >= number;
    }
    return array(vec, this->rows, this->cols);
}

array array::operator<=(double number) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] <= number;
    }
    return array(vec, this->rows, this->cols);
}

array array::operator==(double number) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] == number;
    }
    return array(vec, this->rows, this->cols);
}

array array::operator!=(double number) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] != number;
    }
    return array(vec, this->rows, this->cols);
}

array array::operator>(const array &matrix) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] > matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator<(const array &matrix) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] < matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator>=(const array &matrix) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] >= matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator<=(const array &matrix) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] <= matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator==(const array &matrix) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] == matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator!=(const array &matrix) const
{
    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] != matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

double array::at(int index) const
{
    return this->vector[index];
}

double array::at(int row, int col) const
{
    return this->vector[cols * row + col];
}

double array::operator()(int index) const
{
    return this->vector[index];
}

double array::operator()(int row, int col) const
{
    return this->vector[cols * row + col];
}

array array::operator()(std::vector<int> row_range, std::vector<int> col_range) const
{
    int rows = (row_range[1] - row_range[0]), cols = (col_range[1] - col_range[0]);
    std::vector<double> vec(rows * cols);

    for (int i = row_range[0]; i < row_range[1]; i++)
    {
        for (int j = col_range[0]; j < col_range[1]; j++)
            vec[cols * (i - row_range[0]) + j] = this->vector[i * cols + j];
    }
    return array(vec, rows, cols);
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

std::vector<double> array::get_vector() const
{
    return this->vector;
}
