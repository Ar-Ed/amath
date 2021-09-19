#include "amatrix.h"

// 'v' for values, rows for the number of rows, cols for the number of cols. Rows wrap according to the dimensions
array::array(std::vector<double> v, int rows, int cols) : vector(v)
{
    if (rows < 1 || cols < 1)
        throw std::invalid_argument("\n\t'rows' and 'cols' arguments in constructor should be positive");
    else if (v.size() != rows * cols)
        throw std::invalid_argument("\n\tSize of the 'v' should be equal to the product of 'rows' and 'cols'");
        
    this->rows = rows;
    this->cols = cols;
    this->size = cols * rows;
}

// 'uniform_val' for initilization with an uniform value, rows for the number of rows, cols for the number of cols. Rows wrap according to the dimensions
array::array(double uniform_val, int rows, int cols)
{
    if (rows < 1 || cols < 1)
        throw std::invalid_argument("\n\t'rows' and 'cols' arguments in constructor should be positive");
    this->vector = std::vector<double>(this->size, uniform_val);
}

// comma-seperated values are read from the 'file_path'. Dimensions are the same as the layout in the file
array::array(std::string file_path)
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

    if (vec.size() < 1)
        throw std::invalid_argument("\n\tThere are no values in the file");

    this->rows = _rows;
    this->vector = vec;
    this->size = this->vector.size();
    this->cols = this->size / _rows;
    f.close();
}

// LU Factorization
std::vector<array> array::LUFactor() const // fix L part
{
    if (this->rows != this->cols)
        throw std::invalid_argument("\n\tLU factorization requires square matrices.");

    array total = diagonal(1, this->rows, this->cols);
    std::vector<array> elimination_matrices;

    std::vector<double> vec{this->vector};

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = i + 1; j < this->rows; j++)
        {
            double constant = vec[this->cols * j + i] / vec[(this->cols + 1) * i];
            array elim = diagonal(1, this->rows, this->cols);

            elim.assign(constant, j, i);
            elimination_matrices.push_back(elim);

            for (int k = 0; k < this->cols; k++)
                vec[this->cols * j + k] -= constant * vec[this->cols * i + k];
        }
    }
    for (int i = elimination_matrices.size() - 1; i >= 0; i--)
        total = elimination_matrices[i] * total;

    return std::vector<array>{total, array(vec, this->rows, this->cols)};
}

/* array array::minorMatrix() const
{
    std::vector<double> minors(this->cols * this->rows);
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            minors[this->cols * minors + this->rows] = det(array())
        }
    }
} */

// returns the sum of diagonal values. Even if the matrix is not square it acts as if the matrix is square with dimensions 'least(rows,cols)'
array array::trace() const
{
    std::vector<double> vec(this->cols);
    double least = fmin(this->cols, this->rows);
    for (int i = 0; i < least; i++)
        vec[i] = this->at(i, i);

    return array(vec, 1, least);
}

// returns an array with diagonal as values. Even if the matrix is not square it acts as if matrix is square with dimensions 'least(rows,cols)'
array diagonal(double value, int rows, int cols)
{
    std::vector<double> vec(cols * rows, 0);
    int least = fmin(rows, cols);
    for (int i = 0; i < least; i++)
        vec[(cols + 1) * i] = 1;

    return array(vec, rows, cols);
}

// returns the inverse matrix of the square matrix
array array::inverse()
{
    if (this->get_rows() != this->get_cols())
        throw std::invalid_argument("\n\tMatrix should be square to have an inverse matrix");

    array minverse = diagonal(1, this->get_rows(), this->get_cols());
    std::vector<double> vec{this->get_vector()};
    std::vector<double> x{minverse.get_vector()};

    for (int i = 0; i < this->get_rows(); i++)
    {
        for (int j = i + 1; j < this->get_rows(); j++)
        {
            double constant = vec[this->get_cols() * j + i] / vec[(this->get_cols() + 1) * i];
            for (int k = 0; k < this->get_cols(); k++)
            {
                vec[this->get_cols() * j + k] -= constant * vec[this->get_cols() * i + k];
                x[this->get_cols() * j + k] -= constant * x[this->get_cols() * i + k];
            }
        }
    }

    for (int i = this->get_rows() - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            double constant = vec[this->get_cols() * j + i] / vec[(this->get_cols() + 1) * i];
            for (int k = 0; k < this->get_cols(); k++)
            {
                vec[this->get_cols() * j + k] = vec[this->get_cols() * j + k] - constant * vec[this->get_cols() * i + k];
                x[this->get_cols() * j + k] = x[this->get_cols() * j + k] - constant * x[this->get_cols() * i + k];
            }
        }
    }

    for (int i = 0; i < this->rows; i++)
        for (int j = 0; j < this->cols; j++)
        {
            x[i * cols + j] /= vec[(this->cols + 1) * i];
        }

    return array(x, this->get_rows(), this->get_cols());
}

// returns an array with dimensions 'rows', 'cols'. Values are distributed from the 'start' value to the 'end'
array random(double start, double _end, int rows, int cols)
{
    // https://stackoverflow.com/questions/21516575/fill-a-vector-with-random-numbers-c
    std::random_device rnd_device;
    std::mt19937 mersenne_engine{rnd_device()};
    std::uniform_real_distribution<double> dist{start, _end};

    auto gen = [&dist, &mersenne_engine]()
    { return dist(mersenne_engine); };

    std::vector<double> vec(cols * rows);
    generate(begin(vec), end(vec), gen);

    return array(vec, rows, cols);
}

// writes comma-seperated values. Layout is set according to the 'rows', 'cols'
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

// comma-seperated values are read from the 'file_path'. Dimensions are the same as the layout in the file
array read_file(std::string file_path)
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

    if (vec.size() < 1)
        throw std::invalid_argument("\n\tThere are no values in the file");

    f.close();

    return array(vec, _rows, vec.size() / _rows);
}

// Least Squares method for polynomial curve fitting. returns 'array' containing the coefficients of a polynomial. 'data' should have 2 columns. If the number of data points is equal to the 'order_of_polynomial' then the outcome is a polynomial crossing all the points.
array leastSquares(array data, int order_of_polynomial)
{
    if (data.get_cols() != 2)
        throw std::invalid_argument("\n\tNumber of cols of the 'data' array should be 2");
    else if (order_of_polynomial < 1)
        throw std::invalid_argument("\n\t 'order_of_polynomial' should be positive");

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

// solves the Ax = b equation for square A.
array array::solveSquare(array matrix)
{
    if (matrix.get_rows() != matrix.get_cols())
        throw std::invalid_argument("Matrix should be square.");

    else if (this->rows != matrix.get_rows())
        throw std::invalid_argument("Number of rows of b should be equal to the number of columns of the matrix");

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

// prints with the delimeter between ever value with the specidifed precision
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

// prints with the specified precision
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

// returns the determinant of the array
double array::det() const
{
    if (this->rows != this->cols)
        throw std::invalid_argument("Matrix should be square to compute the determinant");

    array upper = this->Utriangular();

    double determinant = 1;
    for (int i = 0; i < this->cols; i++)
        determinant *= upper(i, i);

    return determinant;
}

// returns an upper triangular matrix
array array::Utriangular() const
{
    if (this->rows != this->cols)
        throw std::invalid_argument("Matrix should be square to find an upper triangular form");

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

// returns an lower triangular matrix
array array::Ltriangular() const
{
    if (this->rows != this->cols)
        throw std::invalid_argument("Matrix should be square to find an lower triangular form");

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

// aggregates the axis with the given function
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
    else if (axis == BOTH)
    {
        std::vector<double> vec(1);
        double res = this->at(0);
        for (int j = 1; j < this->size; j++)
        {
            res = function(res, this->at(j));
        }
        vec[0] = res;
        return array(vec, 1, 1);
    }
    else
        throw std::invalid_argument("axis should be either COL, ROW or BOTH");
}

// returns the transpose of the array
array array::transpose() const
{
    std::vector<double> vec{this->vector};
    for (int i = 0; i < this->cols; i++)
    {
        for (int j = 0; j < this->rows; j++)
            vec[j * this->cols + i] = this->vector[i * this->cols + j];
    }
    return array(vec, this->cols, this->rows);
}

// returns an array with values of value by value application of 'function' to the (*this)
array array::apply(double (*function)(double)) const
{
    std::vector<double> res;
    for (int i = 0; i < this->size; i++)
        res.push_back(function(this->vector[i]));
    return array(res, this->rows, this->cols);
}

// returns an array. Pairwise operations between (*this) and the 'second_array'
array array::pairWise(double (*function)(double, double), array second_array) const
{
    if (second_array.get_rows() != this->rows || second_array.get_cols() != this->cols)
        throw std::invalid_argument("\n\tNumber of cols and rows of the matrices should be equal for pairwise functionality");

    std::vector<double> res(this->size);
    for (int i = 0; i < this->size; i++)
        res[i] = function(this->at(i), second_array.at(i));

    return array(res, this->rows, this->cols);
}

// row and col indexing
array array::operator()(std::vector<int> row_range, std::vector<int> col_range) const
{
    if (row_range.size() != 2 || col_range.size() != 2)
        throw std::invalid_argument('row_range and col_range should be vectors with 2 values');
    else if (row_range[0] < 0 || row_range[1] < 0 || col_range[0] < 0 || col_range[1] < 0)
        throw std::invalid_argument('row_range and col_range should be vectors with poisitive values');

    int rows = (row_range[1] - row_range[0]), cols = (col_range[1] - col_range[0]);
    std::vector<double> vec(rows * cols);

    for (int i = row_range[0]; i < row_range[1]; i++)
    {
        for (int j = col_range[0]; j < col_range[1]; j++)
            vec[cols * (i - row_range[0]) + j] = this->vector[i * cols + j];
    }
    return array(vec, rows, cols);
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
        throw std::invalid_argument("number of cols of the first one should be equal to the number of rows of the second one.");

    std::vector<double> res;
    double product = 0;

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            for (int k = 0; k < this->cols; k++)
                product += (*this)(i, k) * matrix(k, j);
            res.push_back(product);
            product = 0;
        }
    }

    return array(res, this->rows, matrix.cols);
}

array array::operator/(const array &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
        throw std::invalid_argument("Number of cols and rows should match");

    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] / matrix.vector[i]);
    return array(res, this->rows, this->cols);
}

array array::operator-(const array &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
        throw std::invalid_argument("Number of cols and rows should match");

    std::vector<double> res;
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] - matrix.vector[i]);
    return array(res, this->rows, this->cols);
}

array array::operator+(const array &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
        throw std::invalid_argument("Number of cols and rows should match");

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
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the arrays should match for comparison");

    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] > matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator<(const array &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the arrays should match for comparison");

    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] < matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator>=(const array &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the arrays should match for comparison");

    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] >= matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator<=(const array &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the arrays should match for comparison");

    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] <= matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator==(const array &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the arrays should match for comparison");

    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] == matrix(i);
    }
    return array(vec, this->rows, this->cols);
}

array array::operator!=(const array &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the arrays should match for comparison");

    std::vector<double> vec(this->size);
    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] != matrix(i);
    }
    return array(vec, this->rows, this->cols);
}