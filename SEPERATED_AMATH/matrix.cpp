#include "amath.h"

// 'v' for values, rows for the number of rows, cols for the number of cols. Rows wrap according to the dimensions
Matrix::Matrix(const std::vector<double> &v, const int rows, const int cols) : vector(v)
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
Matrix::Matrix(const double uniform_val, const int rows, const int cols) : cols(cols), rows(rows), size(rows * cols)
{
    if (rows < 1 || cols < 1)
        throw std::invalid_argument("\n\t'rows' and 'cols' arguments in constructor should be positive");

    this->vector = std::vector<double>(size, uniform_val);
}

// comma-seperated values are read from the 'file_path'. Dimensions are the same as the layout in the file
Matrix::Matrix(const std::string &file_path)
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

const int Matrix::argmax() const
{
    return std::distance(this->vector.begin(), std::max_element(this->vector.begin(), this->vector.end()));
}

const int Matrix::argmin() const
{
    return std::distance(this->vector.begin(), std::min_element(this->vector.begin(), this->vector.end()));
}

const double Matrix::max() const
{
    return *std::max_element(this->vector.begin(), this->vector.end());
}

const double Matrix::min() const
{
    return *std::min_element(this->vector.begin(), this->vector.end());
}

// returns an Matrix with 1 rows with 'partition_count' values ranging from 'start' to 'end'
Matrix linspace(const double start, const double end, const int partition_count)
{
    std::vector<double> vec(partition_count);
    double dx = (end - start) / (partition_count - 1);

    for (int i = 0; i < partition_count; i++)
    {
        vec[i] = start + i * dx;
    }

    return Matrix(vec, 1, partition_count);
}

// returns an Matrix with 1 rows with 'partition_count' values ranging from 'start' to 'end' logarithmically
Matrix logspace(double start, double end, int partition_count)
{
    std::vector<double> vec(partition_count);
    double logdy = (log10(end) - log10(start)) / (partition_count - 1);

    for (int i = 0; i < partition_count; i++)
    {
        vec[i] = start + pow(10, i * logdy) - 1;
    }

    return Matrix(vec, 1, partition_count);
}

// returns an Matrix with 1 rows with the values distanced 'step_size' ranging from 'start' to 'end'
Matrix arange(const double start, const double end, const int step_size)
{
    std::vector<double> vec((end - start) / step_size + 1);

    vec[0] = start;
    for (int i = 0; vec[i - 1] < end; i++)
    {
        vec[i] = start + i * step_size;
    }

    return Matrix(vec, 1, (end - start) / step_size + 1);
}

// Returns true if the Matrix is a square matrix
bool Matrix::isSquare() const
{
    if (this->cols == this->rows)
        return true;
    return false;
}

// Returns true if the Matrix is skew-symmetric matrix
bool Matrix::isSkewSymmetric() const
{
    if (!this->isSquare())
        return false;
    for (int i = 0; i < this->cols; i++)
    {
        for (int j = i; j < this->cols; j++)
        {
            if (i == j)
                continue;
            if (this->at(i, j) != -this->at(j, i))
                return false;
        }
    }

    return true;
}

// Returns true if the Matrix is symmetric matrix
bool Matrix::isSymmetric() const
{
    if (!this->isSquare())
        return false;
    for (int i = 0; i < this->cols; i++)
    {
        for (int j = i; j < this->cols; j++)
        {
            if (i == j)
                continue;
            if (this->at(i, j) != this->at(j, i))
                return false;
        }
    }

    return true;
}

// Returns true if the Matrix contains only 0 or 1
bool Matrix::isBoolean() const
{
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            if (this->at(i, j) != 0 && this->at(i, j) != 1)
                return false;
        }
    }

    return true;
}

// Returns true if the Matrix contains non-zero values only on the diagonal
bool Matrix::isDiagonal() const
{
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            if (i == j)
                continue;
            if (this->at(i, j) != 0)
                return false;
        }
    }

    return true;
}

// LU Factorization
std::vector<Matrix> Matrix::LUFactor() const
{
    if (this->rows != this->cols)
        throw std::invalid_argument("\n\tLU factorization requires square a matrix.");

    Matrix LMatrix = diagonal(1, this->rows, this->cols);
    std::vector<double> vec{this->vector};

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = i + 1; j < this->rows; j++)
        {
            double constant = vec[this->cols * j + i] / vec[(this->cols + 1) * i];
            LMatrix.assign(constant * LMatrix(this->cols * i + i), j, i);

            for (int k = 0; k < this->cols; k++)
                vec[this->cols * j + k] -= constant * vec[this->cols * i + k];
        }
    }

    return std::vector<Matrix>{LMatrix, Matrix(vec, this->rows, this->cols)};
}

// returns the cofactor C_{ij}
double Matrix::cofactor(int i, int j) const
{
    if (i < 0 || j < 0)
        throw std::invalid_argument("i and j should be nonnegative for the 'minor' function. (ERROR IN 'minor' FUNCTION)");
    else if (i >= this->rows || j >= this->cols)
        throw std::invalid_argument("i and j should be less then the number of cols and rows for the 'minor' function. (ERROR IN 'minor' FUNCTION)");

    std::vector<double> vals;

    for (int x = 0; x < this->rows; x++)
    {
        if (x == i)
            continue;
        for (int y = 0; y < this->cols; y++)
        {
            if (y == j)
                continue;
            vals.push_back(this->at(x, y));
        }
    }

    if ((i + j) % 2 == 0)
        return Matrix(vals, this->rows - 1, this->rows - 1).det();
    else
        return -Matrix(vals, this->rows - 1, this->rows - 1).det();
}

// returns the minor M_{ij}
double Matrix::minor(int i, int j) const
{
    if (i < 0 || j < 0)
        throw std::invalid_argument("i and j should be nonnegative for the 'minor' function. (ERROR IN 'minor' FUNCTION)");
    else if (i >= this->rows || j >= this->cols)
        throw std::invalid_argument("i and j should be less then the number of cols and rows for the 'minor' function. (ERROR IN 'minor' FUNCTION)");

    std::vector<double> vals;
    for (int x = 0; x < this->rows; x++)
    {
        if (x == i)
            continue;
        for (int y = 0; y < this->cols; y++)
        {
            if (y == j)
                continue;
            vals.push_back(this->at(x, y));
        }
    }
    return Matrix(vals, this->rows - 1, this->rows - 1).det();
}

// returns the matrix of minors
Matrix Matrix::minorMatrix() const
{
    std::vector<double> minors(this->cols * this->rows);
    std::vector<double> dets((this->cols - 1) * (this->rows - 1));

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            minors[this->cols * i + j] = this->minor(i, j);
        }
    }

    return Matrix(minors, this->rows, this->cols);
}

// returns the matrix of cofactors
Matrix Matrix::cofactorMatrix() const
{
    std::vector<double> minors(this->cols * this->rows);
    std::vector<double> dets((this->cols - 1) * (this->rows - 1));

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            minors[this->cols * i + j] = this->cofactor(i, j);
        }
    }

    return Matrix(minors, this->rows, this->cols);
}

// returns the sum of diagonal values. Even if the matrix is not square it acts as if the matrix is square with dimensions 'least(rows,cols)'
Matrix Matrix::trace() const
{
    std::vector<double> vec(this->cols);
    double least = fmin(this->cols, this->rows);
    for (int i = 0; i < least; i++)
        vec[i] = this->at(i, i);

    return Matrix(vec, 1, least);
}

// returns an Matrix with diagonal as values. Even if the matrix is not square it acts as if matrix is square with dimensions 'least(rows,cols)'
Matrix diagonal(double value, int rows, int cols)
{
    std::vector<double> vec(cols * rows, 0);
    int least = fmin(rows, cols);
    for (int i = 0; i < least; i++)
        vec[(cols + 1) * i] = 1;

    return Matrix(vec, rows, cols);
}

// returns the inverse matrix of the square matrix
Matrix Matrix::inverse() const
{
    if (this->get_rows() != this->get_cols())
        throw std::invalid_argument("\n\tMatrix should be square to have an inverse matrix");

    Matrix minverse = diagonal(1, this->get_rows(), this->get_cols());
    std::vector<double> vec{this->vector};
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

    return Matrix(x, this->get_rows(), this->get_cols());
}

// returns an Matrix with dimensions 'rows', 'cols'. Values are distributed from the 'start' value to the 'end'
Matrix random(double start, double _end, int rows, int cols)
{
    // https://stackoverflow.com/questions/21516575/fill-a-vector-with-random-numbers-c
    std::random_device rnd_device;
    std::mt19937 mersenne_engine{rnd_device()};
    std::uniform_real_distribution<double> dist{start, _end};

    auto gen = [&dist, &mersenne_engine]()
    { return dist(mersenne_engine); };

    std::vector<double> vec(cols * rows);
    generate(begin(vec), end(vec), gen);

    return Matrix(vec, rows, cols);
}

// writes comma-seperated values. Layout is set according to the 'rows', 'cols'
void Matrix::write_file(const std::string &file_path) const
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
Matrix read_file(const std::string &file_path)
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

    return Matrix(vec, _rows, vec.size() / _rows);
}

// Least Squares method for polynomial curve fitting. returns 'Matrix' containing the coefficients of a polynomial. 'data' should have 2 columns. If the number of data points is equal to the 'order_of_polynomial' then the outcome is a polynomial crossing all the points.
Matrix leastSquares(const Matrix &data, const int order_of_polynomial)
{
    if (data.get_cols() != 2)
        throw std::invalid_argument("\n\tNumber of cols of the 'data' Matrix should be 2");
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

    Matrix b(_b, order_of_polynomial + 1, 1);
    Matrix A(a, order_of_polynomial + 1, order_of_polynomial + 1);

    return b.solveSquare(A);
}

// solves the Ax = b equation for square A.
Matrix Matrix::solveSquare(const Matrix &matrix) const
{
    if (matrix.get_rows() != matrix.get_cols())
        throw std::invalid_argument("Matrix should be square.");

    else if (this->rows != matrix.get_rows())
        throw std::invalid_argument("Number of rows of b should be equal to the number of columns of the matrix");

    std::vector<double> vec{matrix.get_vector()};
    std::vector<double> x{this->vector};

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

    return Matrix(x, this->get_rows(), 1);
}

// prints with the delimeter between ever value with the specidifed precision
void Matrix::print(std::string delimeter, int precision) const
{
    for (int i = 0; i < this->rows; i++)
    {
        std::cout.precision(precision);
        for (int j = 0; j < this->cols; j++)
            std::cout << (*this)(i, j) << delimeter;
        std::cout << "\n";
    }
}

std::ostream &operator<<(std::ostream &out, const Matrix &arr)
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
void Matrix::print(int precision) const
{
    for (int i = 0; i < this->rows; i++)
    {
        std::cout.precision(precision);
        for (int j = 0; j < this->cols; j++)
            std::cout << (*this)(i, j) << " ";
        std::cout << "\n";
    }
}

// returns the determinant of the Matrix
double Matrix::det() const
{
    if (this->rows != this->cols)
        throw std::invalid_argument("Matrix should be square to compute the determinant");

    Matrix upper = this->Utriangular();

    double determinant = 1;
    for (int i = 0; i < this->cols; i++)
        determinant *= upper(i, i);

    return determinant;
}

// returns an upper triangular matrix
Matrix Matrix::Utriangular() const
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
    return Matrix(vec, this->rows, this->cols);
}

// returns an lower triangular matrix
Matrix Matrix::Ltriangular() const
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
    return Matrix(vec, this->rows, this->cols);
}

// aggregates the axis with the given function
Matrix Matrix::aggregate(int axis, double (*function)(double, double)) const
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
        return Matrix(vec, this->rows, 1);
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
        return Matrix(vec, 1, this->cols);
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
        return Matrix(vec, 1, 1);
    }
    else
        throw std::invalid_argument("axis should be either COL, ROW or BOTH");
}

// returns the transpose of the Matrix
Matrix Matrix::transpose() const
{
    std::vector<double> vec{this->vector};
    for (int i = 0; i < this->cols; i++)
    {
        for (int j = 0; j < this->rows; j++)
            vec[j * this->cols + i] = this->vector[i * this->cols + j];
    }
    return Matrix(vec, this->cols, this->rows);
}

// returns an Matrix with values of value by value application of 'function' to the (*this)
Matrix Matrix::apply(double (*function)(double)) const
{
    std::vector<double> res;
    for (int i = 0; i < this->size; i++)
        res.push_back(function(this->vector[i]));
    return Matrix(res, this->rows, this->cols);
}

// returns an Matrix. Pairwise operations between (*this) and the 'second_Matrix'
Matrix Matrix::pairWise(double (*function)(double, double), const Matrix &second_Matrix) const
{
    if (second_Matrix.get_rows() != this->rows || second_Matrix.get_cols() != this->cols)
        throw std::invalid_argument("\n\tNumber of cols and rows of the matrices should be equal for pairwise functionality");

    std::vector<double> res(this->size);
    for (int i = 0; i < this->size; i++)
        res[i] = function(this->at(i), second_Matrix.at(i));

    return Matrix(res, this->rows, this->cols);
}

// row and col indexing
Matrix Matrix::operator()(const std::vector<int> &row_range, const std::vector<int> &col_range) const
{
    if (row_range.size() != 2 || col_range.size() != 2)
        throw std::invalid_argument("row_range and col_range should be vectors with 2 values");
    else if (row_range[0] < 0 || row_range[1] < 0 || col_range[0] < 0 || col_range[1] < 0)
        throw std::invalid_argument("row_range and col_range should be vectors with poisitive values");

    int rows = (row_range[1] - row_range[0]), cols = (col_range[1] - col_range[0]);

    std::vector<double> vec;
    vec.reserve(rows * cols);

    for (int i = row_range[0]; i < row_range[1]; i++)
    {
        for (int j = col_range[0]; j < col_range[1]; j++)
            vec.emplace_back(this->at(i, j));
    }
    return Matrix(vec, rows, cols);
}

Matrix Matrix::operator*(const double &number) const
{
    std::vector<double> res;
    res.reserve(size);
    for (int i = 0; i < size; i++)
        res.emplace_back(this->vector[i] * number);
    return Matrix(res, this->rows, this->cols);
}

Matrix Matrix::operator/(const double &number) const
{
    std::vector<double> res;
    res.reserve(size);
    for (int i = 0; i < size; i++)
        res.emplace_back(this->vector[i] / number);
    return Matrix(res, this->rows, this->cols);
}

Matrix Matrix::operator-(const double &number) const
{
    std::vector<double> res;
    res.reserve(size);
    for (int i = 0; i < size; i++)
        res.emplace_back(this->vector[i] - number);
    return Matrix(res, this->rows, this->cols);
}

Matrix Matrix::operator+(const double &number) const
{
    std::vector<double> res;
    res.reserve(size);
    for (int i = 0; i < size; i++)
        res.push_back(this->vector[i] + number);
    return Matrix(res, this->rows, this->cols);
}

// no threading multiplication, add checks if the element is 0
Matrix Matrix::operator*(const Matrix &matrix) const
{
    if (this->cols != matrix.rows)
        throw std::invalid_argument("number of cols of the first one should be equal to the number of rows of the second one.");

    std::vector<double> res;
    res.reserve(matrix.cols * this->rows);
    double product = 0;

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            for (int k = 0; k < this->cols; k++)
                product += (*this)(i, k) * matrix(k, j);
            res.emplace_back(product);
            product = 0;
        }
    }

    return Matrix(res, this->rows, matrix.cols);
}

// Tried not to access 'res' and locking at the end to copy still slow (used nice example with 2 threads)
/* void Matrix::dotProduct(int i, const Matrix *_this, const Matrix *matrix, std::vector<double> *res)
{
    //std::vector<double> results(matrix->get_cols());
    for (int j = 0; j < matrix->get_cols(); j++)
    {
        double sum = 0.;
        for (int x = 0; x < _this->get_cols(); x++)
        {
            sum += _this->at(i, x) * matrix->at(x, j);
        }
        (*res)[matrix->get_cols() * i + j] = sum;
    }
}

//better algorithm for block multiplications find the sweet spot for threading
Matrix Matrix::operator*(const Matrix &matrix, int thread_count) const
{
    if (this->cols != matrix.get_rows())
        throw std::invalid_argument("for matrix multiplication number of cols of the first one should be equal to the number of rows of the second one.");

    std::vector<double> res(this->rows * matrix.get_cols());
    std::vector<std::thread> threads(this->rows);

    for (int i = 0; i < this->rows; i++)
    {
        threads[i] = std::thread(Matrix::dotProduct, i, this, &matrix, &res);
    }

    for (int i = 0; i < this->rows; i++)
    {
        threads[i].join();
    }

    return Matrix(res, this->rows, matrix.get_cols());
}  */

Matrix Matrix::operator/(const Matrix &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
        throw std::invalid_argument("Number of cols and rows should match");

    std::vector<double> res;
    res.reserve(this->size);

    for (int i = 0; i < this->size; i++)
        res.emplace_back(this->vector[i] / matrix.vector[i]);
    return Matrix(res, this->rows, this->cols);
}

Matrix Matrix::operator-(const Matrix &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
        throw std::invalid_argument("Number of cols and rows should match");

    std::vector<double> res;
    res.reserve(this->size);

    for (int i = 0; i < size; i++)
        res.emplace_back(this->vector[i] - matrix.vector[i]);
    return Matrix(res, this->rows, this->cols);
}

Matrix Matrix::operator+(const Matrix &matrix) const
{
    if (this->rows != matrix.rows || this->cols != matrix.cols)
        throw std::invalid_argument("Number of cols and rows should match");

    std::vector<double> res;
    res.reserve(this->size);

    for (int i = 0; i < size; i++)
        res.emplace_back(this->vector[i] + matrix.vector[i]);
    return Matrix(res, this->rows, this->cols);
}

Matrix Matrix::operator>(const double &number) const
{
    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] > number);
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator<(const double &number) const
{
    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] < number);
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator>=(const double &number) const
{
    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec[i] = this->vector[i] >= number;
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator<=(const double &number) const
{
    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] <= number);
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator==(const double &number) const
{
    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] == number);
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator!=(const double &number) const
{
    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] != number);
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator>(const Matrix &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the Matrixs should match for comparison");

    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] > matrix(i));
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator<(const Matrix &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the Matrixs should match for comparison");

    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] < matrix(i));
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator>=(const Matrix &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the Matrixs should match for comparison");

    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] >= matrix(i));
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator<=(const Matrix &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the Matrixs should match for comparison");

    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] <= matrix(i));
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator==(const Matrix &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the Matrixs should match for comparison");

    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] == matrix(i));
    }
    return Matrix(vec, this->rows, this->cols);
}

Matrix Matrix::operator!=(const Matrix &matrix) const
{
    if (this->rows * this->cols != matrix.rows * matrix.cols)
        throw std::invalid_argument("Sizes of the Matrixs should match for comparison");

    std::vector<double> vec;
    vec.reserve(this->size);

    for (int i = 0; i < this->size; i++)
    {
        vec.emplace_back(this->vector[i] != matrix(i));
    }
    return Matrix(vec, this->rows, this->cols);
}

void plot(const Matrix &Matrix, const std::vector<double> &xlim, const std::vector<double> &ylim, const std::string &file_name)
{
    if (Matrix.get_cols() != 2)
        throw std::invalid_argument("To plot an Matrix, Matrix must have 2 cols");
    std::vector<double> x{Matrix({0, Matrix.get_rows()}, {0, 1}).get_vector()};
    std::vector<double> y{Matrix({0, Matrix.get_rows()}, {1, 2}).get_vector()};

    matplotlibcpp::figure_size(600, 400);

    matplotlibcpp::plot(x, y, "ro");

    matplotlibcpp::xlim(xlim[0], xlim[1]);
    matplotlibcpp::ylim(ylim[0], ylim[1]);

    if (!file_name.empty())
    {
        matplotlibcpp::save(file_name);
    }

    matplotlibcpp::show();
}

void plotModel(const Matrix &model, const Matrix &data, const Matrix &range, const std::vector<double> &xlim, const std::vector<double> &ylim, const std::string &file_name)
{
    if (model.get_cols() != 1)
        throw std::runtime_error("plotModel Matrix must have 1 columns");

    std::vector<double> vec(range.get_size());

    for (int j = 0; j < range.get_size(); j++)
    {
        double y_val = 0;
        for (int i = 0; i < model.get_size(); i++)
        {
            y_val += model(model.get_size() - i - 1) * pow(range(j), i);
        }
        vec[j] = y_val;
    }

    matplotlibcpp::figure_size(600, 400);

    matplotlibcpp::plot(range.get_vector(), vec, "b");
    matplotlibcpp::scatter(data({0, data.get_rows()}, {0, 1}).get_vector(), data({0, data.get_rows()}, {1, 2}).get_vector(), 10.);

    matplotlibcpp::xlim(xlim[0], xlim[1]);
    matplotlibcpp::ylim(ylim[0], ylim[1]);

    if (!file_name.empty())
    {
        matplotlibcpp::save(file_name);
    }

    matplotlibcpp::show();
}