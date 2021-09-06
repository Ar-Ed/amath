#include <iostream>
#include <vector>
#include <string>

struct array
{
    array(std::vector<double> v, int rows, int cols) : vector(v), cols(cols), rows(rows), size(cols * rows) {}
    void print(std::string delimeter = " ")
    {
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
                std::cout << this->at(i, j) << delimeter;
            std::cout << "\n";
        }
    }

    array transpose()
    {
        return array(this->vector, this->cols, this->rows);
    }

    array operator*(double number)
    {
        std::vector<double> res;
        for (int i = 0; i < size; i++)
            res.push_back(this->vector[i] * number);
        return array(res, this->rows, this->cols);
    }

    array operator/(double number)
    {
        std::vector<double> res;
        for (int i = 0; i < size; i++)
            res.push_back(this->vector[i] / number);
        return array(res, this->rows, this->cols);
    }

    array operator-(double number)
    {
        std::vector<double> res;
        for (int i = 0; i < size; i++)
            res.push_back(this->vector[i] - number);
        return array(res, this->rows, this->cols);
    }

    array operator+(double number)
    {
        std::vector<double> res;
        for (int i = 0; i < size; i++)
            res.push_back(this->vector[i] + number);
        return array(res, this->rows, this->cols);
    }

    array operator*(const array &matrix)
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

    array operator/(const array &matrix)
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

    array operator-(const array &matrix)
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

    array operator+(const array &matrix)
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
    array copy() const
    {
        return array(this->vector, this->rows, this->cols);
    }

    void reShape(int rows, int cols)
    {
        this->rows = rows;
        this->cols = cols;
    }
    double at(int row, int col) const
    {
        return this->vector[cols * row + col];
    }
    int get_size()
    {
        return this->size;
    }
    int get_cols()
    {
        return this->cols;
    }
    int get_rows()
    {
        return this->rows;
    }

private:
    std::vector<double> vector;
    int cols;
    int rows;
    int size;
};

int main()
{
    array m1({1, 2, 3, 5, 5, 6, 3, 2, 32, 4}, 2, 5);

    array m2 = m1 + m1;
    m2.print(" ");

    std::cout << "\n";

    m1 = m2.copy();
    m1.print(" ");

    std::cout << "\n";

    (m1 * 3 - m1 / 3).print(" ");

    std::cout << "\n";

    array m3 = m1 * m1.transpose();
    m3.print(" ");

    //array m2({1, 2, 3, 5, 5, 6, 3, 2, 32, 4}, 2, 5);

    //array m3 = m1 + m2;

    //std::cout << m1.get_cols() << " " << m1.get_rows() << " " << m1.get_size();
    return 0;
}
