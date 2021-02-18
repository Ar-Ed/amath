#include "matrix.h"

int col_count = 0;
int row_count = 0;
vector<vector<double> > matrix; // valarray is an option

Matrix::Matrix(int cc, int rc, const vector<double>& array):col_count(cc), row_count(rc),matrix(rc,vector<double>(cc, 0)){
    for(int i = 0; i < row_count;i++){
        for(int j = 0; j < col_count; j++){
            matrix[i][j] = array[i*col_count + j];
        }
    }
} 
    
double Matrix::size()const{return (row_count*col_count);}

Matrix Matrix::transpose()const{
    vector<double> output;
    for(int i = 0; i < col_count; i++){
        for(int j = 0 ; j < row_count; j++){
            output.push_back(matrix[j][i]);
        }
    }
    return Matrix(row_count, col_count, output);
}

Matrix Matrix::reShape(int col, int row) const{
    vector<double> output;
    for(int i = 0; i < row_count; i++){
        for(int j = 0 ; j < col_count; j++){
            output.push_back(matrix[i][j]);
        }
    }
    return Matrix(col,row, output);
}

Matrix Matrix::operator*(const Matrix& m)const{
    vector<double> output(row_count*m.col_count, 0);
    for(int i= 0; i < row_count; i++){
        for(int j = 0; j < m.col_count; j++){
            for(int k = 0; k < m.row_count; k++){
                output[i*m.col_count+j] += (matrix[i][k] * m.matrix[k][j]);
            }
        }
    }
    return Matrix(row_count, m.col_count, output);
}

Matrix Matrix::operator*(double num) const{
    vector<double> output;
    for(int i = 0; i < row_count; i++){
        for(int j = 0 ; j < col_count; j++){
            output.push_back(matrix[i][j]*num);
        }
    }
    return Matrix(col_count, row_count, output);
}

Matrix Matrix::operator+(const Matrix& m)const{
    vector<double> output;
    for(int i = 0; i < row_count; i++){
        for(int j = 0 ; j < col_count; j++){
            output.push_back(matrix[i][j] + m.matrix[i][j]);
        }
    }
    return Matrix(col_count, row_count, output);
}

Matrix Matrix::operator+(double num) const{
    vector<double> output;
    for(int i = 0; i < row_count; i++){
        for(int j = 0 ; j < col_count; j++){
            output.push_back(matrix[i][j] + num);
        }
    }
    return Matrix(col_count, row_count, output);
}

Matrix Matrix::operator-(const Matrix& m)const{
    vector<double> output;
    for(int i = 0; i < row_count; i++){
        for(int j = 0 ; j < col_count; j++){
            output.push_back(matrix[i][j] - m.matrix[i][j]);
        }
    }
    return Matrix(col_count, row_count, output);
}

Matrix Matrix::operator-(double num) const{
    vector<double> output;
    for(int i = 0; i < row_count; i++){
        for(int j = 0 ; j < col_count; j++){
            output.push_back(matrix[i][j] - num);
        }
    }
    return Matrix(col_count, row_count, output);
}

ostream& operator<<(ostream& out,const Matrix& m) {
    for(int i = 0; i < m.row_count; i++){
        for(int j = 0; j < m.col_count; j++){
            out << " " << m.matrix[i][j];
        }out << "\n";
    }
    return out;
}