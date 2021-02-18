#pragma once
#include <iostream>
#include <vector>

using namespace std;

class Matrix{
    public:
        Matrix(int cc, int rc, const vector<double>& array);
         
        double size()const;
        
        Matrix transpose()const;
        
        Matrix reShape(int col, int row) const;

        Matrix operator*(const Matrix& m)const;
        
        Matrix operator*(double num) const;
        
        Matrix operator+(const Matrix& m)const;
        
        Matrix operator+(double num) const;

        Matrix operator-(const Matrix& m)const;

        Matrix operator-(double num) const;
        
        friend ostream& operator<<(ostream& out,const Matrix& m);
        
    private:
        int col_count = 0;
        int row_count = 0;
        vector<vector<double> > matrix; // valarray is an option
};