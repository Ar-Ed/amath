#pragma once
#include <vector>
#include <iostream>
#include "point.h"

using namespace std;

class Poly{
    public:
        Poly(vector<vector<double> > array);
        void setC(int term, double constant);
        void setO(int term, double order);

        double getO(int term)const;
        double getC(int term)const;
        
        int termCount()const;
        
        Point evaluatePoint(double x_coordinate)const;

        double definiteIntegral(double l_bound, double u_bound)const;
        
        double definiteDerivative(double x_coordinate)const;
        
        Poly integral();
        
        Poly derivative();

        Poly operator+(const Poly& p)const;

        Poly operator-(const Poly& p)const;
        
        friend ostream& operator<<(ostream& out, const Poly& p);
        
    private:
        vector<vector<double> > terms;
        int term_count;
};