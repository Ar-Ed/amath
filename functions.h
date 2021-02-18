#pragma once
#include <cmath>

double numericIntegral(double coefficient,double (*function)( double ), double l_bound, double u_bound, int partition);

double numericDerivative(double coefficient, double function(double), double x_coordinate, int precision);

double special_func(double x);

double exp_func(double x);

double gauss_integrand(double x);