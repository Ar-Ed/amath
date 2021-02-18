#include "functions.h"

double numericIntegral(double coefficient,double (*function)( double ), double l_bound, double u_bound, int partition){
    double sum = 0;
    for(int index = 0; index < partition; index++){
        sum+= (u_bound-l_bound) * ((coefficient*(function)(l_bound + (u_bound-l_bound)*index/partition)) + (coefficient*(function)(l_bound + (u_bound-l_bound)*(index+1)/partition)))/partition/2;
    }
    return sum;
}

double numericDerivative(double coefficient, double function(double), double x_coordinate, int precision){
    double step_size = 1;
    while(abs(function(x_coordinate - step_size) - function(x_coordinate + step_size)) > pow(10, -1 * precision)){
        step_size /= 10;
    }
    return (function(x_coordinate + step_size) - function(x_coordinate-step_size))/(2*step_size);
}
// can define preferred functions and use them with numeric integral and derivatives.
double special_func(double x){ // an example of preferred function
    return cos(x)*exp(x);
}

double exp_func(double x){
    return pow(x,2);
}

double gauss_integrand(double x){
    return pow(exp(1), pow(x,2));
}