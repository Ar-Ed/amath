#pragma once
#include <cmath>
using namespace std;

class Shape{// longest_line formula for side_count % 2 != 0 is wrong
    public:
        Shape(int side_count, double side_length);
        
        Shape(int side_count, double side_length, double height);
        Shape(Shape& sh, double height);
        
        double Area();
        
        double Volume();
        
        double Angle();
        
        double longestLine();
        
    private:
        int side_count;
        double side_length;
        double height;
        double area;
        double volume;
        double longest_line;
        static constexpr double PI = 3.141592;
};