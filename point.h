#pragma once
#include <vector>
#include <iostream>

using namespace std;

class Point{
    public:
        Point(const vector<double>&points);
        ~Point();
        
        double getAxis(int index)const;
        
        int axisCount()const;
        
        void setAxis(int index, double value);
        
        static int pointCount();

        Point operator+(const Point &p)const;
        
        Point operator-(const Point&p)const;
        
        void operator++();
        void operator++(int);
        void operator--(int);
        void operator--();
        
        Point rotate(int axis1 , int axis2, double degree);
        
        double distance(const Point& p)const;
        
        friend ostream& operator<<(ostream& out, const Point& p);
        
        friend istream& operator>>(istream& in, const Point& p);

    private:
        vector<double>points;
        int axis_count;
        static int number_of_points;
        static constexpr double PI = 3.141592;
};