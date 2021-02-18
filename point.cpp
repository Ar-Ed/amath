//operator+ lesson to learn -> (const *** & p) and (*** p) both has place and need overlaod 
//SOLUTION NO NEED TO OVERLOAD use (const type& variable_name)
#include "point.h"

vector<double>points;
int axis_count = 0;
int Point::number_of_points = 0;
static constexpr double PI = 3.141592;

Point::Point(const vector<double>&points):points(points){
    number_of_points++; 
    axis_count =  points.size();  
} // vector of points
Point::~Point(){number_of_points--;}

double Point::getAxis(int index)const{return points[index];}

int Point::axisCount()const{return axis_count;}

void Point::setAxis(int index, double value){points[index] = value;}

int Point::pointCount(){return number_of_points;}

Point Point::operator+(const Point &p)const{ // IF YOU TAKE POINT BY REFERENCE IT WON'T EXCEPT NON INITILIAZED POINTS MIGHT NEED TO OVERLOAD
    vector<double> output;
    if(p.axis_count > axis_count){
        for(int i = 0; i <axis_count; i++){
            output.push_back(p.points[i]+ points[i]);
        }for(int i = axis_count; i < p.axis_count;i++){
            output.push_back(p.points[i]);
        }
    }else{
        for(int i = 0; i <p.axis_count; i++){
            output.push_back(p.points[i]+ points[i]);
        }for(int i = p.axis_count; i < axis_count;i++){
            output.push_back(points[i]);
        }
    }
    return Point(output);
}

Point Point::operator-(const Point&p)const{
    vector<double> output;
    if(p.axis_count > axis_count){
        for(int i = 0; i <axis_count; i++){
            output.push_back(-p.points[i]+ points[i]);
        }for(int i = axis_count; i < p.axis_count;i++){
            output.push_back(-p.points[i]);
        }
    }else{
        for(int i = 0; i <p.axis_count; i++){
            output.push_back(-p.points[i]+ points[i]);
        }for(int i = p.axis_count; i < axis_count;i++){
            output.push_back(points[i]);
        }
    }
    return Point(output);
}

void Point::operator++(){
    for(int i = 0; i < axis_count; i++)++points[i];
}
void Point::operator++(int){
    for(int i = 0; i < axis_count; i++)points[i]++;
}
void Point::operator--(int){
    for(int i = 0; i < axis_count; i++)points[i]--;
}
void Point::operator--(){
    for(int i = 0; i < axis_count; i++)--points[i];
}

Point Point::rotate(int axis1 , int axis2, double degree){
    vector<double> output = points;
    double temp = points[axis1];
    output[axis1] = -1*points[axis2]*sin(degree*PI/180) + points[axis1]*cos(degree*PI/180);
    output[axis2] = temp*sin(degree*PI/180) + points[axis2]*cos(degree*PI/180);
    return Point(output);
}

double Point::distance(const Point& p)const{ 
    double sum = 0;
    if(p.axis_count > axis_count){
        for(int i = 0; i <axis_count; i++){
            sum += pow(p.points[i] - points[i],2);
        }
        for(int i = axis_count; i < p.axis_count;i++){
            sum += pow(p.points[i],2);
        }
    }else{
        for(int i = 0; i <p.axis_count; i++){
            sum += pow(p.points[i] - points[i],2);
        }
        for(int i = p.axis_count; i < axis_count;i++){
            sum += pow(points[i],2);
        }
    }
    return sqrt(sum);
}

ostream& operator<<(ostream& out, const Point& p) {
    out << "( ";
    for(double vals : p.points)out << vals << ", ";
    out << ")" << endl;
    return out;
}

istream& operator>>(istream& in, const Point& p){
    for(double vals : p.points)in >> vals;
    return in;
}