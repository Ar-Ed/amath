#include "shape.h"

int side_count = 0;
double side_length = 0;
double height = 0;
double area = 0;
double volume = 0;
double longest_line = 0;
static constexpr double PI = 3.141592;

Shape::Shape(int side_count, double side_length):side_count(side_count), side_length(side_length){
    area = pow(side_length,2)*side_count/4/tan(PI/side_count);
    if(side_count % 2 == 0){
        longest_line = 2 * sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2));   
    }else longest_line = sqrt(pow(side_length/2/tan(PI/side_count),2)+ pow(side_length/2,2)) * sqrt(2*(1-cos(PI*(side_count+1)/side_count)));
    
}

Shape::Shape(int side_count, double side_length, double height):side_count(side_count), side_length(side_length),height(height){
    volume = pow(side_length,2)*side_count/4/tan(PI/side_count) * height;
    area = volume * 2 /height + side_count * side_length * height;
    if(side_count % 2 == 0){
        longest_line = sqrt(pow(2 * sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2)), 2) + pow(height,2));    
    }else longest_line = sqrt(pow(side_length/2/tan(PI/side_count),2)+ pow(side_length/2,2)) * sqrt(2*(1-cos(PI*(side_count+1)/side_count)));
}
Shape::Shape(Shape& sh, double height):side_count(sh.side_count), side_length(sh.side_length), height(height){
    volume = pow(side_length,2)*side_count/4/tan(PI/side_count) * height;
    area = volume * 2 /height + side_count * side_length * height;
    if(side_count % 2 == 0){
        longest_line = sqrt(pow(2 * sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2)), 2) + pow(height,2));    
    }else longest_line = sqrt(pow(sqrt(pow((side_length/2/tan(PI/side_count)),2)+ pow(side_length/2,2)) * sqrt(2*(1-cos(PI*(side_count+1)/side_count))),2) + pow(height,2));
}

double Shape::Area(){
    return area;
}

double Shape::Volume(){
    if(volume){return volume;}
    return 0;
}

double Shape::Angle(){
    return (PI*(1-2/side_count));
}

double Shape::longestLine(){
    return longest_line;
}