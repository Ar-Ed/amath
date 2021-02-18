//OUTDATED OUTDATED OUTDATED OUTDATED
//OUTDATED OUTDATED OUTDATED OUTDATED
//OUTDATED OUTDATED OUTDATED OUTDATED
//OUTDATED OUTDATED OUTDATED OUTDATED
//OUTDATED OUTDATED OUTDATED OUTDATED
//OUTDATED OUTDATED OUTDATED OUTDATED
//LOOK AT OTHER FILES

//determinant, inverse solution, add dimensions to points

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Shape{// longest_line formula for side_count % 2 != 0 is wrong
    public:
        Shape(int side_count, double side_length):side_count(side_count), side_length(side_length){
            area = pow(side_length,2)*side_count/4/tan(PI/side_count);
            if(side_count % 2 == 0){
                longest_line = 2 * sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2));   
            }else longest_line = sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2)) + side_length/2/tan(PI/side_count);
            
        }
        
        Shape(int side_count, double side_length, double height):side_count(side_count), side_length(side_length),height(height){
            volume = pow(side_length,2)*side_count/4/tan(PI/side_count) * height;
            area = volume * 2 /height + side_count * side_length * height;
            if(side_count % 2 == 0){
                longest_line = sqrt(pow(2 * sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2)), 2) + pow(height,2));    
            }else longest_line = sqrt(pow((sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2)) + side_length/2/tan(PI/side_count)),2) + pow(height,2));
        }
        Shape(Shape& sh, double height):side_count(sh.side_count), side_length(sh.side_length), height(height){
            volume = pow(side_length,2)*side_count/4/tan(PI/side_count) * height;
            area = volume * 2 /height + side_count * side_length * height;
            if(side_count % 2 == 0){
                longest_line = sqrt(pow(2 * sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2)), 2) + pow(height,2));    
            }else longest_line = sqrt(pow((sqrt(pow(side_length/2,2)+pow(side_length/2/tan(PI/side_count), 2)) + side_length/2/tan(PI/side_count)),2) + pow(height,2));
        }
        
        double getArea(){
            return area;
        }
        
        double getVolume(){
            if(volume){return volume;}
            return 0;
        }
        
        double getAngle(){
            return (PI-2*PI/side_count);
        }
        
        double longestLine(){
            return longest_line;
        }
        
    private:
        int side_count = 0;
        double side_length = 0;
        double height = 0;
        double area = 0;
        double volume = 0;
        double longest_line = 0;
        static constexpr double PI = 3.141592;
};

class Point{
    public:
        Point(vector<double>points):points(points){
            number_of_points++; 
            axis_count =  points.size();
            
        } // vector of points
        ~Point(){number_of_points--;}
        
        double getAxis(int index)const{return points[index];}
        
        int axisCount()const{return axis_count;}
        
        void setAxis(int index, double value){points[index] = value;}
        
        static int pointCount(){return number_of_points;}

        Point operator+(Point&p){
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
        
        Point operator-(Point&p){
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
        
        void operator++(){
            for(int i = 0; i < axis_count; i++)++points[i];
        }
        void operator++(int){
            for(int i = 0; i < axis_count; i++)points[i]++;
        }
        void operator--(int){
            for(int i = 0; i < axis_count; i++)points[i]--;
        }
        void operator--(){
            for(int i = 0; i < axis_count; i++)--points[i];
        }
        
        Point rotate(int axis1 , int axis2, double degree){
            vector<double> output = points;
            double temp = points[axis1];
            output[axis1] = -1*points[axis2]*sin(degree*PI/180) + points[axis1]*cos(degree*PI/180);
            output[axis2] = temp*sin(degree*PI/180) + points[axis2]*cos(degree*PI/180);
            return Point(output);
        }
        
        double distance(Point& p){ 
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
        
        friend ostream& operator<<(ostream& out, const Point& p) {
            out << "( ";
            for(double vals : p.points)
                out << vals << ", ";
            out << ")" << endl;
            return out;
        }
        
        friend istream& operator>>(istream& in, Point& p){
            for(double vals : p.points)
                in >> vals;
            return in;
        }
    private:
        vector<double>points;
        int axis_count = 0;
        static int number_of_points;
        static constexpr double PI = 3.141592;
};
int Point::number_of_points = 0;

class Matrix{
    public:
        Matrix(int cc, int rc, vector<double>array):col_count(cc), row_count(rc),matrix(rc,vector<double>(cc, 0)){
            for(int i = 0; i < row_count;i++){
                for(int j = 0; j < col_count; j++){
                    matrix[i][j] = array[i*col_count + j];
                }
            }
        } 
         
        double size()const{return (row_count*col_count);}
        
        Matrix transpose()const{
            vector<double> output;
            for(int i = 0; i < col_count; i++){
                for(int j = 0 ; j < row_count; j++){
                    output.push_back(matrix[j][i]);
                }
            }
            return Matrix(row_count, col_count, output);
        }
        
        Matrix reShape(int col, int row) const{
            vector<double> output;
            for(int i = 0; i < row_count; i++){
                for(int j = 0 ; j < col_count; j++){
                    output.push_back(matrix[i][j]);
                }
            }
            return Matrix(col,row, output);
        }
        
        Matrix operator*(Matrix& m)const{
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
        
        Matrix operator*(double num) const{
            vector<double> output;
            for(int i = 0; i < row_count; i++){
                for(int j = 0 ; j < col_count; j++){
                    output.push_back(matrix[i][j]*num);
                }
            }
            return Matrix(col_count, row_count, output);
        }
        
        Matrix operator+(const Matrix& m)const{
            vector<double> output;
            for(int i = 0; i < row_count; i++){
                for(int j = 0 ; j < col_count; j++){
                    output.push_back(matrix[i][j] + m.matrix[i][j]);
                }
            }
            return Matrix(col_count, row_count, output);
        }
        
        Matrix operator+(double num) const{
            vector<double> output;
            for(int i = 0; i < row_count; i++){
                for(int j = 0 ; j < col_count; j++){
                    output.push_back(matrix[i][j] + num);
                }
            }
            return Matrix(col_count, row_count, output);
        }

        Matrix operator-(const Matrix& m)const{
            vector<double> output;
            for(int i = 0; i < row_count; i++){
                for(int j = 0 ; j < col_count; j++){
                    output.push_back(matrix[i][j] - m.matrix[i][j]);
                }
            }
            return Matrix(col_count, row_count, output);
        }

        Matrix operator-(double num) const{
            vector<double> output;
            for(int i = 0; i < row_count; i++){
                for(int j = 0 ; j < col_count; j++){
                    output.push_back(matrix[i][j] - num);
                }
            }
            return Matrix(col_count, row_count, output);
        }
        
        friend ostream& operator<<(ostream& out,const Matrix& m) {
            for(int i = 0; i < m.row_count; i++){
                for(int j = 0; j < m.col_count; j++){
                    out << " " << m.matrix[i][j];
                }out << "\n";
            }
            return out;
        }
        
    private:
        int col_count = 0;
        int row_count = 0;
        vector<vector<double> > matrix; // valarray is an option
};

class Poly{ // vektorlerle genellestirilebilir // vektorlerle polinomlarlar olusturma ve Point argumaniyla egri beya dogru cizme
    public:
        Poly(double co = 1, double power = 1):coefficient(co),power(power){}
        void setC(double c){coefficient = c;}
        void setP(double power){this->power = power;}

        double getP()const{return power;}
        double getC()const{return coefficient;}

        double numericIntegral(double l_bound, double u_bound, int partition)const{
            double lower_sum = 0, upper_sum = 0;
            for(int index = 0; index<partition; index++){
                lower_sum += coefficient * pow(l_bound + (u_bound-l_bound)*index/partition, power) * (u_bound-l_bound)/partition;
                upper_sum += coefficient * pow(l_bound + (u_bound-l_bound)*(index+1)/partition, power) * (u_bound-l_bound)/partition;
            }
            return (upper_sum+lower_sum)/2;
        }
        
        double definiteIntegral(double l_bound, double u_bound)const{
            return coefficient*(pow(u_bound,power+1)-pow(l_bound,power+1))/(power+1);
        }
        
        double derivative(double x)const{
            return coefficient * power * pow(x,power-1);
        }
        
    private:
        double coefficient;
        double power;
};
            
double numericIntegral(double coefficient,double (*function)( double ), double l_bound, double u_bound, int partition){
    double sum = 0;
    for(int index = 0; index < partition; index++){
        sum+= (u_bound-l_bound) * ((coefficient*(function)(l_bound + (u_bound-l_bound)*index/partition)) + (coefficient*(function)(l_bound + (u_bound-l_bound)*(index+1)/partition)))/partition/2;
    }
    return sum;
}

double numericDerivative(double coefficient, double function(double), double x_coordinate, long precision){
    double step_size = 1;
    while(abs(function(x_coordinate - step_size) - function(x_coordinate + step_size)) > pow(10, -1 * precision)){
        step_size /= 10;
    }
    return (function(x_coordinate + step_size) - function(x_coordinate))/(2*step_size);
}

double special_func(double x){
    return cos(x)*exp(x);
}

double exp_func(double x){
    return pow(x,2);
}

double gauss_func(double x){
    return pow(exp(1), pow(x,2));
}

int main()
{
//    Point p1({1,2,3,4});
//    Point p4({1,2,3});
//    p1++;
//    Point p3 = p1.rotate(1,2,90);
//    cout << p1 << endl; 
//    cout << p3 <<endl;
//    cout << p1.distance(p3) << endl;
//    cout << p1 + p4 << endl;
//    cout << p1.distance(p4);
    
    // Poly poly(5,3);
    // cout << "\n" <<poly.numericIntegral(0,1,10);
    // cout << "\n" <<poly.definiteIntegral(0,1);
    // cout << "\n" <<poly.derivative(2);
    
    // cout << numericIntegral(1, special_func, 0, 1, 5) << endl;
    // cout << numericDerivative(1, special_func, 1, 6) << endl;
    // cout << numericIntegral(1,cos,0,1,10) << endl;
    // cout << numericIntegral(1,exp_func,0,1,10) << endl;
    // cout << numericIntegral(1,gauss_func,0,1,10) << endl;
    // cout << numericDerivative(1,gauss_func,1,5) << endl;

    
    Matrix m1(2,3,{
        1, 2,
        3, 4,
        5, 6
    });
    
    Matrix m2 = m1.transpose();
    
    Matrix m3 = m1 * m2;
    Matrix m4 = m1 + m1;
    Matrix m5 = m4 * 4;
    Matrix m6 = m5.reShape(6, 1);

    cout << "m2 size : \n" << m2.size() << endl;
    cout << "m1 : \n" << m1 << endl;
    cout << "m2 : \n" << m2 << endl;
    cout << "m3 : \n" << m3 << endl;
    cout << "m4 : \n" << m4 << endl;
    cout << "m5 : \n" << m5 << endl;
    cout << "m6 : \n" << m6 << endl;

    cout << m3.size();

    Shape rectangle(4,1);
    Shape sh(rectangle, 2);
    cout << rectangle.getArea() << endl;
    cout << rectangle.getAngle() << endl;
    cout << rectangle.longestLine() <<endl;
    
    cout << sh.getArea() << endl;
    cout << sh.getAngle() << endl;
    cout << sh.longestLine() <<endl;
    cout << sh.getVolume() << endl;
    return 0;
}
