#include "poly.h"

vector<vector<double> > terms;
int term_count;

Poly::Poly(vector<vector<double> > array):terms(array){
    for(int i = 0 ; i < terms.size(); i++){
        if(terms[i][0] == 0)terms.erase(terms.begin() + i);
    }
    term_count = terms.size();
    int index = term_count;
    while(index--){
        for(int i = 0; i < index; i++){
            if(terms[i][1] < terms[i+1][1]){
                vector<double> temp = terms[i];
                terms[i]= terms[i+1];
                terms[i+1] = temp;
            }
        }
    }
}
void Poly::setC(int term, double constant){terms[term][0] = constant;}
void Poly::setO(int term, double order){terms[term][1] = order;}

double Poly::getO(int term)const{return terms[term][1];}
double Poly::getC(int term)const{return terms[term][0];}

int Poly::termCount()const{return term_count;}

Point Poly::evaluatePoint(double x_coordinate)const{
    vector<double>coordinate(2,0);
    for(vector<double> term : terms){
        coordinate[1] += term[0] * pow(x_coordinate, term[1]);
    }coordinate[0] = x_coordinate;
    return Point(coordinate);
}

double Poly::definiteIntegral(double l_bound, double u_bound)const{
    double sum = 0;
    for(vector<double> term : terms){
        sum += term[0]*(pow(u_bound,term[1]+1)-pow(l_bound,term[1]+1))/(term[1]+1);
    }
    return sum;
}

double Poly::definiteDerivative(double x_coordinate)const{
    double sum = 0;
    for(vector<double> term: terms){
        sum += term[0] * term[1] * pow(x_coordinate,term[1]-1);    
    }
    return sum;
}

Poly Poly::integral(){
    vector<vector<double> > array;
    for(vector<double> term :terms){
        array.push_back({term[0]/(1 + term[1]),term[1]+1});
    }
    return Poly(array);
}

Poly Poly::derivative(){
    vector<vector<double> > array;
    for(vector<double> term :terms){
        array.push_back({term[0]*term[1],term[1]-1});
    }
    return Poly(array);
}
//polys might not be in the same order, you should sum compatible terms (sorting in earlier steps might be a solution)
Poly Poly::operator+(const Poly& p)const{
    vector<vector<double> > array;
    if(p.term_count > term_count){
        for(int i= 0; i < term_count; i++){
            array.push_back({p.terms[i][0] + terms[i][0],terms[i][1]});
        }for(int i = term_count; i< p.term_count; i++)array.push_back(p.terms[i]);
    }else{
        for(int i= 0; i < p.term_count; i++){
            array.push_back({p.terms[i][0] + terms[i][0],terms[i][1]});
        }for(int i = p.term_count; i< term_count; i++)array.push_back(terms[i]);
    }
    return Poly(array);
}
Poly Poly::operator-(const Poly& p)const{
    vector<vector<double> > array;
    if(p.term_count > term_count){
        for(int i= 0; i < term_count; i++){
            array.push_back({p.terms[i][0] - terms[i][0],terms[i][1]});
        }for(int i = term_count; i< p.term_count; i++)array.push_back(p.terms[i]);
    }else{
        for(int i= 0; i < p.term_count; i++){
            array.push_back({p.terms[i][0] - terms[i][0],terms[i][1]});
        }for(int i = p.term_count; i< term_count; i++)array.push_back(terms[i]);
    }
    return Poly(array);
}

ostream& operator<<(ostream& out, const Poly& p){
    for(vector<double> term : p.terms){
        out << term[0] << "x^(" << term[1] << ") ";
    }
    return out;
}