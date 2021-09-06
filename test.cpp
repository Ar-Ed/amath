#include "lib.h"
#include "cmath"

double myfucntion(double number){
    return cos(number) * log10(number) / 100;
}

int main()
{
    array m1({1, 2, 3, 5, 5, 6, 3, 2, 32, 4}, 2, 5);

    array m2 = m1 + m1;
    m2.print(" ");

    std::cout << "\n";

    m1 = m2.copy();
    m1.print(" ");

    std::cout << "\n";

    (m1 * 3 - m1 / 3).print(" ");

    std::cout << "\n";

    array m3 = m1 * m1.transpose();
    m3.print(" ");

    std::cout << "\n";
    (m3.apply(myfucntion)).print();

    //array m2({1, 2, 3, 5, 5, 6, 3, 2, 32, 4}, 2, 5);

    //array m3 = m1 + m2;

    //std::cout << m1.get_cols() << " " << m1.get_rows() << " " << m1.get_size();
    return 0;
}