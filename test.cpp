#include "array.h"
#include <cmath>

int main()
{
    array m1({
            1, 2, 3, 5, 
            5, 6, 3, 2, 
            32, 4, 1, 5, 
            10, 20, 40, 0}, 4, 4);
    m1.print();
    std::cout << "\n";
    m1.Utriangular().print();
    std::cout << m1.det();

    std::cout << "\n\n";
    m1.aggregate(COL, [](auto x, auto y) {return x + y;}).print();

    std::cout << m1 << "\n" << m1(0, 0) << "\n" << (m1 == m1) << "\n" << m1({1,3}, {0,4});

    return 0;
}