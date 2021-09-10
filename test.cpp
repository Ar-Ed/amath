#include "amatrix.h"
#include <cmath>

int main()
{
    array m1({1, 2, 3, 5,
              5, 6, 3, 2,
              32, 4, 1, 5,
              10, 20, 40, 0},
             4, 4);
    m1.print();
    std::cout << "\n";
    m1.Utriangular().print();
    std::cout << m1.det();

    std::cout << "\n\n";
    m1.aggregate(COL, [](auto x, auto y)
                 { return x + y; })
        .print();

    std::cout << m1 << "\n"
              << m1(0, 0) << "\n"
              << (m1 == m1) << "\n"
              << m1({1, 3}, {0, 4}) << "\n";
    /*     m1.write_file("mycsv.txt");
    m1.write_file("mycsv.csv"); */
    /*     array m2("mycsv.csv");
    std::cout << m2; */
    std::cout << array({123, 2, 3, 4, 5, 6, 7, 11, 17}, 3, 3) << "\n"
              << array({123, 2, 3, 4, 5, 6, 7, 11, 17}, 3, 3).Ltriangular();

    array b({4, 7}, 2, 1);
    array A({3, 2, 5, 3}, 2, 2);

    std::cout << "\n"
              << b.solveSquare(A);

    array b1({1, 1,
              2, 3,
              3, 4},
             3, 2);
    std::cout << "\n\n"
              << leastSquares(b1, 2);

    return 0;
}