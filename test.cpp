#include "amatrix.h"

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

    std::cout << array({123, 2, 3, 4, 5, 6, 7, 11, 17}, 3, 3) << "\n"
              << array({123, 2, 3, 4, 5, 6, 7, 11, 17}, 3, 3).Ltriangular();

    array b({4, 7}, 2, 1);
    array A({3, 2, 5, 3}, 2, 2);

    std::cout << "fdsa\n"
              << b.solveSquare(A);

    array b1({1, 1,
              2, 3,
              3, 4},
             3, 2);

    std::cout << "\n\n"
              << leastSquares(b1, 1);

    m1.write_file("m1.csv");
    array m2("m1.csv");
    m1 = read_file("m1.csv");

    std::cout << "\n\n"
              << m1 << "\n"
              << m2 << m2.get_rows() << " " << m2.get_cols() << " " << m2.get_size();

    std::cout << "\n"
              << m2.transpose() << " \n"
              << m2.aggregate(BOTH, [](double y, double x)
                              { return x + y; });

    std::cout << array(1, 3, 3);

    std::cout << "\n"
              << array(2, 2, 1, 2, 3, 4);

    std::cout << "\n"
              << random(0, 1, 3, 4) * 10 << "\n";

    std::cout << "\n"
              << random(0, 1, 3, 3).pairWise([](double x, double y)
                                             { return x + y; },
                                             random(0, 1, 3, 3));

    std::cout << "\n"
              << array(2, 2, 1, 2, 3, 4).inverse();

    std::cout << "\n"
              << array(2, 2, 1, 2, 3, 4).trace();

    std::vector<array> lu = array(4, 4, 1, 2, 3, 4, 7, 8, 2, 1, 9, 5, 3, 5, 11, 13, 20, 17).LUFactor();
    std::cout << "\n"
              << lu[0] << "\n"
              << lu[1] << "\n"
              << lu[0] * lu[1] << "\n"
              << array(4, 4, 1, 2, 3, 4, 7, 8, 2, 1, 9, 5, 3, 5, 11, 13, 20, 17).Ltriangular() << "\n"
              << array(4, 4, 1, 2, 3, 4, 7, 8, 2, 1, 9, 5, 3, 5, 11, 13, 20, 17).Utriangular();

    std::cout << "\n"
              << array(3, 3, 1, 2, 5, 6, 7, 1, 2, 3, 4).minorMatrix() << "\n"
              << array(3, 3, 1, 2, 5, 6, 7, 1, 2, 3, 4).cofactorMatrix();

    array(3, 3, 1, 2, 5, 6, 7, 1, 2, 3, 4).cofactorMatrix().print(",", 6);
    return 0;
}