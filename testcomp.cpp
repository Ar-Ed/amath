#include "amatrix.h"

int main()
{

    complex n1(1, -2);
    complex n2(9, 5);
    std::cout << n1 * n2 << " \n"
              << n1 / n2 << " \n"
              << ~((n1 + n2) * 4 / 2 + 5) << "\n"
              << complex(60 * PI180);

    return 0;
}