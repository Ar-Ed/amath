#include "amatrix.h"

double func(double x)
{
    return log10(x) * sin(x);
}

int main()
{
    std::cout.precision(15);
    std::cout << approx::derivative(2, func, 10) << "\n"
              << approx::simpInt(0, 1, 0.1, sin) << "\n"
              << approx::simpInt(1e-10, 1, 1000, log10) << "\n"; // Weighted approach is needed in such cases


    sp::BinomialDistribution dist(10, 0.4);
    
    std::cout.precision(15);
    std::cout << dist.pmf(4) << "\n" << dist.pdf_range(0, 4);
    return 0;
}