#include "HEADER_ONLY_AMATH/amath.h"

class Timer
{
    std::chrono::high_resolution_clock::time_point start, end;
    int timer;
    static int timer_count;

public:
    Timer()
    {
        start = std::chrono::high_resolution_clock::now();
        timer = timer_count;
        timer_count++;
    }

    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "\nTimer Id: " << timer << " Finished in: " << duration.count() << "\n";
    }
};

int Timer::timer_count = 0;

int main()
{
    Matrix m1({1, 2, 3, 5,
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

    std::cout << Matrix({123, 2, 3, 4, 5, 6, 7, 11, 17}, 3, 3) << "\n"
              << Matrix({123, 2, 3, 4, 5, 6, 7, 11, 17}, 3, 3).Ltriangular();

    Matrix b({4, 7}, 2, 1);
    Matrix A({3, 2, 5, 3}, 2, 2);

    std::cout << "fdsa\n"
              << b.solveSquare(A);

    Matrix b1({1, 1,
               2, 3,
               3, 4},
              3, 2);

    std::cout << "\n\n"
              << leastSquares(b1, 1);

    m1.write_file("example_data/m1.csv");
    Matrix m2("example_data/m1.csv");
    m1 = read_file("example_data/m1.csv");

    std::cout << "\n\n"
              << m1 << "\n"
              << m2 << m2.get_rows() << " " << m2.get_cols() << " " << m2.get_size();

    std::cout << "\n"
              << m2.transpose() << " \n"
              << m2.aggregate(BOTH, [](double y, double x)
                              { return x + y; });

    std::cout << Matrix(1, 3, 3);

    std::cout << "\n"
              << Matrix(2, 2, 1, 2, 3, 4);

    std::cout << "\n"
              << random(0, 1, 3, 4) * 10 << "\n";

    std::cout << "\n"
              << random(0, 1, 3, 3).pairWise([](double x, double y)
                                             { return x + y; },
                                             random(0, 1, 3, 3));

    std::cout << "\n"
              << Matrix(2, 2, 1, 2, 3, 4).inverse();

    std::cout << "\n"
              << Matrix(2, 2, 1, 2, 3, 4).trace();

    std::vector<Matrix> lu = Matrix(4, 4, 1, 2, 3, 4, 7, 8, 2, 1, 9, 5, 3, 5, 11, 13, 20, 17).LUFactor();
    std::cout << "\n"
              << lu[0] << "\n"
              << lu[1] << "\n"
              << lu[0] * lu[1] << "\n"
              << Matrix(4, 4, 1, 2, 3, 4, 7, 8, 2, 1, 9, 5, 3, 5, 11, 13, 20, 17).Ltriangular() << "\n"
              << Matrix(4, 4, 1, 2, 3, 4, 7, 8, 2, 1, 9, 5, 3, 5, 11, 13, 20, 17).Utriangular();

    std::cout << "\n"
              << Matrix(3, 3, 1, 2, 5, 6, 7, 1, 2, 3, 4).minorMatrix() << "\n"
              << Matrix(3, 3, 1, 2, 5, 6, 7, 1, 2, 3, 4).cofactorMatrix();

    Matrix(3, 3, 1, 2, 5, 6, 7, 1, 2, 3, 4).cofactorMatrix().print(",", 6);

    //plot({0, 7}, {0, 10}, Matrix(4, 2, 3, 4, 4, 2, 5, 4, 6, 7), "myfile.png");

    auto start = std::chrono::high_resolution_clock::now();
    /* for (int i = 0; i < 10; i++)
        Matrix(2., 5, 100) * Matrix(1., 100, 100000); */
    {
        Timer timer;
        Matrix(2., 2, 100) * Matrix(1., 100, 100000);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\n"
              << (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count() / 10000;

    // something wrong with this one ( arange is not reliable)
    //std::cout << logspace(1, 1e9, 10) << "\n"<< linspace(1, 100, 10) << "\n"<< arange(1, 100, 10).get_size();

    /*for (int i = 1; i < 10; i++)
    {
        Matrix coeff = leastSquares(b1, i);
        std::cout << b1 << "\n"
                  << coeff;

        plotModel({0, 4}, {0, 5}, coeff, linspace(0, 4, 100), b1, "save.pdf");
    }

    Matrix coeff = leastSquares(b1, 1);
    plotModel({0, 4}, {0, 5}, coeff, linspace(0, 4, 100), b1, "save.pdf"); */

    Matrix mydata = read_file("example_data/linear.csv");
    Matrix model = leastSquares(mydata, 1);

    //plotModel(model, mydata, linspace(-2, 12, 10), {-2, 12}, {-10, 20});
    std::cout << mydata.max() << "\n"
              << mydata.min() << "\n"
              << mydata.argmax() << "\n"
              << mydata.argmin() << "\n\n"
              << leastSquares(mydata, 1);

    Matrix matrixx(10, 10000, 100000);

    {
        Timer timer1;
        Matrix mm{matrixx > matrixx};
    }


    return 0;
}