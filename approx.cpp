#include "amatrix.h"

unsigned nChoosek(unsigned n, unsigned k)
{
    if (k > n)
        return 0;
    if (k * 2 > n)
        k = n - k;
    if (k == 0)
        return 1;

    unsigned result = n;
    for (unsigned i = 2; i <= k; ++i)
    {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

sp::BinomialDistribution::BinomialDistribution(int trial_count, double success_probabiliity) : trial_count(trial_count), success_probabiliity(success_probabiliity) { ; };

double sp::BinomialDistribution::pdf_range(int min_suc, int max_suc)
{
    double result = 0;
    for (int i = min_suc; i <= max_suc; i++)
    {
        result += nChoosek(this->trial_count, i) * pow(this->success_probabiliity, i) * pow((1 - this->success_probabiliity), this->trial_count - i);
    }
    return result;
}

double sp::BinomialDistribution::pmf(int number_of_success)
{
    return nChoosek(this->trial_count, number_of_success) * pow(this->success_probabiliity, number_of_success) * pow((1 - this->success_probabiliity), this->trial_count - number_of_success);
}

double approx::simpInt(double start, double end, double dx, double (*function)(double))
{
    double sum = function(start) + function(end);
    int partition_count = (int)((end - start) / dx);

    if (partition_count % 2 == 1)
        partition_count++;
    dx = (end - start) / partition_count;

    char cons4 = 1;
    double cur;
    for (double j = start + dx; j < end - dx; j += dx)
    {
        if (std::isnan(cur = function(j)) || std::isinf(cur))
            continue;

        if (cons4)
        {
            sum += 4 * cur;
            cons4 = 0;
        }
        else
        {
            sum += 2 * cur;
            cons4 = 1;
        }
    }
    return dx * sum / 3;
}

double approx::simpInt(double start, double end, int partitions, double (*function)(double))
{
    double dx, sum = function(start) + function(end);
    dx = (end - start) / partitions;

    double cur;
    for (int j = 1; j < partitions; j++)
    {
        if (std::isnan(cur = function(start + dx * j)) || std::isinf(cur))
            continue;
        if (j % 2 == 1)
            sum += 4 * cur;
        else
            sum += 2 * cur;
    }
    return dx * sum / 3;
}

// Be careful about precision > 15
//double nthderivative(double x_val, int derivative_order, double (*function)(double), int precision)
//{
/*     double deltaX = 1;
    std::vector<double> results;

    while (std::isnan(function(x_val + deltaX)) || std::isnan(function(x_val - deltaX)))
    {
        deltaX /= 10;
        results.push_back((function(x_val + deltaX) - function(x_val - deltaX)) / deltaX / 2);
    }
    do
    {
        last_result = cur_result;
        deltaX /= 10;
        cur_result = (function(x_val + deltaX) - function(x_val - deltaX)) / deltaX / 2;
    } while (pow(10, precision) * (cur_result - last_result) > 1); */
//}

// Works probabilistic, doesn't give the righ precision every time
double approx::derivative(double x_val, double (*function)(double), int precision)
{
    double deltaX = 1;
    double last_result, cur_result;

    // might be nan due to deltaX and the starting point(x_val)
    // might not terminate add && deltaX > some number
    // similar for the do/while part

    while (std::isnan(function(x_val + deltaX)) || std::isnan(function(x_val - deltaX)))
    {
        deltaX /= 10;
        last_result = cur_result = (function(x_val + deltaX) - function(x_val - deltaX)) / deltaX / 2;
    }
    do
    {
        last_result = cur_result;
        deltaX /= 10;
        cur_result = (function(x_val + deltaX) - function(x_val - deltaX)) / deltaX / 2;
    } while (pow(10, precision) * (cur_result - last_result) > 1);

    return cur_result;
}