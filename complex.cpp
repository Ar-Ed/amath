#include "amatrix.h"

complex::complex(double real, double imag) : real(real), imag(imag)
{
    this->angle = atan(imag / real);
}

complex::complex(double angle) : angle(angle)
{
    std::vector<double> pair{angle_to_complex(angle)};
    this->imag = pair[1];
    this->real = pair[0];
}

complex complex::operator~() const
{
    return complex(this->real, (-1) * this->imag);
}

complex complex::operator+(double constant) const
{
    return complex(this->real + constant, this->imag);
}

complex complex::operator-(double constant) const
{
    return complex(this->real - constant, this->imag);
}

complex complex::operator*(double constant) const
{
    return complex(this->real * constant, this->imag);
}

complex complex::operator/(double constant) const
{
    if (constant == 0)
    {
        //error
    }
    return complex(this->real / constant, this->imag);
}

complex complex::operator+(const complex &number) const
{
    return complex(this->real + number.get_real(), this->imag + number.get_imag());
}

complex complex::operator-(const complex &number) const
{
    return complex(this->real - number.get_real(), this->imag - number.get_imag());
}
complex complex::operator*(const complex &number) const
{
    return complex(this->real * number.get_real() - this->imag * number.get_imag(), this->imag * number.get_real() + this->real * number.get_imag());
}

complex complex::operator/(const complex &number) const
{
    double scale = pow(number.get_imag(), 2) * pow(number.get_real(), 2);
    if (scale == 0)
    {
        //error
    }
    return complex((this->real * number.get_real() + this->imag * number.get_imag()) / scale, (this->imag * number.get_real() - this->real * number.get_imag()) / scale);
}

double complex::get_real() const
{
    return this->real;
}

double complex::get_imag() const
{
    return this->imag;
}

double complex::get_angle() const
{
    return this->angle;
}

std::vector<double> complex::angle_to_complex(double angle)
{
    std::vector<double> pair = {cos(angle), sin(angle)};
    return pair;
}

std::ostream &operator<<(std::ostream &out, const complex &number)
{
    out << number.get_real() << " + " << number.get_imag() << "i";
    return out;
}