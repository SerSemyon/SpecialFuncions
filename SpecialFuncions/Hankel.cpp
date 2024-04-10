#include "CPUfunctions.h"

#define _USE_MATH_DEFINES // для C++
#include <cmath>
#include <math.h>

void H1(const int v, double* x, double* Re, double* Im, const unsigned int size)
{
	J(v, x, Re, size);
	Neumann(v, x, Im, size, Re);
}

std::complex<double> H1(const int v, double x)
{
    if (v != 0) {
        double Re = J(v, x);
        double Im = Neumann(v, x, Re);
        return std::complex<double>(Re, Im);
    }
    else {
        double Re = J_0(x);
        double Im = Y_0(x, Re);
        return std::complex<double>(Re, Im);
    }
}

std::complex<double> H_negative(int nu, std::complex<double> H)
{
    std::complex<double> e = std::complex<double>(cos(nu * M_PI), sin(nu * M_PI));
    return e * H;
}

void H1(const int v, double* x, std::complex<double>* res, const unsigned int size)
{
    for (int i = 0; i < size; i++) {
        double Re = J(v, x[i]);
        double Im = Neumann(v, x[i], Re);
        res[i] = std::complex<double>(Re, Im);
    }
}