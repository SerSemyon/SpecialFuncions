#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "CPUfunctions.h"
#include "GPUfunctions.h"
#include "Test.h"
#include "log_duration.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <complex>

class CylindricalBody
{
    double a_0 = 1;
    //std::complex<double> J_nk_1p, J_nk

public:
    std::complex<double> B_n(int n, double k_1, double k_2,
        std::complex<double> J_nk_1, std::complex<double> J_nk_2, std::complex<double> dJ_nk_1, std::complex<double> dJ_nk_2,
        std::complex<double> H_n, std::complex<double> dH_n);

    std::complex<double> C_n(int n, double k_1, double k_2,
        std::complex<double> J_nk_1, std::complex<double> J_nk_2, std::complex<double> dJ_nk_1, std::complex<double> dJ_nk_2,
        std::complex<double> H_n, std::complex<double> dH_n);

    std::complex<double> find_B_0(double k_1, double k_2, double R);

    std::complex<double> find_C_0(double k_1, double k_2, double R);

    std::complex<double> find_B_n(int N, double k_1, double k_2, double R);

    std::complex<double> find_C_n(int N, double k_1, double k_2, double R);

    std::complex<double>** E_1(double k_1, double k_2, double* r, size_t count, size_t number_division, double R);

    std::complex<double>** E_2(double k_1, double k_2, double* r, size_t count, size_t number_division, double R);
};

