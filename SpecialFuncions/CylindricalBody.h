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
    double _k_1;
    double _k_2;
    double _R;
    double _a_0;
    // Поля ниже использовались для избавления от повторного вычисления. В этом варианте программы упрощённый вариант с вызовом функций вычисления напрямую, чтобы легче было искать ошибки.
    /*std::complex<double> J_nk_1p, J_nk_1pp, J_nk_1pm,
        J_nk_1n, J_nk_1np, J_nk_1nm;
    std::complex<double> J_nk_2p, J_nk_2pp, J_nk_2pm,
        J_nk_2n, J_nk_2np, J_nk_2nm;
    std::complex<double> H_np, H_npp, H_npm,
        H_nn, H_nnp, H_nnm;
    std::complex<double> dJ_nk_1p, dJ_nk_1n;
    std::complex<double> dJ_nk_2p, dJ_nk_2n;
    std::complex<double> dH_np, dH_nn;*/

public:
    CylindricalBody(double k_1, double k_2, double radiusCircle, double a_0);

    std::complex<double> B_n(int n,
        std::complex<double> J_nk_1, std::complex<double> J_nk_2, std::complex<double> dJ_nk_1, std::complex<double> dJ_nk_2,
        std::complex<double> H_n, std::complex<double> dH_n);

    std::complex<double> C_n(int n,
        std::complex<double> J_nk_1, std::complex<double> J_nk_2, std::complex<double> dJ_nk_1, std::complex<double> dJ_nk_2,
        std::complex<double> H_n, std::complex<double> dH_n);

    std::complex<double> find_B_0();

    std::complex<double> find_C_0();

    std::complex<double> find_B_n(int N);

    std::complex<double> find_C_n(int N);

    std::complex<double>** E_1(double* r, size_t count, size_t number_division);

    std::complex<double>** E_2(double* r, size_t count, size_t number_division);
};

