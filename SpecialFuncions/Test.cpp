#pragma once
#include <ctime>
#include "Test.h"
#include "GPUfunctions.h"
#include "CPUfunctions.h"
#include <iostream>
#include "log_duration.h"

#include <iostream> 
#include <iomanip> 
#include <chrono> 

double epsilon = 1E-7;
int nJ0 = 1000000;
int nJ1 = 1000000;
int nJ = 1000000;
int nY0 = 10;
int nY1 = 1000;
double b0 = 8;
double b1 = 8;
double bJ = 8;
double h0 = b0 / nJ0;
double h1 = b1 / nJ1;
double hJ = bJ / nJ;
double hY0 = b0/ nY0;
double hY1 = b1 / nY1;

size_t count_experiments = 100;
double time_first_function;
double time_second_function;

/* TODO - Встроенная реализация даёт низкую точность и, чем дальше от нуля, тем выше ошибка.
Поэтому для проверки значений нужно будет использовать таблицы с более точными результатами,иначе тест всегда будет проваливаться. 
Текущая реализация теста говорит, что тест провален, даже если значения верны. */
void TestBesselCPU()
{
    std::cout << "TestBesselCPU started" << std::endl;
    int v = 1;
    int n = 1000;
    bool successfully = true;
    double* res1 = new double[n];
    double* res2 = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = i * 0.01;
        res1[i] = std::cyl_bessel_i(v, x[i]);//__std_smf_cyl_bessel_i(v, x[i]);
    }
    J(v, x, res2, n);
    for (int i = 0; i < n; i++)
    {
        if (abs(res1[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestBesselCPU failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestBesselCPU OK" << std::endl << std::endl;
}

void TestJ0()
{
    std::cout << "TestJ0 started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* res1 = new double[nJ0];
    double* res2 = new double[nJ0];
    double* x = new double[nJ0];
    for (int i = 0; i < nJ0; i++)
    {
        x[i] = i * h0;
    }
    time_first_function = 0;
    for (size_t i = 0; i < count_experiments; i++) 
    {
        auto begin = std::chrono::steady_clock::now();

        J(v, x, res1, nJ0);

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_first_function += elapsed_ms.count();
    }
    std::cout << "J " << time_first_function / count_experiments << std::endl;
    time_second_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        for (int i = 0; i < nJ0; i++)
        {
            res2[i] = J_0(x[i]);
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_second_function += elapsed_ms.count();
    }
    std::cout << "J0 " << time_second_function / count_experiments << std::endl;
    for (int i = 0; i < nJ0; i++)
    {
        if (abs(res1[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ0 failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ0 OK" << std::endl << std::endl;
}

void TestJ1()
{
    std::cout << "TestJ1 started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* res1 = new double[nJ1];
    double* res2 = new double[nJ1];
    double* x = new double[nJ1];
    for (int i = 0; i < nJ1; i++)
    {
        x[i] = i * h1;
    }
    {
        LOG_DURATION("J");
        J(v, x, res1, nJ1);
    }
    {
        LOG_DURATION("J1");
        for (int i = 0; i < nJ1; i++)
        {
            res2[i] = J_1(x[i]);
        }
    }
    for (int i = 0; i < nJ1; i++)
    {
        if (abs(res1[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ1 failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ1 OK" << std::endl << std::endl;
}

void TestNeumann_one_point()
{
    std::cout << "TestNeumann_one_point started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* res1 = new double[nY0];
    double* res2 = new double[nY0];
    double* x = new double[nY0];
    for (int i = 0; i < nY0; i++)
    {
        x[i] = i * hY0;
        res1[i] = __std_smf_cyl_neumann(v, x[i]);
    }
    double* J_pos = new double[nY0];
    double* J_neg = new double[nY0];
    J(v, x, J_pos, nY0);
    {
        LOG_DURATION("Neumann_one_point");
        for (int i = 0; i < nY0; i++) {
            res2[i] = Neumann(v, x[i], J_pos[i]);
        }
    }
    for (int i = 0; i < nY0; i++)
    {
        std::cout << x[i] << " " << std::fixed << std::setprecision(10) << res1[i] << " " << res2[i] << std::endl;
        if (abs(res1[i] - res2[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestNeumann_one_point failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestNeumann_one_point OK" << std::endl << std::endl;
}

void TestNeumannCPU()
{
    std::cout << "TestNeumannCPU started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* res1 = new double[nY0];
    double* res2 = new double[nY0];
    double* x = new double[nY0];
    for (int i = 0; i < nY0; i++)
    {
        x[i] = i * hY0;
        res1[i] = __std_smf_cyl_neumann(v, x[i]);
    }
    double* J_pos = new double[nY0];
    double* J_neg = new double[nY0];
    J(v, x, J_pos, nY0);
    J_negative(v, J_neg, nY0, J_pos);
    {
        LOG_DURATION("Neumann");
        Neumann(v, x, res2,nY0, J_pos, J_neg);
    }
    for (int i = 0; i < nY0; i++)
    {
        std::cout << x[i] << " " << std::fixed << std::setprecision(10) << res1[i] << " " << res2[i] << std::endl;
        if (abs(res1[i] - res2[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestNeumannCPU failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestNeumannCPU OK" << std::endl << std::endl;
}

void TestNeumann_CUDA()
{
    std::cout << "TestNeumann_CUDA started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* res1 = new double[nY0];
    double* res2 = new double[nY0];
    double* x = new double[nY0];
    for (int i = 0; i < nY0; i++)
    {
        x[i] = (i + 1) * hY0;
    }
    double* J_pos = new double[nY0];
    double* J_neg = new double[nY0];
    J(v, x, J_pos, nY0);
    J_negative(v, J_neg, nY0, J_pos);
    {
        LOG_DURATION("CPU");
        Neumann(v, x, res1, nY0, J_pos, J_neg);
    }
    {
        LOG_DURATION("GPU");
        Neumann_CUDA(v, x, res2, nY0, J_pos, J_neg);
    }
    for (int i = 0; i < nY0; i++)
    {
        //std::cout << x[i] << " " << res1[i] << " " << res2[i] << std::endl;
        if (abs(res1[i] - res2[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestNeumann_CUDA failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestNeumann_CUDA OK" << std::endl << std::endl;
}

void TestY0()
{
    std::cout << "TestY0 started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* res1 = new double[nY0];
    double* res2 = new double[nY0];
    double* res3 = new double[nY0];
    double* x = new double[nY0];
    for (int i = 0; i < nY0; i++)
    {
        x[i] = i * hY0;
    }
    double* Js = new double[nY0];
    J(v, x, Js, nY0);
    {
        LOG_DURATION("Y_0");
        for (int i = 0; i < nY0; i++)
        {
            res1[i] = Y_0(x[i], Js[i]);
        }
    }
    {
        LOG_DURATION("Neumann");
        Neumann(v, x, res2, nY0, Js);
    }
    {
        LOG_DURATION("Y_0");
        Y_0(x, res3, nY0, Js);
    }
    for (int i = 0; i < nY0; i++)
    {
        std::cout << x[i] << " " << res1[i] << " " << res2[i] << " " << res3[i] << std::endl;
        if (abs(res1[i] - res2[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestY0 failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
        if (abs(res1[i] - res3[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestY0 failed!" << x[i] << " " << res1[i] << " " << res3[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    delete[] res3;
    if (successfully)
        std::cout << "TestY0 OK" << std::endl << std::endl;
}

void TestY1()
{
    std::cout << "TestY1 started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* res1 = new double[nY1];
    double* res2 = new double[nY1];
    double* res3 = new double[nY1];
    double* x = new double[nY1];
    for (int i = 0; i < nY1; i++)
    {
        x[i] = i * hY1;
    }
    double* Js = new double[nY1];
    J(v, x, Js, nY1);
    {
        LOG_DURATION("Y_1");
        for (int i = 0; i < nY1; i++)
        {
            res1[i] = Y_1(x[i], Js[i]);
        }
    }
    {
        LOG_DURATION("Neumann");
        Neumann(v, x, res2, nY1, Js);
    }
    {
        LOG_DURATION("Y_1");
        Y_1(x, res3, nY1, Js);
    }
    for (int i = 0; i < nY1; i++)
    {
        if (abs(res1[i] - res2[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestY1 failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
        if (abs(res1[i] - res3[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestY1 failed!" << x[i] << " " << res1[i] << " " << res3[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    delete[] res3;
    if (successfully)
        std::cout << "TestY1 OK" << std::endl << std::endl;
}

void TestBessel_CUDA()
{
    std::cout << "TestBesselCuda started" << std::endl;
    int v = 0;
    bool successfully = true;
    double h = bJ / nJ;
    double* x = new double[nJ];
    double* resGPU = new double[nJ];
    double* resCPU = new double[nJ];
    for (int i = 0; i < nJ; i++)
    {
        x[i] = h * i;
    }

    {
        LOG_DURATION("CPU clock");
        J(v, x, resCPU, nJ);
    }

    {
        LOG_DURATION("GPU clock");
        BesselWithCuda(v, x, resGPU, nJ);
    }

    for (int i = 0; i < nJ; i++)
    {
        if (abs(resGPU[i] - resCPU[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestBesselCuda failed! point:" << x[i] << " |resGPU-resCPU|=" << abs(resGPU[i] - resCPU[i])  << std::endl << std::endl;
            break;
        }
    }
    delete[] x;
    delete[] resGPU;
    delete[] resCPU;
    if (successfully)
        std::cout << "TestBesselCuda OK" << std::endl << std::endl;
}

void TestJ0_CUDA()
{
    std::cout << "TestJ0_CUDA started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* x = new double[nJ0];
    double* res1 = new double[nJ0];
    double* res2 = new double[nJ0];
    for (int i = 0; i < nJ0; i++)
    {
        x[i] = i * h0;
    }

    {
        LOG_DURATION("Bessel");
        BesselWithCuda(v, x, res1, nJ0);
    }

    {
        LOG_DURATION("J0");
        J0_CUDA(x, res2, nJ0);
    }

    for (int i = 0; i < nJ0; i++)
    {
        if (abs(res1[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ0_CUDA failed! point:" << x[i] << " |res1-resCPU|=" << abs(res1[i] - res2[i]) << std::endl << std::endl;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ0_CUDA OK" << std::endl << std::endl;
}

void TestJ1_CUDA()
{
    std::cout << "TestJ1_CUDA started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* x = new double[nJ1];
    double* res1 = new double[nJ1];
    double* res2 = new double[nJ1];
    for (int i = 0; i < nJ1; i++)
    {
        x[i] = i * h1;
    }

    {
        LOG_DURATION("Bessel");
        BesselWithCuda(v, x, res1, nJ1);
    }

    {
        LOG_DURATION("J1");
        J1_CUDA(x, res2, nJ1);
    }

    for (int i = 0; i < nJ1; i++)
    {
        if (abs(res1[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ1_CUDA failed! point:" << x[i] << " |res1-resCPU|=" << abs(res1[i] - res2[i]) << std::endl << std::endl;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ1_CUDA OK" << std::endl << std::endl;
}

void TestY0_CUDA()
{
    std::cout << "TestY0_CUDA started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* res1 = new double[nY0];
    double* res2 = new double[nY0];
    double* x = new double[nY0];
    for (int i = 0; i < nY0; i++)
    {
        x[i] = (i+1) * hY0;
    }
    double* Js = new double[nY0];
    J(v, x, Js, nY0);
    {
        LOG_DURATION("CPU");
        Y_0(x, res1, nY0, Js);
        //Neumann(v, x, res1, n, Js);
    }
    {
        LOG_DURATION("GPU");
        Y0_CUDA(x, res2, nY0, Js);
    }
    for (int i = 0; i < nY0; i++)
    {
        if (abs(res1[i] - res2[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestY0_CUDA failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestY0_CUDA OK" << std::endl << std::endl;
}

void TestY1_CUDA()
{
    std::cout << "TestY1_CUDA started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* res1 = new double[nY1];
    double* res2 = new double[nY1];
    double* x = new double[nY1];
    for (int i = 0; i < nY1; i++)
    {
        x[i] = (i + 1) * hY1;
    }
    double* Js = new double[nY1];
    J(v, x, Js, nY1);
    {
        LOG_DURATION("CPU");
        Y_1(x, res1, nY1, Js);
        //Neumann(v, x, res1, n, Js);
    }
    {
        LOG_DURATION("GPU");
        Y1_CUDA(x, res2, nY1, Js);
    }
    for (int i = 0; i < nY1; i++)
    {
        if (abs(res1[i] - res2[i]) > 1E-4)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestY1_CUDA failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestY1_CUDA OK" << std::endl << std::endl;
}

double T_recursively(int n, double x)
{
    if (n == 0)
        return 1;
    if (n == 1)
        return x;
    return 2 * x * T_recursively(n - 1, x) - T_recursively(n - 2, x);
}

void TestChebyshevPolynomials()
{
    LOG_DURATION("TestChebyshevPolynomials");
    std::cout << "TestChebyshevPolynomials started" << std::endl;
    double t1, t2;
    bool successfully = true;
    for (int i = 0; i < 11; i++)
    {
        t1 = T_recursively(i, 0.2);
        t2 = T(i, 0.2);
        if ((t1 - t2) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestChebyshevPolynomials failed! order " << i << "T_recursively - T = " << t1 - t2 << std::endl << std::endl;
            successfully = false;
        }
    }
    if (successfully)
        std::cout << "TestChebyshevPolynomials OK" << std::endl << std::endl;
}

void TestJnew()
{
    std::cout << "TestJnew started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* x = new double[nJ];
    double* res1 = new double[nJ];
    double* res2 = new double[nJ];
    for (int i = 0; i < nJ; i++)
    {
        x[i] = hJ * i;
    }

    {
        LOG_DURATION("J");
        J(v, x, res1, nJ);
    }

    {
        LOG_DURATION("Jnew");
        Jnew(v, x, res2, nJ);
    }

    for (int i = 0; i < nJ; i++)
    {
        if (abs(res1[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJnew failed! point:" << x[i] << " |res1-res2|=" << abs(res1[i] - res2[i]) << std::endl << std::endl;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJnew OK" << std::endl << std::endl;
}

void TestBesselOrderedSet()
{
    std::cout << "TestBesselOrderedSet started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* x = new double[nJ];
    double* res1 = new double[nJ];
    double* res2 = new double[nJ];
    for (int i = 0; i < nJ; i++)
    {
        x[i] = hJ * i;
    }

    {
        LOG_DURATION("J");
        J(v, x, res1, nJ);
    }

    {
        LOG_DURATION("BesselOrderedSet");
        BesselOrderedSet(v, x, res2, nJ);
    }

    for (int i = 0; i < nJ; i++)
    {
        if (abs(res2[i] - res1[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestBesselOrderedSet failed! point:" << x[i] << " |res1-res2|=" << abs(res2[i] - res1[i]) << std::endl << std::endl;
            break;
        }
    }
    delete[] x;
    delete[] res2;
    delete[] res1;
    if (successfully)
        std::cout << "TestBesselOrderedSet OK" << std::endl << std::endl;
}

void TestZ_vNext()
{
    std::cout << "TestZ_vNext started" << std::endl;
    int v = 0;
    int n = 49;
    bool successfully = true;
    double* res0 = new double[n];
    double* res1 = new double[n];
    double* res2 = new double[n];
    double* resZ = new double[n];
    double* resN = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = (i+1) * 0.1;
    }
    double* Js = new double[n];
    {
        LOG_DURATION("J_0");
        J(v, x, res0, n);
    }
    v = 1;
    {
        LOG_DURATION("J_1");
        J(v, x, res1, n);
    }
    v = 2;
    {
        LOG_DURATION("Neumann");
        J(v, x, res2, n);
    }
    {
        LOG_DURATION("Z_2");
        for (int i = 0; i < n; i++)
        {
            std::cout << x[i] << " " << res0[i] << " " << res1[i] << std::endl;
            resZ[i] = cyl_next_order(v-1, x[i], res1[i], res0[i]);
        }
    }
    {
        LOG_DURATION("Neumann");
        cyl_next_order(v-1, x, resN, n, res1, res0);
    }
    for (int i = 0; i < n; i++)
    {
        if (abs(resZ[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestZ_vNext failed! " << i << " " << x[i] << " " << resZ[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
        if (abs(resZ[i] - resN[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestZ_vNext failed! " << i << " " << x[i] << " " << resZ[i] << " " << resN[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res0;
    delete[] res1;
    delete[] res2;
    delete[] resZ;
    delete[] resN;
    if (successfully)
        std::cout << "TestZ_vNext OK" << std::endl << std::endl;
}

void TestCyl_next_order_CUDA()
{
    std::cout << "TestCyl_next_order_CUDA started" << std::endl;
    int v = 0;
    int n = 49;
    bool successfully = true;
    double* res0 = new double[n];
    double* res1 = new double[n];
    double* resCPU = new double[n];
    double* resGPU = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = (i + 1) * 0.1;
    }
    double* Js = new double[n];
    {
        LOG_DURATION("J_0");
        J(v, x, res0, n);
    }
    v = 1;
    {
        LOG_DURATION("J_1");
        J(v, x, res1, n);
    }
    v = 2;
    {
        LOG_DURATION("CPU");
        cyl_next_order(v - 1, x, resCPU, n, res1, res0);
    }
    {
        LOG_DURATION("GPU");
        cyl_next_order_CUDA(v - 1, x, resGPU, n, res1, res0);
    }
    for (int i = 0; i < n; i++)
    {
        if (abs(resCPU[i] - resGPU[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestCyl_next_order_CUDA failed! " << i << " " << x[i] << " " << resCPU[i] << " " << resGPU[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res0;
    delete[] res1;
    delete[] resCPU;
    delete[] resGPU;
    if (successfully)
        std::cout << "TestCyl_next_order_CUDA OK" << std::endl << std::endl;
}

void TestJ_asymptotic()
{
    std::cout << "TestJ_asymptotic started" << std::endl;
    int v = 10;
    bool successfully = true;
    int n = 100;
    double* res1 = new double[n];
    double* res2 = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = i * 0.1;
    }
    {
        LOG_DURATION("J");
        J(v, x, res1, n);
    }
    {
        LOG_DURATION("J_asymptotic");
        for (int i = 0; i < n; i++)
        {
            res2[i] = J_asymptotic(v,x[i]);
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (abs(res1[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ_asymptotic failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ_asymptotic OK" << std::endl << std::endl;
}

void TestY_asymptotic()
{
    std::cout << "TestY_asymptotic started" << std::endl;
    int v = 1;
    bool successfully = true;
    int n = 1;
    double* res0 = new double[n];
    double* res1 = new double[n];
    double* res2 = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = (i+1) * 10;
    }
    {
        LOG_DURATION("J");
        J(v, x, res0, n);
    }
    {
        LOG_DURATION("Y");
        Neumann(v, x, res1, n, res0);
    }
    {
        LOG_DURATION("Y_asymptotic");
        for (int i = 0; i < n; i++)
        {
            res2[i] = Y_asymptotic(v, x[i]);
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (abs(res1[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestY_asymptotic failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res0;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestY_asymptotic OK" << std::endl << std::endl;
}

void Test_dZ()
{
    std::cout << "Test_dZ started" << std::endl;
    int v = 0;
    bool successfully = true;
    int n = 15000000;
    double* res0 = new double[n];
    double* res1 = new double[n];
    double* res2 = new double[n];
    double* resGPU = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = (i + 1) * 0.00000001;
    }
    {
        LOG_DURATION("J_0");
        J(v, x, res0, n);
    }
    v = 2;
    {
        LOG_DURATION("J_2");
        J(v, x, res1, n);
    }
    v = 1;
    {
        LOG_DURATION("dJ_1");
        dZ(v, x, res2, n, res0, res1);
    }
    {
        LOG_DURATION("dJ_1_CUDA");
        dZ_CUDA(x, resGPU, n, res0, res1);
    }
    for (int i = 0; i < n; i++)
    {
        if (abs(resGPU[i] - res2[i]) > epsilon)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "Test_dZ failed!" << x[i] << " " << resGPU[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res0;
    delete[] res1;
    delete[] res2;
    delete[] resGPU;
    if (successfully)
        std::cout << "Test_dZ OK" << std::endl << std::endl;
}

void Testrec_CUDA()
{
    std::cout << "Testrec_CUDA started" << std::endl;
    int v = 2;
    bool successfully = true;
    double* res1 = new double[nJ1];
    double* res2 = new double[nJ1];
    double* x = new double[nJ1];
    for (int i = 0; i < nJ1; i++)
    {
        x[i] = i * h1;
    }
    {
        LOG_DURATION("J");
        BesselWithCuda(v, x, res1, nJ1, 2);
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ1 OK" << std::endl << std::endl;
}

void TestJ_0_T()
{
    std::cout << "TestJ_0_T started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* res1 = new double[nJ0];
    double* res2 = new double[nJ0];
    double* x = new double[nJ0];
    for (int i = 0; i < nJ0; i++)
    {
        x[i] = i * h0;
    }

    time_first_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        for (int i = 0; i < nJ0; i++)
        {
            res1[i] = J_0(x[i]);
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_first_function += elapsed_ms.count();
    }
    std::cout << "J_0 " << time_first_function / count_experiments << std::endl;
    time_second_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        for (int i = 0; i < nJ0; i++)
        {
            res2[i] = J_0_T(x[i]);
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_second_function += elapsed_ms.count();
    }
    std::cout << "J_0_T " << time_second_function / count_experiments << std::endl;

    for (int i = 0; i < nJ0; i++)
    {
        if (abs(res1[i] - res2[i]) > 1E-5)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ_0_T failed!" << x[i] << " " << res1[i] << " " << res2[i] << " " << std::cyl_bessel_i(0, x[i]) << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ_0_T OK" << std::endl << std::endl;
}

void TestJ_1_T()
{
    std::cout << "TestJ_1_T started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* res1 = new double[nJ1];
    double* res2 = new double[nJ1];
    double* x = new double[nJ1];
    for (int i = 0; i < nJ1; i++)
    {
        x[i] = i * h1;
    }

    time_first_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        for (int i = 0; i < nJ1; i++)
        {
            res1[i] = J_1(x[i]);
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_first_function += elapsed_ms.count();
    }
    std::cout << "J_1 " << time_first_function / count_experiments << std::endl;
    time_second_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        for (int i = 0; i < nJ1; i++)
        {
            res2[i] = J_1_T(x[i]);
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_second_function += elapsed_ms.count();
    }
    std::cout << "J_1_T " << time_second_function / count_experiments << std::endl;

    for (int i = 0; i < nJ1; i++)
    {
        if (abs(res1[i] - res2[i]) > 1E-6)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ_1_T failed!" << x[i] << " " << res1[i] << " " << res2[i] << std::endl << std::endl;
            successfully = false;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ_1_T OK" << std::endl << std::endl;
}

void TestJ_0_T_CUDA()
{
    std::cout << "TestJ_0_T_CUDA started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* x = new double[nJ0];
    double* res1 = new double[nJ0];
    double* res2 = new double[nJ0];
    for (int i = 0; i < nJ0; i++)
    {
        x[i] = i * h0;
    }

    time_first_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        BesselWithCuda(v, x, res1, nJ0);

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_first_function += elapsed_ms.count();
    }
    std::cout << "BesselWithCuda " << time_first_function / count_experiments << std::endl;
    time_second_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        J_0_T_CUDA(x, res2, nJ0);

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_second_function += elapsed_ms.count();
    }
    std::cout << "J_0_T_CUDA " << time_second_function / count_experiments << std::endl;

    for (int i = 0; i < nJ0; i++)
    {
        if (abs(res1[i] - res2[i]) > 1E-5)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ0_CUDA failed! point:" << x[i] << " |res1-resCPU|=" << abs(res1[i] - res2[i]) << std::endl << std::endl;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ_0_T_CUDA OK" << std::endl << std::endl;
}

void TestJ_1_T_CUDA()
{
    std::cout << "TestJ_1_T_CUDA started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* x = new double[nJ1];
    double* res1 = new double[nJ1];
    double* res2 = new double[nJ1];
    for (int i = 0; i < nJ1; i++)
    {
        x[i] = i * h1;
    }

    time_first_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        BesselWithCuda(v, x, res1, nJ1);

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_first_function += elapsed_ms.count();
    }
    std::cout << "BesselWithCuda " << time_first_function / count_experiments << std::endl;
    time_second_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        J_1_T_CUDA(x, res2, nJ1);

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_second_function += elapsed_ms.count();
    }
    std::cout << "J_1_T_CUDA " << time_second_function / count_experiments << std::endl;

    for (int i = 0; i < nJ1; i++)
    {
        if (abs(res1[i] - res2[i]) > 1E-6)
        {
            std::cout << "WARNING!!!" << std::endl;
            std::cout << "TestJ_1_T_CUDA failed! point:" << x[i] << " |res1-resCPU|=" << abs(res1[i] - res2[i]) << std::endl << std::endl;
            break;
        }
    }
    delete[] x;
    delete[] res1;
    delete[] res2;
    if (successfully)
        std::cout << "TestJ_1_T_CUDA OK" << std::endl << std::endl;
}

void Measure_J0_Time()
{
    std::cout << "Measure_J0_Time started" << std::endl;
    int v = 0;
    bool successfully = true;
    double* res1 = new double[nJ0];
    double* res2 = new double[nJ0];
    double* x = new double[nJ0];
    for (int i = 0; i < nJ0; i++)
    {
        x[i] = i * h0;
    }
    {
        J(v, x, res1, nJ0);
    }
    std::cout << "J ended" << std::endl;
    time_second_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        for (int i = 0; i < nJ0; i++)
        {
            res2[i] = J_0(x[i]);
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_second_function += elapsed_ms.count();
    }
    std::cout << "J0 " << time_second_function / count_experiments << std::endl;
    double accuracy = 0;
    double diff;
    for (int i = 0; i < nJ0; i++)
    {
        diff = std::abs(res1[i] - res2[i]);
        if (diff > accuracy) {
            accuracy = diff;
        }
    }
    std::cout << "Accuracy " << accuracy << std::endl;
    delete[] x;
    delete[] res1;
    delete[] res2;
}

void Measure_J1_Time()
{
    std::cout << "Measure_J1_Time started" << std::endl;
    int v = 1;
    bool successfully = true;
    double* res1 = new double[nJ1];
    double* res2 = new double[nJ1];
    double* x = new double[nJ1];
    for (int i = 0; i < nJ1; i++)
    {
        x[i] = i * h1;
    }
    {
        J(v, x, res1, nJ1);
    }
    std::cout << "J ended" << std::endl;
    time_second_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        for (int i = 0; i < nJ1; i++)
        {
            res2[i] = J_1(x[i]);
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_second_function += elapsed_ms.count();
    }
    std::cout << "J1 " << time_second_function / count_experiments << std::endl;
    double accuracy = 0;
    double diff;
    for (int i = 0; i < nJ1; i++)
    {
        diff = std::abs(res1[i] - res2[i]);
        if (diff > accuracy) {
            accuracy = diff;
        }
    }
    std::cout << "Accuracy " << accuracy << std::endl;
    delete[] x;
    delete[] res1;
    delete[] res2;
}

void Measure_Y0_Time()
{
    std::cout << "Measure_Y0_Time started" << std::endl;
    int v = 0;
    double* res1 = new double[nY0];
    double* res2 = new double[nY0];
    double* x = new double[nY0];
    for (int i = 0; i < nY0; i++)
    {
        x[i] = i * hY0;
    }
    double* Js = new double[nY0];
    J(v, x, Js, nY0);
    {
        Neumann(v, x, res1, nY0, Js);
        for (int i = 0; i < nY0; i++) 
        {
            std::cout << x[i] << " " << res1[i] << std::endl;
        }
    }
    std::cout << "Neumann ended" << std::endl;
    time_second_function = 0;
    for (size_t i = 0; i < count_experiments; i++)
    {
        auto begin = std::chrono::steady_clock::now();

        for (int i = 0; i < nY0; i++)
        {
            res2[i] = Y_0(x[i], Js[i]);
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        time_second_function += elapsed_ms.count();
    }
    std::cout << "Y_0 " << time_second_function / count_experiments << std::endl;
    double accuracy = 0;
    double diff;
    for (int i = 0; i < nY0; i++)
    {
        diff = std::abs(res1[i] - res2[i]);
        if (diff > accuracy) {
            accuracy = diff;
        }
    }
    std::cout << "Accuracy " << accuracy << std::endl;
    delete[] x;
    delete[] res1;
    delete[] res2;
}