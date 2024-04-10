#include "CylindricalBody.h"

std::complex<double> CylindricalBody::B_n(int n, double k_1, double k_2,
    std::complex<double> J_nk_1, std::complex<double> J_nk_2, std::complex<double> dJ_nk_1, std::complex<double> dJ_nk_2,
    std::complex<double> H_n, std::complex<double> dH_n)
{
    std::complex<double> res = a_0 * (J_nk_1 * dJ_nk_2 * k_2 - J_nk_2 * dJ_nk_1 * k_1)
        / (J_nk_2 * dH_n * k_1 - H_n * dJ_nk_2 * k_2);
    if (n == 0)
        return res;
    if (n % 2 != 0)
        return res * std::complex<double>(0, 1);
    return -res;
}

std::complex<double> CylindricalBody::C_n(int n, double k_1, double k_2,
    std::complex<double> J_nk_1, std::complex<double> J_nk_2, std::complex<double> dJ_nk_1, std::complex<double> dJ_nk_2,
    std::complex<double> H_n, std::complex<double> dH_n)
{
    std::complex<double> res = a_0 * (H_n * dJ_nk_1 * k_1 - J_nk_1 * dH_n * k_1)
        / (H_n * dJ_nk_2 * k_2 - J_nk_2 * dH_n * k_1);
    if (n == 0)
        return res;
    if (n % 2 != 0)
        return res * std::complex<double>(0, 1);
    return -res;
}

std::complex<double> CylindricalBody::find_B_0(double k_1, double k_2, double R)
{
    std::complex<double> prevJk_1;
    std::complex<double> prevJk_2;
    std::complex<double> prevH;
    std::complex<double> Jk_1;
    std::complex<double> Jk_2;
    std::complex<double> H;
    std::complex<double> nextJk_1;
    std::complex<double> nextJk_2;
    std::complex<double> nextH;
    std::complex<double> dJk_1;
    std::complex<double> dJk_2;
    std::complex<double> dH;
    double x_1 = k_1 * R;
    double x_2 = k_2 * R;
    Jk_1 = J_0(x_1);
    Jk_2 = J_0(x_2);
    H = std::complex<double>(Jk_1.real(), Y_0(x_1, Jk_1.real()));
    nextJk_1 = J(x_1, 1);
    nextJk_2 = J(x_2, 1);
    nextH = H1(1, x_1);
    prevJk_1 = -nextJk_1;
    prevJk_2 = -nextJk_2;
    prevH = H_negative(1, nextH);
    dJk_1 = dZ(0, prevJk_1, nextJk_1);
    dJk_2 = dZ(0, prevJk_2, nextJk_2);
    dH = dZ(0, prevH, nextH);
    return B_n(0, k_1, k_2, Jk_1, Jk_2, dJk_1, dJk_2, H, dH);
}

std::complex<double> CylindricalBody::find_C_0(double k_1, double k_2, double R)
{
    std::complex<double> prevJk_1;
    std::complex<double> prevJk_2;
    std::complex<double> prevH;
    std::complex<double> Jk_1;
    std::complex<double> Jk_2;
    std::complex<double> H;
    std::complex<double> nextJk_1;
    std::complex<double> nextJk_2;
    std::complex<double> nextH;
    std::complex<double> dJk_1;
    std::complex<double> dJk_2;
    std::complex<double> dH;
    double x_1 = k_1 * R;
    double x_2 = k_2 * R;
    Jk_1 = J_0(x_1);
    Jk_2 = J_0(x_2);
    H = std::complex<double>(Jk_1.real(), Y_0(x_1, Jk_1.real()));
    nextJk_1 = J(x_1, 1);
    nextJk_2 = J(x_2, 1);
    nextH = H1(1, x_1);
    prevJk_1 = -nextJk_1;
    prevJk_2 = -nextJk_2;
    prevH = H_negative(1, nextH);
    dJk_1 = dZ(0, prevJk_1, nextJk_1);
    dJk_2 = dZ(0, prevJk_2, nextJk_2);
    dH = dZ(0, prevH, nextH);
    return C_n(0, k_1, k_2, Jk_1, Jk_2, dJk_1, dJk_2, H, dH);
}

std::complex<double> CylindricalBody::find_B_n(int N, double k_1, double k_2, double R)
{
    std::complex<double> prevJk_1;
    std::complex<double> prevJk_2;
    std::complex<double> prevH;
    std::complex<double> Jk_1;
    std::complex<double> Jk_2;
    std::complex<double> H;
    std::complex<double> nextJk_1;
    std::complex<double> nextJk_2;
    std::complex<double> nextH;
    std::complex<double> dJk_1;
    std::complex<double> dJk_2;
    std::complex<double> dH;
    bool positive = (N > 0);
    if (!positive)
        N = -N;
    double x_1 = k_1 * R;
    double x_2 = k_2 * R;
    Jk_1 = J(x_1, N);
    Jk_2 = J(x_2, N);
    H = std::complex<double>(Jk_1.real(), Neumann(N, x_1, Jk_1.real()));
    nextJk_1 = J(x_1, N + 1);
    nextJk_2 = J(x_2, N + 1);
    nextH = H1(N + 1, x_1);
    prevJk_1 = J(x_1, N - 1);
    prevJk_2 = J(x_2, N - 1);
    prevH = H1(N - 1, x_1);
    if (!positive) {
        Jk_1 = J_negative(Jk_1, N);
        Jk_2 = J_negative(Jk_2, N);
        H = H_negative(N, H);
        nextJk_1 = J_negative(nextJk_1, N + 1);
        nextJk_2 = J_negative(nextJk_2, N + 1);
        nextH = H_negative(N + 1, nextH);
        prevJk_1 = J_negative(prevJk_1, N - 1);
        prevJk_2 = J_negative(prevJk_2, N - 1);
        prevH = H_negative(N - 1, prevH);
    }
    dJk_1 = dZ(N, prevJk_1, nextJk_1);
    dJk_2 = dZ(N, prevJk_2, nextJk_2);
    dH = dZ(N, prevH, nextH);
    return B_n(N, k_1, k_2, Jk_1, Jk_2, dJk_1, dJk_2, H, dH);
}

std::complex<double> CylindricalBody::find_C_n(int N, double k_1, double k_2, double R)
{
    std::complex<double> prevJk_1;
    std::complex<double> prevJk_2;
    std::complex<double> prevH;
    std::complex<double> Jk_1;
    std::complex<double> Jk_2;
    std::complex<double> H;
    std::complex<double> nextJk_1;
    std::complex<double> nextJk_2;
    std::complex<double> nextH;
    std::complex<double> dJk_1;
    std::complex<double> dJk_2;
    std::complex<double> dH;
    bool positive = (N > 0);
    if (!positive)
        N = -N;
    double x_1 = k_1 * R;
    double x_2 = k_2 * R;
    Jk_1 = J(x_1, N);
    Jk_2 = J(x_2, N);
    H = std::complex<double>(Jk_1.real(), Neumann(N, x_1, Jk_1.real()));
    nextJk_1 = J(x_1, N + 1);
    nextJk_2 = J(x_2, N + 1);
    nextH = H1(N + 1, x_1);
    prevJk_1 = J(x_1, N - 1);
    prevJk_2 = J(x_2, N - 1);
    prevH = H1(N - 1, x_1);
    if (!positive) {
        Jk_1 = J_negative(Jk_1, N);
        Jk_2 = J_negative(Jk_2, N);
        H = H_negative(N, H);
        nextJk_1 = J_negative(nextJk_1, N + 1);
        nextJk_2 = J_negative(nextJk_2, N + 1);
        nextH = H_negative(N + 1, nextH);
        prevJk_1 = J_negative(prevJk_1, N - 1);
        prevJk_2 = J_negative(prevJk_2, N - 1);
        prevH = H_negative(N - 1, prevH);
    }
    dJk_1 = dZ(N, prevJk_1, nextJk_1);
    dJk_2 = dZ(N, prevJk_2, nextJk_2);
    dH = dZ(N, prevH, nextH);
    return C_n(N, k_1, k_2, Jk_1, Jk_2, dJk_1, dJk_2, H, dH);
}

std::complex<double>** CylindricalBody::E_1(double k_1, double k_2, double* r, size_t count, size_t number_division, double R) {
    double* k_1r = new double[count];
    for (size_t i = 0; i < count; i++) {
        k_1r[i] = k_1 * r[i];
    }
    std::complex<double>* H_values = new std::complex<double>[count];
    std::complex<double>* H_negative_values = new std::complex<double>[count];
    size_t number_iterations = 4;
    std::complex<double>** E_z = new std::complex<double>*[count];
    for (int i = 0; i < count; i++) {
        E_z[i] = new std::complex<double>[number_division];
    }
    for (int i = 0; i < count; i++) {
        H_values[i] = H1(0, k_1r[i]);
    }
    std::complex<double> b_0 = find_B_0(k_1, k_2, R);
    double alpha = 2 * M_PI / number_division;
    std::complex<double> exp_in_pow;
    for (size_t i = 0; i < count; i++) {
        for (size_t j = 0; j < number_division; j++) {
            exp_in_pow = std::complex<double>(cos(j * alpha), sin(j * alpha));
            E_z[i][j] = b_0 * H_values[i] * exp_in_pow;
        }
    }
    std::complex<double> b_n;
    for (int i = 1; i < number_iterations; i++) {
        for (int j = 0; j < count; j++) {
            H_values[j] = H1(i, k_1r[j]);
        }
        b_n = find_B_n(i, k_1, k_2, R);
        for (size_t j = 0; j < count; j++) {
            for (size_t k = 0; k < number_division; k++) {
                exp_in_pow = std::complex<double>(cos(i * k * alpha), sin(i * k * alpha));
                E_z[j][k] += b_n * H_values[j] * exp_in_pow;
            }
        }
        for (int j = 0; j < count; j++) {
            H_negative_values[j] = H_negative(i, H_values[j]);
        }
        b_n = find_B_n(-i, k_1, k_2, R);
        for (size_t j = 0; j < count; j++) {
            for (size_t k = 0; k < number_division; k++) {
                exp_in_pow = std::complex<double>(cos(i * k * alpha), sin(i * k * alpha));
                E_z[j][k] += b_n * H_negative_values[j] * exp_in_pow;
            }
        }
    }
    return E_z;
}

std::complex<double>** CylindricalBody::E_2(double k_1, double k_2, double* r, size_t count, size_t number_division, double R) {
    double* k_2r = new double[count];
    for (size_t i = 0; i < count; i++) {
        k_2r[i] = k_2 * r[i];
    }
    std::complex<double>* J_values = new std::complex<double>[count];
    std::complex<double>* J_negative_values = new std::complex<double>[count];
    size_t number_iterations = 4;
    std::complex<double>** E_z = new std::complex<double>*[count];
    for (int i = 0; i < count; i++) {
        E_z[i] = new std::complex<double>[number_division];
    }
    BesselOrderedSet(k_2r, 0, J_values, count);
    std::complex<double> c_0 = find_C_0(k_1, k_2, R);
    double alpha = 2 * M_PI / number_division;
    std::complex<double> exp_in_pow;
    for (size_t i = 0; i < count; i++) {
        for (size_t j = 0; j < number_division; j++) {
            exp_in_pow = std::complex<double>(cos(j * alpha), sin(j * alpha));
            E_z[i][j] = c_0 * J_values[i] * exp_in_pow;
        }
    }
    std::complex<double> c_n;
    for (int i = 1; i < number_iterations; i++) {
        BesselOrderedSet(k_2r, i, J_values, count);
        c_n = find_C_n(i, k_1, k_2, R);
        for (size_t j = 0; j < count; j++) {
            for (size_t k = 0; k < number_division; k++) {
                exp_in_pow = std::complex<double>(cos(i * k * alpha), sin(i * k * alpha));
                E_z[j][k] += c_n * J_values[j] * exp_in_pow;
            }
        }
        J_negative(J_values, i, J_negative_values, count);
        c_n = find_C_n(-i, k_1, k_2, R);
        for (size_t j = 0; j < count; j++) {
            for (size_t k = 0; k < number_division; k++) {
                exp_in_pow = std::complex<double>(cos(i * k * alpha), sin(i * k * alpha));
                E_z[j][k] += c_n * J_negative_values[j] * exp_in_pow;
            }
        }
    }
    return E_z;
}