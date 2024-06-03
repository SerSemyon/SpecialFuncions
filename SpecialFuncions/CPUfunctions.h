#pragma once
#include <complex>

long long Fact(int x);

/// <summary>
/// Вычисление гамма-функции при x>-3
/// </summary>
double Gamma(double x);

/// <summary>
/// Вычисляет значения любой цилиндрической фукнции порядка v+1 используя значения известные значения порядка v, v-1.
/// </summary>
/// <param name="v"> Порядок функции </param>
/// <param name="x"> Значение параметра </param>
/// <param name="Z_v"> Значение фукнции порядка v </param>
/// <param name="Z_vPrev"> Значение функции порядка v-1</param>
/// <returns></returns>
double cyl_next_order(double v, double x, double value_v, double value_v_minus_1);

void cyl_next_order(double v, double* x, double* result, unsigned int size, double* value_v, double* value_v_minus_1);

/// <summary>
/// Вычисление функции Бесселя на CPU для одной точки
/// </summary>
/// <param name="v"> порядок функции </param>
/// <param name="x"> значения параметра </param>
double J(double v, double x);

double J_asymptotic(const double v, const double x);

double Y_asymptotic(const double v, const double x);

/// <summary>
/// Вычисление функции Бесселя на CPU
/// </summary>
/// <param name="v"> порядок функции </param>
/// <param name="x"> значения параметра </param>
/// <param name="result"> полученные значения </param>
/// <param name="size"> количество точек </param>
void J(const double v, const double* x, double* result, const unsigned int size);

/// <summary>
/// Вычисление функции Бесселя нулевого порядка на отрезке [-8;8] 
/// </summary>
/// <param name="x"> Значение параметра </param>
double J_0(double x);

/// <summary>
/// Вычисление функции Бесселя нулевого порядка на отрезке [-8;8] 
/// </summary>
/// <param name="x"> Значение параметра </param>
double J_1(double x);

/// <summary>
/// Вычисление полинома Чебышёва первого рода
/// </summary>
/// <param name="n"> порядок полинома </param>
/// <param name="x"> значения параметра </param>
double T(int n, double x);

/// <summary>
/// Вычисление функции Неймана целого порядка
/// </summary>
/// <param name="v"> Порядок функции </param>
/// <param name="x"> Значение параметра </param>
/// <param name="res"> результат вычислений </param>
/// <param name="n"> количество точек </param>
/// <param name="Jpositive"> значения функции Бесселя порядка v </param>
void Neumann(int v, double* x, double* res, int n, double* Jpositive);
double Neumann(int v, double x, double Jpositive);

/// <summary>
/// Вычисление функции Неймана
/// </summary>
/// <param name="v"> Порядок функции </param>
/// <param name="x"> Значение параметра </param>
/// <param name="res"> результат вычислений </param>
/// <param name="n"> количество точек </param>
/// <param name="Jpositive"> значения функции Бесселя порядка v </param>
/// <param name="Jnegative"> значения функции Бесселя порядка -v </param>
void Neumann(double v, double* x, double* res, int n, double* Jpositive, double* Jnegative);

/// <summary>
/// Вычисление функции Неймана нулевого порядка на (0;8]
/// </summary>
/// <param name="x"> Значение параметра </param>
/// <param name="J0"> значения функции Бесселя нулевого порядка </param>
double Y_0(double x, double J0);

/// <summary>
/// Вычисление функции Неймана первого порядка на (0;8]
/// </summary>
/// <param name="x"> Значение параметра </param>
/// <param name="J1"> Значение функции Бесселя первого порядка </param>
double Y_1(double x, double J1);

/// <summary>
/// Вычисление функции Неймана нулевого порядка на (0;8]
/// </summary>
/// <param name="x"> Значения параметра </param>
/// <param name="res"> Результат вычислений </param>
/// <param name="n"> Количество точек </param>
/// <param name="J0"> Значения функции Бесселя нулевого порядка </param>
void Y_0(double* x, double* res, int n, double* J0);

/// <summary>
/// Вычисление функции Неймана первого порядка на (0;8]
/// </summary>
/// <param name="x"> Значения параметра </param>
/// <param name="res"> Результат вычислений </param>
/// <param name="n"> Количество точек </param>
/// <param name="J1"> Значения функции Бесселя первого порядка </param>
void Y_1(double* x, double* res, int n, double* J1);

/// <summary>
/// Мой вариант вычисления функции Бесселя для упорядоченного по возрастанию набора точек
/// </summary>
/// <param name="v"> Порядок функции </param>
/// <param name="x"> Значение параметра </param>
/// <param name="res"> Указатель на результат </param>
/// <param name="n"> Количество точек </param>
void BesselOrderedSet(double v, double* x, double* res, int n);
void BesselOrderedSet(double* x, double v, std::complex<double>* res, int n);

/// <summary>
/// Мой вариант вычисления функции Бесселя
/// </summary>
/// <param name="v"> Порядок функции </param>
/// <param name="x"> Значение параметра </param>
/// <param name="res"> Указатель на результат </param>
/// <param name="n"> Количество точек </param>
void Jnew(double v, double* x, double* res, int n);

/// <summary>
/// Вычисление производной любой цилиндрической функции через известные значения функции порядка v-1 и v+1.
/// </summary>
/// <param name="v"> Порядок функции </param>
/// <param name="Z_vPrev"> Значение функции прядка v-1 </param>
/// <param name="Z_vNext"> Значение функции прядка v+1 </param>
/// <returns></returns>
double dZ(double v, double Z_vPrev, double Z_vNext);

std::complex<double> dZ(double v, std::complex<double> Z_vPrev, std::complex<double> Z_vNext);

void H1(const double v, const double* const x, double* Re, double* Im, const unsigned int size);

std::complex<double> H1(const int v, double x);

void dZ(double v, double* x, double* result, unsigned int size, double* Z_vPrev, double* Z_vNext);

void J_negative(const int n, double* result, const unsigned int size, const double* const J_positive);

/// <summary>
/// Вычисление функции Бесселя отрицательного целочисленного порядка, через уже вычисленное значение положительного порядка
/// </summary>
/// <param name="J_positive"> Вычисленные значения функции Бесселя порядка n </param>
/// <param name="n"> -Порядок функции </param>
/// <param name="res"> Указатель на результат, в который будет сохранены значения функции Бесселя порядка -n </param>
/// <param name="size"> Количество точек </param>
void J_negative(std::complex<double>* J_positive, int n, std::complex<double>* res, int size);

std::complex<double> J_negative(std::complex<double> J_positive, int n);

std::complex<double> H_negative(int nu, std::complex<double> H);

void H1(const int v, double* x, std::complex<double>* res, const unsigned int size);

/// <summary>
/// Вычисление функции Бесселя с ожидаемой точностью 10^-6
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double J_0_T(double x);

/// <summary>
/// Вычисление функции Бесселя с ожидаемой точностью 10^-7
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double J_1_T(double x);