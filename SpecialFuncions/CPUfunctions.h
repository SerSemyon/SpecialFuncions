#pragma once
#include <complex>

long long Fact(int x);

/// <summary>
/// ���������� �����-������� ��� x>-3
/// </summary>
double Gamma(double x);

/// <summary>
/// ��������� �������� ����� �������������� ������� ������� v+1 ��������� �������� ��������� �������� ������� v, v-1.
/// </summary>
/// <param name="v"> ������� ������� </param>
/// <param name="x"> �������� ��������� </param>
/// <param name="Z_v"> �������� ������� ������� v </param>
/// <param name="Z_vPrev"> �������� ������� ������� v-1</param>
/// <returns></returns>
double cyl_next_order(double v, double x, double value_v, double value_v_minus_1);

void cyl_next_order(double v, double* x, double* result, unsigned int size, double* value_v, double* value_v_minus_1);

/// <summary>
/// ���������� ������� ������� �� CPU ��� ����� �����
/// </summary>
/// <param name="v"> ������� ������� </param>
/// <param name="x"> �������� ��������� </param>
double J(double v, double x);

double J_asymptotic(const double v, const double x);

double Y_asymptotic(const double v, const double x);

/// <summary>
/// ���������� ������� ������� �� CPU
/// </summary>
/// <param name="v"> ������� ������� </param>
/// <param name="x"> �������� ��������� </param>
/// <param name="result"> ���������� �������� </param>
/// <param name="size"> ���������� ����� </param>
void J(const double v, const double* x, double* result, const unsigned int size);

/// <summary>
/// ���������� ������� ������� �������� ������� �� ������� [-8;8] 
/// </summary>
/// <param name="x"> �������� ��������� </param>
double J_0(double x);

/// <summary>
/// ���������� ������� ������� �������� ������� �� ������� [-8;8] 
/// </summary>
/// <param name="x"> �������� ��������� </param>
double J_1(double x);

/// <summary>
/// ���������� �������� �������� ������� ����
/// </summary>
/// <param name="n"> ������� �������� </param>
/// <param name="x"> �������� ��������� </param>
double T(int n, double x);

/// <summary>
/// ���������� ������� ������� ������ �������
/// </summary>
/// <param name="v"> ������� ������� </param>
/// <param name="x"> �������� ��������� </param>
/// <param name="res"> ��������� ���������� </param>
/// <param name="n"> ���������� ����� </param>
/// <param name="Jpositive"> �������� ������� ������� ������� v </param>
void Neumann(int v, double* x, double* res, int n, double* Jpositive);
double Neumann(int v, double x, double Jpositive);

/// <summary>
/// ���������� ������� �������
/// </summary>
/// <param name="v"> ������� ������� </param>
/// <param name="x"> �������� ��������� </param>
/// <param name="res"> ��������� ���������� </param>
/// <param name="n"> ���������� ����� </param>
/// <param name="Jpositive"> �������� ������� ������� ������� v </param>
/// <param name="Jnegative"> �������� ������� ������� ������� -v </param>
void Neumann(double v, double* x, double* res, int n, double* Jpositive, double* Jnegative);

/// <summary>
/// ���������� ������� ������� �������� ������� �� (0;8]
/// </summary>
/// <param name="x"> �������� ��������� </param>
/// <param name="J0"> �������� ������� ������� �������� ������� </param>
double Y_0(double x, double J0);

/// <summary>
/// ���������� ������� ������� ������� ������� �� (0;8]
/// </summary>
/// <param name="x"> �������� ��������� </param>
/// <param name="J1"> �������� ������� ������� ������� ������� </param>
double Y_1(double x, double J1);

/// <summary>
/// ���������� ������� ������� �������� ������� �� (0;8]
/// </summary>
/// <param name="x"> �������� ��������� </param>
/// <param name="res"> ��������� ���������� </param>
/// <param name="n"> ���������� ����� </param>
/// <param name="J0"> �������� ������� ������� �������� ������� </param>
void Y_0(double* x, double* res, int n, double* J0);

/// <summary>
/// ���������� ������� ������� ������� ������� �� (0;8]
/// </summary>
/// <param name="x"> �������� ��������� </param>
/// <param name="res"> ��������� ���������� </param>
/// <param name="n"> ���������� ����� </param>
/// <param name="J1"> �������� ������� ������� ������� ������� </param>
void Y_1(double* x, double* res, int n, double* J1);

/// <summary>
/// ��� ������� ���������� ������� ������� ��� �������������� �� ����������� ������ �����
/// </summary>
/// <param name="v"> ������� ������� </param>
/// <param name="x"> �������� ��������� </param>
/// <param name="res"> ��������� �� ��������� </param>
/// <param name="n"> ���������� ����� </param>
void BesselOrderedSet(double v, double* x, double* res, int n);
void BesselOrderedSet(double* x, double v, std::complex<double>* res, int n);

/// <summary>
/// ��� ������� ���������� ������� �������
/// </summary>
/// <param name="v"> ������� ������� </param>
/// <param name="x"> �������� ��������� </param>
/// <param name="res"> ��������� �� ��������� </param>
/// <param name="n"> ���������� ����� </param>
void Jnew(double v, double* x, double* res, int n);

/// <summary>
/// ���������� ����������� ����� �������������� ������� ����� ��������� �������� ������� ������� v-1 � v+1.
/// </summary>
/// <param name="v"> ������� ������� </param>
/// <param name="Z_vPrev"> �������� ������� ������ v-1 </param>
/// <param name="Z_vNext"> �������� ������� ������ v+1 </param>
/// <returns></returns>
double dZ(double v, double Z_vPrev, double Z_vNext);

std::complex<double> dZ(double v, std::complex<double> Z_vPrev, std::complex<double> Z_vNext);

void H1(const double v, const double* const x, double* Re, double* Im, const unsigned int size);

std::complex<double> H1(const int v, double x);

void dZ(double v, double* x, double* result, unsigned int size, double* Z_vPrev, double* Z_vNext);

void J_negative(const int n, double* result, const unsigned int size, const double* const J_positive);

/// <summary>
/// ���������� ������� ������� �������������� �������������� �������, ����� ��� ����������� �������� �������������� �������
/// </summary>
/// <param name="J_positive"> ����������� �������� ������� ������� ������� n </param>
/// <param name="n"> -������� ������� </param>
/// <param name="res"> ��������� �� ���������, � ������� ����� ��������� �������� ������� ������� ������� -n </param>
/// <param name="size"> ���������� ����� </param>
void J_negative(std::complex<double>* J_positive, int n, std::complex<double>* res, int size);

std::complex<double> J_negative(std::complex<double> J_positive, int n);

std::complex<double> H_negative(int nu, std::complex<double> H);

void H1(const int v, double* x, std::complex<double>* res, const unsigned int size);

/// <summary>
/// ���������� ������� ������� � ��������� ��������� 10^-6
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double J_0_T(double x);

/// <summary>
/// ���������� ������� ������� � ��������� ��������� 10^-7
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double J_1_T(double x);