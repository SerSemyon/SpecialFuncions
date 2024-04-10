#define _USE_MATH_DEFINES
#include <math.h>
#include "CPUfunctions.h"
#include <cmath>

const double epsilon = 1E-12;
const double C = 0.5772156;

int max_iter = 10000;

void Neumann(int v, double* x, double* res, int size, double* Jpositive)
{
	unsigned int factV_minus_one = Fact(v - 1);
	unsigned int factV = factV_minus_one * v;
	unsigned int factKN = factV;
	double S_2 = 0;
	unsigned int factK = 1;
	int sign = 1;
	double M = 1.0 / factV;
	for (int i = 1; i <= v; i++)
	{
		S_2 += 1.0 / i;
	}
	M *= S_2;
	double* sumWithPsi = new double[size];
	for (int i = 0; i < size; i++)
	{
		sumWithPsi[i] = M;
	}
	for (int k = 1; k < 15; k++)
	{
		sign = -sign;
		factK *= k;
		factKN *= (v + k);
		S_2 += 1.0 / k + 1.0 / (k + v);
		M = sign * S_2 / (factK * factKN);
		for (int i = 0; i < size; i++)
		{
			sumWithPsi[i] += std::pow(0.5 * x[i], 2 * k) * M;
		}
	}
	for (int i = 0; i < size; i++)
	{
		sumWithPsi[i] *= std::pow(0.5 * x[i], v);
	}
	double* S_1 = new double[size];
	double f_1 = factV_minus_one;
	for (int i = 0; i < size; i++)
	{
		S_1[i] = f_1;
	}
	for (int k = 1; k < v; k++)
	{
		f_1 *= (v - k - 1);
		f_1 /= k;
		for (int i = 0; i < size; i++)
		{
			S_1[i] += f_1 * std::pow(0.5 * x[i], 2 * k);
		}
	}
	for (int i = 0; i < size; i++)
	{
		res[i] = 2 * Jpositive[i] * (std::log(0.5 * x[i]) + C);
		res[i] -= sumWithPsi[i];
		if (v > 0)
			res[i] -= std::pow(0.5 * x[i], -v) * S_1[i];
		res[i] /= M_PI;
	}
}

double Neumann(int v, double x, double Jpositive)
{
	long long factV_minus_one = Fact(v - 1);
	long long factV = Fact(v);
	double S_2 = 0;
	int sign = 1;
	double M = 1.0 / factV;
	for (int i = 1; i <= v; i++)
	{
		S_2 += 1.0 / i;
	}
	M *= S_2;
	double sumWithPsi = M;
	double prev;
	double diff;
	double cur = M;
	double mul = 1.0/ factV;
	int k = 1;
	do {
		prev = cur;
		sign = -sign;
		mul /= k * (v + k);
		S_2 += 1.0 / k + 1.0 / (k + v);
		M = sign * S_2 * mul;
		cur = std::pow(0.5 * x, 2 * k) * M;
		sumWithPsi += cur;
		k++;
		if (k > max_iter)
			break;
		diff = abs(cur - prev);
	} while (diff > epsilon);
	sumWithPsi *= std::pow(0.5 * x, v);
	double f_1 = factV_minus_one;
	double S_1 = f_1;
	for (int k = 1; k < v; k++)
	{
		double multiply = v - k - 1;
		if (multiply != 0)
			f_1 *= multiply;
		f_1 /= k;
		S_1 += f_1 * std::pow(0.5 * x, 2 * k);
	}
	double res = 2 * Jpositive * (std::log(0.5 * x) + C);
	res -= sumWithPsi;
	if (v > 0)
		res -= std::pow(0.5 * x, -v) * S_1;
	res /= M_PI;
	return res;
}

void Neumann(double v, double* x, double* res, int n, double* Jpositive, double* Jnegative)
{
	double arg = v * M_PI;
	if (v != (long long)v) // Если не целое
	{
		for (int i = 0; i < n; i++)
		{
			res[i] = (Jpositive[i] * cos(arg) - Jnegative[i]) / sin(arg);
		}
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			res[i] = Neumann(v, x[i], Jpositive[i]);
		}
		//Neumann(v, x, res, n, Jpositive);
	}
}

double Y_asymptotic(const double v, const double x)
{
	return -1.0 / M_PI * pow(x / 2, -v) * Gamma(v) * exp(0.25 * x * x / v);
}

const double eps = 1E-12;
const double b0[] = {
		-0.02150'51114'49657'55061,
		-0.27511'81330'43518'79146,
		0.19860'56347'02554'15556,
		0.23425'27461'09021'80210,
		-0.16563'59817'13650'41312,
		0.04462'13795'40669'28217,
		-0.00693'22862'91523'18829,
		0.00071'91174'03752'30309,
		-0.00005'39250'79722'93939,
		0.00000'30764'93288'10848,
		-0.00000'01384'57181'23009,
		0.00000'00050'51054'36909,
		-0.00000'00001'52582'85043,
		0.00000'00000'03882'86747,
		-0.00000'00000'00084'42875,
		0.00000'00000'00001'58748,
		-0.00000'00000'00000'02608,
		0.00000'00000'00000'00038
};

double Y_0(double x, double J0) 
{
	double T2 = x / 8.0;
	T2 = 2.0 * T2 * T2 - 1.0;
	double T_previous = 1.0; double T_current = T2;
	double T;
	double sum = b0[0] * T_previous + b0[1] * T_current;
	for (int n = 2; n <= 17; n++) {
		T = 2.0 * T2 * T_current - T_previous;
		sum += b0[n] * T;
		T_previous = T_current; T_current = T;
	};
	sum += (log(x / 2.0) + C) * J0 * 2.0 / M_PI;
	return sum;
};

const double b1[] = {
	-0.04017'29465'44414'07579,
	-0.44444'71476'30558'06261,
	-0.02271'92444'28417'73587,
	0.20664'45410'17490'51976,
	-0.08667'16970'56948'52366,
	0.01763'67030'03163'13441,
	-0.00223'56192'94485'09524,
	0.00019'70623'02701'54078,
	-0.00001'28858'53299'24086,
	0.00000'06528'47952'35852,
	-0.00000'00264'50737'17479,
	0.00000'00008'78030'11712,
	-0.00000'00000'24343'27870,
	0.00000'00000'00572'61216,
	-0.00000'00000'00011'57794,
	0.00000'00000'00000'20347,
	-0.00000'00000'00000'00314,
	0.00000'00000'00000'00004
};

double Y_1(double x, double J1) 
{
	double z = x / 8.0;
	double T_previous = 1.0, T_current = z;
	double T;
	double s = b1[0] * T_current;
	for (int n = 1; n <= 17; n++) {
		T = 2.0 * z * T_current - T_previous;
		T_previous = T_current; T_current = T;
		T = 2.0 * z * T_current - T_previous;
		s += b1[n] * T;
		T_previous = T_current; T_current = T;
	};
	s += (C + log(x / 2.0)) * J1 * 2.0 / M_PI - 2.0 / (M_PI * x);
	return s;
};


void Y_0(double* x, double* res, int n, double* J0)
{
	for (int i = 0; i < n; i++)
	{
		double T2 = x[i] / 8.0;
		T2 = 2.0 * T2 * T2 - 1.0;
		double T_previous = 1.0; double T_current = T2;
		double T;
		double sum = b0[0] * T_previous + b0[1] * T_current;
		for (int n = 2; n <= 17; n++) {
			T = 2.0 * T2 * T_current - T_previous;
			sum += b0[n] * T;
			T_previous = T_current; T_current = T;
		};
		sum += (log(x[i] / 2.0) + C) * J0[i] * 2.0 / M_PI;
		res[i] = sum;
	}
}

void Y_1(double* x, double* res, int n, double* J1) 
{
	for (int i = 0; i < n; i++)
	{
		double z = x[i] / 8.0;
		double T_previous = 1.0, T_current = z;
		double T;
		double s = b1[0] * T_current;
		for (int n = 1; n <= 17; n++) {
			T = 2.0 * z * T_current - T_previous;
			T_previous = T_current; T_current = T;
			T = 2.0 * z * T_current - T_previous;
			s += b1[n] * T;
			T_previous = T_current; T_current = T;
		};
		s += (C + log(x[i] / 2.0)) * J1[i] * 2.0 / M_PI - 2.0 / (M_PI * x[i]);
		res[i] = s;
	}
}