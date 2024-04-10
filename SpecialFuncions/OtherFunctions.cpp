#include <math.h>
#include "CPUfunctions.h"
#include <cmath>

double cyl_next_order(double v, double x, double value_v, double value_v_minus_1)
{
	return 2 * v * value_v / x - value_v_minus_1;
}

void cyl_next_order(double v, double* x, double* result, unsigned int size, double* value_v, double* value_v_minus_1)
{
	for (int i = 0; i < size; i++)
	{
		result[i] = 2 * v * value_v[i] / x[i] - value_v_minus_1[i];
	}
}

void dZ(double v, double* x, double* result, unsigned int size, double* Z_vPrev, double* Z_vNext)
{
	for (int i = 0; i < size; i++)
	{
		result[i] = 0.5 * (Z_vPrev[i] - Z_vNext[i]);
	}
}

double dZ(double v, double Z_vPrev, double Z_vNext)
{
	return 0.5 * (Z_vPrev - Z_vNext);
}

std::complex<double> dZ(double v, std::complex<double> Z_vPrev, std::complex<double> Z_vNext)
{
	return 0.5 * (Z_vPrev - Z_vNext);
}