#include "CPUfunctions.h"

double T(int n, double x)
{
	if (n == 0)
		return 1;
	if (n == 1)
		return x;
	double T_n2 = 1; // T_{n-2}
	double T_n1 = x; // T_{n-1}
	double T_n;
	for (int i = 2; i <= n; i++)
	{
		T_n = 2 * x * T_n1 - T_n2;
		T_n2 = T_n1;
		T_n1 = T_n;
	}
	return T_n;
}