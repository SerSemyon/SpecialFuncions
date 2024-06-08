#include <math.h>
#include "CPUfunctions.h"
#include <cmath>

const double eps = 1E-12;

double J(double x, double v) {
	int k = 0;
	double aprev = 1 / Gamma(v + 1);
	double aNext;
	double summ = aprev;
	double diff;
	do {
		aNext = -x * x * aprev / ((k + 1) * (v + k + 1) * 4);
		summ += aNext;
		diff = abs(aprev - aNext);
		aprev = aNext;
		k++;
	} while (diff > eps);
	return summ * pow(x * 0.5, v);
}

void J(const double v, const double* x, double* result, const unsigned int size) {
	int max_iterations = 100;
    double aNext;
    double diff;
    for (unsigned int i = 0; i < size; i++) {
        int k = 0;
        double aprev = 1 / Gamma(v + 1);
        double summ = aprev;
        do {
            aNext = -x[i] * x[i] * aprev / ((k + 1) * (v + k + 1) * 4);
            summ += aNext;
            diff = abs(aprev - aNext);
            aprev = aNext;
            k++;
        } while (diff > eps);
        result[i] = summ * pow(x[i] * 0.5, v);
    }
}

double J_asymptotic (const double v, const double x) {
	return pow(x / 2, v) / Gamma(v + 1) * exp(-x * x / (4 * (v + 1)));
}

const double a0[18] =
{
	0.15772'79714'74890'11956,
	-0.00872'34423'52852'22129,
	0.26517'86132'03336'80987,
	-0.37009'49938'72649'77903,
	0.15806'71023'32097'26128,
	-0.03489'37694'11408'88516,
	0.00481'91800'69467'60450,
	-0.00046'06261'66206'27505,
	0.00003'24603'28821'00508,
	-0.00000'17619'46907'76215,
	0.00000'00760'81635'92419,
	-0.00000'00026'79253'53056,
	0.00000'00000'78486'96314,
	-0.00000'00000'01943'83469,
	0.00000'00000'00041'25321,
	-0.00000'00000'00000'75885,
	0.00000'00000'00000'01222,
	-0.00000'00000'00000'00017
};

double J_0(double x) {
	double z = x / 8.0; z = 2.0 * z * z - 1.0;
	double T_previous = 1.0, T_current = z;
	double T;
	double s = a0[0] * T_previous + a0[1] * T_current;
	for (int n = 2; n <= 17; n++) {
		T = 2.0 * z * T_current - T_previous;
		s += a0[n] * T;
		T_previous = T_current; T_current = T;
	};
	return s;
};

const double a1[18] =
{
	0.05245'81903'34656'48458,
	0.04809'64691'58230'37394,
	0.31327'50823'61567'18380,
	-0.24186'74084'47407'48475,
	0.07426'67962'16787'03781,
	-0.01296'76273'11735'17510,
	0.00148'99128'96667'63839,
	-0.00012'22786'85054'32427,
	0.00000'75626'30229'69605,
	-0.00000'03661'30855'23363,
	0.00000'00142'77324'38731,
	-0.00000'00004'58570'03076,
	0.00000'00000'12351'74811,
	-0.00000'00000'00283'17735,
	0.00000'00000'00005'59509,
	-0.00000'00000'00000'09629,
	0.00000'00000'00000'00146,
	-0.00000'00000'00000'00002
};

double J_1(double x) {
	double z = x / 8.0;
	double T_previous = z, T_current = 4 * z * z * z - 3 * z;
	double T;
	double s = a1[0] * T_previous + a1[1] * T_current;
	for (int n = 2; n <= 17; n++) {
		T = (4.0 * z * z - 2) * T_current - T_previous;
		s += a1[n] * T;
		T_previous = T_current; T_current = T;
	};
	return s;
};

void Jnew(double v, double* x, double* res, int n) {
	double eps = 1E-12;
	double aNext;
	double diff;
	int max_iterations = 100;
	bool* converged = new bool[n];
	int notConverged = n;
	double* aPrev = new double[n];
	double a_0 = 1 / Gamma(v + 1);
	for (int i = 0; i < n; i++)
	{
		aPrev[i] = a_0;
		converged[i] = false;
		res[i] = a_0;
	}
	double p_k;
	for (int k = 0; k < max_iterations; k++)
	{
		if (notConverged == 0)
			break;
		p_k = -1 / (4 * (v + k + 1) * (k + 1));
		for (int i = 0; i < n; i++)
		{
			if (!converged[i])
			{
				aNext = p_k * aPrev[i] * x[i] * x[i];
				res[i] += aNext;
				diff = abs(aPrev[i] - aNext);
				aPrev[i] = aNext;
				if (diff < eps)
				{
					converged[i] = true;
					notConverged--;
				}
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		res[i] *= pow(x[i] * 0.5, v);
	}
	delete[] converged;
	delete[] aPrev;
}

void BesselOrderedSet(double v, double* x, double* res, int n) {
	double aNext;
	double diff;
	int max_iterations = 100;
	int ind_begin_set = 0;
	double* aPrev = new double[n];
	double a_0 = 1 / Gamma(v + 1);
	for (int i = 0; i < n; i++)
	{
		aPrev[i] = a_0;
		res[i] = a_0;
	}
	double p_k;
	for (int k = 0; k < max_iterations; k++)
	{
		if (ind_begin_set > n - 2)
			break;
		p_k = -1 / (4 * (v + k + 1) * (k + 1));
		for (int i = ind_begin_set; i < n; i++)
		{
			aNext = p_k * aPrev[i] * x[i] * x[i];
			res[i] += aNext;
			diff = abs(aPrev[i] - aNext);
			aPrev[i] = aNext;
			if (diff < eps)
			{
				ind_begin_set++;
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		res[i] *= pow(x[i] * 0.5, v);
	}
	delete[] aPrev;
}

void BesselOrderedSet(double* x, double v, std::complex<double>* res, int n) {
	double eps = 1E-12;
	double aNext;
	double diff;
	int max_iterations = 100;
	int ind_begin_set = 0;
	double* aPrev = new double[n];
	double a_0 = 1 / Gamma(v + 1);
	for (int i = 0; i < n; i++)
	{
		aPrev[i] = a_0;
		res[i] = a_0;
	}
	double p_k;
	for (int k = 0; k < max_iterations; k++)
	{
		if (ind_begin_set > n - 2)
			break;
		p_k = -1 / (4 * (v + k + 1) * (k + 1));
		for (int i = ind_begin_set; i < n; i++)
		{
			aNext = p_k * aPrev[i] * x[i] * x[i];
			res[i] += aNext;
			diff = abs(aPrev[i] - aNext);
			aPrev[i] = aNext;
			if (diff < eps)
			{
				ind_begin_set++;
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		res[i] *= pow(x[i] * 0.5, v);
	}
}

void J_negative(const int n, double* result, const unsigned int size, const double* const J_positive)
{
	if (n % 2 == 0)
	{
		for (unsigned int i = 0; i < size; i++)
		{
			result[i] = J_positive[i];
		}
	}
	else
	{
		for (unsigned int i = 0; i < size; i++)
		{
			result[i] = -J_positive[i];
		}
	}
}

void J_negative(std::complex<double>* J_positive, int n, std::complex<double>* res, int size)
{
	if (n % 2 != 0)
	{
		for (int i = 0; i < size; i++)
		{
			res[i] = -J_positive[i];
		}
	}
	else
	{
		for (int i = 0; i < size; i++)
		{
			res[i] = J_positive[i];
		}
	}
}

std::complex<double> J_negative(std::complex<double> J_positive, int n)
{
	if (n % 2 != 0)
	{
		return -J_positive;
	}
	else
	{
		return J_positive;
	}
}

const double c_T_0[9] = {
	0.1577'2797'1474'8901'2,
	-0.008723442352852221,
	0.2651786132033368,
	-0.37009499387264977,
	0.15806710233209725,
	-0.034893769411408884,
	0.004819180069467605,
	-0.00046062616620627504,
	0.00003246032882100508
};

double J_0_T(double x) {
	double z = x / 8.0;
	double T_previous = 1.0, T_current = 2 * z * z - 1;
	double T;
	double s = c_T_0[0] * T_previous + c_T_0[1] * T_current;
	for (int n = 2; n < 9; n++) {
		T = (4.0 * z * z - 2) * T_current - T_previous;
		s += c_T_0[n] * T;
		T_previous = T_current; T_current = T;
	};
	return s;
}

const double c_T_1[9] = {
	 0.05245819033465648458,
	 0.04809646915823037394,
	 0.31327508236156718380,
	-0.24186740844740748475,
	 0.07426679621678703781,
	-0.01296762731173517510,
	 0.00148991289666763839,
	-0.00012227868505432427,
	 0.00000756263022969605
};

double J_1_T(double x) {
	double z = x / 8.0;
	double T_previous = z, T_current = 4 * z * z * z - 3 * z;
	double T;
	double s = c_T_1[0] * T_previous + c_T_1[1] * T_current;
	for (int n = 2; n < 9; n++) {
		T = (4.0 * z * z - 2) * T_current - T_previous;
		s += c_T_1[n] * T;
		T_previous = T_current; T_current = T;
	};
	return s;
}