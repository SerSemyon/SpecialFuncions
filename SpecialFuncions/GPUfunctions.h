#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

/// <summary>
/// Вычисление функции Бесселя на видеокарте NVidia
/// </summary>
/// <param name="v"> порядок функции </param>
/// <param name="x"> значения параметра </param>
/// <param name="result"> полученные значения </param>
/// <param name="size"> количество точек </param>
void BesselWithCuda(const double v, const double* const x, double* result, const unsigned int size);
//cudaError_t BesselWithCudaNew(const double v, const double* x, double* result, const unsigned int size);

/// <summary>
/// Вычисление функции Бесселя нулевого порядка на отрезке [-8;8] 
/// </summary>
/// <param name="x"> значения параметра </param>
/// <param name="result"> полученные значения </param>
/// <param name="size"> количество точек </param>
void J0_CUDA(const double* const x, double* result, const unsigned int size);

/// <summary>
/// Вычисление функции Бесселя нулевого порядка на отрезке [-8;8] 
/// </summary>
/// <param name="x"> значения параметра </param>
/// <param name="result"> полученные значения </param>
/// <param name="size"> количество точек </param>
void J1_CUDA(const double* const x, double* result, const unsigned int size);

void Neumann_CUDA(double v, const double* const x, double* result, const unsigned int size, const double* const Jpositive, const double* const Jnegative);

/// <summary>
/// Вычисление функции Неймана нулевого порядка на (0;8]
/// </summary>
/// <param name="x"> значения параметра </param>
/// <param name="result"> значения параметра </param>
/// <param name="size"> количество точек </param>
/// <param name="J0"> значения функции Бесселя нулевого порядка </param>
void Y0_CUDA(const double* const x, double* result, const unsigned int size, const double* const J0);

/// <summary>
/// Вычисление функции Неймана нулевого порядка на (0;8]
/// </summary>
/// <param name="x"> значения параметра </param>
/// <param name="result"> значения параметра </param>
/// <param name="size"> количество точек </param>
/// <param name="J1"> значения функции Бесселя первого порядка </param>
void Y1_CUDA(const double* const x, double* result, const unsigned int size, const double* const J1);

void H1_CUDA(const double v, const double* const x, double* Re, double* Im, const unsigned int size);

void dZ_CUDA(double* x, double* result, unsigned int size, double* Z_vPrev, double* Z_vNext);

void cyl_next_order_CUDA(double v, double* x, double* result, unsigned int size, double* value_v, double* value_v_minus_1);

void BesselWithCuda(const double v, const double* const x, double* result, const unsigned int size_x, const unsigned int size_v);

__global__ void cyl_next_order_OneThread(double v, double* x, double* result, unsigned int size, double* value_v, double* value_v_minus_1);

void J_0_T_CUDA(const double* const x, double* result, const unsigned int size);
void J_1_T_CUDA(const double* const x, double* result, const unsigned int size);