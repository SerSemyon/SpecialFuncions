#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#include "CPUfunctions.h"
#include "GPUfunctions.h"
#include <iostream>
#include "log_duration.h"

__global__ void dZ_OneThread(double* x, double* result, unsigned int size, double* Z_vPrev, double* Z_vNext)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < size)
    {
        result[i] = 0.5 * (Z_vPrev[i] - Z_vNext[i]);
        i += blockDim.x * gridDim.x;
    }
}

void dZ_CUDA(double* x, double* result, unsigned int size, double* Z_vPrev, double* Z_vNext)
{
    double* dev_x = 0;
    double* dev_res = 0;
    double* dev_Z_vPrev = 0;
    double* dev_Z_vNext = 0;

    cudaMalloc((void**)&dev_x, size * sizeof(double));
    cudaMemcpy(dev_x, x, size * sizeof(double), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&dev_res, size * sizeof(double));

    cudaMalloc((void**)&dev_Z_vPrev, size * sizeof(double));
    cudaMemcpy(dev_Z_vPrev, Z_vPrev, size * sizeof(double), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&dev_Z_vNext, size * sizeof(double));
    cudaMemcpy(dev_Z_vNext, Z_vNext, size * sizeof(double), cudaMemcpyHostToDevice);

    {
        LOG_DURATION("GPU without data transfers");
        dZ_OneThread << <(size + 127) / 128, 128 >> > (dev_x, dev_res, size, dev_Z_vPrev, dev_Z_vNext);

        cudaGetLastError();
        cudaDeviceSynchronize();
    }

    cudaGetLastError();
    cudaDeviceSynchronize();

    cudaMemcpy(result, dev_res, size * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(dev_res);
    cudaFree(dev_x);
    cudaFree(dev_Z_vPrev);
    cudaFree(dev_Z_vNext);
}

__global__ void cyl_next_order_OneThread(double v, double* x, double* result, unsigned int size, double* value_v, double* value_v_minus_1)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < size)
    {
        result[i] = 2 * v * value_v[i] / x[i] - value_v_minus_1[i];
        i += blockDim.x * gridDim.x;
    }
}

void cyl_next_order_CUDA(double v, double* x, double* result, unsigned int size, double* value_v, double* value_v_minus_1)
{
    double* dev_x = 0;
    double* dev_res = 0;
    double* dev_value_v = 0;
    double* dev_value_v_minus_1 = 0;

    cudaMalloc((void**)&dev_x, size * sizeof(double));
    cudaMemcpy(dev_x, x, size * sizeof(double), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&dev_res, size * sizeof(double));

    cudaMalloc((void**)&dev_value_v, size * sizeof(double));
    cudaMemcpy(dev_value_v, value_v, size * sizeof(double), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&dev_value_v_minus_1, size * sizeof(double));
    cudaMemcpy(dev_value_v_minus_1, value_v_minus_1, size * sizeof(double), cudaMemcpyHostToDevice);

    {
        LOG_DURATION("GPU without data transfers");
        cyl_next_order_OneThread << <(size + 127) / 128, 128 >> > (v, dev_x, dev_res, size, dev_value_v, dev_value_v_minus_1);

        cudaGetLastError();
        cudaDeviceSynchronize();
    }

    cudaGetLastError();
    cudaDeviceSynchronize();

    cudaMemcpy(result, dev_res, size * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(dev_res);
    cudaFree(dev_x);
    cudaFree(dev_value_v);
    cudaFree(dev_value_v_minus_1);
}