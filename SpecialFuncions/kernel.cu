#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "CPUfunctions.h"
#include "GPUfunctions.h"
#include "Test.h"
#include "log_duration.h"

void tests()
{
    TestBessel_CUDA();
    TestBessel_CUDA();
    TestJnew();
    TestBesselOrderedSet();
    TestJ0();
    TestJ0_CUDA();
    TestJ1();
    TestJ1_CUDA();
    ////TestBesselNew();
    TestNeumannCPU();
    TestY0();
    TestY1();
    TestY0_CUDA();
    TestY1_CUDA();
    TestZ_vNext();
    TestCyl_next_order_CUDA();
    TestChebyshevPolynomials();
    Test_dZ();
    TestJ_asymptotic();
    TestY_asymptotic();
    TestNeumann_CUDA();

    Testrec_CUDA();
}