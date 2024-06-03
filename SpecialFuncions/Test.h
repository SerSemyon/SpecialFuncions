#pragma once

unsigned int FindExecutionTime(void method());
void TestBesselCPU();
//void TestBesselNew();
void TestNeumannCPU();
void TestJ0(); 
void TestJ1();
void TestY0();
void TestY1();
void TestBessel_CUDA();
void TestJ0_CUDA();
void TestJ1_CUDA();
void TestNeumann_CUDA();
void TestY0_CUDA();
void TestY1_CUDA();
void TestChebyshevPolynomials();
void TestJnew();
void TestBesselOrderedSet();
void TestZ_vNext();
void TestJ_asymptotic();
void TestY_asymptotic();
void Test_dZ();
void TestCyl_next_order_CUDA();
void Testrec_CUDA();
void TestJ_0_T();
void TestJ_1_T();