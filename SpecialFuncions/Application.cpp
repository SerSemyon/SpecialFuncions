
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "CPUfunctions.h"
#include "GPUfunctions.h"
#include "Test.h"
#include "log_duration.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <complex>
#include "CylindricalBody.h"

int main(void)
{
    GLFWwindow* window;

    if (!glfwInit())
        return -1;

    window = glfwCreateWindow(1000, 1000, "Model", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    size_t n = 500;
    size_t angles_count = 2000;
    double* r = new double[n];
    double* externalR = new double[n];
    
    double k_1 = 1;
    double k_2 = 0.5;
    double radiusCircle = 1;
    
    double stepRadius = radiusCircle / n;
    for (size_t i = 0; i < n; i++) {
        r[i] = stepRadius * i;
        externalR[i] = radiusCircle + r[i];
    }
    //std::complex<double>** res = E_2(1, 0.1, r, n, angles_count, 0.1);
    // 
    //std::complex<double>** res = E_2(1, 0.5, r, n, angles_count, 10);
    //std::complex<double>** res = E_2(1, 0.5, r, n, angles_count, 10);
    //std::complex<double>** res = E_2(0.5, 0.5, r, n, angles_count, 0.1);
    // 
    //std::complex<double>** E_2_values = E_2(0.1, 1, r, n, angles_count, 0.1);
    //std::complex<double>** res = E_2(0.1, 1, r, n, angles_count, 10);
    //std::complex<double>** res = E_2(0.1, 1, r, n, angles_count, 1);
    //std::complex<double>** res = E_2(1, 1, r, n, angles_count, 1);
    
    CylindricalBody cb;
    std::complex<double>** E_1_values = cb.E_1(k_1, k_2, externalR, n, angles_count, radiusCircle);
    std::complex<double>** E_2_values = cb.E_2(k_1, k_2, r, n, angles_count, radiusCircle);
    
    double max = 0;
    for (int i = 0; i < n; i++) {
        if (abs(E_1_values[i][0]) > max) {
            max = abs(E_1_values[i][0]);
        }
        if (abs(E_2_values[i][0]) > max) {
            max = abs(E_2_values[i][0]);
        }
    }

    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT);
        
        glBegin(GL_POINTS);
        for (int i = 0; i < n; i++) {
            //glUniform4f(uColorLocation, 0.0f, 0.0f, 1.0f, 1.0f);
            for (int j = 0; j < angles_count; j++) {
                //glColor3f(0.0, res[i][j].real() / max, res[i][j].imag() / max);
                glColor3f(abs(E_2_values[i][j]) / max, E_2_values[i][j].real() / max, E_2_values[i][j].imag() / max);
                //glColor3f(abs(res[i][j]), res[i][j].real(), res[i][j].imag());
                double angle = 2 * j * M_PI / angles_count;
                double radius = r[i] / radiusCircle / 2;
                glVertex2f(radius * cos(angle), radius * sin(angle));
        
                glColor3f(abs(E_1_values[i][j]) / max, E_1_values[i][j].real() / max, E_1_values[i][j].imag() / max);
                radius = externalR[i] / radiusCircle / 2;
                glVertex2f(radius * cos(angle), radius * sin(angle));
            }
        }
        glEnd;
        
        glfwSwapBuffers(window);
        
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
