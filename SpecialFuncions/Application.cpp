
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

#include "Test.h"

int ShowWindowWithDiagram() {
    GLFWwindow* window;

    if (!glfwInit())
        return -1;

    window = glfwCreateWindow(1000, 1000, "Model cylindrical body", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    // число делений по радиусу
    size_t n = 500;
    // число делений по углу
    size_t angles_count = 4000;
    // значения функций внутри цилиндра на отрезке [0,R)
    double* r = new double[n];
    // значения функций внутри цилиндра на отрезке [R, 2R)
    double* externalR = new double[n];

    // Параметры цилиндра и среды
    double a_0 = 1;
    double k_1 = 1;
    double k_2 = 0.5;
    double radiusCircle = 1;

    double stepRadius = radiusCircle / n;
    for (size_t i = 0; i < n; i++) {
        r[i] = stepRadius * i;
        externalR[i] = radiusCircle + r[i];
    }

    // Создание объекта цилиндра с параметрами
    CylindricalBody cb(k_1, k_2, radiusCircle, a_0);
    std::complex<double>** E_1_values = cb.E_1(externalR, n, angles_count);
    std::complex<double>** E_2_values = cb.E_2(r, n, angles_count);

    // Нахождение максимального значения модуля для нормирования значений по нему при отображении
    double max = 0;
    for (int i = 0; i < n; i++) {
        if (abs(E_1_values[i][0]) > max) {
            max = abs(E_1_values[i][0]);
        }
        if (abs(E_2_values[i][0]) > max) {
            max = abs(E_2_values[i][0]);
        }
    }

    glClearColor(1, 1, 1, 1);

    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT);

        glBegin(GL_POINTS);
        for (int i = 0; i < n; i++) {
            //glUniform4f(uColorLocation, 0.0f, 0.0f, 1.0f, 1.0f);
            for (int j = 0; j < angles_count; j++) {
                //glColor3f(0.0, res[i][j].real() / max, res[i][j].imag() / max);

                // красный - модуль, зелёный - вещественная, синий - мнимая
                glColor3d(abs(E_2_values[i][j]) / max, E_2_values[i][j].real() / max, E_2_values[i][j].imag() / max);
                //glColor3f(abs(res[i][j]), res[i][j].real(), res[i][j].imag());
                double angle = 2 * j * M_PI / angles_count;
                double radius = r[i] / radiusCircle / 2;
                glVertex2d(radius * cos(angle), radius * sin(angle));

                glColor3d(abs(E_1_values[i][j]) / max, E_1_values[i][j].real() / max, E_1_values[i][j].imag() / max);
                radius = externalR[i] / radiusCircle / 2;
                glVertex2d(radius * cos(angle), radius * sin(angle));
            }
        }
        glEnd;

        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

int main(void)
{
    //ShowWindowWithDiagram();
    TestJ_0_T();
    TestJ_1_T();
    return 0;
}
