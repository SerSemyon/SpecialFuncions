
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

void write_abs_to_csv(std::string name, std::complex<double>** values, double r, size_t angles_count, size_t n) {
    std::ofstream out;
    out.open(name + ".csv");
    if (out.is_open())
    {
        out << ", ";
        for (size_t i = 0; i < angles_count; i++)
        {
            out << round( 2 * i * M_PI * 100 / angles_count) / 100 << ',';
        }
        out << '\n';
        for (size_t i = 0; i < n; i++) {
            out << i << ',';
            for (size_t j = 0; j < angles_count - 1; j++) {
                out << round(std::abs(values[i][j]) * 100) / 100 << ',';
            }
            out << round(std::abs(values[i][angles_count - 1]) * 100) / 100 << '\n';
        }
    }
    out.close();
}

void write_to_csv(std::string name, std::complex<double>** values, double r, size_t angles_count, size_t n) {
    std::ofstream out;
    out.open(name + ".csv");
    if (out.is_open())
    {
        out << ", ";
        for (size_t i = 0; i < angles_count; i++)
        {
            out << round(2 * i * M_PI * 100 / angles_count) / 100 << ", , ";
        }
        out << '\n';
        for (size_t i = 0; i < n; i++) {
            out << i << ',';
            for (size_t j = 0; j < angles_count - 1; j++) {
                out << values[i][j] << ',';
            }
            out << values[i][angles_count - 1] << '\n';
        }
    }
    out.close();
}

int ShowWindowWithDiagram() {
    GLFWwindow* window;

    if (!glfwInit())
        return -1;

    window = glfwCreateWindow(1500, 750, "Model cylindrical body", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    // число делений по радиусу
    size_t n = 500;
    //size_t n = 10;
    // число делений по углу
    size_t angles_count = 2000;
    //size_t angles_count = 8;
    // значения функций внутри цилиндра на отрезке [0,R)
    double* r = new double[n];
    // значения функций внутри цилиндра на отрезке [R, 2R)
    double* externalR = new double[n];

    // Параметры цилиндра и среды
    double a_0 = 1;
    double k_1 = 1;
    double k_2 = 3;
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
    double min = 1000;
    for (int i = 0; i < n; i++) {
        if (std::abs(E_1_values[i][0]) > max) {
            max = std::abs(E_1_values[i][0]);
        }
        if (std::abs(E_2_values[i][0]) > max) {
            max = std::abs(E_2_values[i][0]);
        }
    }
    std::cout << max << std::endl;

    /*write_to_csv("E1", E_1_values, radiusCircle, angles_count, n);
    write_to_csv("E2", E_2_values, radiusCircle, angles_count, n);

    write_abs_to_csv("E1_abs", E_1_values, radiusCircle, angles_count, n);
    write_abs_to_csv("E2_abs", E_2_values, radiusCircle, angles_count, n);*/

    glClearColor(1, 1, 1, 1);
    //glClearColor(0, 0, 0, 0);

    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT);

        glBegin(GL_POINTS);
        for (int i = 0; i < n; i++) {
            //glUniform4f(uColorLocation, 0.0f, 0.0f, 1.0f, 1.0f);
            for (int j = 0; j < angles_count; j++) {
                //glColor3f(0.0, res[i][j].real() / max, res[i][j].imag() / max);

                // модуль
                glColor3d(std::abs(E_2_values[i][j]) / max, std::abs(E_2_values[i][j]) / max, std::abs(E_2_values[i][j]) / max);

                double angle = 2 * j * M_PI / angles_count;
                double radius = r[i] / radiusCircle / 2;
                glVertex2d(radius / 3 * cos(angle) - 0.5, radius * sin(angle) / 1.5);


                // Только положительные
                //glColor3d(0, E_2_values[i][j].real() / max, E_2_values[i][j].imag() / max);
                // Абсолютные значения вещественной и мнимой частей
                glColor3d(0, std::abs(E_2_values[i][j].real()) / max, std::abs(E_2_values[i][j].imag()) / max);
                glVertex2d(radius / 3 * cos(angle) + 0.5, radius * sin(angle) / 1.5);

                glColor3d(std::abs(E_1_values[i][j]) / max, std::abs(E_1_values[i][j]) / max, std::abs(E_1_values[i][j]) / max);
                radius = externalR[i] / radiusCircle / 2;
                glVertex2d(radius / 3 * cos(angle) - 0.5, radius * sin(angle) / 1.5);

                //glColor3d(0, E_1_values[i][j].real() / max, E_1_values[i][j].imag() / max);
                glColor3d(0, std::abs(E_1_values[i][j].real()) / max, std::abs(E_1_values[i][j].imag()) / max);
                glVertex2d(radius / 3 * cos(angle) + 0.5, radius * sin(angle) / 1.5);
            }
        }
        glEnd();

        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

int main(void)
{
    //Measure_J0_Time();
    //Measure_J1_Time();
    //TestY0();
    //Measure_Y0_Time();
    //TestJ0();
    ShowWindowWithDiagram();
    /*TestNeumann_one_point();
    TestNeumannCPU();
    TestJ1_CUDA();*/
    //TestNeumann_CUDA();
    /*TestJ_0_T();
    TestJ_1_T(); 
    TestJ_0_T_CUDA();
    TestJ_1_T_CUDA();*/
    return 0;
}
