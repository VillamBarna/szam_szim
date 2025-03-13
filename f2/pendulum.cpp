#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <Python.h>

using namespace std;

#include "vector.hpp"        // vectors with components of type double
#include "odeint.hpp"        // ODE integration routines, Runge-Kutta ...
using namespace cpl;

const double pi = 4 * atan(1.0);

const double g = 9.8;        // acceleration of gravity

double L = 1.0;              // length of pendulum
double L2 = 1.0; 
double q = 0.5;              // damping coefficient
double Omega_D = 2.0/3.0;    // frequency of driving force
double F_D = 0.9;            // amplitude of driving force
bool nonlinear;              // linear if false
char* algorithm;
double theta, omega, tMax;
double theta2, omega2;

void pendulum();
void double_pendulum();
double flipover();

static PyObject* py_simulation_pendulum(PyObject* self, PyObject* args) {
    if (!PyArg_ParseTuple(args, "spddddddd", &algorithm, &nonlinear, &L, &q, &Omega_D, &F_D, &theta, &omega, &tMax)) {
        return NULL;
    }
    pendulum();
    Py_RETURN_NONE;
}

static PyObject* py_simulation_double_pendulum(PyObject* self, PyObject* args) {
    if (!PyArg_ParseTuple(args, "spddddddddd", &algorithm,&nonlinear,  &L, &q, &Omega_D, &F_D, &theta, &theta2, &omega, &omega2, &tMax)) {
        return NULL;
    }
    double_pendulum();
    Py_RETURN_NONE;
}

static PyObject* py_simulation_flipover(PyObject* self, PyObject* args) {
    if (!PyArg_ParseTuple(args, "dddddd", &L, &theta, &omega, &theta2, &omega2, &tMax)) {
        return NULL;
    }
    double result = flipover();
    return PyFloat_FromDouble(result);
}

static PyMethodDef SimeMethods[] = {
    {"run_pendulum", (PyCFunction) py_simulation_pendulum, METH_VARARGS, "Run simulation"},
    {"run_double_pendulum", (PyCFunction) py_simulation_double_pendulum, METH_VARARGS, "Run simulation"},
    {"flipover", (PyCFunction) py_simulation_flipover, METH_VARARGS, "Run simulation"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef simmodule = {
    PyModuleDef_HEAD_INIT,
    "pendulum",
    NULL,
    -1,
    SimeMethods
};

PyMODINIT_FUNC 
PyInit_pendulum(void) {
    return PyModule_Create(&simmodule);
}

Vector f(const Vector& x) {  // extended derivative vector
    double t = x[0];
    double omega = x[2];
    double theta = x[1];
    Vector f(3);             // Vector with 3 components
    f[0] = 1;

    f[1] = omega;
    if (nonlinear)
        f[2] = - (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    else
        f[2] = - (g/L) * theta - q * omega + F_D * sin(Omega_D * t);
    return f;
}

Vector f_double(const Vector& x) {
    double theta = x[1];
    double omega = x[2];
    double theta2 = x[3];
    double omega2 = x[4];
    double beta = theta - theta2;
    Vector f(5);
    f[0] = 1;
    f[1] = omega;
    f[2] = (- g * 3 * sin(theta)
		  - g * sin(theta-2*theta2)
		  - 2 * omega2*omega2 * L * sin(beta)
		  - omega*omega * L * sin(2*beta))/
            L/(3-cos(2*beta));
    f[3] = omega2;
    f[4] = 2*sin(beta)*(2*omega*omega * L
		  + g * 2 * cos(theta)
		  + omega2*omega2*L*cos(beta))/
            L/(3-cos(2*beta));
    return f;
}

void Euler(Vector& x, double dt, Vector derivs(const Vector&)) {
	x += dt * derivs(x);
}

void EulerCromer(Vector& x, double dt, Vector derivs(const Vector&)) {
    Vector d=derivs(x);
    d[2] += d[3]*dt;
    x += d*dt;
}

void pendulum() {
    double dt = 0.01;
    double accuracy = 1e-6;
    ofstream dataFile("pendulum.data");

    double t = 0;
    Vector x(3);
    x[0] = t;
    x[1] = theta;
    x[2] = omega;

    while (t < tMax) {
        if (strcmp(algorithm, "euler") == 0) {
            Euler(x, dt, f);
        }
        else if (strcmp(algorithm, "euler_cromer") == 0) {
            EulerCromer(x, dt, f);
        }
        else if (strcmp(algorithm, "RK45_adaptive") == 0) {
            adaptiveRK4Step(x, dt, accuracy, f);
        }
        else if (strcmp(algorithm, "RK45") == 0) {
            RK4Step(x, dt, f);
        }
        else {
            cout << "Available algorithms: euler, euler_cromer, rk45, rk45_adaptive" << endl;
            return;
        }
        t = x[0], theta = x[1], omega = x[2];
        if (nonlinear) {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
        }
        dataFile << t << '\t' << theta << '\t' << omega << '\n';
    }

    cout << " Output data to file pendulum.data" << endl;
    dataFile.close();
}

void double_pendulum() {
    double dt = 0.05;
    double accuracy = 1e-6;
    ofstream dataFile("pendulum.data");

    double t = 0;
    Vector x(5);
    x[0] = t;
    x[1] = theta;
    x[2] = omega;
    x[3] = theta2;
    x[4] = omega2;

    while (t < tMax) {
	   adaptiveRK4Step(x, dt, accuracy, f_double);
        t = x[0], theta = x[1], omega = x[2], theta2 = x[3], omega2 = x[4];
        if (nonlinear) {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
            /* while (theta2 >= pi) theta2 -= 2 * pi; */
            /* while (theta2 < -pi) theta2 += 2 * pi; */
        }
        dataFile << t << '\t' << theta << '\t' << omega << '\t' << theta2 << '\t' << omega2 << '\n';
    }

    cout << " Output data to file pendulum.data" << endl;
    dataFile.close();
}

double flipover() {
    double dt = 0.05;
    double accuracy = 1e-6;

    double t = 0;
    Vector x(5);
    x[0] = t;
    x[1] = theta;
    x[2] = omega;
    x[3] = theta2;
    x[4] = omega2;

    while (x[0] < tMax) {
	   adaptiveRK4Step(x, dt, accuracy, f_double);
	   if ((abs(x[1]) >= 2*pi) || (abs(x[3]) >= 2*pi)) {
		   return x[0];
	   }
    }
    return -1;
}

