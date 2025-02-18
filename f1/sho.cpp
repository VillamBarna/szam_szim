#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

double omega;          // the natural frequency
double x, v;           // position and velocity at time t
int periods;           // number of periods to integrate
int stepsPerPeriod;    // number of time steps dt per period
string fileName;       // name of output file
bool eulerCromer;

void Euler(double dt);
void EulerCromer(double dt);     // takes an Euler-Cromer step
double energy();                 // computes the energy
void simulation();

static PyObject* py_simulation(PyObject* self, PyObject* args) {
    if (!PyArg_ParseTuple(args, "dddiisp", &omega, &x, &v, &periods, &stepsPerPeriod, &fileName, &eulerCromer)) {
        return NULL;
    }
    simulation();
    Py_RETURN_NONE;
}

static PyMethodDef SimeMethods[] = {
    {"run_simulation", (PyCFunction) py_simulation, METH_VARARGS, "Run simulation"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef simmodule = {
    PyModuleDef_HEAD_INIT,
    "sho",
    NULL,
    -1,
    SimeMethods
};

PyMODINIT_FUNC 
PyInit_sho(void) {
    return PyModule_Create(&simmodule);
}

void EulerCromer (double dt) {
    double a = - omega * omega * x;
    v += a * dt;
    x += v * dt;
}

void Euler (double dt) {
    double a = - omega * omega *x;
    x += v * dt;
    v += a * dt;
}

double energy ( ) {
    return 0.5 * (v * v + omega * omega * x * x);
}

void simulation ( ) {
    ofstream file(fileName.c_str());
    if (!file) {
        cerr << "Cannot open " << fileName << "\nExiting ...\n";
        return ;
    }
    const double pi = 4 * atan(1.0);
    double T = 2 * pi / omega;
    double dt = T / stepsPerPeriod;
    double t = 0;
    file << t << '\t' << x << '\t' << v << '\n';
    for (int p = 1; p <= periods; p++) {
        for (int s = 0; s < stepsPerPeriod; s++) {
            if (eulerCromer) {
                EulerCromer(dt);
            }
            else {
                Euler(dt);
            }
            t += dt;
            file << t << '\t' << x << '\t' << v << '\n';
        }
        cout << "Period = " << p << "\tt = " << t
             << "\tx = " << x << "\tv = " << v
             << "\tenergy = " << energy() << endl;
    }
    file.close();
}

