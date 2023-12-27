#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h> 
#include <pybind11/chrono.h>
#include "helpers.h"
#include "rungekutta.h"

namespace py = pybind11;

PYBIND11_MODULE(integrators, m) {
    m.doc() = "A pybind11 plugin for TDSE integrators";
    py::class_<RungeKutta>(m, "RungeKutta")
        .def(py::init<>())
        .def_static("runPropagator", &RungeKutta::runPropagator);
    m.def("diagonalise", &diagonalise, "A function that diagonalises a matrix");
    m.def("calc_expt", &calc_expt, ""); 
    m.def("cwiseExpcd", &cwiseExpcd, "");
    m.def("testmul", &testmul, "");
}