#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include "../integrators/helpers.h"
#include "../integrators/rungekutta.h"

namespace py = pybind11;
using namespace pybind11::literals;

// Export the integrators

void export_integrators(py::module &m)
{
    py::class_<RungeKutta>(m, "RungeKutta")
        .def(py::init<>())
        .def_static("runPropagator", &RungeKutta::runPropagator);
    m.def("diagonalise", &diagonalise, "A function that diagonalises a matrix");
    m.def("calc_expt", &calc_expt, "");
    m.def("cwiseExpcd", &cwiseExpcd, "");
}
