// rungekutta.h
#pragma once

#include "helpers.h"
// I/O and fmt headers for formatted output.
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <fmt/core.h>
#include <fmt/printf.h>
//  std library headers
#include <vector>
#include <complex>
#include <string>
// Eigen headers
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>


class RungeKutta {
public:
    RungeKutta();
    static void runPropagator(
        RowMatrixXd A,
        Eigen::VectorXcd y0,
        std::tuple<double, double, double, int> time_params,
        MatrixList ops_list,
        std::string headers,
        std::string outfilename);
};
