// rungekutta.cc
#include "rungekutta.h"

RungeKutta::RungeKutta() {} 
void RungeKutta::runPropagator(
    RowMatrixXd A, 
    Eigen::VectorXcd y0,
    std::tuple<double, double, double, int> time_params, 
    MatrixList ops_list,
    std::string headers,
    std::string outfilename) {
    // Reading time and field params    
    double ti, tf, dt;
    int print_nstep;
    double fs_to_au = 41.341374575751;
    std::tie (ti, tf, dt, print_nstep) = time_params;
    // declaring RK4 vetors
    Eigen::VectorXcd yi, k1, k2, k3, k4;
    Eigen::VectorXd expt_vals;
    int nops=ops_list.size(); 
    // writing headers to output file.
    FILE * fout = fopen(outfilename.c_str(), "w"); // output file object fout
    fmt::fprintf(fout, "%20s\t", "time");
    fmt::fprintf(fout, "%20s\t", "norm");
    fmt::fprintf(fout, "%20s\t", "autocorr");
    fmt::fprintf(fout, "%20s\t", headers);
    fmt::fprintf(fout, "\n");
    // writing expectation values at initial time ti.
    expt_vals = calc_expt(y0, nops, ops_list, y0);
    fmt::fprintf(fout, "%20.16f\t", ti/fs_to_au);
    for (int i = 0; i < 2 + nops; i++) {
        fmt::fprintf(fout,"%20.16f\t", expt_vals[i]);
    };
    fmt::fprintf(fout, "\n");
    // starting the propagation loop 
    while (ti <= tf) {
        for (int i = 0; i <= print_nstep; i++) {
            // Adding relevant dipoles to the hamiltonian
            // Runge Kutta 4 method
            k1.noalias() = A*yi;
            k2.noalias() = A*(yi + (dt*0.5)*k1);
            k3.noalias() = A*(yi + (dt*0.5)*k2);
            k4.noalias() = A*(yi + (dt*1.0)*k3);
            yi.noalias() += (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            ti += dt;
	        if (ti >= tf) {
		        std::cout<<"breaking loop"<<std::endl;
		        break;
	        };
        };
        // writing expectation values at ti
        fmt::fprintf(fout, "%20.16f\t", ti/fs_to_au);
        for (int i = 0; i < 2 + nops; i++) {
            fmt::fprintf(fout, "%20.16f\t", expt_vals[i]);
        };
        fmt::fprintf(fout,"\n", "");
    };
    fclose(fout);
    return;
};
