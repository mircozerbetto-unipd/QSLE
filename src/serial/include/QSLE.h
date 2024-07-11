/***********************************************************************************
 * QSLE v1.0 - Program to simulate the Quantum-Stochastic Liouville Equation       *
 * Copyright (C) 2024  Riccardo Cortivo, Mirco Zerbetto                            *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or any later version.                                           *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 ***********************************************************************************
 * Authors: Riccardo Cortivo, Mirco Zerbetto                                       *
 * Dipartimento di Scienze Chimiche - Universita' di Padova - Italy                *
 * E-mail: mirco.zerbetto@unipd.it                                                 *
 ***********************************************************************************/

#ifndef QSLE_H
#define QSLE_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include "spline.h"
#include "armadillo"
#include "omp.h"
#include "DiTe2.h"

using namespace std;
using namespace arma;

class QSLE {
public:
    explicit QSLE(const string&);
    virtual ~QSLE();
    void recover_data();
    void set_grid();
    void set_potential_vectors();
    void set_nonadiabatic_coeff();
    void set_g0();
    void set_friction();
    //void compute_parabolic_decoherence_time();
    void compute_decoherence_time();
    //double parabolic_fitting();
    void compute_transition_coefficients();
    void find_c();
    void compute_QSLE_matrix();
    void compute_LRBF_matrix();
    void compute_FD_matrix();
    void set_initial_condition();
    void time_evolution();
    
    void print_m();
    void print_g0();
    //void print_forces();
    //void print_testM();
    //void print_dg0();
    void print_restart();
    bool restart_check();
    
    long long int duration_dite2, duration_m, duration_evo, duration_FP;

private:
    stringstream restart;
    string base_name, PES_file, pdb_file, d01_file, m01_file, m10_file, distribution, transition, dite_input, reff, C, viscosity, method;
    long long int angle_points, mom_points, number_of_points, normalization_step, restart_step, total_steps;
    int close_q, close_p;
    double csi, tau_dec, ts, mom_max, max_Ek, temperature, initial_time, dx, dy, start_angle, start_mom, normalization, cq, cp;
    vector<double> angle_grid, mom_grid, supp_PES, MS, GS, ES, FG, FE, deltaE, deltaE_01, supp_d01, d01, supp_g0, g0, dg0, m01, m10;
    //vector<double> testM, tempo;
    vector<int> dihedral;
    sp_mat GammaM;
    mat Phi_q, Phi_p;
    vec prob_dens, X, Y, integral_X1, integral_X2, integral_Y1, integral_Y2;
    //vec k1, k2, k3, k4,;

};


#endif // QSLE_H
