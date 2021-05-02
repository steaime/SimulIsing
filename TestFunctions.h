#pragma once

#include "stdafx.h"
#include "Simulation.h"
#include "IsingParameters.h"
#include "Constants.h"
#include "InitialConditions.h"

#ifndef  TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H

int check_indexing(IsingParameters my_params, int test_site);
void check_neighbors(Simulation* my_sim, int* sites, int num_sites);
void check_noise(Noise *my_noise);
void check_meanfield_gamma(IsingParameters my_params, double start_rate = 1e-5, double end_rate = 1, double ppd_rate = 10);
void export_meanfield_gamma(IsingParameters my_params, double* rates, int num_rates, std::string out_file);
void check_meanfield_rate(Simulation* my_sim, double start_gamma = 1e-3, double end_gamma = 1, double ppd_gamma = 10);
void check_norm_mfparams(Simulation* my_sim, std::ofstream* fout);
void check_real_roots();
bool test_roots_P3(double R1, double R2, double R3, double Pref=1.0);

#endif // ! TEST_FUNCTIONS_H
