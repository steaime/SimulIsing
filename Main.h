#pragma once
#include "stdafx.h"
#include "Simulation.h"
#include "IsingParameters.h"
#include "Constants.h"
#include "InitialConditions.h"
#include "TestFunctions.h"
#include "ConfigReader.h"

// Global variable: list of simulation objects
// Multiple simulations can be run in parallel, e.g. to improve statistics
Simulation* simList = NULL;

void InitializeSimulations(const ConfigParams &config_reader);

void SetupSimulation(Simulation* my_sim, bool print_output, int num_sites, double Tau0, double omega, double K, double alpha, 
				double beta, double n, bool UsePBC, double temp = DEF_SIM_TEMP, bool force_meanfield = false);
void SetupSimulation(Simulation* my_sim, Noise AlphaNoise, Noise InitConditions, bool print_output, int num_sites, double Tau0,
				double omega, double K, double beta, double n, bool UsePBC, double temp, bool force_meanfield = false);
void SetupSimulation(Simulation* my_sim, const IsingParameters &Parameters, const Noise &InitConditions, const MCParams &tstep_params, 
				bool force_meanfield = false, bool print_output = false);
void SetupSimulation(Simulation* my_sim, std::string sParamFile);
void SetupSimulation(Simulation* my_sim, const ConfigParams &params);

void RunStrainSweep(int sim_number, double start_gamma, double end_gamma, int ppd_gamma, bool preshear, bool ascdesc,
	bool print_std, bool print_nloops, bool print_fastsites, std::string out_file_name, std::string* out_file_stats = NULL,
	std::string* out_file_histo = NULL, std::string* out_file_eq = NULL, std::string* out_file_corr = NULL, std::string* out_raw = NULL, 
	std::string* out_raw_final = NULL, std::string out_file_histalpha = "", bool mute_out=false);
void RunStrainSweep(int sim_number, double* gamma_list, int num_gamma, bool preshear, bool ascdesc, bool print_std, bool print_nloops,
	bool print_fastsites, std::string out_file_name, std::string* out_file_stats = NULL,
	std::string* out_file_histo = NULL, std::string* out_file_eq = NULL, std::string* out_file_corr = NULL, std::string* out_raw = NULL, 
	std::string* out_raw_final = NULL, std::string out_file_histalpha = "", bool mute_out = false);
void StrainSweepCore(int sim_number, double* gamma_list, int gamma_num, bool preshear, bool ascdesc, bool print_std, bool print_nloops, 
	bool print_fastsites, std::string out_file_name, double **data, double*** all_data, double**** dyn_data, bool**** bln_fast_data, bool mute_out = false);

void print_sim_params(std::string out_file_name, int sim_number);


void SweepAndSave(const ConfigParams &config_reader);
double SweepAndCompare(const ConfigParams &config_reader, std::string sRefData);

void FreeMem(const ConfigParams &config_reader);

int main(int argc, char* argv[]);