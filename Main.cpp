#include "Main.h"

// uncomment this to debug SweepAndCompare function
//#define DEBUG_COMPARE "C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\config\\ExpLudox.txt"


int main(int argc, char* argv[]) {

	std::string current_exec_name = argv[0]; // Name of the current exec program
	std::string first_arg;					 // If specified, path of INI config file
	std::vector<std::string> all_args;
	if (argc > 1) {
		first_arg = argv[1];
		all_args.assign(argv + 1, argv + argc);
	}
	else {
		std::string exe_dir = GetExeDir();
		first_arg = JoinPath(exe_dir, DEF_INI_FILE);
	}

	int iCompare = std::distance(all_args.begin(), std::find(all_args.begin(), all_args.end(), "-COMPARE"));
	bool bSilentRun = (std::find(all_args.begin(), all_args.end(), "-SILENT") != all_args.end());
	bool bPauseAfter = (std::find(all_args.begin(), all_args.end(), "-PAUSE") != all_args.end());

	std::time_t rand_seed = time(NULL);
	srand(rand_seed);

#if DEBUG_MODE
	check_real_roots();
#endif

#ifndef DEBUG_COMPARE:
	if (!CheckFileExists(first_arg)) {
		std::cout << "ERROR: configuration file '" << first_arg << "' not found. Execution aborted" << std::endl;
		return 1;
	}
	ConfigParams config_reader(first_arg);
	config_reader.iRandSeed = rand_seed;
	CopyFileToFolder(first_arg, config_reader.out_folder);
	InitializeSimulations(config_reader);
	if (all_args.size() > 0 && iCompare < all_args.size()-1) {
		SweepAndCompare(config_reader, all_args[iCompare + 1]);
#else
	ConfigParams config_reader("C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\config\\simParams_0.ini");
	InitializeSimulations(config_reader);
	if (true){
		SweepAndCompare(config_reader, DEBUG_COMPARE);
#endif
	}
	else {
		SweepAndSave(config_reader);
	}
	FreeMem(config_reader);

	if (!bSilentRun) {
		std::cout << "...Program terminated!" << std::endl;
	}
#ifndef DEBUG_COMPARE:
	if (bPauseAfter) 
#endif
	{
		std::cin.ignore();
	}
	return 0;

}



double SweepAndCompare(const ConfigParams &config_reader, std::string sRefData) {

	// Load reference data
	int iNumGamma = FileCountLines(sRefData);
	int iNumCols = config_reader.iRefFileColNum;
	double **ppdExpData = new double*[iNumGamma];
	for (int i = 0; i < iNumGamma; i++) {
		ppdExpData[i] = new double[iNumCols];
	}
	FileLoadMultiColumn(sRefData, iNumCols, ppdExpData);
	double* pdGammaList = new double[iNumGamma];
	for (int i = 0; i < iNumGamma; i++) {
		pdGammaList[i] = ppdExpData[i][config_reader.iRefFileColIdxStrain];
	}

	// Initialize memory to store simulation results
	double** ppdAvgRates = new double*[config_reader.iNumSimulations];
	double*** pppdAllRates = new double**[config_reader.iNumSimulations];
	int iNumGammaRes = iNumGamma;
	if (config_reader.bAscDesc) {
		iNumGammaRes = 2 * iNumGamma;
	}
	for (int i = 0; i < config_reader.iNumSimulations; i++) {
		ppdAvgRates[i] = new double[iNumGammaRes];
		pppdAllRates[i] = new double*[iNumGammaRes];
		for (int j = 0; j < iNumGammaRes; j++) {
			pppdAllRates[i][j] = new double[::simList[0].GetNumSites()];
		}
	}

	// Run strain sweep, don't save output, store rates in memory
	StrainSweepCore(config_reader.iNumSimulations, pdGammaList, iNumGamma, config_reader.bPreshear, config_reader.bAscDesc, false,
		false, false, "", ppdAvgRates, pppdAllRates, NULL, NULL, true);
	
	// Evaluate parameters to compare simulation and reference
	double *pdSimGamma = new double[iNumGammaRes];
	double *pdSimChi = new double[iNumGammaRes];
	CalcColumnAverage(ppdAvgRates, config_reader.iNumSimulations, iNumGammaRes, pdSimGamma);
	for (int i = 0; i < iNumGamma; i++) {
		int cur_nFast = 0;
		for (int j = 0; j < config_reader.iNumSimulations; j++) {
			cur_nFast += CountGreaterThan(pppdAllRates[j][i], ::simList[j].GetNumSites(), config_reader.dFastSitesRelThr * 1.0 / config_reader.Tau0);
		}
		pdSimChi[i] = 1 - (double)cur_nFast / (1.0 * ::simList[0].GetNumSites() * config_reader.iNumSimulations);
	}

	// Evaluate difference between simulation and reference
	double dDiffChi = 0;
	double dDiffRate = 0;
	for (int i = 0; i < iNumGamma; i++) {
		dDiffChi += pow(pdSimChi[i] - ppdExpData[i][config_reader.iRefFileColIdxChi], 2) / iNumGamma;
		double cur_rate_diff;
		if (config_reader.bDiffRelRate) {
			cur_rate_diff = log(pdSimGamma[i]/ppdExpData[i][config_reader.iRefFileColIdxRate]);
		}
		else {
			cur_rate_diff = pdSimGamma[i] - ppdExpData[i][config_reader.iRefFileColIdxRate];
		}
		dDiffRate += pow(cur_rate_diff, 2) / iNumGamma;
#if DEBUG_MODE:
		std::cout << pdGammaList[i] << "\t" << pdSimGamma[i] << "\t" << ppdExpData[i][config_reader.iRefFileColIdxRate] << 
			"\t" << pow(cur_rate_diff, 2) << "\t" << pdSimChi[i] << "\t" << ppdExpData[i][config_reader.iRefFileColIdxChi] << 
			"\t" << pow(pdSimChi[i] - ppdExpData[i][config_reader.iRefFileColIdxChi], 2) << std::endl;
#endif
	}
	double dDiff = config_reader.dDiffWeightChi * dDiffChi + config_reader.dDiffWeightRate * dDiffRate;
	if (config_reader.bDiffCombine) {
		std::cout << dDiff << std::endl;
	}
	else {
		std::cout << dDiffRate << " " << dDiffChi << std::endl;
	}

	// Free memory
	for (int i = 0; i < config_reader.iNumSimulations; i++) {
		for (int j = 0; j < iNumGammaRes; j++) {
			delete[] pppdAllRates[i][j];
			pppdAllRates[i][j] = NULL;
		}
		delete[] ppdAvgRates[i];
		delete[] pppdAllRates[i];
		delete[] ppdExpData[i];
		ppdAvgRates[i] = NULL;
		pppdAllRates[i] = NULL;
		ppdExpData[i] = NULL;
	}
	delete[] ppdAvgRates;
	delete[] pdGammaList;
	delete[] pppdAllRates;
	delete[] ppdExpData;
	delete[] pdSimGamma;
	delete[] pdSimChi;
	ppdAvgRates = NULL;
	pdGammaList = NULL;
	pppdAllRates = NULL;
	ppdExpData = NULL;
	pdSimGamma = NULL;
	pdSimChi = NULL;

	return dDiff;
}





void SweepAndSave(const ConfigParams &config_reader) {

	// Initialize file names (eventually)
	std::string* hist_names = NULL;
	std::string* corr_names = NULL;
	std::string* outraw_names = NULL;
	std::string out_hist_alpha = "";
	std::string out_names = config_reader.out_folder + config_reader.fout_prefix;
	std::string stats_names = config_reader.out_folder + config_reader.fstat_prefix;
	std::string fname_prefix = "";
	if (config_reader.bSaveHist) {
		hist_names = new std::string[config_reader.iNumSimulations];
	}
	if (config_reader.bSaveCorr) {
		corr_names = new std::string[config_reader.iNumSimulations];
	}
	if (config_reader.bSaveRaw) {
		outraw_names = new std::string[config_reader.iNumSimulations];
	}
	if (config_reader.bSaveHistAlpha) {
		if (config_reader.out_folder.back() != '\\') {
			out_hist_alpha = config_reader.out_folder + "\\" + config_reader.histalpha_fname + config_reader.txt_ext;
		}
		else {
			out_hist_alpha = config_reader.out_folder + config_reader.histalpha_fname + config_reader.txt_ext;
		}
	}
	if (config_reader.bAscDesc) {
		fname_prefix = fname_prefix + config_reader.ascdesc_prefix;
	}
	if (config_reader.bPreshear) {
		fname_prefix = fname_prefix + config_reader.preshear_prefix;
	}
	out_names = out_names + fname_prefix + "all" + config_reader.txt_ext;
	stats_names = stats_names + fname_prefix + "all" + config_reader.txt_ext;

	for (int i = 0; i < config_reader.iNumSimulations; i++) {
		if (config_reader.bSaveHist) {
			if (config_reader.out_folder.back() != '\\') {
				hist_names[i] = config_reader.out_folder + "\\" + config_reader.hist_prefix + std::to_string(i) + config_reader.txt_ext;
			}
			else {
				hist_names[i] = config_reader.out_folder + config_reader.hist_prefix + std::to_string(i) + config_reader.txt_ext;
			}
		}
		if (config_reader.bSaveCorr) {
			if (config_reader.out_folder.back() != '\\') {
				corr_names[i] = config_reader.out_folder + "\\" + config_reader.corr_prefix + std::to_string(i) + config_reader.txt_ext;
			}
			else {
				corr_names[i] = config_reader.out_folder + config_reader.corr_prefix + std::to_string(i) + config_reader.txt_ext;
			}
		}
		if (config_reader.bSaveRaw) {
			outraw_names[i] = config_reader.out_folder;
			if (config_reader.out_folder.back() != '\\') {
				outraw_names[i]  = outraw_names[i] + "\\";
			}
			outraw_names[i] = outraw_names[i] + config_reader.rawfinal_prefix + fname_prefix + std::to_string(i) + config_reader.raw_ext;
		}
	}

	// Get list of gamma values
	int iNumGamma = config_reader.GetNumGamma();
	double* pdGammaList = new double[iNumGamma];
	config_reader.GetListGamma(pdGammaList);

	// Run simulations
	RunStrainSweep(config_reader.iNumSimulations, pdGammaList, iNumGamma, config_reader.bPreshear, config_reader.bAscDesc, config_reader.bPrintStd, config_reader.bPrintNloops,
		config_reader.bPrintFastSites, out_names, &stats_names, hist_names, NULL, corr_names, NULL, outraw_names, out_hist_alpha, config_reader.bMuteOut);

	// Free memory
	delete[] hist_names;
	delete[] corr_names;
	delete[] outraw_names;
	delete[] pdGammaList;
	hist_names = NULL;
	corr_names = NULL;
	outraw_names = NULL;
	pdGammaList = NULL;

}



void InitializeSimulations(const ConfigParams &config_reader) {

	// Initialize simulation list
	::simList = new Simulation[config_reader.iNumSimulations];

	// Initialize simulations
	for (int i = 0; i < config_reader.iNumSimulations; i++) {
		SetupSimulation(&::simList[i], config_reader);
	}

	// Print simulation parameters
	print_sim_params(config_reader.out_folder + config_reader.fparams, config_reader.iNumSimulations);

}

void FreeMem(const ConfigParams &config_reader) {
	for (int i = 0; i < config_reader.iNumSimulations; i++) {
		::simList[i].Clear();
	}
	delete[] ::simList;
	::simList = NULL;
}

void SetupSimulation(Simulation* simList, bool print_output, int num_sites, double Tau0, double omega, double K,
	double alpha, double beta, double n, bool UsePBC, double temp, bool force_meanfield) {
	Noise alpha_noise(NoiseTypes::NONOISE, alpha);
	Noise init_conditions(NoiseTypes::NONOISE, 1.0 / Tau0);
	SetupSimulation(simList, alpha_noise, init_conditions, print_output, num_sites, Tau0, omega, K, beta, n, UsePBC, temp, force_meanfield);
}

void SetupSimulation(Simulation* simList, Noise AlphaNoise, Noise InitConditions, bool print_output, int num_sites, double Tau0,
	double omega, double K, double beta, double n, bool UsePBC, double temp, bool force_meanfield) {
	int grid_size[N_DIM];
	for (int i = 0; i < N_DIM; i++) {
		grid_size[i] = num_sites;
	}
#if DEBUG_MODE:
	std::cout << "   ... initializing simulation parameters" << std::endl;
#endif
	IsingParameters my_params(grid_size, n, K, Tau0, omega, beta, AlphaNoise, UsePBC);
	if (print_output) {
		check_real_roots();
		check_indexing(my_params, my_params.NumSites() / 7);
		check_meanfield_gamma(my_params, 1e-4, 0.01, 10);
	}
#if DEBUG_MODE:
	std::cout << "   ... setting up simulation" << std::endl;
#endif
	MCParams tstep_params(temp);
	SetupSimulation(simList, my_params, InitConditions, tstep_params, force_meanfield, print_output);
#if DEBUG_MODE:
	std::cout << "       setup_simul() ended correctly!" << std::endl;
#endif
}

void SetupSimulation(Simulation* my_sim, std::string sParamFile) {
	ConfigParams params(sParamFile);
	SetupSimulation(my_sim, params);
}

void SetupSimulation(Simulation * my_sim, const ConfigParams & params)
{
	IsingParameters myParams(params);
#if 0
	double *_pdInitParams = new double[NOISE_MAXPARAMS];		// This pointer is destroyed during Noise initialization.
																// Not the cleanest approach...
	for (int i = 0; i < NOISE_MAXPARAMS; i++) {
		_pdInitParams[i] = params.pdInitParams[i];
	}
	Noise init_conditions(params.iInitType, _pdInitParams, params.dInitAvg);
#else
	Noise init_conditions(params.iInitType, params.pdInitParams, params.dInitAvg);
#endif
	MCParams tstep_params(params);
	SetupSimulation(my_sim, myParams, init_conditions, tstep_params, params.bForceMeanField);
	if (params.bSaveAlphaSnapshot) {
		int _iNumAlpha = my_sim->CountAlphas();
		double* _pdAlphas = new double[_iNumAlpha];
		my_sim->FlattenAlphas(_pdAlphas);
		std::string out_alpha_fname = params.out_folder + params.sAlphaSnapshotPrefix + std::to_string(my_sim->ID()) + params.raw_ext;
		ExportDataToRawBinary(_pdAlphas, &_iNumAlpha, 1, out_alpha_fname);
		delete[] _pdAlphas;
		_pdAlphas = NULL;
	}
	if (params.bSaveInitSnapshot) {
		// The current rates have just been initialized: they are a genuine initial condition
		int _iNumSites = my_sim->GetNumSites();
		std::string out_init_fname = params.out_folder + params.sInitSnapshotPrefix + std::to_string(my_sim->ID()) + params.raw_ext;
		ExportDataToRawBinary(my_sim->GetRates(), &_iNumSites, 1, out_init_fname);
	}
}

void SetupSimulation(Simulation* my_sim, const IsingParameters &Parameters, const Noise &InitConditions, 
	                 const MCParams &tstep_params, bool force_meanfield, bool print_output) {
	my_sim->Setup(Parameters, InitConditions, tstep_params);
#if DEBUG_MODE:
	std::cout << "my_sim->Setup() ended correctly!" << std::endl;
#endif
	if (print_output) {
		std::cout << "\n\nSimulation initialized. Average alpha=" << my_sim->GetAverageAlpha()
			<< "; average initial 1/tau0=" << my_sim->GetAverageRate() << std::endl;
		if (my_sim->GetNumSites() >= 3) {
			int sites[3] = { 0, 1, my_sim->GetNumSites() / 3 };
			check_neighbors(my_sim, sites, 3);
		}
		check_meanfield_rate(my_sim, 0.01, 1, 5);
	}
}



/*
gamma_list:		list of strains (in strain units) to impose
gamma_num:		number of elements in gamma_list
preshear:		if true: sample is reinitialized between each point
				if false: point n+1 uses the output of point n as initial condition
ascdesc:		if true: after the sweep gamma_list is reversed and the sweep is run again
				(only makes sense if preshear==false, otherwise it's trivial)
out_file_name:	file to save average relaxation rates for all simulations (eventually asc/desc),
				with error bars based on data dispersion
				Set to empty string ("") not to write output
avg_data:		pointer for exporting average rates.
				avg_data[i][j] will contain the average rate of i-th simulation at j-th gamma
				if ascdesc==true, the list of gamma will be actually 2*gamma_num long
				default value is NULL, in which case nothing is exported
all_data:		pointer for exporting individual rates.
				all_data[i][j][k] will be the rate of k-th site in i-th simulation at j-th gamma
				default value is NULL, in which case nothing is exported
dyn_data:		pointer for exporting equilibration data.
				dyn_data[i][j][k][l] will be quantity (e.g. [avRate, avChange, deviation])
				relative to the k-th equilbration step of j-th gamma in i-th simulation
				Computed quantities are:
				- k=0 average rate
				- k=1 average rate change for each site
				- k=2 deviation from target configuration
				- k=3 standard deviation of rates
				- k=4 fraction of "fast" sites
bln_fast_data:	pointer for exporting fast sites.
				bln_fast_data[i][j][k][l] will be 1(0) if l-th site of i-th simulation under j-th gamma
				at k-th simulation step during equilibration was fast(slow) according to the default threshold
				(see Simulation.h for threshold definition)
*/
void StrainSweepCore(int sim_number, double* gamma_list, int gamma_num, bool preshear, bool ascdesc, bool print_std,
	bool print_nloops, bool print_fastsites, std::string out_file_name, double** avg_data, double*** all_data,
	double**** dyn_data, bool**** bln_fast_data, bool mute_out) {

	bool sim_init = true;
	double **rate_asc = new double*[sim_number];
	double **rate_desc = new double*[sim_number];
	double **ratestd_asc = new double*[sim_number];
	double **ratestd_desc = new double*[sim_number];
	double **fastsites_asc = new double*[sim_number];
	double **fastsites_desc = new double*[sim_number];
	int **nloop_asc = new int*[sim_number];
	int **nloop_desc = new int*[sim_number];
	for (int j = 0; j < sim_number; j++) {
		rate_asc[j] = new double[gamma_num];
		rate_desc[j] = new double[gamma_num];
		ratestd_asc[j] = new double[gamma_num];
		ratestd_desc[j] = new double[gamma_num];
		fastsites_asc[j] = new double[gamma_num];
		fastsites_desc[j] = new double[gamma_num];
		nloop_asc[j] = new int[gamma_num];
		nloop_desc[j] = new int[gamma_num];
	}
	double *list_unique_n = new double[sim_number];

	if (!mute_out) {
		std::cout << "MULTI STRAIN SWEEP ON " << sim_number << " SIMULATIONS" << std::endl;
		std::cout << gamma_num << " gammas, from " << gamma_list[0] << " to " << gamma_list[gamma_num - 1] << ": [";
		for (int i = 0; i < gamma_num; i++) {
			std::cout << gamma_list[i] << "; ";
			if (i > 4) {
				std::cout << "...";
				break;
			}
		}
		std::cout << "]" << std::endl;
	}
	for (int j = 0; j < sim_number; j++) {
		if (::simList[j].GetStatus() == SimStatus::NOT_INITIALIZED) {
			sim_init = false;
			break;
		}
	}
	if (sim_init = false) {
		std::cout << "ERROR: at least one simulation was not initialized correctly" << std::endl;
		return;
	}

	list_unique_n[0] = ::simList[0].GetParam_n();
	int num_unique_n = 1;
	for (int i = 1; i < sim_number; i++) {
		double cur_n = ::simList[i].GetParam_n();
		bool duplicate = false;
		for (int j = 0; j < i; j++) {
			if (cur_n == list_unique_n[j]) {
				duplicate = true;
				break;
			}
		}
		if (duplicate == false) {
			list_unique_n[num_unique_n] = cur_n;
			num_unique_n++;
		}
	}

	// Write file header
	std::ofstream fout;
	if (out_file_name != "") {
		fout.open(out_file_name);
		fout << "gamma0";
		for (int j = 0; j < sim_number; j++) {
			if (ascdesc) {
				fout << "\tAvgRate_" << j << "_ASC[1/s]";
				if (print_std) {
					fout << "\tStdRate_" << j << "_ASC[1/s]";
				}
				if (print_nloops) {
					fout << "\tnLoops_" << j << "_ASC";
				}
				if (print_fastsites) {
					fout << "\tfastSites_" << j << "_ASC";
				}
				fout << "\tAvgRate_" << j << "_DESC[1/s]";
				if (print_std) {
					fout << "\tStdRate_" << j << "_DESC[1/s]";
				}
				if (print_nloops) {
					fout << "\tnLoops_" << j << "_DESC";
				}
				if (print_fastsites) {
					fout << "\tfastSites_" << j << "_DESC";
				}
			}
			else {
				fout << "\tAvgRate_" << j << "[1/s]";
				if (print_std) {
					fout << "\tStdRate_" << j << "[1/s]";
				}
				if (print_nloops) {
					fout << "\tnLoops_" << j;
				}
				if (print_fastsites) {
					fout << "\tfastSites_" << j;
				}
			}
		}
		fout << "\t1/gamma0";
		for (int j = 0; j < num_unique_n; j++) {
			fout << "\t1/gamma0^n=" << list_unique_n[j];
		}
		fout << std::endl;
	}

	// Run strain sweep
	for (int i = 0; i < gamma_num; i++) {
		if (!mute_out) {
			std::cout << "Current gamma=" << gamma_list[i] << std::endl;
		}
		for (int j = 0; j < sim_number; j++) {
			if (preshear) {
				::simList[j].GenInitialConditions();
			}
			::simList[j].SetGamma(gamma_list[i], false);
			double** dyn_data_ptr = NULL;
			bool** fast_data_ptr = NULL;
			if (dyn_data == NULL) {
				dyn_data_ptr = NULL;
			}
			else {
				dyn_data_ptr = dyn_data[j][i];
			}
			if (bln_fast_data == NULL) {
				fast_data_ptr = NULL;
			}
			else {
				fast_data_ptr = bln_fast_data[j][i];
			}
			nloop_asc[j][i] = ::simList[j].AdjustLattice(dyn_data_ptr, fast_data_ptr);
			rate_asc[j][i] = ::simList[j].GetAverageRate();
			ratestd_asc[j][i] = ::simList[j].GetRateStd();
			fastsites_asc[j][i] = ::simList[j].GetFractionFastSites();
			if (avg_data != NULL) {
				avg_data[j][i] = rate_asc[j][i];
			}
			if (all_data != NULL) {
				for (int k = 0; k < ::simList[j].GetNumSites(); k++) {
					all_data[j][i][k] = ::simList[j].GetRate(k);
				}
			}
		}
	}
	if (ascdesc) {
		for (int i = gamma_num - 1; i >= 0; i--) {
			if (!mute_out) {
				std::cout << "Current gamma=" << gamma_list[i] << std::endl;
			}
			for (int j = 0; j < sim_number; j++) {
				if (preshear) {
					::simList[j].GenInitialConditions();
				}
				::simList[j].SetGamma(gamma_list[i], false);
				if (dyn_data != NULL) {
					nloop_desc[j][i] = ::simList[j].AdjustLattice(dyn_data[j][2 * gamma_num - 1 - i]);
				}
				else {
					nloop_desc[j][i] = ::simList[j].AdjustLattice();
				}
				rate_desc[j][i] = ::simList[j].GetAverageRate();
				ratestd_desc[j][i] = ::simList[j].GetRateStd();
				fastsites_desc[j][i] = ::simList[j].GetFractionFastSites();
				if (avg_data != NULL) {
					avg_data[j][2 * gamma_num - 1 - i] = rate_asc[j][i];
				}
				if (all_data != NULL) {
					for (int k = 0; k < ::simList[j].GetNumSites(); k++) {
						all_data[j][2 * gamma_num - 1 - i][k] = ::simList[j].GetRate(k);
					}
				}
			}
		}
	}

	// Write data to file
	if (out_file_name != "") {
		for (int i = 0; i < gamma_num; i++) {
			fout << std::setprecision(10) << gamma_list[i];
			for (int j = 0; j < sim_number; j++) {
				fout << "\t" << rate_asc[j][i];
				if (print_std) {
					fout << "\t" << ratestd_asc[j][i];
				}
				if (print_nloops) {
					fout << "\t" << nloop_asc[j][i];
				}
				if (print_fastsites) {
					fout << "\t" << fastsites_asc[j][i];
				}
				if (ascdesc) {
					fout << "\t" << rate_desc[j][i];
					if (print_std) {
						fout << "\t" << ratestd_desc[j][i];
					}
					if (print_nloops) {
						fout << "\t" << nloop_desc[j][i];
					}
					if (print_fastsites) {
						fout << "\t" << fastsites_desc[j][i];
					}
				}
			}
			fout << "\t" << 1.0 / gamma_list[i];
			for (int j = 0; j < num_unique_n; j++) {
				fout << "\t" << 1.0 / pow(gamma_list[i], list_unique_n[j]) << std::endl;
			}
		}
		fout.close();
	}


	for (int j = 0; j < sim_number; j++) {
		delete[] rate_asc[j];
		delete[] rate_desc[j];
		delete[] ratestd_asc[j];
		delete[] ratestd_desc[j];
		delete[] fastsites_asc[j];
		delete[] fastsites_desc[j];
		delete[] nloop_asc[j];
		delete[] nloop_desc[j];
		rate_asc[j] = NULL;
		rate_desc[j] = NULL;
		ratestd_asc[j] = NULL;
		ratestd_desc[j] = NULL;
		fastsites_asc[j] = NULL;
		fastsites_desc[j] = NULL;
		nloop_asc[j] = NULL;
		nloop_desc[j] = NULL;
	}
	delete[] rate_asc;
	delete[] rate_desc;
	delete[] ratestd_asc;
	delete[] ratestd_desc;
	delete[] fastsites_asc;
	delete[] fastsites_desc;
	delete[] nloop_asc;
	delete[] nloop_desc;
	rate_asc = NULL;
	rate_desc = NULL;
	ratestd_asc = NULL;
	ratestd_desc = NULL;
	fastsites_asc = NULL;
	fastsites_desc = NULL;
	nloop_asc = NULL;
	nloop_desc = NULL;

	delete[] list_unique_n;
	list_unique_n = NULL;
}

void RunStrainSweep(int sim_number, double start_gamma, double end_gamma, int ppd_gamma,
	bool preshear, bool ascdesc, bool print_std, bool print_nloops, bool print_fastsites,
	std::string out_file_name, std::string* out_file_stats, std::string* out_file_histo,
	std::string* out_file_eq, std::string* out_file_corr, std::string* out_raw,
	std::string* out_raw_final, std::string out_file_histalpha, bool mute_out) {

	int num_gamma = LogSpace_CalcNum(start_gamma, end_gamma, ppd_gamma);
	double* gamma_list = new double[num_gamma];
	LogSpaceNum(start_gamma, end_gamma, num_gamma, gamma_list);
#if VERBOSE:
	std::cout << num_gamma << " gammas, from " << start_gamma << " to " << end_gamma << ": [";
	for (int i = 0; i < num_gamma; i++) {
		std::cout << gamma_list[i] << "; ";
		if (i > 4) {
			std::cout << "...";
			break;
		}
	}
	std::cout << "]" << std::endl;
#endif

	RunStrainSweep(sim_number, gamma_list, num_gamma, preshear, ascdesc, print_std, print_nloops, print_fastsites,
		out_file_name, out_file_stats, out_file_histo, out_file_eq, out_file_corr, out_raw, out_raw_final, out_file_histalpha, mute_out);

	delete[] gamma_list;
	gamma_list = NULL;
}

void RunStrainSweep(int sim_number, double* gamma_list, int num_gamma,
	bool preshear, bool ascdesc, bool print_std, bool print_nloops, bool print_fastsites,
	std::string out_file_name, std::string* out_file_stats, std::string* out_file_histo,
	std::string* out_file_eq, std::string* out_file_corr, std::string* out_raw,
	std::string* out_raw_final, std::string out_file_histalpha, bool mute_out) {

#if DEBUG_MODE:
	std::cout << "Strain sweep setup: preparing memory." << std::endl;
#endif

	int len_res = num_gamma;
	if (ascdesc) {
		len_res = 2 * num_gamma;
	}

	double** data = NULL;
	double*** all_data = NULL;
	double**** eq_data = NULL;
	bool**** bln_fast_data = NULL;

	if (out_file_stats != NULL || out_file_histo != NULL || out_file_corr != NULL) {
		data = new double*[sim_number];
		for (int i = 0; i < sim_number; i++) {
			data[i] = new double[len_res];
		}
#if DEBUG_MODE:
		std::cout << "   [DEBUG]: double** data initialized [" << sim_number << "," << len_res << "]" << std::endl;
#endif
	}
	else {
		data = NULL;
	}
	if (out_file_histo != NULL || out_file_corr != NULL || out_raw_final != NULL) {
		all_data = new double**[sim_number];
		for (int i = 0; i < sim_number; i++) {
			all_data[i] = new double*[len_res];
			for (int j = 0; j < len_res; j++) {
				all_data[i][j] = new double[::simList[i].GetNumSites()];
#if DEBUG_MODE:
				for (int k = 0; k < ::simList[i].GetNumSites(); k++) {
					all_data[i][j][k] = BAD_VALUE;
				}
#endif
			}
		}
#if DEBUG_MODE:
		std::cout << "   [DEBUG]: double*** all_data initialized [" << sim_number << "," 
			      << len_res << "," << ::simList[0].GetNumSites() << "]" << std::endl;
#endif
	}
	else {
		all_data = NULL;
	}
	if (out_raw != NULL) {
		bln_fast_data = new bool***[sim_number];
		for (int i = 0; i < sim_number; i++) {
			bln_fast_data[i] = new bool**[len_res];
			for (int j = 0; j < len_res; j++) {
				bln_fast_data[i][j] = new bool*[::simList[i].GetMaxStepNum()];
				for (int k = 0; k < ::simList[i].GetMaxStepNum(); k++) {
					bln_fast_data[i][j][k] = new bool[::simList[i].GetNumSites()];
					for (int l = 0; l < ::simList[i].GetNumSites(); l++) {
						bln_fast_data[i][j][k][l] = false;
					}
				}
			}
		}
#if DEBUG_MODE:
		std::cout << "   [DEBUG]: bool**** bln_fast_data initialized [" << sim_number << "," << len_res << "," << ::simList[0].GetMaxStepNum() 
			      << "," << ::simList[0].GetNumSites() << "]" << std::endl;
#endif
	}
	else {
		bln_fast_data = NULL;
	}
	if (out_file_eq != NULL) {
		eq_data = new double***[sim_number];
		for (int i = 0; i < sim_number; i++) {
			eq_data[i] = new double**[len_res];
			for (int j = 0; j < len_res; j++) {
				eq_data[i][j] = new double*[NUM_EQ_DETAILS];
				for (int k = 0; k < NUM_EQ_DETAILS; k++) {
					eq_data[i][j][k] = new double[::simList[i].GetMaxStepNum()];
					for (int l = 0; l < ::simList[i].GetMaxStepNum(); l++) {
						eq_data[i][j][k][l] = BAD_VALUE;
					}
				}
			}
		}
#if DEBUG_MODE:
		std::cout << "   [DEBUG]: double**** eq_data initialized [" << sim_number << "," << len_res << "," 
			      << NUM_EQ_DETAILS << "," << ::simList[0].GetMaxStepNum() << "]" << std::endl;
#endif
	}
	else {
		eq_data = NULL;
	}

#if VERBOSE:
	std::cout << "Starting strain sweep." << std::endl;
#endif

	////////////////////////////////////////////////////////////////////////////////////////////
	StrainSweepCore(sim_number, gamma_list, num_gamma, preshear, ascdesc, print_std, print_nloops, print_fastsites,
		out_file_name, data, all_data, eq_data, bln_fast_data);
	////////////////////////////////////////////////////////////////////////////////////////////

#if VERBOSE:
	std::cout << "Strain sweep completed. Now saving output files." << std::endl;
#endif

	if (out_file_stats != NULL) {

#if VERBOSE:
		std::cout << "... Now saving statistics to file " << *out_file_stats << std::endl;
#endif

		std::ofstream fout(*out_file_stats);
		fout << "gamma0\tAvgRate[1/s]\tStdRate[1/s]" << std::endl;
		double *vals = new double[sim_number];
		for (int i = 0; i < num_gamma; i++) {
			for (int j = 0; j < sim_number; j++) {
				vals[j] = data[j][i];
			}
			double cur_avg = CalcAverage(vals, sim_number);
			double cur_std = CalcStd(vals, sim_number);
			//std::cout << gamma_list[i] << std::endl;
			fout << std::setprecision(10) << gamma_list[i] << "\t" << cur_avg << "\t";
			if (sim_number == 1) {
				fout << "-" << std::endl;;
			}
			else {
				fout << cur_std << std::endl;;
			}
		}
		if (ascdesc) {
			for (int i = num_gamma - 1; i >= 0; i--) {
				for (int j = 0; j < sim_number; j++) {
					vals[j] = data[j][2 * num_gamma - 1 - i];
				}
				double cur_avg = CalcAverage(vals, sim_number);
				double cur_std = CalcStd(vals, sim_number);
				fout << gamma_list[i] << "\t" << cur_avg << "\t" << cur_std << std::endl;
			}
		}
		delete[] vals;
		vals = NULL;
		fout.close();
	}

	if (out_file_histo != NULL || out_file_histalpha != "") {

#if VERBOSE:
		std::cout << "... Now saving histograms..." << std::endl;
#endif

		double** alphas = NULL;
		double** logalphas = NULL;
		double** glassiness = NULL;
		double** logdistfromgmax = NULL;
		double** jumpstrains = NULL;
		double** logjumpstrains = NULL;
		int* numalphas = NULL;
		double** logavgdata = NULL;
		double*** logdata = NULL;
		int* len2 = NULL;
		int** len3 = NULL;
		double g_max = 0, g_min = 0, la_min = 0, la_max = 0, a_max = 0, a_min = 0, gl_max = 0, gl_min = 0;
		double nd_max = 0, nd_min = 0, js_max = 0, js_min = 0, ljs_max = 0, ljs_min = 0;

		if (out_file_histalpha != "") {
			alphas = new double*[sim_number];
			logalphas = new double* [sim_number];
			glassiness = new double* [sim_number];
			logdistfromgmax = new double* [sim_number];
			jumpstrains = new double* [sim_number];
			logjumpstrains = new double* [sim_number];
			numalphas = new int[sim_number];
			for (int i = 0; i < sim_number; i++) {
				numalphas[i] = ::simList[i].CountAlphas();
				alphas[i] = new double[numalphas[i]];
				logalphas[i] = new double[numalphas[i]];
				glassiness[i] = new double[numalphas[i]];
				logdistfromgmax[i] = new double[numalphas[i]];
				jumpstrains[i] = new double[numalphas[i]];
				logjumpstrains[i] = new double[numalphas[i]];
				::simList[i].FlattenAlphas(alphas[i]);
				double cur_gmax = ::simList[i].GetGMax();
				for (int j = 0; j < numalphas[i]; j++) {
					if (alphas[i][j] <= 0) {
						alphas[i][j] = BAD_VALUE_2;
						logalphas[i][j] = BAD_VALUE_2;
						glassiness[i][j] = BAD_VALUE_2;
						logdistfromgmax[i][j] = BAD_VALUE_2;
						jumpstrains[i][j] = -1;
						logjumpstrains[i][j] = BAD_VALUE_2;
					}
					else {
						logalphas[i][j] = log10(alphas[i][j]);
						glassiness[i][j] = ::simList[i].GetGlassinessFromAlpha(alphas[i][j]);
						if (glassiness[i][j] < cur_gmax) {
							logdistfromgmax[i][j] = log10(::simList[i].GetDistFromGmax(alphas[i][j]));
							jumpstrains[i][j] = ::simList[i].GetYieldStrainFromAlpha(alphas[i][j]);
							if (jumpstrains[i][j] > 0) logjumpstrains[i][j] = log10(jumpstrains[i][j]);
							else {
								logjumpstrains[i][j] = BAD_VALUE_2;
							}
						}
						else {
							logdistfromgmax[i][j] = BAD_VALUE_2;
							jumpstrains[i][j] = -1;
							logjumpstrains[i][j] = BAD_VALUE_2;
						}
					}
				}
			}
			//Log10Array(alphas, sim_number, numalphas, logalphas);
			CalcGlobalMaxMin(alphas, sim_number, numalphas, a_max, a_min, BAD_VALUE_2);
			CalcGlobalMaxMin(logalphas, sim_number, numalphas, la_max, la_min, BAD_VALUE_2);
			CalcGlobalMaxMin(glassiness, sim_number, numalphas, gl_max, gl_min, BAD_VALUE_2);
			CalcGlobalMaxMin(logdistfromgmax, sim_number, numalphas, nd_max, nd_min, BAD_VALUE_2);
			CalcGlobalMaxMin(jumpstrains, sim_number, numalphas, js_max, js_min, -1);
			CalcGlobalMaxMin(logjumpstrains, sim_number, numalphas, ljs_max, ljs_min, BAD_VALUE_2);
		}
		if (out_file_histo != NULL) {
			logavgdata = new double*[sim_number];
			logdata = new double**[sim_number];
			len2 = new int[sim_number];
			len3 = new int*[sim_number];
			for (int i = 0; i < sim_number; i++) {
				len2[i] = len_res;
				len3[i] = new int[len_res];
				logavgdata[i] = new double[len_res];
				logdata[i] = new double*[len_res];
				for (int j = 0; j < len_res; j++) {
					len3[i][j] = ::simList[i].GetNumSites();
					logdata[i][j] = new double[::simList[i].GetNumSites()];
				}
			}
			Log10Array(data, sim_number, len_res, logavgdata);
			Log10Array(all_data, sim_number, len2, len3, logdata);
			CalcGlobalMaxMin(logdata, sim_number, len2, len3, g_max, g_min);
		}

		for (int j = 0; j < sim_number; j++) {

#if DEBUG_MODE:
			std::cout << "    - Building histograms for simulation " << j + 1 << "/" << sim_number << "..." << std::endl;
#endif

			//int nbins = (int)pow(::simList[j].GetNumSites(), 1.0 / 3) + 1;
			int nbins = (int)pow(::simList[j].GetNumSites(), 1.0 / 2) + 1;
#if DEBUG_MODE:
			std::cout << "      ... current histogram will have " << nbins << " bins, lin spaced from " << g_min << " to " << g_max << std::endl;
#endif
			if (out_file_histo != NULL) {
				double* bins = new double[nbins];
				LinSpaceNum(g_min, g_max + (g_max - g_min) / nbins, nbins, bins);
#if DEBUG_MODE:
				std::cout << "      ... bin sample: [";
				for (int b = 0; b < nbins; b++) {
					if (b > 0) {
						std::cout << ";";
					}
					std::cout << bins[b];
					if (b > 4) {
						std::cout << ";...";
						break;
					}
				}
				std::cout << "]" << std::endl;
#endif

				double** cols = new double* [len_res];
				for (int k = 0; k < len_res; k++) {
#if DEBUG_MODE:
					std::cout << "       - Processing gamma " << k + 1 << "/" << len_res << "..." << std::endl;
					std::cout << "         initializing cols[" << k << "] with " << nbins << "-double array" << std::endl;
#endif
					cols[k] = new double[nbins];
					int new_nbins = nbins;
#if DEBUG_MODE:
					std::cout << "          cols[" << k << "] initialized" << std::endl;
					for (int t = 0; t < nbins; t++) {
						cols[k][t] = 0;
					}
					std::cout << "          sample logdata[" << j << "][" << k << "] (len=" << ::simList[j].GetNumSites() << "): [";
					for (int t = 0; t < ::simList[j].GetNumSites(); t++) {
						if (t > 0) {
							std::cout << ";";
						}
						std::cout << logdata[j][k][t];
						if (t > 4) {
							std::cout << ";...";
							break;
						}
					}
					std::cout << "]" << std::endl;
					std::cout << "          len3[" << j << "][" << k << "] =  " << len3[j][k] << std::endl;
					std::cout << "          sample bins before (len=" << nbins << "): [";
					for (int t = 0; t < nbins; t++) {
						if (t > 0) {
							std::cout << ";";
						}
						std::cout << bins[t];
						if (t > 4) {
							std::cout << ";...;" << bins[nbins - 1];
							break;
						}
					}
					std::cout << "]" << std::endl;
#endif
					BuildHistogram(logdata[j][k], len3[j][k], bins, cols[k], new_nbins, true, false, false, false);
					if (nbins != new_nbins) {
						std::cout << "\n\n\nWARNING! BuildHistogram() function modified number of bins from " << nbins << " to " 
							      << new_nbins << ". Possible data corruption...\n\n\n" << std::endl;
					}
#if DEBUG_MODE:
					std::cout << "          sample bins after (len=" << nbins << "): [";
					for (int t = 0; t < nbins; t++) {
						if (t > 0) {
							std::cout << ";";
						}
						std::cout << bins[t];
						if (t > 4) {
							std::cout << ";...;" << bins[nbins - 1];
							break;
						}
					}
					std::cout << "]" << std::endl;
					std::cout << "          sample cols[" << k << "] after (len=" << nbins << "): [";
					for (int t = 0; t < nbins; t++) {
						if (t > 0) {
							std::cout << ";";
						}
						std::cout << cols[k][t];
						if (t > 4) {
							std::cout << ";...;" << cols[k][nbins - 1];
							break;
						}
					}
					std::cout << "]" << std::endl;

#endif
				}


#if DEBUG_MODE:
				std::cout << "      ... Saving histograms to file " << out_file_histo[j] << std::endl;
#endif
				std::ofstream fout(out_file_histo[j]);
				fout << "rate[1/s]";
				for (int i = 0; i < num_gamma; i++) {
					fout << "\tcount_g=" << gamma_list[i];
				}
				if (ascdesc) {
					for (int i = num_gamma - 1; i >= 0; i--) {
						fout << "\tcount_g=" << gamma_list[i];
					}
				}
				fout << std::endl;
#if DEBUG_MODE:
				if ((num_gamma > len_res) || (ascdesc && (2 * num_gamma > len_res))) {
					std::cout << "\nWARNING! Inconsistend len_res=" << len_res << " with " << num_gamma 
						<< " gamma values (ascdesc=" << ascdesc << ") found in building histogram!\n" << std::endl;
				}
#endif
				for (int k = 0; k < nbins; k++) {
					fout << pow(10, bins[k]);
					for (int i = 0; i < num_gamma; i++) {
						fout << "\t" << cols[i][k];
					}
					if (ascdesc) {
						for (int i = num_gamma - 1; i >= 0; i--) {
							fout << "\t" << cols[2 * num_gamma - 1 - i][k];
						}
					}
					fout << std::endl;
				}
				fout.close();

#if DEBUG_MODE:
				std::cout << "          ... ok! Releasing memory for current histogram." << std::endl;
#endif
				for (int i = 0; i < len_res; i++) {
					delete[] cols[i];
					cols[i] = NULL;
				}
				delete[] cols;
				delete[] bins;
				cols = NULL;
				bins = NULL;
			}

			if (out_file_histalpha != "") {
				double* alphabins = new double[nbins];
				double** alphacols = new double* [sim_number];
				double* logalphabins = new double[nbins];
				double** logalphacols = new double* [sim_number];
				double* glassbins = new double[nbins];
				double** glasscols = new double* [sim_number];
				double* normglassbins = new double[nbins];
				double** lognormglasscols = new double* [sim_number];
				double* jumpstrainbins = new double[nbins];
				double** jumpstraincols = new double* [sim_number];
				double* logjumpstrainbins = new double[nbins];
				double** logjumpstraincols = new double* [sim_number];
				LinSpaceNum(a_min, a_max + (a_max - a_min) / nbins, nbins, alphabins);
				LinSpaceNum(la_min, la_max + (la_max - la_min) / nbins, nbins, logalphabins);
				LinSpaceNum(gl_min, gl_max + (gl_max - gl_min) / nbins, nbins, glassbins);
				LinSpaceNum(nd_min, nd_max + (nd_max - nd_min) / nbins, nbins, normglassbins);
				LinSpaceNum(js_min, js_max + (js_max - js_min) / nbins, nbins, jumpstrainbins);
				LinSpaceNum(ljs_min, ljs_max + (ljs_max - ljs_min) / nbins, nbins, logjumpstrainbins);
				for (int j = 0; j < sim_number; j++) {
					alphacols[j] = new double[nbins];
					logalphacols[j] = new double[nbins];
					glasscols[j] = new double[nbins];
					lognormglasscols[j] = new double[nbins];
					jumpstraincols[j] = new double[nbins];
					logjumpstraincols[j] = new double[nbins];
					int new_nbins = nbins;
					int bindiff = 0;
					BuildHistogram(alphas[j], numalphas[j], alphabins, alphacols[j], new_nbins, true, false, false, false);
					bindiff += new_nbins - nbins;
					BuildHistogram(logalphas[j], numalphas[j], logalphabins, logalphacols[j], new_nbins, true, false, false, false);
					bindiff += new_nbins - nbins;
					BuildHistogram(glassiness[j], numalphas[j], glassbins, glasscols[j], new_nbins, true, false, false, false);
					bindiff += new_nbins - nbins;
					BuildHistogram(logdistfromgmax[j], numalphas[j], normglassbins, lognormglasscols[j], new_nbins, true, false, false, false);
					bindiff += new_nbins - nbins;
					BuildHistogram(jumpstrains[j], numalphas[j], jumpstrainbins, jumpstraincols[j], new_nbins, true, false, false, false);
					bindiff += new_nbins - nbins;
					BuildHistogram(logjumpstrains[j], numalphas[j], logjumpstrainbins, logjumpstraincols[j], new_nbins, true, false, false, false);
					bindiff += new_nbins - nbins;
					if (bindiff != 0) {
						std::cout << "\n\n\nWARNING! BuildHistogram() function modified number of bins from " << nbins << " to " 
								  << new_nbins << ". Possible data corruption...\n\n\n" << std::endl;
					}
				}
				std::ofstream fout(out_file_histalpha);
				fout << "alpha(lin)";
				for (int i = 0; i < sim_number; i++) {
					fout << "\tP(alpha)_" << i;
					fout << "\tC(alpha)_" << i;
				}
				fout << "\talpha(log)";
				for (int i = 0; i < sim_number; i++) {
					fout << "\tP(alpha)_" << i;
					fout << "\tC(alpha)_" << i;
				}
				fout << "\tg";
				for (int i = 0; i < sim_number; i++) {
					fout << "\tP(g)_" << i;
					fout << "\tC(g)_" << i;
				}
				fout << "\t1-g/gmax";
				for (int i = 0; i < sim_number; i++) {
					fout << "\tP(1-g/gmax)_" << i;
					fout << "\tC(1-g/gmax)_" << i;
				}
				fout << "\tgamma_y";
				for (int i = 0; i < sim_number; i++) {
					fout << "\tP(gamma_y)_" << i;
					fout << "\tC(gamma_y)_" << i;
				}
				fout << "\tgamma_y(log)";
				for (int i = 0; i < sim_number; i++) {
					fout << "\tP(gamma_y)_" << i;
					fout << "\tC(gamma_y)_" << i;
				}
				fout << std::endl;
				double* c_alpha = new double[sim_number];
				double* c_logalpha = new double[sim_number];
				double* c_glass = new double[sim_number];
				double* c_lnglass = new double[sim_number];
				double* c_jumpstrain = new double[sim_number];
				double* c_logjumpstrain = new double[sim_number];
				for (int i = 0; i < sim_number; i++) {
					c_alpha[i] = 0;
					c_logalpha[i] = 0;
					c_glass[i] = 0;
					c_lnglass[i] = 0;
					c_jumpstrain[i] = 0;
					c_logjumpstrain[i] = 0;
				}
				for (int k = 0; k < nbins; k++) {
					fout << alphabins[k];
					for (int i = 0; i < sim_number; i++) {
						fout << "\t" << alphacols[i][k];
						c_alpha[i] += alphacols[i][k];
						fout << "\t" << c_alpha[i];
					}
					fout << "\t" << pow(10, logalphabins[k]);
					for (int i = 0; i < sim_number; i++) {
						fout << "\t" << logalphacols[i][k];
						c_logalpha[i] += logalphacols[i][k];
						fout << "\t" << c_logalpha[i];
					}
					fout << "\t" << glassbins[k];
					for (int i = 0; i < sim_number; i++) {
						fout << "\t" << glasscols[i][k];
						c_glass[i] += glasscols[i][k];
						fout << "\t" << c_glass[i];
					}
					fout << "\t" << pow(10, normglassbins[k]);
					for (int i = 0; i < sim_number; i++) {
						fout << "\t" << lognormglasscols[i][k];
						c_lnglass[i] += lognormglasscols[i][k];
						fout << "\t" << c_lnglass[i];
					}
					fout << "\t" << jumpstrainbins[k];
					for (int i = 0; i < sim_number; i++) {
						fout << "\t" << jumpstraincols[i][k];
						c_jumpstrain[i] += jumpstraincols[i][k];
						fout << "\t" << c_jumpstrain[i];
					}
					fout << "\t" << pow(10, logjumpstrainbins[k]);
					for (int i = 0; i < sim_number; i++) {
						fout << "\t" << logjumpstraincols[i][k];
						c_logjumpstrain[i] += logjumpstraincols[i][k];
						fout << "\t" << c_logjumpstrain[i];
					}
					fout << std::endl;
				}
				fout.close();
				for (int j = 0; j < sim_number; j++) {
					delete[] alphacols[j];
					delete[] logalphacols[j];
					delete[] glasscols[j];
					delete[] lognormglasscols[j];
					delete[] jumpstraincols[j];
					delete[] logjumpstraincols[j];
					alphacols[j] = NULL;
					logalphacols[j] = NULL;
					glasscols[j] = NULL;
					lognormglasscols[j] = NULL;
					jumpstraincols[j] = NULL;
					logjumpstraincols[j] = NULL;
				}
				delete[] alphabins;
				delete[] alphacols;
				delete[] logalphabins;
				delete[] logalphacols;
				delete[] glassbins;
				delete[] glasscols;
				delete[] normglassbins;
				delete[] lognormglasscols;
				delete[] jumpstrainbins;
				delete[] jumpstraincols;
				delete[] logjumpstrainbins;
				delete[] logjumpstraincols;
				delete[] c_alpha;
				delete[] c_logalpha;
				delete[] c_glass;
				delete[] c_lnglass;
				delete[] c_jumpstrain;
				delete[] c_logjumpstrain;
				alphabins = NULL;
				alphacols = NULL;
				logalphabins = NULL;
				logalphacols = NULL;
				glassbins = NULL;
				glasscols = NULL;
				normglassbins = NULL;
				lognormglasscols = NULL;
				jumpstrainbins = NULL;
				jumpstraincols = NULL;
				logjumpstrainbins = NULL;
				logjumpstraincols = NULL;
				c_alpha = NULL;
				c_logalpha = NULL;
				c_glass = NULL;
				c_lnglass = NULL;
				c_jumpstrain = NULL;
				c_logjumpstrain = NULL;
			}

		}
#if DEBUG_MODE:
		std::cout << "... Now releasing memory used to compute all histograms" << std::endl;
#endif
		if (out_file_histalpha != "") {
			for (int i = 0; i < sim_number; i++) {
				delete[] alphas[i];
				delete[] logalphas[i];
				delete[] glassiness[i];
				delete[] logdistfromgmax[i];
				delete[] jumpstrains[i];
				delete[] logjumpstrains[i];
				alphas[i] = NULL;
				logalphas[i] = NULL;
				glassiness[i] = NULL;
				logdistfromgmax[i] = NULL;
				jumpstrains[i] = NULL;
				logjumpstrains[i] = NULL;
			}
			delete[] alphas;
			delete[] logalphas;
			delete[] numalphas;
			delete[] glassiness;
			delete[] logdistfromgmax;
			delete[] jumpstrains;
			delete[] logjumpstrains;
			alphas = NULL;
			logalphas = NULL;
			numalphas = NULL;
			glassiness = NULL;
			logdistfromgmax = NULL;
			jumpstrains = NULL;
			logjumpstrains = NULL;
		}
		if (out_file_histo != NULL) {
			for (int i = 0; i < sim_number; i++) {
				for (int j = 0; j < len_res; j++) {
					delete[] logdata[i][j];
					logdata[i][j] = NULL;
				}
				delete[] logdata[i];
				delete[] logavgdata[i];
				delete[] len3[i];
				logdata[i] = NULL;
				logavgdata[i] = NULL;
				len3[i] = NULL;
			}
			delete[] logavgdata;
			delete[] logdata;
			delete[] len3;
			delete[] len2;
			logdata = NULL;
			logavgdata = NULL;
			len3 = NULL;
			len2 = NULL;
		}
#if DEBUG_MODE:
		std::cout << "    ... done!" << std::endl;
#endif
	}

	if (out_file_corr != NULL) {

#if VERBOSE:
		std::cout << "... Now saving correlation functions..." << std::endl;
#endif

		int piLatticeShape[N_DIM];
		int* len_corr = new int[sim_number];
		int* len_corr_list = new int[sim_number];
		int** r_list = new int*[sim_number];
		double*** corr_list = new double**[sim_number];
		for (int j = 0; j < sim_number; j++) {
			::simList[j].CopyLatticeShape(piLatticeShape);
			if (::simList[j].GetParam_UsePBC()) {
				len_corr[j] = piLatticeShape[0] / 2;
			}
			else {
				len_corr[j] = piLatticeShape[0];
			}
			r_list[j] = new int[len_corr[j]];
			if (ascdesc) {
				len_corr_list[j] = 2 * num_gamma;
			}
			else {
				len_corr_list[j] = num_gamma;
			}
			corr_list[j] = new double*[len_corr_list[j]];
			for (int i = 0; i < num_gamma; i++) {
				//all_data[i][j][k] will be the rate of k-th site in i-th simulation at j-th gamma
				corr_list[j][i] = new double[len_corr[j]];
#if CORR_BINARY
				bool* bln_data = new bool[::simList[j].GetNumSites()];
				ArrayGreaterThan(all_data[j][i], ::simList[j].GetNumSites(), MAX_SOLID_RATE_RATIO / ::simList[j].GetParam_Tau0(), bln_data);
				len_corr[j] = AutoCorrelation(bln_data, piLatticeShape, corr_list[j][i], r_list[j], -1, ::simList[j].GetParam_UsePBC(), true, true);
				delete[] bln_data;
				bln_data = NULL;
#else
				len_corr[j] = AutoCorrelation(all_data[j][i], ::simList[j].GetParams().LatticeShape, corr_list[j][i], r_list[j], -1, ::simList[j].GetParams().UsePBC, true, true);
#endif
			}
			if (ascdesc) {
				for (int i = num_gamma - 1; i >= 0; i++) {
					corr_list[j][2 * num_gamma - 1 - i] = new double[len_corr[j]];
					len_corr[j] = AutoCorrelation(all_data[j][2 * num_gamma - 1 - i], piLatticeShape, corr_list[j][2 * num_gamma - 1 - i], r_list[j], -1, ::simList[j].GetParam_UsePBC());
				}
			}
		}

		for (int j = 0; j < sim_number; j++) {

#if VERBOSE:
			std::cout << "    - Correlation functions for simulation " << j + 1 << "/" << sim_number << " to file " << out_file_corr[j] << std::endl;
#endif

			std::ofstream fout(out_file_corr[j]);
			fout << "r[sites]";
			for (int i = 0; i < num_gamma; i++) {
				fout << "\tcorr_g=" << gamma_list[i];
			}
			if (ascdesc) {
				for (int i = num_gamma - 1; i >= 0; i++) {
					fout << "\tcount_g=" << gamma_list[i];
				}
			}
			fout << std::endl;
			for (int k = 0; k < len_corr[j]; k++) {
				fout << r_list[j][k];
				for (int i = 0; i < num_gamma; i++) {
#if CORR_NORM_X0
					fout << "\t" << (corr_list[j][i][k] / corr_list[j][i][0]);
#else
					fout << "\t" << corr_list[j][i][k];
#endif
				}
				if (ascdesc) {
					for (int i = num_gamma - 1; i >= 0; i++) {
#if CORR_NORM_X0
						fout << "\t" << (corr_list[j][2 * num_gamma - 1 - i][k] / corr_list[j][2 * num_gamma - 1 - i][0]);
#else
						fout << "\t" << corr_list[j][2 * num_gamma - 1 - i][k];
#endif
					}
				}
				fout << std::endl;
			}
			fout.close();
		}
		for (int j = 0; j < sim_number; j++) {
			for (int k = 0; k < len_corr_list[j]; k++) {
				delete[] corr_list[j][k];
				corr_list[j][k] = NULL;
			}
			delete[] corr_list[j];
			delete[] r_list[j];
			corr_list[j] = NULL;
			r_list[j] = NULL;
		}
		delete[] len_corr;
		delete[] len_corr_list;
		delete[] r_list;
		delete[] corr_list;
		len_corr = NULL;
		len_corr_list = NULL;
		r_list = NULL;
		corr_list = NULL;
	}

	if (out_raw != NULL) {

#if VERBOSE:
		std::cout << "... Now saving raw images..." << std::endl;
#endif

		for (int j = 0; j < sim_number; j++) {
			if (::simList[j].GetNumSites() % 8 != 0) {
				std::cout << "WARNING: site number (" << ::simList[j].GetNumSites() << ") for simulation " << j
					<< " is not multiple of 8. Raw output may be trimmed." << std::endl;
			}
			int nImgs;
			if (ascdesc) {
				nImgs = 2 * num_gamma;
			}
			else {
				nImgs = num_gamma;
			}
			const int exp_dim = 3;
			int* exp_shape = new int[exp_dim];
			exp_shape[0] = nImgs;
			exp_shape[1] = ::simList[j].GetMaxStepNum();
			exp_shape[2] = ::simList[j].GetNumSites();
			bool* blnExport = new bool[CalcProduct(exp_shape, exp_dim)];
			for (int ii = 0; ii < exp_shape[0]; ii++) {
				for (int jj = 0; jj < exp_shape[1]; jj++) {
					for (int kk = 0; kk < exp_shape[2]; kk++) {
						blnExport[ii*exp_shape[1] * exp_shape[2] + jj * exp_shape[2] + kk] = bln_fast_data[j][ii][jj][kk];
					}
				}
			}
			ExportDataToRawBinary(blnExport, exp_shape, exp_dim, out_raw[j]);
			delete[] exp_shape;
			delete[] blnExport;
			exp_shape = NULL;
			blnExport = NULL;
		}
	}

	if (out_raw_final != NULL) {
		int tot_num = 0;
		int* len2 = new int[sim_number];
		int** len3 = new int*[sim_number];
		int l2;
		if (ascdesc) {
			l2 = 2 * num_gamma;
		}
		else {
			l2 = num_gamma;
		}
		for (int i = 0; i < sim_number; i++) {
			len2[i] = l2;
			len3[i] = new int[l2];
			for (int j = 0; j < l2; j++) {
				len3[i][j] = ::simList[i].GetNumSites();
				tot_num += ::simList[i].GetNumSites();
			}
		}
		double* flat_all_data;
		flat_all_data = new double[tot_num];
		Flatten3DArray(all_data, sim_number, len2, len3, flat_all_data);
		int exp_shape[2];
		exp_shape[0] = sim_number;
		exp_shape[1] = l2 * ::simList[0].GetNumSites();
		ExportDataToRawBinary(flat_all_data, exp_shape, 2, out_raw_final[0]);
		delete[] flat_all_data;
		for (int i = 0; i < sim_number; i++) {
			delete[] len3[i];
		}
		delete[] len2;
		delete[] len3;
	}

	if (out_file_eq != NULL) {

#if VERBOSE:
		std::cout << "... Now saving equilibration data..." << std::endl;
#endif

		for (int j = 0; j < sim_number; j++) {

#if VERBOSE:
			std::cout << "    - Simulation " << j + 1 << "/" << sim_number << " to file " << out_file_corr[j] << std::endl;
#endif

			std::ofstream fout(out_file_eq[j]);
			fout << "sim_step";
			for (int k = 0; k < NUM_EQ_DETAILS; k++) {
				for (int i = 0; i < len_res; i++) {
					fout << "\t" << strEqLabels[k] << "_g=" << std::setprecision(10) << gamma_list[i];
				}
			}
			fout << std::endl;
			for (int l = 0; l < ::simList[j].GetMaxStepNum(); l++) {
				fout << l;
				for (int k = 0; k < NUM_EQ_DETAILS; k++) {
					for (int i = 0; i < len_res; i++) {
						fout << "\t";
						if (eq_data[j][i][k][l] != BAD_VALUE) {
							fout << eq_data[j][i][k][l];
						}
					}
				}
				fout << std::endl;
			}
			fout.close();
		}
	}
#if DEBUG_MODE:
	std::cout << "...All output saved succesfully! Now releasing memory..." << std::endl;
#endif
	if (data != NULL) {
		for (int i = 0; i < sim_number; i++) {
			delete[] data[i];
			data[i] = NULL;
		}
#if DEBUG_MODE:
		std::cout << "   [DEBUG]: double** data released" << std::endl;
#endif
	}
	if (all_data != NULL) {
		for (int i = 0; i < sim_number; i++) {
			for (int j = 0; j < len_res; j++) {
				delete[] all_data[i][j];
				all_data[i][j] = NULL;
			}
			delete[] all_data[i];
			all_data[i] = NULL;
		}
#if DEBUG_MODE:
		std::cout << "   [DEBUG]: double*** all_data released" << std::endl;
#endif
	}
	if (bln_fast_data != NULL) {
		for (int i = 0; i < sim_number; i++) {
			for (int j = 0; j < len_res; j++) {
				for (int k = 0; k < ::simList[i].GetMaxStepNum(); k++) {
					delete[] bln_fast_data[i][j][k];
					bln_fast_data[i][j][k] = NULL;
				}
				delete[] bln_fast_data[i][j];
				bln_fast_data[i][j] = NULL;
			}
			delete[] bln_fast_data[i];
			bln_fast_data[i] = NULL;
		}
#if DEBUG_MODE:
		std::cout << "   [DEBUG]: bool**** bln_fast_data released" << std::endl;
#endif
	}
	if (eq_data != NULL) {
		for (int i = 0; i < sim_number; i++) {
			for (int j = 0; j < len_res; j++) {
				for (int k = 0; k < NUM_EQ_DETAILS; k++) {
					delete[] eq_data[i][j][k];
					eq_data[i][j][k] = NULL;
				}
				delete[] eq_data[i][j];
				eq_data[i][j] = NULL;
			}
			delete[] eq_data[i];
			eq_data[i] = NULL;
		}
#if DEBUG_MODE:
		std::cout << "   [DEBUG]: double**** eq_data released" << std::endl;
#endif
	}
	delete[] data;
	delete[] all_data;
	delete[] bln_fast_data;
	delete[] eq_data;
	data = NULL;
	all_data = NULL;
	bln_fast_data = NULL;
	eq_data = NULL;
	if (!mute_out) {
		std::cout << "...done!" << std::endl;
	}
}

void print_sim_params(std::string out_file_name, int sim_number) {
	double _gK = -1;
	double _gn = -1;
	double _gTau0 = -1;
	double _gOmega = -1;
	double _gBeta = -1;
	bool _common_params = true;

	std::ofstream fout(out_file_name);
	fout << "#\tID\tNumSites\tOmega\tTau0\tK\tBeta\tn\tAlphaAvg\tAlphaVar\tAlphaNormVar\tAlphaAvg_nom\tAlphaNoise\tAlphaParams" 
		 << "\tInitAvg\tInitNoise\tInitParams\tTempProtocol\tTempParams\tForceMeanField\tGamma0\tGammac\tgammac\tKc"
		 << "\tKr=psi\tg\t1-g/gmax\tAlphaMax_MF\tcount_unphysical" << std::endl;
	for (int i = 0; i < sim_number; i++) {
		IsingParameters cur_params;
		Noise init_noise;
		::simList[i].CopyParams(&cur_params);
		::simList[i].CopyInitNoise(&init_noise);
		fout << i << "\t" << ::simList[i].ID() << "\t" << cur_params.NumSites() << "\t" << cur_params.Omega 
			<< "\t" << cur_params.Tau0 << "\t" << cur_params.K;
		if (cur_params.UseNewEqn) {
			fout << "\t[new_eqn]";
		}
		else {
			fout << "\t" << cur_params.Beta;
		}
		double cur_alphaAvg = ::simList[i].GetAverageAlpha();
		double cur_alphaVar = ::simList[i].GetAlphaVariance();
		fout << "\t" << cur_params.n << "\t" << cur_alphaAvg << "\t" << cur_alphaVar << "\t" << cur_alphaVar / (cur_alphaAvg * cur_alphaAvg);
		fout << "\t" << cur_params.AlphaAverage
			<< "\t" << cur_params.AlphaNoise.GetType() << "\t";
		if (cur_params.AlphaNoise.NumParams() == 1) {
			fout << cur_params.AlphaNoise.GetParameter(0);
		}
		else if (cur_params.AlphaNoise.NumParams() > 1) {
			fout << "{";
			int num_print_param = cur_params.AlphaNoise.NumParams();
			if (num_print_param > MAX_PRINT_PARAMS) num_print_param = MAX_PRINT_PARAMS;
			for (int i = 0; i < num_print_param; i++) {
				fout << cur_params.AlphaNoise.GetParameter(i) << ";";
			}
			if (num_print_param < cur_params.AlphaNoise.NumParams()) {
				fout << "...";
				for (int i = cur_params.AlphaNoise.NumParams()-3; i < cur_params.AlphaNoise.NumParams(); i++) {
					fout << ";" << cur_params.AlphaNoise.GetParameter(i);
				}
			}
			fout << "}";
		}
		fout << "\t" << init_noise.GetAverage() << "\t" << init_noise.GetType() << "\t";
		if (init_noise.NumParams() == 1) {
			fout << init_noise.GetParameter(0);
		}
		else if (init_noise.NumParams() > 1) {
			fout << "{";
			int num_print_pinit = init_noise.NumParams();
			if (num_print_pinit > MAX_PRINT_PARAMS) num_print_pinit = MAX_PRINT_PARAMS;
			for (int i = 0; i < num_print_pinit; i++) {
				fout << init_noise.GetParameter(i) << ";";
			}
			if (num_print_pinit < init_noise.NumParams()) {
				fout << "...";
				for (int i = init_noise.NumParams() - 3; i < init_noise.NumParams(); i++) {
					fout << ";" << init_noise.GetParameter(i);
				}
			}
			fout << "}";
		}
		//double* _temp_params = ::simList[i].GetEquilibrationParameters();
		fout << "\t" << ::simList[i].GetEquilibrationProtocol() << "\t";
		if (::simList[i].GetEquilibrationParameterNumber() == 1) {
			fout << ::simList[i].GetEquilibrationParameter(0);
		}
		else {
			fout << "{";
			for (int j = 0; j < ::simList[i].GetEquilibrationParameterNumber(); j++) {
				if (j > 0) {
					fout << ";";
				}
				fout << ::simList[i].GetEquilibrationParameter(j);
			}
			fout << "}";
		}
		fout << "\t" << ::simList[i].GetForceMeanField();
		fout << "\t" << cur_params.Gamma0 << "\t" << cur_params.Gammac << "\t" << cur_params.gammac 
			 << "\t" << cur_params.Kc << "\t" << cur_params.psi << "\t" << cur_params.Glassiness_MF 
			 << "\t" << cur_params.NormDistFromGmax_MF << "\t" << cur_params.CalcAlphaFromNormGlassiness(1)
			 << "\t" << ::simList[i].CountUnphysicalSites() << std::endl;

		if (i == 0) {
			_gK = cur_params.K;
			_gn = cur_params.n;
			_gTau0 = cur_params.Tau0;
			_gOmega = cur_params.Omega;
			_gBeta = cur_params.Beta;
		}
		else {
			if (_gK != cur_params.K || _gn != cur_params.n || _gTau0 != cur_params.Tau0 || _gOmega != cur_params.Omega || _gBeta != cur_params.Beta) {
				_common_params = false;
			}
		}
	}

	fout << "\n\n-------------\n\n";
	if (_common_params == false) {
		fout << "WARNING : parameters are not shared between simulations. What follows only refers to the first one\n\n";
	}

	check_norm_mfparams(&(::simList[0]), &fout);

	fout.close();
}