#include "TestFunctions.h"


int check_indexing(IsingParameters my_params, int test_site) {
	std::cout << "\n\n" << N_DIM << "D Ising lattice created, lattice shape [";
	for (int i = 0; i < N_DIM; i++) {
		std::cout << my_params.LatticeShape[i] << ";";
	}
	std::cout << "], Total number of sites: " << my_params.NumSites() << std::endl;
	std::cout << "Ising parameters: K=" << my_params.K << "; Omega=" << my_params.Omega << "; Beta="
		<< my_params.Beta << "; n=" << my_params.n << "; tau0=" << my_params.Tau0
		<< "; alpha(avg)=" << my_params.AlphaAverage << "; usePBC=" << my_params.UsePBC << std::endl;
	int coords[N_DIM];
	my_params.GetSiteCoordinates(test_site, coords);
	int site_id = my_params.SiteID(coords);
	int area = my_params.ProjectionArea(0);
	std::cout << "Test indexing (site " << test_site << ")... Coordinates: [";
	for (int i = 0; i < N_DIM; i++) {
		std::cout << coords[i] << ";";
	}
	std::cout << "]. Site ID: " << site_id << ". Section area for axis 0: " << area << std::endl;
	return site_id;
}

void check_neighbors(Simulation* my_sim, int* sites, int num_sites) {
	std::cout << "Check nearest neighbors:" << std::endl;
	for (int i = 0; i < num_sites; i++) {
		std::cout << "#" << sites[i] << " - coords [";
		int coords[N_DIM];
		my_sim->GetSiteCoordinates(sites[i], coords);
		for (int j = 0; j < N_DIM; j++) {
			std::cout << coords[j] << ";";
		}
		int nnnum = my_sim->GetNNNumber(sites[i]);
		int** nnadd = new int*[nnnum];
		my_sim->GetNNAddress(sites[i], nnadd, &nnnum);
		std::cout << "] " << nnnum << " neighbors: ";
		for (int j = 0; j < nnnum; j++) {
			std::cout << "[";
			for (int k = 0; k < N_DIM; k++) {
				std::cout << nnadd[j][k] << ";";
			}
			std::cout << "]; ";
		}
		std::cout << std::endl;
		for (int j = 0; j < nnnum; j++) {
			delete[] nnadd[j];
			nnadd[j] = NULL;
		}
		delete[] nnadd;
		nnadd = NULL;
	}
}

void check_noise(Noise *my_noise) {

	std::cout << "\n\nNOISE DETAILS:\n - average: " << my_noise->GetAverage() << "\n - type: ";
	if (my_noise->GetType() == NoiseTypes::NONOISE) {
		std::cout << "NoiseTypes::NONOISE";
	}
	else if (my_noise->GetType() == NoiseTypes::WHITE_FLAT) {
		std::cout << "NoiseTypes::WHITE (FLAT P)\n - relative variance: " << my_noise->GetParameter(0);
	}
	else if (my_noise->GetType() == NoiseTypes::WHITE_GAUSS) {
		std::cout << "NoiseTypes::WHITE (GAUSS P)\n - relative variance: " << my_noise->GetParameter(0);
	}
	const int len_sample = 100;
	double sample[len_sample];
	my_noise->Sample(len_sample, sample);
	std::cout << "\n - " << len_sample << " elements sample: [";
	for (int i = 0; i < len_sample; i++) {
		std::cout << sample[i] << "; ";
		if (i > 4) {
			std::cout << "...";
			break;
		}
	}
	std::cout << "]\n   - average: " << CalcAverage(sample, len_sample)
		<< "\n   - variance: " << CalcVariance(sample, len_sample)
		<< "\n   - rel variance: " << CalcVariance(sample, len_sample) / pow(CalcAverage(sample, len_sample), 2)
		<< "\n   - stdev: " << CalcStd(sample, len_sample) << std::endl;

	//delete[] sample;
}

void check_meanfield_gamma(IsingParameters my_params, double start_rate, double end_rate, double ppd_rate) {
	std::cout << "\n\nCHECK MEAN FIELD EQ.\n1/tau0 - 1/gamma0" << std::endl;
	double cur_rate = start_rate;
	if (cur_rate < 1.0 / my_params.Tau0) {
		cur_rate = 1.0 / my_params.Tau0;
	}
	while (cur_rate <= end_rate) {
		double cur_gamma = my_params.GetMeanFieldGamma(cur_rate);
		std::cout << cur_rate << " - " << 1.0 / cur_gamma << std::endl;
		cur_rate *= pow(10, 1.0 / ppd_rate);
	}
}

void export_meanfield_gamma(IsingParameters my_params, double* rates, int num_rates, std::string out_file) {
	std::ofstream fout(out_file);
	fout << "gamma0\tGamma[1/s]\tGamma/Gamma0\tGamma/Gammac\t1/gamma0\tgammac/gamma0\t(gammac/gamma0)^n" << std::endl;
	double cur_gamma;
	for (int i = 0; i < num_rates; i++) {
		if (rates[i] > 1.0 / my_params.Tau0) {
			cur_gamma = my_params.GetMeanFieldGamma(rates[i]);
			fout << cur_gamma << "\t" << rates[i] << "\t" << rates[i] * my_params.Tau0 << "\t" << rates[i] / my_params.Gammac 
				 << "\t" << 1.0 / cur_gamma << "\t" << my_params.gammac / cur_gamma << "\t" << pow(my_params.gammac / cur_gamma, my_params.n) << std::endl;
		}
	}
}

void check_meanfield_rate(Simulation* my_sim, double start_gamma, double end_gamma, double ppd_gamma) {
	std::cout << "\n\nSTRAIN SWEEP - MEAN FIELD\nGamma0 - AvgRate" << std::endl;
	double cur_gamma = start_gamma;
	double cur_rate = 0;
	while (cur_gamma < end_gamma) {
		cur_rate = my_sim->GetMeanFieldRate(cur_gamma, cur_rate);
		std::cout << cur_gamma << " - " << cur_rate << std::endl;
		cur_gamma *= pow(10, 1.0 / ppd_gamma);
	}
	while (cur_gamma >= start_gamma) {
		cur_rate = my_sim->GetMeanFieldRate(cur_gamma, cur_rate);
		std::cout << cur_gamma << " - " << cur_rate << std::endl;
		cur_gamma /= pow(10, 1.0 / ppd_gamma);
	}
}

void check_norm_mfparams(Simulation* my_sim, std::ofstream* fout)
{

	double xmin = 1e-8;
	double xmax = 1e3;
	int numpts = LogSpace_CalcNum(xmin, xmax, 20);
	double* dst_from_gnorm = new double[numpts];
	LogSpaceNum(xmin, xmax, numpts, dst_from_gnorm);

	IsingParameters my_param;
	my_sim->CopyParams(&my_param);
	my_param.UpdateDerivedParams();

	*fout << "COSTANT BASE PARAMETERS:\nK\t= " << my_param.K << "\nn\t= " << my_param.n << "\nTau0\t= " << my_param.Tau0 << "\nOmega\t= "
		<< my_param.Omega << "\nBeta\t= " << my_param.Beta << "\nCONSTANT DERIVED PARAMETERS:\nGamma0\t= " << my_param.Gamma0
		<< "\nGammac\t= " << my_param.GetCriticalGamma() << "\npsi_min\t= " << my_param.GetPsiMin() << "\ng_max\t= " << my_param.GetGMax() << std::endl;
#if DEBUG_MODE
	*fout << "\nVARIABLE (ALPHA-DEPENDENT) PARAMETERS:\n1-g/gmax\t1-gnorm(alpha)\tg\tg(alpha)\tpsi(g)\tpsi(alpha)\talpha(psi)\talpha(g)\tKc\tgammac\tgammay(alpha)"
		<< "\tdiff_1-gnorm\tdiff_g\tdiff_psi\tdiff_alpha" << std::endl;
	double tot_sqrdiff = 0;
#else
	*fout << "\nVARIABLE (ALPHA-DEPENDENT) PARAMETERS:\n1-g/gmax\tg\tpsi\talpha\tKc\tgammac\tgammay(alpha)" << std::endl;
#endif
	for (int i = 0; i < numpts; i++) {
		double curx = dst_from_gnorm[i];
		double curg = (1 - curx) * my_param.GetGMax();
		double curp = my_param.GetPsiFromNormGlassiness(1 - curx);
		double cura = my_param.GetAlphaFromNormGlassiness(1 - curx);
#if DEBUG_MODE
		double curx_test = my_param.GetDistFromGmax(cura);
		double curg_test = my_param.GetGlassinessFromAlpha(cura);
		double curp_test = my_param.GetPsiFromAlpha(cura);
		double cura_test = my_param.GetAlphaFromPsi(curp);
		*fout << curx << "\t" << curx_test << "\t" << curg << "\t" << curg_test << "\t" << curp << "\t" << curp_test << "\t" << cura_test << "\t" << cura
			<< "\t" << my_param.GetCriticalKFromAlpha(cura) << "\t" << my_param.GetCriticalStrainFromAlpha(cura) << "\t" << my_param.GetYieldStrainFromAlpha(cura) 
			<< "\t" << curx - curx_test << "\t" << curg - curg_test << "\t" << curp - curp_test << "\t" << cura - cura_test << std::endl;
		tot_sqrdiff += pow(curx - curx_test, 2) + pow(curg - curg_test, 2) + pow(curp - curp_test, 2) + pow(cura - cura_test, 2);
#else
		*fout << curx << "\t" << curg << "\t" << curp << "\t" << cura << "\t" << my_param.GetCriticalKFromAlpha(cura) << "\t" << 
			my_param.GetCriticalStrainFromAlpha(cura) << "\t" << my_param.GetYieldStrainFromAlpha(cura) << std::endl;
#endif
	}
#if DEBUG_MODE
	*fout << "\nTotal squared diff : " << tot_sqrdiff << std::endl;
#endif

	delete[] dst_from_gnorm;
	dst_from_gnorm = NULL;

#if 0
	if (out_fname != "") {
		double dfng_min = 1e-5, dfng_max = 1e2;
		int _ppd = 50;
		int num_data = LogSpace_CalcNum(dfng_min, dfng_max, _ppd);
		double* dfng_arr = new double(num_data+1);
		num_data = LogSpace(dfng_min, dfng_max, _ppd, dfng_arr);
		std::ofstream fout;
		fout.open(out_fname, 'w');
		fout << "CONSTANT PARAMETERS:" << std::endl;
		fout << ""
		fout << "#\t1-gnorm\tgnorm\tg\tpsi\talpha(g)\t" << std::endl;

		fout.close();
		delete[] dfng_arr;
		dfng_arr = NULL;
	}
#endif
}

bool test_roots_P3(double R1, double R2, double R3, double Pref) {
	double coeff[4];
	coeff[0] = Pref;
	coeff[1] = - Pref * (R1 + R2 + R3);
	coeff[2] = Pref * (R1 * R2 + R2 * R3 + R3 * R1);
	coeff[3] = -Pref * R1 * R2 * R3;
	double res[3];
	int nres = RealRoots_P3(coeff[0], coeff[1], coeff[2], coeff[3], res);
	std::cout << "Computing roots of: " << Pref << "(x-" << R1 << ")(x-" << R2 << ")(x-" << R3 << ") ::: ";
	std::cout << "Coefficients " << coeff[0] << "," << coeff[1] << "," << coeff[2] << "," << coeff[3] <<
		" ::: " << nres << " real roots: ";
	for (int i = 0; i < nres; i++) {
		std::cout << res[i] << ",";
	}
	std::cout << std::endl;
	bool check = true;
	if (nres != 3) check = false;
	else {
		for (int i = 0; i < nres; i++) {
			if (res[i] != R1 && res[i] != R2 && res[i] != R3) check = false;
		}
	}
	if (!check) std::cout << "   >>> ERROR! <<<" << std::endl;
	return check;
}

void check_real_roots()
{
	std::cout << "\nCHECK REAL ROOTS:" << std::endl;
	test_roots_P3(1, 1, 1);
	test_roots_P3(1, -1, 4);
	test_roots_P3(1e-2, 1e4, 0);
	test_roots_P3(3.14, 3.141, 3.1415);
	test_roots_P3(-1, 0, -5);
}
