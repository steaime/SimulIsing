#include "Simulation.h"

Simulation::Simulation()
{
	_ID = rand();
	_status = SimStatus::NOT_INITIALIZED;
	_alpha = NULL;
	_rates = NULL;
	_target_rates = NULL;
	_rate_step = NULL;
	_nn_address = NULL;
	_nn_number = NULL;
}


Simulation::Simulation(int n_sites, bool use_pbc, double omega, double K, double tau0,
				double alpha, double beta, int alpha_noise_type, double* alpha_noise_params,
				int init_noise_type, double* init_noise_params, double effective_T)
{
	_ID = rand();
	_status = SimStatus::NOT_INITIALIZED;

	IsingParameters tmp_params;
	for (int i = 0; i < N_DIM; i++) {
		tmp_params.LatticeShape[i] = n_sites;
	}
	tmp_params.UsePBC = use_pbc;
	tmp_params.K = K;
	tmp_params.Tau0 = tau0;
	tmp_params.AlphaAverage = alpha;
	tmp_params.Beta = beta;
	tmp_params.Omega = omega;
	Noise alpha_noise(alpha_noise_type, alpha_noise_params, alpha);
	tmp_params.AlphaNoise = alpha_noise;

	Noise init_noise(init_noise_type, init_noise_params);
	MCParams tstep_params(effective_T);

	Setup(tmp_params, init_noise, tstep_params);
}

Simulation::Simulation(const IsingParameters &model_parameters, const Noise &init_conditions, const MCParams &tstep_params)
{
	_ID = rand();
	_status = SimStatus::NOT_INITIALIZED;
	Setup(model_parameters, init_conditions, tstep_params);
}

void Simulation::Setup(const IsingParameters &model_parameters, const Noise &init_conditions, const MCParams &tstep_params, bool force_meanfield) {
	_force_meanfield = force_meanfield;
	_init_noise.CopyFrom(init_conditions, true);/// true ###
	_mc_params.CopyFrom(tstep_params);
	UpdateParams(model_parameters);
#if DEBUG_MODE:
	std::cout << "   ... initializing simulation (Simulation::Initialize())" << std::endl;
#endif
	Initialize();
#if DEBUG_MODE:
	std::cout << "       Simulation::Initialize() ended correctly!" << std::endl;
#endif
}

Simulation::~Simulation()
{
	if (_status == SimStatus::INITIALIZED) {
		Clear();
	}
}

void Simulation::Clear() {
	for (int i = 0; i < _params.NumSites(); i++) {
		delete[] _alpha[i];
		delete[] _nn_address[i];
		_alpha[i] = NULL;
		_nn_address[i] = NULL;
		if (_normglassiness != NULL) {
			delete[] _normglassiness[i];
			_normglassiness[i] = NULL;
		}
	}
	delete[] _alpha;
	delete[] _rates;
	delete[] _target_rates;
	delete[] _rate_step;
	delete[] _nn_address;
	delete[] _nn_number;
	_alpha = NULL;
	_rates = NULL;
	_target_rates = NULL;
	_rate_step = NULL;
	_nn_address = NULL;
	_nn_number = NULL;
	if (_normglassiness != NULL) {
		delete[] _normglassiness;
		_normglassiness = NULL;
	}
	_init_noise.Reset();
	_params.AlphaNoise.Reset();
	_status = SimStatus::NOT_INITIALIZED;
}

bool Simulation::GetParam_UsePBC() const {
	return _params.UsePBC;
}

double Simulation::GetParam_K() const
{
	return _params.K;
}

double Simulation::GetParam_n() const
{
	return _params.n;
}

double Simulation::GetParam_Tau0() const
{
	return _params.Tau0;
}

double Simulation::GetParam_Omega() const
{
	return _params.Omega;
}

double Simulation::GetParam_Beta() const
{
	return _params.Beta;
}

double Simulation::GetParam_AlphaAvg() const
{
	return _params.AlphaAverage;
}

int Simulation::GetParam_NumSites(bool force_calc) const
{
	return _params.NumSites(force_calc);
}

double Simulation::GetDerivedParam_CritGamma() const
{
	return _params.CalcCriticalGamma();
}

double Simulation::GetDerivedParam_PsiMin() const
{
	return _params.CalcPsiMin();
}

double Simulation::GetDerivedParam_GMax() const
{
	return _params.CalcGMax();
}

double Simulation::GetDerivedParam_PsiFromGNorm(double norm_g) const
{
	return _params.CalcPsiFromNormGlassiness(norm_g);
}

double Simulation::GetDerivedParam_AlphaFromGNorm(double norm_g) const
{
	return _params.CalcAlphaFromNormGlassiness(norm_g);
}

double Simulation::GetDerivedParam_CritKFromAlpha(double alpha) const
{
	return _params.CalcCriticalKFromAlpha(alpha);
}

double Simulation::GetDerivedParam_CritStrainFromAlpha(double alpha) const
{
	return _params.CalcCriticalStrainFromAlpha(alpha);
}

double Simulation::GetDerivedParam_YieldStrainFromAlpha(double alpha) const
{
	return _params.CalcYieldStrainFromAlpha(alpha);
}

void Simulation::Initialize()
{
#if DEBUG_MODE:
	std::cout << "   ... updating nearest neighbors (Simulation::UpdateNearestNeighbors())" << std::endl;
#endif
	UpdateNearestNeighbors();
#if DEBUG_MODE:
	std::cout << "   ... initializing alpha_ij (Simulation::InitializeAlpha())" << std::endl;
#endif
	InitializeAlpha();
#if DEBUG_MODE:
	std::cout << "   ... generating initial conditions (Simulation::GenInitialConditions())" << std::endl;
#endif
	GenInitialConditions();
#if DEBUG_MODE:
	std::cout << "   ... initializing auxiliary variables" << std::endl;
#endif
	_count_iter = 0;
	_rate_step = new double[_params.NumSites()];
	_status = SimStatus::INITIALIZED;
#if DEBUG_MODE:
	std::cout << "   ... ok" << std::endl;
#endif
}

void Simulation::UpdateNearestNeighbors() {
	UpdateNNNumber();
	UpdateNNAddress();
}

void Simulation::UpdateNNNumber() {
	int num_sites = _params.NumSites();
	_nn_number = new int[num_sites];
	for (int i = 0; i < num_sites; i++) {
		GetNNNumber(i, true);
	}
}

void Simulation::UpdateNNAddress() {
	int num_sites = _params.NumSites();
	_nn_address = new int*[num_sites];
	for (int i = 0; i < num_sites; i++) {
		int cur_nnn = GetNNNumber(i);
		_nn_address[i] = new int[cur_nnn];
		GetNNAddress(i, &_nn_address[i], &cur_nnn, true);
	}
}

void Simulation::InitializeAlpha() {
	double _min_pos_alpha = -1;
	double max_alpha = BAD_VALUE;
	_alpha = new double*[_params.NumSites()];
	_params.AlphaNoise.Reset();
	if (_params.AlphaNoise.GetType() == NoiseTypes::SORTED_LIST && _params.SiteBasedNoise) {
		std::cout << "WARNING: AlphaNoise of type SORTED_LIST is not compatible with site-based noise!" << std::endl;
	}
	if (_params.NoiseOnGlassiness) {
		_normglassiness = new double* [_params.NumSites()];
	}
	for (int i = 0; i < _params.NumSites(); i++) {
		_alpha[i] = new double[_nn_number[i]];
		if (_params.NoiseOnGlassiness) {
			_normglassiness[i] = new double[_nn_number[i]];
			if (_params.SiteBasedNoise) {
				_params.AlphaNoise.Sample(1, _normglassiness[i]);
				for (int j = 1; j < _nn_number[i]; j++) {
					_normglassiness[i][j] = _normglassiness[i][0];
				}
			}
			else {
				_params.AlphaNoise.Sample(_nn_number[i], _normglassiness[i]);
			}
			for (int j = 0; j < _nn_number[i]; j++) {
				// What we sampled was actually 1-g/gmax
				_normglassiness[i][j] = 1 - _normglassiness[i][j];
				// This piece of code is in the wrong place: validation moved to ValidateAlphas() function
#if 0
				if (_normglassiness[i][j] >= 1 - GLASSY_EPS) {
					double dummy = -1;
					for (int iatt = 0; iatt < MAX_RANDOM_ATTEMPTS; iatt++) {
						_params.AlphaNoise.Sample(1, &dummy);
						if (dummy > GLASSY_EPS) {
							_normglassiness[i][j] = 1 - dummy;
							break;
						}
					}
				}
#endif
				_alpha[i][j] = _params.CalcAlphaFromNormGlassiness(_normglassiness[i][j]);
#if CLIP_ALPHA_POSMIN
				if (_min_pos_alpha < 0 || (_alpha[i][j] > 0 && _alpha[i][j] < _min_pos_alpha)) {
					_min_pos_alpha = _alpha[i][j];
				}
#endif
			}
		}
		else {
			if (_params.SiteBasedNoise) {
				_params.AlphaNoise.Sample(1, _alpha[i]);
				for (int j = 1; j < _nn_number[i]; j++) {
					_alpha[i][j] = _alpha[i][0];
				}
			}
			else {
				_params.AlphaNoise.Sample(_nn_number[i], _alpha[i]);
			}
		}
	}
#if DEBUG_MODE
	std::cout << "Generated random Alpha distribution with average " << GetAverageAlpha() << " and variance " << GetAlphaVariance() << std::endl;
	std::cout << "(target average " << _params.AlphaAverage << " and variance ";
	if (_params.AlphaNoise.GetType() == NoiseTypes::WHITE_GAUSS || _params.AlphaNoise.GetType() == NoiseTypes::WHITE_LOGNORM) {
		std::cout << _params.AlphaAverage * _params.AlphaAverage * _params.AlphaNoise.GetParameter(0) << ")" << std::endl;
	}
	else std::cout << "ND" << ")" << std::endl;
	std::cout << "Now proceeding with validation step..." << std::endl;
#endif
	if (_params.ForbidUnphysicalAlpha) {
		max_alpha = _params.CalcAlphaFromNormGlassiness(1.0 - GLASSY_EPS);
	}
	if (_min_pos_alpha < 0) {
		_min_pos_alpha = ALPHA_EPS;
	}
	int nbad = ValidateAlphas(_min_pos_alpha, max_alpha);
#if DEBUG_MODE
	std::cout << "After validation Alpha distribution has average " << GetAverageAlpha() << " and variance " << GetAlphaVariance() 
		<< " (" << nbad << " values clipped)" << std::endl;
#endif
}

void Simulation::GenInitialConditions() {
	_rates = new double[_params.NumSites()];
	_init_noise.Reset();
	_init_noise.Sample(_params.NumSites(), _rates);
	ValidateRates();
	GetAverageRate(true);
	UpdateTargetRates();
}

int Simulation::ValidateRates() {
	int res = 0;
	for (int i = 0; i < _params.NumSites(); i++) {
		if (_rates[i] <= 0) {
			_rates[i] = RATE_EPS;
			res++;
		}
	}
	return res;
}

int Simulation::ValidateAlphas(double min_val, double max_val) {
	// Make sure that all alphas are greater than min_val and that alpha_ij == alpha_ji
	int redrawn = 0;
	int clipped = 0;
	for (int i = 0; i < _params.NumSites(); i++) {
		for (int j = 0; j < _nn_number[i]; j++) {
			int _cur_nnadd = _nn_address[i][j];
			if (_cur_nnadd < i) {
				int _idx = -1;
				for (int k = 0; k < _nn_number[_cur_nnadd]; k++) {
					if (_nn_address[_cur_nnadd][k] == i) {
						_idx = k;
						break;
					}
				}
				if (_params.SiteBasedNoise) {
					double cur_avg = CalcAverage(_alpha[i][j], _alpha[_cur_nnadd][_idx], _params.SBNoiseAvgType);
					if (cur_avg != cur_avg) { // check for NaNs (may occur in harmonic or geometric averages...)
						cur_avg = _alpha[i][j];
					}
					_alpha[i][j] = cur_avg;
					_alpha[_cur_nnadd][_idx] = cur_avg;
				}
				else {
					_alpha[i][j] = _alpha[_cur_nnadd][_idx];
				}

				if (_alpha[i][j] < min_val || (max_val != BAD_VALUE && _alpha[i][j] >= max_val)) {
					redrawn++;
					for (int iatt = 0; iatt < MAX_RANDOM_ATTEMPTS; iatt++) {
						_params.AlphaNoise.Sample(1, &_alpha[i][j]);
						if (_alpha[i][j] >= 0 && _alpha[i][j] < max_val) {
							break;
						}
					}
					if (max_val != BAD_VALUE && _alpha[i][j] >= max_val) {
						_alpha[i][j] = max_val;
						clipped++;
					}
					else if (_alpha[i][j] < min_val) {
						_alpha[i][j] = min_val;
						clipped++;
					}
					_alpha[_cur_nnadd][_idx] = _alpha[i][j];
				}
			}
		}
	}
#if DEBUG_MODE
	int err = 0;
	int nosymm = 0;
	int nan = 0;
	for (int i = 0; i < _params.NumSites(); i++) {
		for (int j = 0; j < _nn_number[i]; j++) {
			int _cur_nnadd_test = _nn_address[i][j];
			int _idx_test = -1;
			for (int k = 0; k < _nn_number[_cur_nnadd_test]; k++) {
				if (_nn_address[_cur_nnadd_test][k] == i) {
					_idx_test = k;
					break;
				}
			}
			if (_alpha[i][j] != _alpha[_cur_nnadd_test][_idx_test]) {
				nosymm++;
			}
			if (_alpha[i][j] <= min_val) {
				_alpha[i][j] = min_val;
				err++;
			}
			else if (max_val != BAD_VALUE && _alpha[i][j] >= max_val) {
				_alpha[i][j] = min_val;
				err++;
			}
			else if (_alpha[i][j] != _alpha[i][j]) {
				nan++;
			}
		}
	}
	std::cout << "Alpha validation step: " << redrawn << " samples redrawn, " << clipped << " clipped (check: " << err/2 << "), " << nosymm 
		<< " asymmetric, " << nan << " NaNs" << std::endl;
#endif
	return redrawn;
}


void Simulation::CopyParams(IsingParameters * copyTo) const
{
	copyTo->CopyFrom(_params);
}

void Simulation::CopyLatticeShape(int * copyTo) const
{
	for (int i = 0; i < N_DIM; i++) {
		copyTo[i] = _params.LatticeShape[i];
	}
}

/*Noise Simulation::GetInitNoise() const  {
	return _init_noise;
}*/

/*IsingParameters Simulation::GetParams() const
{
	IsingParameters res;
	res.CopyFrom(_params);
	return res;
}*/

int Simulation::GetNumSites(bool force_calc) const {
	return _params.NumSites(force_calc);
}

int Simulation::GetNumSlowSites(double thr_rate_ratio) const {
	int res = 0;
	for (int i = 0; i < _params.NumSites(); i++) {
		if (IsSlowSite(i, thr_rate_ratio)) {
			res++;
		}
	}
	return res;
}

bool Simulation::IsFastSite(int site, double thr_rate_ratio) const {
	return !IsSlowSite(site, thr_rate_ratio);
}

bool Simulation::IsSlowSite(int site, double thr_rate_ratio) const {
	return (_rates[site] < thr_rate_ratio / _params.Tau0);
}

int Simulation::GetNumFastSites(double thr_rate_ratio) const {
	return _params.NumSites() - GetNumSlowSites(thr_rate_ratio);
}

double Simulation::GetFractionSlowSites(double thr_rate_ratio) const {
	double res = GetNumSlowSites(thr_rate_ratio) * 1.0 / _params.NumSites();
	return res;
}

double Simulation::GetFractionFastSites(double thr_rate_ratio) const {
	double res = 1.0 - GetFractionSlowSites(thr_rate_ratio);
	return res;
}


void Simulation::GetSiteCoordinates(int site, int* res) const {
	_params.GetSiteCoordinates(site, res);
}

int Simulation::GetRateCorrelation(double* res, int* x, int num_x)
{
#if N_DIM==1:
	return Correlate1D(_rates, _params.NumSites(), res, x, num_x, _params.UsePBC);
#elif N_DIM==2:
	return Correlate2D(_rates, _params.LatticeShape, res, x, num_x, _params.UsePBC);
#elif N_DIM==3:
	return Correlate3D(_rates, _params.NumSites(), res, x, num_x, _params.UsePBC);
#else:
	return -1;
#endif
}

int Simulation::GetStatus() {
	return _status;
}

int Simulation::ID() {
	return _ID;
}

int Simulation::StepCount() {
	return _count_iter;
}

void Simulation::CopyInitNoise(Noise * copyTo) const
{
	copyTo->CopyFrom(_init_noise);
}

int Simulation::UpdateParams(const IsingParameters &new_params)
{
	if (ValidateParams(new_params)) {
		_params.CopyFrom(new_params); //// true for deep copy ###
		return 1;
	}
	else {
		return -1;
	}
}

bool Simulation::ValidateParams(const IsingParameters &params) const
{
	bool bln_good = true;
	for (int i = 0; i < N_DIM; i++) {
		if (params.LatticeShape[i] <= 0) {
			bln_good = false;
		}
	}
	if (params.AlphaAverage < 0) {
		bln_good = false;
	}
	if (params.Beta <= 0) {
		bln_good = false;
	}
	if (params.K < 0) {
		bln_good = false;
	}
	if (params.Omega <= 0) {
		bln_good = false;
	}
	if (params.Tau0 <= 0) {
		bln_good = false;
	}
	return bln_good;
}

double Simulation::GetGamma() const {
	return _gamma;
}

int Simulation::SetGamma(double new_gamma, bool adjust_lattice) {
	if (new_gamma >= 0) {
		_gamma = new_gamma;
		UpdateTargetRates();
		if (adjust_lattice) {
			return AdjustLattice();
		}
		else 
		{
			return -1;
		}
	}
	else
	{
		return -2;
	}
}

double Simulation::GetGMax()
{
	return _params.GetGMax();
}

double Simulation::GetAlphaFromNormGlassiness(double norm_g)
{
	return _params.CalcAlphaFromNormGlassiness(norm_g);
}

double Simulation::GetMeanFieldGamma(double rate) const
{
	return _params.GetMeanFieldGamma(rate);
}

double Simulation::GetGlassinessFromAlpha(double alpha)
{
	return _params.GetGlassinessFromAlpha(alpha);
}

double Simulation::GetNormGlassinessFromAlpha(double alpha)
{
	return _params.GetNormGlassinessFromAlpha(alpha);
}

double Simulation::GetDistFromGmax(double alpha)
{
	return _params.GetDistFromGmax(alpha);
}

double Simulation::GetYieldStrainFromAlpha(double alpha)
{
	return _params.GetYieldStrainFromAlpha(alpha);
}

int Simulation::CountUnphysicalSites()
{
	int res = 0;
#if DEBUG_MODE
	int check_res = 0;
#endif
	double max_alpha = _params.CalcAlphaFromNormGlassiness(1);
	for (int i = 0; i < _params.NumSites(); i++) {
		bool is_unphysical = false;
		for (int j = 0; j < _nn_number[i]; j++) {
			if (_alpha[i][j] > max_alpha) {
				res++;
#if DEBUG_MODE==0
				break;
#endif
			}
#if DEBUG_MODE
			if (_params.GetNormGlassinessFromAlpha(_alpha[i][j]) > 1) {
				check_res++;
				break;
			}
		}
	}
	if (check_res != res) {
		std::cout << "\n\n\n\n-----------\nWARNING: Simulation.CountUnphysicalSites() got inconsistent results:"
			<< res << " (alpha output) vs " << check_res << " (g output). Check equations used!!!\n-----------\n\n\n\n" << std::endl;
	}
#else
		}
	}
#endif
	return res;
}

int Simulation::IsSiteOnEdge(int site) const {
	if (site >= 0 && site < _params.NumSites()) {
		if (_params.NumSites() == 1) {
			return N_DIM;
		}
		else {
			int res = 0;
			int site_coord[N_DIM];
			_params.GetSiteCoordinates(site, site_coord);
			for (int i = 0; i < N_DIM; i++) {
				if (site_coord[i] == 0) {
					res++;
				}
				if (site_coord[i] == _params.LatticeShape[i] - 1) {
					res++;
				}
			}
			return res;
		}
	}
	else {
		return -1;
	}
}

void Simulation::FindEdgesSite(int site, int* axes) const {
	if (site >= 0 && site < _params.NumSites()) {
		if (_params.NumSites() == 1) {
			for (int i = 0; i < N_DIM; i++) {
				axes[i] = 2;
			}
		}
		else {
			int res = 0;
			int site_coord[N_DIM];
			_params.GetSiteCoordinates(site, site_coord);
			for (int i = 0; i < N_DIM; i++) {
				axes[i] = 0;
				if (site_coord[i] == 0) {
					axes[i]++;
				}
				if (site_coord[i] == _params.LatticeShape[i] - 1) {
					axes[i]++;
				}
			}
		}
	}
	else {
		for (int i = 0; i < N_DIM; i++) {
			axes[i] = -1;
		}
	}
}


int Simulation::GetNNNumber(int site, bool force_calc) const {
	if (site >= 0 && site < _params.NumSites()) {
		if (force_calc) {
			int res = 0;
			int edges[N_DIM];
			if (_params.UsePBC == false) {
				FindEdgesSite(site, edges);
			}
			for (int i = 0; i < N_DIM; i++) {
				if (_params.LatticeShape[i] > 2) {
					if (_params.UsePBC) {
						res += 2;
					}
					else {
						res += 2 - edges[i];
					}
				}
				else if (_params.LatticeShape[i] > 1) {
					res += 1;
				}
				else {
					res += 0;
				}
			}
			_nn_number[site] = res;
			return res;
		}
		else
		{
			return _nn_number[site];
		}
	}
	else {
		return -1;
	}
}

void Simulation::GetNNAddress(int site, int** address, int* number, bool force_calc) {
	if (site >= 0 && site < _params.NumSites()) {
		int num_addresses = GetNNNumber(site);
		if (number != NULL) {
			*number = GetNNNumber(site);
		}
		if (force_calc) {
			if (num_addresses > 0) {
				//address = new int[num_addresses];
#if COORDS_IN_HEAP
				int* site_coord = new int[N_DIM];
#else
				int site_coord[N_DIM];
#endif
				int cur_idx = 0;
				_params.GetSiteCoordinates(site, site_coord);
				for (int i = 0; i < N_DIM; i++) {
					int new_coord_p[N_DIM];
					int new_coord_m[N_DIM];
					for (int j = 0; j < N_DIM; j++) {
						if (i == j) {
							new_coord_m[j] = site_coord[j] - 1;
							new_coord_p[j] = site_coord[j] + 1;
						}
						else {
							new_coord_m[j] = site_coord[j];
							new_coord_p[j] = site_coord[j];
						}
					}
					if (new_coord_m[i] < 0 && _params.UsePBC && _params.LatticeShape[i] > 2) {
						new_coord_m[i] = _params.LatticeShape[i] - 1;
					}
					if (new_coord_m[i] >= 0) {
						(*address)[cur_idx] = _params.SiteID(new_coord_m);
						cur_idx++;
					}
					if (new_coord_p[i] >= _params.LatticeShape[i] && _params.UsePBC && _params.LatticeShape[i] > 2) {
						new_coord_p[i] = 0;
					}
					if (new_coord_p[i] < _params.LatticeShape[i]) {
						(*address)[cur_idx] = _params.SiteID(new_coord_p);
						cur_idx++;
					}
				}
				_nn_address[site] = *address;
#if COORDS_IN_HEAP
				delete[] site_coord;
				site_coord = NULL;
#endif
			}
		}
		else
		{
			*address = _nn_address[site];
		}
	}
	else
	{
		if (number != NULL) {
			*number = -1;
		}
	}
}

double** Simulation::GetAlphas() const {
	return _alpha;
}

int Simulation::CountAlphas() const
{
	int iCount = 0;
	for (int i = 0; i < _params.NumSites(); i++) {
		for (int j = 0; j < _nn_number[i]; j++) {
			iCount++;
		}
	}
	return iCount;
}

int Simulation::FlattenAlphas(double* res) const
{
	return Flatten2DArray(_alpha, _params.NumSites(), _nn_number, res);
}

double* Simulation::GetRates() const {
	return _rates;
}

double* Simulation::GetTargetRates(bool update) {
	if (update) {
		UpdateTargetRates();
	}
	return _target_rates;
}

void Simulation::UpdateTargetRates() {
	if (_status == SimStatus::NOT_INITIALIZED) {
		_target_rates = new double[_params.NumSites()];
	}
	for (int i = 0; i < _params.NumSites(); i++) {
		UpdateTargetRate(i);
	}
}

void Simulation::UpdateTargetRate(int site) {
	if (_gamma > 0) {
		double dblCoupl = 0;
		if (_nn_number[site] == 0) {
			dblCoupl = _params.AlphaAverage * pow(_params.Omega / _avg_rate, _params.Beta);
		}
		else if (_force_meanfield) {
#if MEANFIELD_USE_LOCALALPHA
			dblCoupl = GetSiteAveragedAlpha(site) * pow(_params.Omega / _avg_rate, _params.Beta);
#else
			dblCoupl = _params.AlphaAverage * pow(_params.Omega / _avg_rate, _params.Beta);
#endif
		}
		else {
			if (_params.UseNewEqn) {
				for (int i = 0; i < _nn_number[site]; i++) {
					dblCoupl += GetAlphaNN(site, i) * pow(_params.Omega, 2.0) / (_rates[site] * _rates[_nn_address[site][i]]); // New equation
				}
			}
			else {
				for (int i = 0; i < _nn_number[site]; i++) {
					dblCoupl += GetAlphaNN(site, i) * pow(_params.Omega / _rates[_nn_address[site][i]], _params.Beta); // Domenico's eq. 19
					// Works as well, but it's much slower:
					// dblCoupl += GetAlpha(site, _nn_address[site][i]) * pow(_params.Omega / _rates[_nn_address[site][i]], _params.Beta); // Domenico's eq. 19
				}
			}
			dblCoupl /= _nn_number[site];
		}
		double dblSingle = 1.0 / pow(_gamma, _params.n);
		_target_rates[site] = 1.0 / _params.Tau0 + _params.Omega * _params.K / (dblSingle + dblCoupl); // Domenico's eq. 19
	}
	else
	{
		_target_rates[site] = 1.0 / _params.Tau0;
	}

}

double Simulation::GetAlpha(int site_i, int site_j) const {
	int _nn_idx = -1;
	for (int i = 0; i < _nn_number[site_i]; i++) {
		if (_nn_address[site_i][i] == site_j) {
			_nn_idx = i;
			break;
		}
	}
#if 0
	int _nn2_idx = -1;
	for (int i = 0; i < _nn_number[site_j]; i++) {
		if (_nn_address[site_j][i] == site_i) {
			_nn2_idx = i;
			break;
		}
	}
	if (_alpha[site_i][_nn_idx] != _alpha[site_j][_nn2_idx]) {
		std::cout << "ERROR: alpha(" << site_i << "," << site_j << ") != alpha (" << site_j << "," << site_i << ")" << std::endl;
	}
#endif
	return _alpha[site_i][_nn_idx];
}

double Simulation::GetSiteAveragedAlpha(int site_index) const
{
	if (_nn_number[site_index] == 0) {
		return _params.AlphaAverage;
	}
	else {
		double alpha_res = 0;
		for (int i = 0; i < _nn_number[site_index]; i++) {
			alpha_res += GetAlphaNN(site_index, i);
		}
		return alpha_res / _nn_number[site_index];
	}
}

double Simulation::GetAlphaNN(int site_i, int neighbor_index) const {
	return _alpha[site_i][neighbor_index];
}

double Simulation::GetRate(int site) const {
	return _rates[site];
}

double Simulation::GetTargetRate(int site, bool force_calc) {
	if (force_calc) {
		UpdateTargetRate(site);
	}
	return _target_rates[site];
}

double Simulation::GetAverageRate(bool force_calc) {
	if (force_calc) {
		_avg_rate = CalcAverage(_rates, _params.NumSites(), _params.RateAvgType);
	}
	return _avg_rate;
}

double Simulation::GetRateStd() const {
	return CalcStd(_rates, _params.NumSites());
}

double Simulation::GetTauStd(bool statistic) {
	if (statistic) {
		return GetRateStd() / (pow(GetAverageRate(), 2));
	}
	else {
		double* taus;
		taus = new double[_params.NumSites()];
		for (int i = 0; i < _params.NumSites(); i++) {
			taus[i] = 1.0 / _rates[i];
		}
		double res = CalcStd(taus, _params.NumSites());
		delete[] taus;
		taus = NULL;
		return res;
	}
}

double Simulation::GetAverageAlpha() const {
	return CalcAverage(_alpha, _params.NumSites(), _nn_number);
}

double Simulation::GetAlphaVariance() const
{
	return CalcVariance(_alpha, _params.NumSites(), _nn_number);
}

double Simulation::GetTemperature() const
{
	return _mc_params.GetTemperature();
}

int Simulation::GetMaxStepNum() const {
	return _mc_params.GetMaxIterations()*_mc_params.GetEquilibrationRunNumber();
}

int Simulation::GetMaxStepPerRun() const {
	return _mc_params.GetMaxIterations();
}

int Simulation::SetEquilibrationRunNumber(int new_val) {
	return _mc_params.SetEquilibrationRunNumber(new_val);
}

bool Simulation::_UseMeanfield(int site) {
	if (site >= 0) {
		return (_force_meanfield || (_nn_number[site] == 0));
	}
	else {
		return (_force_meanfield || (_params.NumSites() <= 1));
	}
}

int Simulation::SimStep(bool ForceUpdateAvRate, double* rate_change) {
	int res;
	_mc_params.RandomJumps(_rate_step, _target_rates, _params.NumSites());
	bool bln_accept = true; // Obsolete variable. In the original version there was
							// an extra Montecarlo accept/reject decision
	if (bln_accept) {
		if (rate_change != NULL) {
			for (int i = 0; i < _params.NumSites(); i++) {
				rate_change[i] = _rate_step[i] - _rates[i];
			}
		}
		for (int i = 0; i < _params.NumSites(); i++) {
			_rates[i] = _rate_step[i];
		}
		if (_UseMeanfield() || ForceUpdateAvRate || DEBUG_MODE) {
			GetAverageRate(true);
		}
		UpdateTargetRates();
		res = 1;
	}
	else {
		if (rate_change != NULL) {
			for (int i = 0; i < _params.NumSites(); i++) {
				rate_change[i] = 0;
			}
		}
		res = 0;
	}
	_count_iter++;
#if DEBUG_MODE
	std::cout << _gamma << " \t| " << _count_iter << " \t| " << GetAverageRate() << " \t| " << GetRateStd() << std::endl;
#endif
	return res;
}

int Simulation::GetEquilibrationProtocol() const {
	return _mc_params.GetTempProtocol();
}

int Simulation::GetEquilibrationParameterNumber() const {
	return _mc_params.TempProtocolNumParams();
}

/*double* Simulation::GetEquilibrationParameters() const {
	return _mc_params.GetTempProtocolParams();
}*/

double Simulation::GetEquilibrationParameter(int paramIdx) const
{
	return _mc_params.GetTempProtocolParam(paramIdx);
}

int Simulation::SetEquilibrationProtocol(int newProtocol, double* protocolParams) {
	return _mc_params.SetTempProtocol(newProtocol, protocolParams);
}


int Simulation::AdjustLattice(double** eqDetails, bool** blnFastSites) {
	double* ratechange;
#if VERBOSE
	std::cout << "Now adjusting lattice to external strain gamma0=" << _gamma << std::endl;
	std::cout << "_rates[0] - _target_rates[0] - deviation" << std::endl;
#endif
	if (eqDetails != NULL) {
		ratechange = new double[_params.NumSites()];
	}
	else {
		ratechange = NULL;
	}
	int cumul_nloops = 0;
	for (int i = 0; i < _mc_params.GetEquilibrationRunNumber(); i++) {
		int nloops = 0;
		while (_mc_params.Equilibrated(_rates, _target_rates, _params.NumSites())==false)
		{
			_mc_params.SetTemperatureStep(nloops);

			if (eqDetails != NULL) {
				eqDetails[EqDetails::AVRATE][cumul_nloops] = GetAverageRate();
				eqDetails[EqDetails::STDRATE][cumul_nloops] = GetRateStd();
				eqDetails[EqDetails::FASTSITES][cumul_nloops] = GetFractionFastSites();
				eqDetails[EqDetails::DEV][cumul_nloops] = _mc_params.GetDeviation(_rates, _target_rates, _params.NumSites());
			}

			/////////
			SimStep(eqDetails != NULL, ratechange);
			/////////

			if (eqDetails != NULL) {
				eqDetails[EqDetails::RATECHANGE][cumul_nloops] = sqrt(CalcAvgPower(ratechange, _params.NumSites(), 2.0));
			}
			if (blnFastSites != NULL) {
				for (int j = 0; j < _params.NumSites(); j++) {
					blnFastSites[cumul_nloops][j] = IsFastSite(j);
				}
			}

			nloops++;
			cumul_nloops++;
#if VERBOSE
			std::cout << nloops << " - " << _rates[0] << " - " << _target_rates[0] << " - " << _mc_params.GetDeviation(_rates, _target_rates, _params.NumSites()) << std::endl;
#endif
			if (nloops >= _mc_params.GetMaxIterations()) {
				break;
			}
		}
	}
	GetAverageRate(true);
	if (eqDetails != NULL) {
		delete[] ratechange;
		ratechange = NULL;
	}
	return cumul_nloops;
}

bool Simulation::GetForceMeanField() const {
	return _force_meanfield;
}
void Simulation::SetForceMeanField(bool val) {
	_force_meanfield = val;
	GetAverageRate(true);
}

double Simulation::GetMeanFieldRate(double gamma, double start_rate, double precision) {
	if (gamma < 0) {
		gamma = _gamma;
	}
	double cur_rate = start_rate;
	if (cur_rate <= 0) {
		cur_rate = GetAverageRate();
	}
	double dblSingle = 1.0 / pow(gamma, _params.n);
	double new_rate;
	double dblCoupl;
	double rel_diff;
	do {
		dblCoupl = _params.AlphaAverage * pow(_params.Omega / cur_rate, _params.Beta);
		new_rate = 1.0 / _params.Tau0 + _params.Omega * _params.K / (dblSingle + dblCoupl);
		rel_diff = (new_rate - cur_rate) / new_rate;
		cur_rate = new_rate;
	} while (abs(rel_diff) > precision);
	return new_rate;
}
