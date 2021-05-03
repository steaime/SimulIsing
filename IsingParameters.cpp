#include "IsingParameters.h"

IsingParameters::IsingParameters()
{
	for (int i = 0; i < N_DIM; i++) {
		LatticeShape[i] = 0;
	}
	UsePBC = true;
	UseNewEqn = true;
	K = 1.0;
	Tau0 = 1.0;
	Omega = 1.0;
	Beta = 1.0;
	AlphaAverage = 0.0;
	// initialize AlphaNoise?
}

IsingParameters::IsingParameters(int* _LatticeShape, double _n, double _K, double _Tau0, double _Omega, double _Beta,
								const Noise &_AlphaNoise, bool _UsePBC, bool _UseNewEqn, int _avgType) {
	_num_sites = 1;
	for (int i = 0; i < N_DIM; i++) {
		LatticeShape[i] = _LatticeShape[i];
		_num_sites *= LatticeShape[i];
	}
	K = _K;
	n = _n;
	Tau0 = _Tau0;
	Omega = _Omega;
	Beta = _Beta;
	UsePBC = _UsePBC;
	UseNewEqn = _UseNewEqn;
	RateAvgType = _avgType;
	NoiseOnGlassiness = false; // WARNING: NoiseOnGlassiness not supported with this constructor!
	AlphaAverage = _AlphaNoise.GetAverage();
	AlphaNoise.CopyFrom(_AlphaNoise, true);
}

IsingParameters::IsingParameters(std::string sParamFile) {
	ConfigParams params(sParamFile);
	Initialize(params);
}

IsingParameters::IsingParameters(const ConfigParams &params) {
	Initialize(params);
}

IsingParameters::~IsingParameters()
{
	if (_cbroots != NULL) {
		delete[] _cbroots;
		_cbroots = NULL;
	}
}

void IsingParameters::Initialize(const ConfigParams &params) {
	int _LatticeShape[INI_MAXNUMDIM];
	_cbroots = new double[4];

	for (int i = 0; i < INI_MAXNUMDIM; i++) {
		_LatticeShape[i] = (int)(params.piLatticeShape[i]); // don't remove this explicit cast.
															// don't know why but it gives weird results without it
	}
	int _nSites = params.iNumSites;
	_num_sites = 1;
	for (int i = 0; i < N_DIM; i++) {
		if (_LatticeShape[i] <= 0) {
			LatticeShape[i] = _nSites;
		}
		else {
			LatticeShape[i] = _LatticeShape[i];
		}
		_num_sites *= LatticeShape[i];
	}

	K = params.K;
	n = params.n;
	Tau0 = params.Tau0;
	Omega = params.Omega;
	Beta = params.Beta;
	UsePBC = params.dUsePBC;
	UseNewEqn = params.bUseNewEqn;
	NoiseOnGlassiness = params.bAlphaFromGlassiness;
	ForbidUnphysicalAlpha = params.bForbidUnphysicalAlpha;
	SiteBasedNoise = params.bSiteBasedNoise;
	SBNoiseAvgType = params.iSBNoiseAvgType;
	RateAvgType = params.iRateAvgType;

	if (NoiseOnGlassiness) {
		AlphaNoise.Initialize(params.iAlphaNoiseType, params.dNormDistFromGmax, params.pdAlphaParams);
		AlphaAverage = CalcAlphaFromNormGlassiness(1 - params.dNormDistFromGmax);
	}
	else {
		AlphaNoise.Initialize(params.iAlphaNoiseType, params.dAlphaAvg, params.pdAlphaParams);
		AlphaAverage = params.dAlphaAvg;
	}

	UpdateDerivedParams();
}

void IsingParameters::UpdateDerivedParams()
{
	Gamma0 = 1.0 / Tau0;
	Gammac = GetCriticalGamma();
	Kc = CalcCriticalKFromAlpha(AlphaAverage);
	gammac = CalcCriticalStrainFromAlpha(AlphaAverage);
	psi_min = GetPsiMin();
	psi = GetPsiFromAlpha(AlphaAverage);
	Glassiness_MF = 1 - psi;
	NormDistFromGmax_MF = (psi - psi_min) / (1 - psi_min);

	_params_initialized = true;
}

void IsingParameters::CopyFrom(const IsingParameters &params, bool deep_copy)
{
	_num_sites = 1;
	for (int i = 0; i < N_DIM; i++) {
		LatticeShape[i] = params.LatticeShape[i];
		_num_sites *= LatticeShape[i];
	}
	UsePBC = params.UsePBC;
	K = params.K;
	n = params.n;
	Tau0 = params.Tau0;
	Omega = params.Omega;
	AlphaAverage = params.AlphaAverage;
	Beta = params.Beta;
	UsePBC = params.UsePBC;
	UseNewEqn = params.UseNewEqn;
	NoiseOnGlassiness = params.NoiseOnGlassiness;
	ForbidUnphysicalAlpha = params.ForbidUnphysicalAlpha;
	SiteBasedNoise = params.SiteBasedNoise;
	SBNoiseAvgType = params.SBNoiseAvgType;
	RateAvgType = params.RateAvgType;
	UpdateDerivedParams();

	if (deep_copy) {
		AlphaNoise.Initialize(params.AlphaNoise.GetType(), params.AlphaNoise.GetAverage(), params.AlphaNoise.GetParameters());
	}
	else {
		AlphaNoise = params.AlphaNoise;
	}
}

double IsingParameters::_psi_normalization()
{
	return 4 * Beta / (Beta * Beta - 1);
}

int IsingParameters::NumSites(bool force_calc) const {
	if (force_calc) {
		int res = 1;
		for (int i = 0; i < N_DIM; i++) {
			res *= LatticeShape[i];
		}
		return res;
	}
	else {
		return _num_sites;
	}
}

void IsingParameters::GetSiteCoordinates(int site, int* res) const {
	for (int i = 0; i < N_DIM; i++) {
		int sect = 1;
		for (int j = i+1; j < N_DIM; j++) {
			sect *= LatticeShape[j];
		}
		int coord = site / sect;
		site -= coord * sect;
		res[i] = coord;
	}
}

int IsingParameters::SiteID(int* site_coords) const {
	int res = 0;
	for (int i = 0; i < N_DIM; i++) {
		int sect = 1;
		for (int j = i + 1; j < N_DIM; j++) {
			sect *= LatticeShape[j];
		}
		res += site_coords[i] * sect;
		}
	return res;
}

int IsingParameters::ProjectionArea(int axis) const {
	if (axis >= 0 && axis < N_DIM) {
		int res = 1;
		for (int i = 0; i < N_DIM; i++) {
			if (i != axis) {
				res *= LatticeShape[i];
			}
		}
		return res;
	}
	else {
		return -1;
	}
}

double IsingParameters::GetMeanFieldGamma(double rate, double custom_alpha) const {
	double use_alpha = custom_alpha;
	if (use_alpha < 0) use_alpha = AlphaAverage;
	if (rate > 1 / Tau0) {
		return 1.0 / pow(K / (rate / Omega - 1.0 / (Omega * Tau0)) - use_alpha * pow(Omega / rate, Beta), 1.0 / n); // Domenico's eq. 18
	}
	else if (rate == 1 / Tau0) {
		return 0;
	}
	else {
		return -1;
	}
}

double IsingParameters::GetGlassinessFromAlpha(double alpha)
{
	if (_params_initialized == false) {
		UpdateDerivedParams();
	}
	if (alpha == BAD_VALUE) {
		return Glassiness_MF;
	}
	else {
		return 1 - (pow((Beta + 1) / (Beta - 1), Beta) * K / (alpha * pow(Omega * Tau0, Beta - 1))) / _psi_normalization();
	}
}

double IsingParameters::GetPsiFromAlpha(double alpha)
{
	return pow((Beta + 1) / (Beta - 1), Beta) * K / (alpha * pow(Omega * Tau0, Beta - 1)) / _psi_normalization();
}

double IsingParameters::CalcPsiFromNormGlassiness(double norm_g) const
{
	return 1 - norm_g * CalcGMax();
}

double IsingParameters::GetNormGlassinessFromAlpha(double alpha)
{
	return GetGlassinessFromAlpha(alpha) / GetGMax();
}

double IsingParameters::GetDistFromGmax(double alpha)
{
	return 1 - GetNormGlassinessFromAlpha(alpha);
}

double IsingParameters::CalcCriticalStrainFromAlpha(double alpha) const
{
	return pow(alpha * pow(Omega * Tau0, Beta) * pow((Beta - 1) / (Beta + 1), Beta + 1), -1.0 / n);
}

double IsingParameters::CalcCriticalKFromAlpha(double alpha) const
{
	return 4* alpha * Beta* pow(Omega, Beta - 1)* pow(Tau0, -2) / (pow(Beta - 1, 2) * pow(CalcCriticalGamma(), Beta + 1));
}

double IsingParameters::GetYieldStrainFromAlpha(double alpha)
{
	double C3 = K * Omega;
	double C2 = -2 * alpha * Omega * Omega;
	double C1 = 4 * alpha * Omega * Omega / Tau0;
	double C0 = -2 * alpha * Omega * Omega / (Tau0 * Tau0);
	if (_cbroots == NULL) {
		_cbroots = new double[4];
	}
	int nroots = RealRoots_P3(C3, C2, C1, C0, _cbroots);
	double jump_rate = -1;
	for (int i = 0; i < nroots; i++) {
		if (_cbroots[i] > 1.0 / Tau0) {
			if (jump_rate < 0 || jump_rate > _cbroots[i]) {
				jump_rate = _cbroots[i];
			}
		}
	}
	if (jump_rate < 0) {
		return -1;
	}
	else {
		return GetMeanFieldGamma(jump_rate, alpha);
	}
}

double IsingParameters::CalcYieldStrainFromAlpha(double alpha) const
{
	double C3 = K * Omega;
	double C2 = -2 * alpha * Omega * Omega;
	double C1 = 4 * alpha * Omega * Omega / Tau0;
	double C0 = -2 * alpha * Omega * Omega / (Tau0 * Tau0);
	double *tmp_cbroots = new double[4];
	int nroots = RealRoots_P3(C3, C2, C1, C0, tmp_cbroots);
	double jump_rate = -1;
	for (int i = 0; i < nroots; i++) {
		if (tmp_cbroots[i] > 1.0 / Tau0) {
			if (jump_rate < 0 || jump_rate > tmp_cbroots[i]) {
				jump_rate = tmp_cbroots[i];
			}
		}
	}
	if (jump_rate < 0) {
		return -1;
	}
	else {
		return GetMeanFieldGamma(jump_rate, alpha);
	}
	delete[] tmp_cbroots;
	tmp_cbroots = NULL;
}

double IsingParameters::GetCriticalGamma()
{
	Gammac = CalcCriticalGamma();
	return Gammac;
}

double IsingParameters::CalcPsiMin() const
{
	return (pow((Beta + 1) / Beta, Beta) / (Beta - 1)) / (4 * Beta / (Beta * Beta - 1));//= _old_psi_min / _old_psi_c;
}

double IsingParameters::CalcGMax() const
{
	return 1 - CalcPsiMin();
}

double IsingParameters::CalcCriticalGamma() const
{
	return (1 + 2 / (Beta - 1)) / Tau0;
}

double IsingParameters::GetPsiMin()
{
	psi_min = CalcPsiMin();
	return psi_min;
}

double IsingParameters::GetGMax()
{
	return 1 - GetPsiMin();
}

double IsingParameters::CalcAlphaFromNormGlassiness(double norm_g) const
{
	return CalcAlphaFromPsi(CalcPsiFromNormGlassiness(norm_g));
}

double IsingParameters::CalcAlphaFromPsi(double psi_val) const
{
	return K * pow((Beta + 1) / (Beta - 1), Beta) / (psi_val * (4 * Beta / (Beta * Beta - 1)) * pow(Omega * Tau0, Beta - 1));
}
