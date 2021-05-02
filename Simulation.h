#pragma once
#ifndef SIMULATION_H
#define SIMULATION_H

#include "stdafx.h"
#include "IsingParameters.h"
#include "InitialConditions.h"
#include "Noise.h"
#include "MontecarloParameters.h"
#include "Constants.h"
#include "SharedFunctions.h"


class Simulation
{
public:
	Simulation();
	Simulation(int n_sites, bool use_pbc, double omega, double K, double tau0,
		double alpha, double beta, int alpha_noise_type, double* alpha_noise_params,
		int init_noise_type, double* init_noise_params, double effective_T = DEF_SIM_TEMP);
	Simulation(const IsingParameters &model_parameters, const Noise &init_conditions, const MCParams &tstep_params);
	~Simulation();
	void Setup(const IsingParameters &model_parameters, const Noise &init_conditions, const MCParams &tstep_params, bool force_meanfield = false);
	void Clear();

	//IsingParameters GetParams() const;
	bool GetParam_UsePBC() const;
	double GetParam_K() const;
	double GetParam_n() const;
	double GetParam_Tau0() const;
	double GetParam_Omega() const;
	double GetParam_Beta() const;
	double GetParam_AlphaAvg() const;
	int GetParam_NumSites(bool force_calc = false) const;

	void CopyParams(IsingParameters *copyTo) const;
	void CopyLatticeShape(int *copyTo) const;
	//Noise GetInitNoise() const;
	void CopyInitNoise(Noise *copyTo) const;
	int UpdateParams(const IsingParameters &new_params);
	bool ValidateParams(const IsingParameters &params) const;

	void Initialize();
	void InitializeAlpha();
	void GenInitialConditions();
	int ValidateRates();
	int ValidateAlphas(double min_val = ALPHA_EPS, double max_val = BAD_VALUE);

	int IsSiteOnEdge(int site) const;
	void FindEdgesSite(int site, int* axes) const;
	int GetNNNumber(int site, bool force_calc = false) const;
	void GetNNAddress(int site, int** address, int* number = NULL, bool force_calc = false);

	void UpdateNearestNeighbors();
	void UpdateNNNumber();
	void UpdateNNAddress();

	double** GetAlphas() const;
	int CountAlphas() const;
	int FlattenAlphas(double* res) const;
	double* GetRates() const;
	double* GetTargetRates(bool update = true);
	void UpdateTargetRates();
	void UpdateTargetRate(int site);
	double GetAlpha(int site_i, int site_j) const;
	double GetSiteAveragedAlpha(int site_index) const;
	double GetAlphaNN(int site_i, int neighbor_index) const;
	double GetRate(int site) const;
	double GetTargetRate(int site, bool force_calc = false);

	double GetMeanFieldRate(double gamma = -1, double start_rate = -1, double precision = DEF_SIM_CONV_PREC);

	double GetAverageRate(bool force_calc = false);
	double GetRateStd() const;
	double GetTauStd(bool statistic = true);
	double GetAverageAlpha() const;
	double GetAlphaVariance() const;

	double GetTemperature() const;
	int GetMaxStepNum() const;
	int GetMaxStepPerRun() const;
	int GetEquilibrationProtocol() const;
	int GetEquilibrationParameterNumber() const;
	//double* GetEquilibrationParameters() const;
	double GetEquilibrationParameter(int paramIdx) const;
	int SetEquilibrationProtocol(int newProtocol, double* protocolParams);
	int SetEquilibrationRunNumber(int new_val);

	double GetGamma() const;
	int SetGamma(double new_gamma, bool adjust_lattice = true);

	double GetGMax();
	double GetAlphaFromNormGlassiness(double norm_g);
	double GetMeanFieldGamma(double rate) const;
	double GetGlassinessFromAlpha(double alpha = BAD_VALUE);
	double GetNormGlassinessFromAlpha(double alpha = BAD_VALUE);
	double GetDistFromGmax(double alpha = BAD_VALUE);
	double GetYieldStrainFromAlpha(double alpha);
	int CountUnphysicalSites();

	int SimStep(bool ForceUpdateAvRate = false, double* rate_change = NULL);
	int AdjustLattice(double** eqDetails = NULL, bool** blnFastSites = NULL);

	bool GetForceMeanField() const;
	void SetForceMeanField(bool val);
	bool IsFastSite(int site, double thr_rate_ratio = MAX_SOLID_RATE_RATIO) const;
	bool IsSlowSite(int site, double thr_rate_ratio = MAX_SOLID_RATE_RATIO) const;

	int GetNumSites(bool force_calc = false) const;
	int GetNumSlowSites(double thr_rate_ratio = MAX_SOLID_RATE_RATIO) const;
	double GetFractionSlowSites(double thr_rate_ratio = MAX_SOLID_RATE_RATIO) const;
	int GetNumFastSites(double thr_rate_ratio = MAX_SOLID_RATE_RATIO) const;
	double GetFractionFastSites(double thr_rate_ratio = MAX_SOLID_RATE_RATIO) const;
	void GetSiteCoordinates(int site, int* res) const;

	int GetRateCorrelation(double* res, int* x = NULL, int num_x = -1);

	int GetStatus();
	int ID();
	int StepCount();

private:
	IsingParameters _params;
	Noise _init_noise;
	MCParams _mc_params;
	double _gamma;
	double** _alpha = NULL;
	double** _normglassiness = NULL; // only instantiated and populated when noise is defined on glassiness
	double* _rates = NULL;
	double* _target_rates = NULL;
	double* _rate_step = NULL;
	double _avg_rate;
	int* _nn_number = NULL;
	int** _nn_address = NULL;
	bool _force_meanfield;
	int _ID = -1;
	int _status = SimStatus::NOT_INITIALIZED;
	int _count_iter;

	bool _UseMeanfield(int site = -1);
};

#endif // !LATTICE_H

