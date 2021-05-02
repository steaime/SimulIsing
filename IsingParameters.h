#pragma once
#ifndef ISING_PARAMETERS_H
#define ISING_PARAMETERS_H

#include "stdafx.h"
#include "Constants.h"
#include "Noise.h"
#include "ConfigReader.h"

class IsingParameters
{
public:
	IsingParameters();
	IsingParameters(int* _LatticeShape, double _n, double _K, double _Tau0, double _Omega, double _Beta, 
					const Noise &_AlphaNoise, bool _UsePBC, bool _UseNewEqn = true, int _avgType = AvgType::ARITMETIC);
	IsingParameters(std::string sParamFile);
	IsingParameters(const ConfigParams &params);
	~IsingParameters();

	int LatticeShape[N_DIM];
	bool UsePBC = true;
	bool UseNewEqn = true;
	double K = NULL;
	double n = NULL;
	double Tau0 = NULL;
	double Omega = NULL;
	double Beta = NULL;
	double AlphaAverage = NULL;
	int RateAvgType = AvgType::ARITMETIC;
	int SBNoiseAvgType = AvgType::GEOMETRIC;
	Noise AlphaNoise;

	bool NoiseOnGlassiness = false;
	bool ForbidUnphysicalAlpha = true;
	bool SiteBasedNoise = true;

	double Gamma0 = NULL;
	double Gammac = NULL;
	double Kc = NULL;
	double gammac = NULL;
	double psi = NULL;
	double psi_min = NULL;
	double Glassiness_MF = NULL;
	double NormDistFromGmax_MF = NULL;

	int NumSites(bool force_calc = false) const;
	int ProjectionArea(int axis) const;
	void GetSiteCoordinates(int site, int* res) const;
	int SiteID(int* site_coords) const;

	double GetPsiMin();
	double GetGMax();
	double GetCriticalGamma();

	double GetPsiFromNormGlassiness(double norm_g);
	double GetAlphaFromNormGlassiness(double norm_g);
	double GetAlphaFromPsi(double psi_val);
	double GetPsiFromAlpha(double alpha);
	double GetGlassinessFromAlpha(double alpha = BAD_VALUE);
	double GetNormGlassinessFromAlpha(double alpha = BAD_VALUE);
	double GetDistFromGmax(double alpha = BAD_VALUE);
	double GetCriticalStrainFromAlpha(double alpha);
	double GetCriticalKFromAlpha(double alpha);
	double GetYieldStrainFromAlpha(double alpha);

	double GetMeanFieldGamma(double rate, double custom_alpha = -1) const;
	
	void Initialize(const ConfigParams &params);
	void UpdateDerivedParams();
	void CopyFrom(const IsingParameters &params, bool deep_copy = false);

private:
	int _num_sites = 0;

	bool _params_initialized = false;
	double* _cbroots = NULL;	 // container for cubic roots used to compute yield strain
	double _psi_normalization(); // it returns what was once psi_c in the old model (=8/3 for Beta=2)
};

#endif // !ISING_PARAMETERS_H
