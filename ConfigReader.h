#pragma once
#include "stdafx.h"
#include "Constants.h"
#include "SharedFunctions.h"

#ifndef CONFIGREADER_H
#define CONFIGREADER_H

class ConfigParams
{
public:

	ConfigParams();
	ConfigParams(std::string sParamFile);
	~ConfigParams();

	double dMinGamma = 0.0;
	double dMaxGamma = 0.0;
	double dInitAvg = 0;
	double K = 0;
	double n = 0;
	double Tau0 = 0;
	double Omega = 0;
	double Beta = 0;
	double dAlphaAvg = 0;
	double dNormDistFromGmax = 0.0;
	double dTemperature = DEF_SIM_TEMP;
	double dTempNoiseCoeff = RANDOM_NOISE_COEFF;
	double dConvPrecision = DEF_SIM_CONV_PREC;
	double dFastSitesRelThr = MAX_SOLID_RATE_RATIO;
	double dDiffWeightRate = 0;
	double dDiffWeightChi = 0;
	//double pdAlphaParams[NOISE_MAXPARAMS];			// TODO: change this to dynamic!
	//double pdInitParams[NOISE_MAXPARAMS];			// TODO: change this to dynamic!
	double *pdAlphaParams = NULL;
	double *pdInitParams = NULL;
	double pdTempProtParams[INI_TPROT_MAXPARAMS];

	int iNumSimulations = 0;
	int iPPDGamma = 0;
	int iInitType = 0;
	int iNumSites = 0;
	int iAlphaNoiseType = NoiseTypes::NONOISE;
	int iRateAvgType = AvgType::ARITMETIC;
	int iSBNoiseAvgType = AvgType::GEOMETRIC;
	int iInitFileHdrLen = 0;
	int iAlphaFileHdrLen = 0;
	int iMaxIter = DEF_SIMSTEP_MAX_ITER;
	int iEqRunNum = DEF_EQ_RUN_NUM;
	int iTempProt = TemperatureProtocol::CONSTANT;
	int iRefFileColNum = 3;
	int iRefFileColIdxStrain = 0;
	int iRefFileColIdxRate = 1;
	int iRefFileColIdxChi = 2;
	int piLatticeShape[INI_MAXNUMDIM];

	std::time_t iRandSeed = -1;

	bool bPreshear = false;
	bool bAscDesc = false;
	bool bSaveHist = false;
	bool bSaveHistAlpha = false;
	bool bSaveCorr = false;
	bool bSaveRaw = false;
	bool bPrintStd = true;
	bool bSaveAlphaSnapshot = false;
	bool bSaveInitSnapshot = false;
	bool bAlphaFileSwapEndian = false;
	bool bInitFileSwapEndian = false;
	bool bPrintNloops = true;
	bool bPrintFastSites = true;
	bool bMuteOut = false;
	bool bGammaFile = false;
	bool dUsePBC = true;
	bool bUseNewEqn = true;
	bool bAlphaFromGlassiness = false;
	bool bSiteBasedNoise = true;
	bool bForbidUnphysicalAlpha = false;
	bool bForceMeanField = false;
	bool bTempNoiseLog = RANDOM_NOISE_LOG;
	bool bLogDev = DEF_DEVIATION_LOG;
	bool bDiffRelRate = true;
	bool bDiffCombine = true;

	std::string sAlphaParamFile = "";
	std::string sInitParamFile = "";
	std::string sAlphaSnapshotPrefix = "";
	std::string sInitSnapshotPrefix = "";
	std::string sGammaListFile = "";
	std::string out_folder = "";
	std::string txt_ext = ".txt";
	std::string raw_ext = ".raw";
	std::string hist_prefix = "hist_";
	std::string histalpha_fname = "hist_alpha";
	std::string corr_prefix = "corr_";
	std::string rawfinal_prefix = "raw_";
	std::string preshear_prefix = "ps_";
	std::string ascdesc_prefix = "ad_";
	std::string fout_prefix = "out_";
	std::string fstat_prefix = "stat_";
	std::string fparams = "sim_params.txt";

	char cInitFileFmt = 'd';
	char cAlphaFileFmt = 'd';

	void ReadINI(std::string sParamFile);

	int GetNumGamma() const;
	void GetListGamma(double* pdOut) const;
	int GetNumSites() const;

	void Clear();

private:

	bool _init;
};

#endif