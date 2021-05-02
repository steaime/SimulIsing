#include "ConfigReader.h"

ConfigParams::ConfigParams()
{
}

ConfigParams::ConfigParams(std::string sParamFile)
{
	ReadINI(sParamFile);
}

ConfigParams::~ConfigParams()
{
}

void ConfigParams::ReadINI(std::string sParamFile)
{
	INIReader simParams(sParamFile);

	iNumSites = simParams.GetInteger("model", "NumSites", 0);
	for (int i = 0; i < INI_MAXNUMDIM; i++) {
		piLatticeShape[i] = simParams.GetInteger("model", "LatticeShape" + std::to_string(i), 0);
	}


	K = simParams.GetReal("model", "K", 0);
	n = simParams.GetReal("model", "n", 0);
	Tau0 = simParams.GetReal("model", "Tau0", 0);
	Omega = simParams.GetReal("model", "Omega", 0);
	Beta = simParams.GetReal("model", "Beta", 2.0);
	dUsePBC = simParams.GetBoolean("model", "UsePBC", true);
	bUseNewEqn = simParams.GetBoolean("model", "UseNewEqn", true);
	bAlphaFromGlassiness = simParams.GetBoolean("model", "NoiseOnGlassiness", false);
	bSiteBasedNoise = simParams.GetBoolean("model", "SiteBasedNoise", true);
	iSBNoiseAvgType = simParams.GetInteger("model", "SBNoiseAvgType", AvgType::GEOMETRIC);
	bForbidUnphysicalAlpha = simParams.GetBoolean("model", "ForbidUnphysicalAlpha", false);
	iRateAvgType = simParams.GetInteger("model", "RateAvgType", AvgType::ARITMETIC);
	iAlphaNoiseType = simParams.GetInteger("model", "AlphaNoiseType", 0);
	if (bAlphaFromGlassiness) {
		dNormDistFromGmax = 1 - simParams.GetReal("model", "NormGlassiness", 0);
		dAlphaAvg = -1.0;
	}
	else {
		dAlphaAvg = simParams.GetReal("model", "AlphaAverage", 0);
		dNormDistFromGmax = -1.0;
	}

	// This for loop shouldn't be necessary, but to avoid problems...
	for (int i = 0; i < NOISE_MAXPARAMS; i++) {
		pdAlphaParams[i] = 0;
		pdInitParams[i] = 0;
	}
	sAlphaParamFile = simParams.Get("model", "AlphaParamFile", "");
	if (sAlphaParamFile == "") {
		for (int i = 0; i < INI_MAXPARAMS; i++) {
			pdAlphaParams[i] = simParams.GetReal("model", "AlphaParam" + std::to_string(i), 0);
		}
	}
	else {
		int64_t num_alphapars = -1;
		ImportDataFromRawBinary(sAlphaParamFile, &(pdAlphaParams[1]), num_alphapars, NOISE_MAXPARAMS);
		pdAlphaParams[0] = (double)num_alphapars;
	}
	sInitParamFile = simParams.Get("initialization", "InitParamFile", "");
	if (sInitParamFile == "") {
		for (int i = 0; i < INI_MAXPARAMS; i++) {
			pdInitParams[i] = simParams.GetReal("initialization", "InitParam" + std::to_string(i), 0);
		}
	}
	else {
		int64_t num_initpars = -1;
		ImportDataFromRawBinary(sInitParamFile, &(pdInitParams[1]), num_initpars, NOISE_MAXPARAMS);
		pdInitParams[0] = (double)num_initpars;
	}

	bForceMeanField = simParams.GetBoolean("simulation", "ForceMeanField", false);
	iNumSimulations = simParams.GetInteger("simulation", "SimNumber", 0);
	bGammaFile = simParams.GetBoolean("simulation", "GammaListFromFile", false);
	sGammaListFile = simParams.Get("simulation", "GammaListFile", "");
	dMinGamma = simParams.GetReal("simulation", "MinGamma", 0);
	dMaxGamma = simParams.GetReal("simulation", "MaxGamma", 0);
	iPPDGamma = simParams.GetInteger("simulation", "PPDGamma", 0);
	bPreshear = simParams.GetBoolean("simulation", "RunPreshear", false);
	bAscDesc = simParams.GetBoolean("simulation", "RunAscDesc", false);

	out_folder = simParams.Get("output", "OutFolder", "");
	if (out_folder.back() != '\\') {
		out_folder.append("\\");
	}
	bMuteOut = !(simParams.GetBoolean("output", "PrintOutput", false));
	bSaveHist = simParams.GetBoolean("output", "SaveRateHistogram", false);
	bSaveHistAlpha = simParams.GetBoolean("output", "SaveAlphaHistogram", false);
	bSaveCorr = simParams.GetBoolean("output", "SaveCorrelations", false);
	bSaveRaw = simParams.GetBoolean("output", "SaveRawSnapshots", false);
	bPrintStd = simParams.GetBoolean("output", "Out_PrintStdRate", true);
	bPrintNloops = simParams.GetBoolean("output", "Out_PrintNloops", true);
	bPrintFastSites = simParams.GetBoolean("output", "Out_PrintFastSites", true);
	bSaveAlphaSnapshot = simParams.GetBoolean("output", "SaveAlphaSnapshot", true);
	bSaveInitSnapshot = simParams.GetBoolean("output", "SaveInitSnapshot", true);

	iInitType = simParams.GetInteger("initialization", "Type", 0);
	dInitAvg = simParams.GetReal("initialization", "Average", 0);
	if (dInitAvg < 0) {
		dInitAvg = 1.0 / Tau0;
	}

	dTemperature = simParams.GetReal("equilibration", "Temperature", DEF_SIM_TEMP);
	dTempNoiseCoeff = simParams.GetReal("equilibration", "NoiseCoefficient", RANDOM_NOISE_COEFF);
	bTempNoiseLog = simParams.GetBoolean("equilibration", "NoiseLog", RANDOM_NOISE_LOG);
	dConvPrecision = simParams.GetReal("equilibration", "ConvergencePrecision", DEF_SIM_CONV_PREC);
	iMaxIter = simParams.GetInteger("equilibration", "MaxIterations", DEF_SIMSTEP_MAX_ITER);
	bLogDev = simParams.GetBoolean("equilibration", "DeviationLog", DEF_DEVIATION_LOG);
	iEqRunNum = simParams.GetInteger("equilibration", "EquilibrationRunNumber", DEF_EQ_RUN_NUM);
	iTempProt = simParams.GetInteger("equilibration", "Protocol", TemperatureProtocol::CONSTANT);
	for (int i = 0; i < INI_TPROT_MAXPARAMS; i++) {
		pdTempProtParams[i] = simParams.GetReal("equilibration", "ProtocolParam" + std::to_string(i), 0);
	}

	txt_ext = simParams.Get("format", "TextExt", ".txt");
	raw_ext = simParams.Get("format", "RawExt", ".raw");
	hist_prefix = simParams.Get("format", "HistogramPrefix", "hist_");
	histalpha_fname = simParams.Get("format", "AlphaHistogramName", "hist_alpha");
	corr_prefix = simParams.Get("format", "CorrelationPrefix", "corr_");
	rawfinal_prefix = simParams.Get("format", "RawFinalPrefix", "raw_");
	preshear_prefix = simParams.Get("format", "PreshearPrefix", "ps_");
	ascdesc_prefix = simParams.Get("format", "AscdescPrefix", "ad_");
	fout_prefix = simParams.Get("format", "OutputFilePrefix", "out_");
	fstat_prefix = simParams.Get("format", "StatsFilePrefix", "stat_");
	fparams = simParams.Get("format", "SimParamsName", "sim_params.txt");
	sAlphaSnapshotPrefix = simParams.Get("format", "AlphaSnapshotPrefix", "alpharaw_");
	sInitSnapshotPrefix = simParams.Get("format", "InitSnapshotPrefix", "initraw_");

	iRefFileColNum = simParams.GetInteger("analysis", "RefFileColNum", 3);
	iRefFileColIdxStrain = simParams.GetInteger("analysis", "RefFileColIdxStrain", 0);
	iRefFileColIdxRate = simParams.GetInteger("analysis", "RefFileColIdxRate", 1);
	iRefFileColIdxChi = simParams.GetInteger("analysis", "RefFileColIdxChi", 2);
	dFastSitesRelThr = simParams.GetReal("analysis", "FastSitesRelThr", MAX_SOLID_RATE_RATIO);
	bDiffRelRate = simParams.GetBoolean("analysis", "DiffRelRate", true);
	dDiffWeightRate = simParams.GetReal("analysis", "DiffWeightRate", 0);
	dDiffWeightChi = simParams.GetReal("analysis", "DiffWeightChi", 0);
	bDiffCombine = simParams.GetBoolean("analysis", "DiffCombine", true);
}

int ConfigParams::GetNumGamma() const
{
	if (bGammaFile) {
		return FileCountLines(sGammaListFile);
	}
	else {
		return LogSpace_CalcNum(dMinGamma, dMaxGamma, iPPDGamma);
	}
}

void ConfigParams::GetListGamma(double* pdOut) const
{
	if (bGammaFile) {
		FileLoadValues(sGammaListFile, pdOut);
	}
	else {
		LogSpaceNum(dMinGamma, dMaxGamma, GetNumGamma(), pdOut);
	}
}
