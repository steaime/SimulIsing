#include "ConfigReader.h"

ConfigParams::ConfigParams()
{
	_init = false;
}

ConfigParams::ConfigParams(std::string sParamFile)
{
	_init = false;
	ReadINI(sParamFile);
}

ConfigParams::~ConfigParams()
{
	if (_init) {
		delete[] pdAlphaParams;
		pdAlphaParams = NULL;
		delete[] pdInitParams;
		pdInitParams = NULL;
	}
}

void ConfigParams::Clear() {

	dMinGamma = 0.0;
	dMaxGamma = 0.0;
	dInitAvg = 0;
	K = 0;
	n = 0;
	Tau0 = 0;
	Omega = 0;
	Beta = 0;
	dAlphaAvg = 0;
	dNormDistFromGmax = 0.0;
	dTemperature = DEF_SIM_TEMP;
	dTempNoiseCoeff = RANDOM_NOISE_COEFF;
	dConvPrecision = DEF_SIM_CONV_PREC;
	dFastSitesRelThr = MAX_SOLID_RATE_RATIO;
	dDiffWeightRate = 0;
	dDiffWeightChi = 0;
	for (int i = 0; i < INI_TPROT_MAXPARAMS; i++) {
		pdTempProtParams[i] = 0;
	}

	iNumSimulations = 0;
	iPPDGamma = 0;
	iInitType = 0;
	iNumSites = 0;
	iAlphaNoiseType = NoiseTypes::NONOISE;
	iRateAvgType = AvgType::ARITMETIC;
	iSBNoiseAvgType = AvgType::GEOMETRIC;
	iMaxIter = DEF_SIMSTEP_MAX_ITER;
	iEqRunNum = DEF_EQ_RUN_NUM;
	iTempProt = TemperatureProtocol::CONSTANT;
	iRefFileColNum = 3;
	iRefFileColIdxStrain = 0;
	iRefFileColIdxRate = 1;
	iRefFileColIdxChi = 2;
	for (int i = 0; i < INI_MAXNUMDIM; i++) {
		piLatticeShape[i] = 0;
	}

	bPreshear = false;
	bAscDesc = false;
	bSaveHist = false;
	bSaveHistAlpha = false;
	bSaveCorr = false;
	bSaveRaw = false;
	bPrintStd = true;
	bSaveAlphaSnapshot = false;
	bSaveInitSnapshot = false;
	bPrintNloops = true;
	bPrintFastSites = true;
	bMuteOut = false;
	bGammaFile = false;
	dUsePBC = true;
	bUseNewEqn = true;
	bAlphaFromGlassiness = false;
	bSiteBasedNoise = true;
	bForbidUnphysicalAlpha = false;
	bForceMeanField = false;
	bTempNoiseLog = RANDOM_NOISE_LOG;
	bLogDev = DEF_DEVIATION_LOG;
	bDiffRelRate = true;
	bDiffCombine = true;

	sAlphaParamFile = "";
	sInitParamFile = "";
	sAlphaSnapshotPrefix = "";
	sInitSnapshotPrefix = "";
	sGammaListFile = "";
	out_folder = "";
	txt_ext = ".txt";
	raw_ext = ".raw";
	hist_prefix = "hist_";
	histalpha_fname = "hist_alpha";
	corr_prefix = "corr_";
	rawfinal_prefix = "raw_";
	preshear_prefix = "ps_";
	ascdesc_prefix = "ad_";
	fout_prefix = "out_";
	fstat_prefix = "stat_";
	fparams = "sim_params.txt";

	if (_init) {
		delete[] pdAlphaParams;
		pdAlphaParams = NULL;
		delete[] pdInitParams;
		pdInitParams = NULL;
	}

	_init = false;
}

void ConfigParams::ReadINI(std::string sParamFile)
{

	if (_init) Clear();

	INIReader simParams(sParamFile);

	out_folder = simParams.Get("output", "OutFolder", "");
	if (out_folder.back() != '\\') {
		out_folder.append("\\");
	}

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
	sAlphaParamFile = simParams.Get("model", "AlphaParamFile", "");
	if (sAlphaParamFile != "") {
		if (CheckPathRelative(sAlphaParamFile)) sAlphaParamFile = JoinPath(out_folder, sAlphaParamFile);
		if (!CheckFileExists(sAlphaParamFile)) {
			std::cout << "WARNING: specified model/AlphaParamFile '" << sAlphaParamFile << "' not found and bypassed" << std::endl;
			sAlphaParamFile = "";
		}
	}
	if (sAlphaParamFile == "") {
		pdAlphaParams = new double[INI_MAXPARAMS];
		for (int i = 0; i < INI_MAXPARAMS; i++) {
			pdAlphaParams[i] = simParams.GetReal("model", "AlphaParam" + std::to_string(i), 0);
		}
	}
	else {
		int32_t nimgs_alpharaw;
		int64_t npx_alpharaw;
		std::fstream alpha_fin = OpenRawBinary(sAlphaParamFile, nimgs_alpharaw, npx_alpharaw);
		int64_t num_alphapars = nimgs_alpharaw * npx_alpharaw;
		if (num_alphapars + 1 > NOISE_MAXPARAMS) num_alphapars = NOISE_MAXPARAMS - 1;
		pdAlphaParams = new double[num_alphapars + 1];
		ImportDataFromRawBinary(&alpha_fin, &(pdAlphaParams[1]), num_alphapars, NOISE_MAXPARAMS - 1);
		pdAlphaParams[0] = (double)num_alphapars;
	}

	bForceMeanField = simParams.GetBoolean("simulation", "ForceMeanField", false);
	iNumSimulations = simParams.GetInteger("simulation", "SimNumber", 0);
	bGammaFile = simParams.GetBoolean("simulation", "GammaListFromFile", false);
	sGammaListFile = simParams.Get("simulation", "GammaListFile", "");
	if (sGammaListFile != "") {
		if (CheckPathRelative(sGammaListFile)) sGammaListFile = JoinPath(out_folder, sGammaListFile);
		if (!CheckFileExists(sGammaListFile)) {
			std::cout << "WARNING: specified simulation/GammaListFile '" << sGammaListFile << "' not found and bypassed" << std::endl;
			sGammaListFile = "";
		}
	}
	dMinGamma = simParams.GetReal("simulation", "MinGamma", 0);
	dMaxGamma = simParams.GetReal("simulation", "MaxGamma", 0);
	iPPDGamma = simParams.GetInteger("simulation", "PPDGamma", 0);
	bPreshear = simParams.GetBoolean("simulation", "RunPreshear", false);
	bAscDesc = simParams.GetBoolean("simulation", "RunAscDesc", false);

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
	sInitParamFile = simParams.Get("initialization", "InitParamFile", "");
	if (sInitParamFile != "") {
		if (CheckPathRelative(sInitParamFile)) sInitParamFile = JoinPath(out_folder, sInitParamFile);
		if (!CheckFileExists(sInitParamFile)) {
			std::cout << "WARNING: specified initialization/sInitParamFile '" << sInitParamFile << "' not found and bypassed" << std::endl;
			sInitParamFile = "";
		}
	}
	iInitFileHdrLen = simParams.GetInteger("initialization", "InitFileHeader", 0);
	std::string sInitFileFmt = simParams.Get("initialization", "InitFileFormat", "d");
	cInitFileFmt = sInitFileFmt.c_str()[0];
	bInitFileSwapEndian = simParams.GetBoolean("initialization", "InitFileSwapEndian", false);
	if (sInitParamFile == "") {
		pdInitParams = new double[INI_MAXPARAMS];
		for (int i = 0; i < INI_MAXPARAMS; i++) {
			pdInitParams[i] = simParams.GetReal("initialization", "InitParam" + std::to_string(i), 0);
		}
	}
	else {
		int32_t nimgs_initraw = 1;
		int64_t npx_initraw = -1;
		std::fstream init_fin = OpenRawBinary(sInitParamFile, nimgs_initraw, npx_initraw, iInitFileHdrLen);
		int64_t num_initpars = nimgs_initraw * npx_initraw;
		if (num_initpars < 0) num_initpars = GetNumSites();
		if (num_initpars + 1 > NOISE_MAXPARAMS) num_initpars = NOISE_MAXPARAMS - 1;
		pdInitParams = new double[num_initpars + 1];
		ImportDataFromRawBinary(&init_fin, &(pdInitParams[1]), num_initpars, NOISE_MAXPARAMS - 1, true, cInitFileFmt, bInitFileSwapEndian);
		pdInitParams[0] = (double)num_initpars;
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

	_init = true;

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

int ConfigParams::GetNumSites() const
{
	int num_sites = 1;
	for (int i = 0; i < N_DIM; i++) {
		if (piLatticeShape[i] > 0) num_sites *= piLatticeShape[i];
		else num_sites *= iNumSites;
	}
	return num_sites;
}
