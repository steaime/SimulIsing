#pragma once

#ifndef MONTECARLO_PARAMETERS_H
#define MONTECARLO_PARAMETERS_H

#include "stdafx.h"
#include "Constants.h"
#include "SharedFunctions.h"
#include "ConfigReader.h"

#define MC_CHECK_NEGATIVES 1

class MCParams
{
public:
	MCParams();
	MCParams(double effT, double RandomNoiseCoeff = RANDOM_NOISE_COEFF, bool RandomNoiseLog = RANDOM_NOISE_LOG, 
			int MaxIter = DEF_SIMSTEP_MAX_ITER, double ConvPrecision = DEF_SIM_CONV_PREC, bool LogDeviation = DEF_DEVIATION_LOG,
			int EqProtocol = TemperatureProtocol::CONSTANT, double* EqProtParams = NULL, int EqRunNumber = DEF_EQ_RUN_NUM);
	MCParams(std::string sParamFile);
	MCParams(const ConfigParams &params);
	~MCParams();
	double RandomJump(double val_target);
	void RandomJumps(double* result, const double* val_target, int num_vals);
	double GetDeviation(const double* val_current, const double* val_target, int num_vals) const;
	bool Equilibrated(const double* val_current, const double* val_target, int num_vals) const;
	double GetTemperature() const;
	double GetNoiseCoeff() const;
	bool GetNoiseLog() const;
	void SetTemperature(double NewT);
	void SetTemperatureStep(int EqStepIndex);
	int SetTempProtocol(int newProtocol, double* ProtocolParams);
	int GetTempProtocol() const;
	int TempProtocolNumParams() const;
	//double* GetTempProtocolParams() const;
	double GetTempProtocolParam(int paramIdx) const;
	void CopyTempProtParams(double *copyTo) const;
	int GetEquilibrationRunNumber() const;
	int SetEquilibrationRunNumber(int new_val);
	double GetConvergencePrecision() const;
	int GetMaxIterations() const;
	bool GetLogDeviation() const;

	void Initialize(const ConfigParams &params);
	void CopyFrom(const MCParams &params);

private:
	double _noise_coeff;
	bool _noise_log;
	double _eff_T;
	int _T_protocol;
	int _T_protocol_numparams;
	double* _T_protocol_params;
	int _eq_run_number;
	double _conv_prec;
	int _max_iter;
	bool _log_dev;
	std::random_device _rd{};
	std::mt19937 _gen{ _rd() };
};


#endif // !MONTECARLO_PARAMETERS_H
