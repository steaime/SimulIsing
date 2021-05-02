#include "MontecarloParameters.h"

MCParams::MCParams()
{
	_eff_T = DEF_SIM_TEMP;
	_noise_coeff = RANDOM_NOISE_COEFF;
	_noise_log = RANDOM_NOISE_LOG;
	_conv_prec = DEF_SIM_CONV_PREC;
	_max_iter = DEF_SIMSTEP_MAX_ITER;
	_log_dev = DEF_DEVIATION_LOG;
	_eq_run_number = DEF_EQ_RUN_NUM;
	_T_protocol = TemperatureProtocol::CONSTANT;
	_T_protocol_params = new double[1];
	_T_protocol_params[0] = _eff_T;
	_T_protocol_numparams = 1;
}


MCParams::MCParams(double effT, double RandomNoiseCoeff, bool RandomNoiseLog, int MaxIter, double ConvPrecision, bool LogDeviation,
					int EqProtocol, double* EqProtParams, int EqRunNumber)
{
	if (effT > 0) {
		_eff_T = effT;
	}
	else {
		_eff_T = 0;
	}
	_noise_coeff = RandomNoiseCoeff;
	_noise_log = RandomNoiseLog;
	_conv_prec = ConvPrecision;
	_max_iter = MaxIter;
	_log_dev = LogDeviation;
	_eq_run_number = EqRunNumber;
	_T_protocol = TemperatureProtocol::CONSTANT;
	_T_protocol_params = new double[1];
	_T_protocol_params[0] = _eff_T;
	_T_protocol_numparams = 1;
	if (EqProtParams != NULL) {
		SetTempProtocol(EqProtocol, EqProtParams);
	}
}

MCParams::MCParams(std::string sParamFile) {
	ConfigParams params(sParamFile);
	Initialize(params);
}

MCParams::MCParams(const ConfigParams & params)
{
	Initialize(params);
}

void MCParams::Initialize(const ConfigParams & params) {
	_eff_T = params.dTemperature;
	_noise_coeff = params.dTempNoiseCoeff;
	_noise_log = params.bTempNoiseLog;
	_conv_prec = params.dConvPrecision;
	_max_iter = params.iMaxIter;
	_log_dev = params.bLogDev;
	_eq_run_number = params.iEqRunNum;
	_T_protocol = params.iTempProt;
	double *tmp_params = new double[INI_TPROT_MAXPARAMS];
	tmp_params[0] = _eff_T;
	for (int i = 1; i < INI_TPROT_MAXPARAMS; i++) {
		tmp_params[i] = params.pdTempProtParams[i];
	}
	SetTempProtocol(_T_protocol, tmp_params);
	delete[] tmp_params;
	tmp_params = NULL;
}

void MCParams::CopyFrom(const MCParams &params) {
	_noise_coeff = params.GetNoiseCoeff();
	_noise_log = params.GetNoiseLog();
	_eff_T = params.GetTemperature();
	_T_protocol = params.GetTempProtocol();
	_T_protocol_numparams = params.TempProtocolNumParams();
	double *tmp_params = new double[TempProtocolNumParams()];
	params.CopyTempProtParams(tmp_params);
	SetTempProtocol(_T_protocol, tmp_params);
	_eq_run_number = params.GetEquilibrationRunNumber();
	_conv_prec = params.GetConvergencePrecision();
	_max_iter = params.GetMaxIterations();
	_log_dev = params.GetLogDeviation();
	delete[] tmp_params;
	tmp_params = NULL;
}

MCParams::~MCParams()
{
}

double MCParams::RandomJump(double val_target) {
	if (_eff_T > 0 && _noise_coeff > 0) {
		double cur_avg;
		if (_noise_log) {
			cur_avg = val_target;
		}
		else {
			cur_avg = log10(val_target);
		}
		std::normal_distribution<> d{ cur_avg, _noise_coeff*_eff_T };
		double cur_value = d(_gen);
		if (_noise_log) {
			return pow(10, cur_value);
		}
		else {
#if MC_CHECK_NEGATIVES
			if (cur_value <= 0) {
				cur_value = RATE_EPS;
			}
#endif
			return cur_value;
		}
	}
	else {
		return val_target;
	}
}

void MCParams::RandomJumps(double* result, const double* val_target, int num_vals) {
	if (_eff_T > 0 && _noise_coeff > 0) {
		std::normal_distribution<> d{ 0, _noise_coeff*_eff_T };
		for (int i = 0; i < num_vals; i++) {
			if (_noise_log) {
				result[i] = pow(10, d(_gen) + log10(val_target[i]));
			}
			else {
				result[i] = d(_gen) + val_target[i];
			}
		}
	}
	else {
		for (int i = 0; i < num_vals; i++) {
			result[i] = val_target[i];
		}
	}
}

double MCParams::GetDeviation(const double* val_current, const double* val_target, int num_vals) const {
	double res = 0;
	for (int i = 0; i < num_vals; i++) {
		if (_log_dev) {
			res += pow((log10(val_current[i]) - log10(val_target[i])), 2)/pow(log10(val_target[i]), 2);
		}
		else {
			res += (val_current[i] - val_target[i])*(val_current[i] - val_target[i])/(val_target[i]* val_target[i]);
		}
	}
	return res / num_vals;
}

bool MCParams::Equilibrated(const double* val_current, const double* val_target, int num_vals) const {
	if (_conv_prec > 0) {
		double dev = GetDeviation(val_current, val_target, num_vals);
		bool res = (dev < _conv_prec);
		return res;
	}
	else {
		return false;
	}
}

int MCParams::GetMaxIterations() const {
	return _max_iter;
}

bool MCParams::GetLogDeviation() const
{
	return _log_dev;
}

double MCParams::GetTemperature() const
{
	return _eff_T;
}

double MCParams::GetNoiseCoeff() const
{
	return _noise_coeff;
}

bool MCParams::GetNoiseLog() const
{
	return _noise_log;
}

double MCParams::GetConvergencePrecision() const
{
	return _conv_prec;
}

void MCParams::SetTemperature(double NewT) {
	if (NewT >= 0) {
		_eff_T = NewT;
	}
}

/*
- TemperatureProtocol::CONSTANT
	constant temperature
	1 parameter required (the effective temperature)
- TemperatureProtocol::LINRAMP
	temperature linearly decreasing from T0 to 0
	2 parameters required
	- ProtocolParams[0] is T0
	  If negative, the current temperature will be used as T0
	- ProtocolParams[1] is the rate (T[i]-T[i+1], where i is the equilibration step index)
	  It MUST be positive
- TemperatureProtocol::EXPRAMP
	temperature exponentially decreasing from T0 to T1, then set to 0
	3 parameters required
	- ProtocolParams[0] is T0
	  If negative, the current temperature will be used as T0
	- ProtocolParams[1] is T1
	  It MUST be >= T0 (or >= to the current temperature, if T0<0)
	- ProtocolParams[2] is the exponential rate (log10(T[i])-log10(T[i+1]), where i is the equilibration step index)
	  It MUST be positive
- TemperatureProtocol::EXPRAMP
	temperature decreasing as a power law from T0 to T1, then set to 0
	3 parameters required
	- ProtocolParams[0] is T0
	If negative, the current temperature will be used as T0
	- ProtocolParams[1] is T1
	It MUST be >= T0 (or >= to the current temperature, if T0<0)
	- ProtocolParams[2] is the exponent (It MUST be negative)
*/
int MCParams::SetTempProtocol(int newProtocol, double* ProtocolParams) {
	if (newProtocol == TemperatureProtocol::CONSTANT) {
		SetTemperature(ProtocolParams[0]);
		_T_protocol = newProtocol;
		_T_protocol_numparams = 1;
	}
	else if (newProtocol == TemperatureProtocol::LINRAMP) {
		if (ProtocolParams[1] > 0) {
			_T_protocol_params = new double[2];
			SetTemperature(ProtocolParams[0]);
			_T_protocol_params[0] = _eff_T;
			_T_protocol_params[1] = ProtocolParams[1];
			_T_protocol = newProtocol;
			_T_protocol_numparams = 2;
		}
	}
	else if (newProtocol == TemperatureProtocol::EXPRAMP) {
		if (ProtocolParams[2] > 0) {
			double tmp_cur_T = _eff_T;
			SetTemperature(ProtocolParams[0]);
			if (_eff_T > ProtocolParams[1]) {
				_T_protocol_params = new double[3];
				_T_protocol_params[0] = _eff_T;
				_T_protocol_params[1] = ProtocolParams[1];
				_T_protocol_params[2] = ProtocolParams[2];
				_T_protocol = newProtocol;
				_T_protocol_numparams = 3;
			}
			else {
				SetTemperature(tmp_cur_T);
			}
		}
	}
	else if (newProtocol == TemperatureProtocol::PWRRAMP) {
		if (ProtocolParams[2] < 0) {
			double tmp_cur_T = _eff_T;
			SetTemperature(ProtocolParams[0]);
			if (_eff_T > ProtocolParams[1]) {
				_T_protocol_params = new double[3];
				_T_protocol_params[0] = _eff_T;
				_T_protocol_params[1] = ProtocolParams[1];
				_T_protocol_params[2] = ProtocolParams[2];
				_T_protocol = newProtocol;
				_T_protocol_numparams = 3;
			}
			else {
				SetTemperature(tmp_cur_T);
			}
		}
	}
	return _T_protocol;
}


void MCParams::SetTemperatureStep(int EqStepIndex) {
	if (_T_protocol == TemperatureProtocol::CONSTANT) {
		SetTemperature(_T_protocol_params[0]);
	}
	else if (_T_protocol == TemperatureProtocol::LINRAMP) {
		double test_temp = _T_protocol_params[0] - EqStepIndex * 1.0 / _T_protocol_params[1];
		if (test_temp > 0) {
			SetTemperature(test_temp);
		}
		else {
			SetTemperature(0);
		}
	}
	else if (_T_protocol == TemperatureProtocol::EXPRAMP) {
		double test_temp = _T_protocol_params[0] * pow(10, -EqStepIndex * 1.0 / _T_protocol_params[2]);
		if (test_temp >= _T_protocol_params[1]) {
			SetTemperature(test_temp);
		}
		else {
			SetTemperature(0);
		}
	}
	else if (_T_protocol == TemperatureProtocol::PWRRAMP) {
		double test_temp = _T_protocol_params[0] * pow(1.0*EqStepIndex, _T_protocol_params[2]);
		if (test_temp >= _T_protocol_params[1]) {
			SetTemperature(test_temp);
		}
		else {
			SetTemperature(0);
		}
	}
}

int MCParams::GetTempProtocol() const {
	return _T_protocol;
}

int MCParams::TempProtocolNumParams() const {
	return _T_protocol_numparams;
}

/*double* MCParams::GetTempProtocolParams() const {
	return _T_protocol_params;
}*/

double MCParams::GetTempProtocolParam(int paramIdx) const
{
	return _T_protocol_params[paramIdx];
}

void MCParams::CopyTempProtParams(double * copyTo) const
{
	if (_T_protocol_params != NULL) {
		for (int i = 0; i < TempProtocolNumParams(); i++) {
			copyTo[i] = _T_protocol_params[i];
		}
	}
}

int MCParams::GetEquilibrationRunNumber() const {
	return _eq_run_number;
}

int MCParams::SetEquilibrationRunNumber(int new_val) {
	if (new_val > 0) {
		_eq_run_number = new_val;
		return new_val;
	}
	else {
		return -1;
	}
}
