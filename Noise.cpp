#include "Noise.h"

Noise::Noise()
{
	Initialize(NoiseTypes::NONOISE);
}

Noise::~Noise()
{
	// REMEMBER: keep this commented! For some reason, deleting _params here
	// messes up memory a lot. Never really understood why, but just keep this commented.
	/*delete[] _params;
	_params = NULL;*/
}

Noise::Noise(int type, double avg) {
	Initialize(type, avg);
}

Noise::Noise(int type, const double* params, double avg) {
	Initialize(type, avg, params);
}

void Noise::Initialize(int type, double avg, const double* params, bool deep_copy) {
	_avg = avg;
	_type = type;
	_num_sampled = 0;
	NumParams(true, params);
#if FORCE_DEEP_COPY==0:
	if (deep_copy) {
#endif
		delete[] _params;
		if (_num_params > 0) {
			_params = new double[NumParams()];
			for (int i = 0; i < NumParams(); i++) {
				_params[i] = params[i];
			}
		}
#if FORCE_DEEP_COPY==0:
	}
	else {
		_params = params;
	}
#endif
}

void Noise::CopyFrom(const Noise &old_noise, bool deep_copy) {
	_avg = old_noise.GetAverage();
	_type = old_noise.GetType();
	_num_sampled = 0;
	NumParams(true);
	if (deep_copy) {
		delete[] _params;
		_params = new double[NumParams()];
		for (int i = 0; i < NumParams(); i++) {
			_params[i] = old_noise.GetParameters()[i];
		}
	}
	else {
		_params = old_noise.GetParameters();
	}
}

double Noise::GetAverage() const {
	return _avg;
}

int Noise::GetType() const {
	return _type;
}

double* Noise::GetParameters() const {
	return _params;
}

int Noise::GetParameters(double* res) const {
	res = _params;
	return _num_params;
}

double Noise::GetParameter(int index) const {
	if (index >= 0 && index < _num_params) {
		return _params[index];
	}
	else {
		return -1;
	}
}

int Noise::NumParams(bool force_calc, const double* params) {
	if (force_calc) {
		if (_type == NoiseTypes::NONOISE) {
			_num_params = 0;
		}
		else if (_type == NoiseTypes::WHITE_FLAT) {
			_num_params = 1;
		}
		else if (_type == NoiseTypes::WHITE_GAUSS) {
			_num_params = 1;
		}
		else if (_type == NoiseTypes::WHITE_LOGNORM) {
			_num_params = 1;
		}
		else if (_type == NoiseTypes::WHITE_CUSTOM) {
			if (params == NULL) {
				_num_params = -1;
			}
			else {
				_num_params = 3 * params[1] + 2;
			}
		}
		else if (_type == NoiseTypes::SORTED_LIST) {
			if (params == NULL) {
				_num_params = -1;
			}
			else {
				_num_params = params[0] + 1;
			}
		}
	}
	return _num_params;
}

void Noise::Reset()
{
	_num_sampled = 0;
}

void Noise::Sample(int num_samples, double* result) {
	//result = new double[num_samples];
	if (_type == NoiseTypes::NONOISE) {
		for (int i = 0; i < num_samples; i++) {
			result[i] = _avg;
		}
	} 
	else if (_type == NoiseTypes::WHITE_FLAT) {
		double variance = _params[0];
		if (_avg != 0) {
			variance = abs(variance * _avg);
		}
		double amplitude = sqrt(3 * variance);
		for (int i = 0; i < num_samples; i++) {
			result[i] = _avg + 2 * amplitude * ((double)rand() / (RAND_MAX)) - amplitude;
		}
	}
	else if (_type == NoiseTypes::WHITE_GAUSS) {
		double variance = _params[0] * _avg * _avg;
		std::random_device _rd{};
		std::mt19937 _gen{ _rd() };
		std::normal_distribution<> d{ _avg, sqrt(variance) };
		for (int i = 0; i < num_samples; i++) {
			result[i] = d(_gen);
		}
	}
	else if (_type == NoiseTypes::WHITE_LOGNORM) {
		/*
		Lognormal distribution: 
		     X = exp(mu+sqrt(sigma2)*Z)
		where:
		- Z is a standard normal variable
		- mu is the average of ln(X)
		- sigma2 is the variance of ln(X)
		Aritmetic moments of the distribution:
		- <X> = exp(mu+0.5*sigma2)
		- Var(X) = exp(2*mu+sigma2)*(exp(sigma2)-1)
		- Var(X)/<X>^2 = exp(sigma2) - 1
		Parameters (mu, sigma2) from distribution moments (<X>=mean, Var(X)/<X>^2=normvar):
		- sigma2 = ln(normvar + 1)
		- mu = ln(mean/sqrt(normvar+1)) = ln(mean) - 0.5 * sigma2
		*/
		double _sigma2 = log(_params[0] + 1.0);
		double _mu = log(_avg) - 0.5 * _sigma2;
		std::random_device _rd{};
		std::mt19937 _gen{ _rd() };
		std::normal_distribution<> Z{ _mu, sqrt(_sigma2) };
		for (int i = 0; i < num_samples; i++) {
			double _exponent = Z(_gen);
			result[i] = exp(_exponent);
		}
	}
	else if (_type == NoiseTypes::WHITE_CUSTOM) {
		/*
		Custom cumulative distribution obtained by linear interpolation of N points
		_params: {Median, N, x0, C0, Cer0, x1, C1, Cer1, ..., xN, CN, CerN}
		To avoid risky extrapolations, make sure that C0=0 and CN=1
		*/
		for (int i = 0; i < num_samples; i++) {
			double _rnd = ((double)rand() / (RAND_MAX));
			int cur_idx = -1;
			for (int j = 1; j < _params[1]; j++) {
				if (_params[3 * j + 3] > _rnd) {
					cur_idx = j;
					break;
				}
			}
			if (cur_idx < 0) cur_idx = _params[1] - 1;
			double _xA = _params[3 * cur_idx - 1], _yA = _params[3 * cur_idx], _yErA = _params[3 * cur_idx + 1];
			double _xB = _params[3 * cur_idx + 2], _yB = _params[3 * cur_idx + 3], _yErB = _params[3 * cur_idx + 4];
			double _xavg;
			if (_yB != _yA) {
				_xavg = _xA + (_xB - _xA) * (_rnd - _yA) / (_yB - _yA);
			}
			else {
				_xavg = 0.5 * (_xA + _xB);
			}
			if (_yErA <= 0 && _yErB <= 0) {
				result[i] = _xavg;
			}
			else {
				double _yerr = _yErA + (_yErB - _yErA) * (_xavg - _xA) / (_xB - _xA);
				double _xerr;
				if (_yB != _yA) {
					_xerr = _yerr * (_xB - _xA) / (_yB - _yA);
				}
				else {
					_xerr = 0.5 * (_xB - _xA);
				}
				std::random_device _rd{};
				std::mt19937 _gen{ _rd() };
				std::normal_distribution<> d{ _xavg, _xerr };
				double val = d(_gen);
				if (val > _params[3 * int(_params[1]) - 1]) {
					int cacca = -1;
				}
				result[i] = val;
			}
		}
	}
	else if (_type == NoiseTypes::SORTED_LIST) {
		for (int i = 0; i < num_samples; i++) {
			result[i] = _params[1 + (_num_sampled + i) % (_num_params - 1)];
		}
	}
	_num_sampled += num_samples;
}

double Noise::SingleSample()
{
	double buf[1];
	Sample(1, buf);
	return buf[0];
}
