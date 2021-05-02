#pragma once
#ifndef NOISE_H
#define NOISE_H

#include "stdafx.h"
#include "Constants.h"
#include "ConfigReader.h"





class Noise
{
public:
	Noise();
	Noise(int type, double avg = 0.0);
	Noise(int type, const double* params, double avg = 0.0);
	~Noise();

	double GetAverage() const;
	int GetType() const;
	double* GetParameters() const;
	int GetParameters(double* res) const;
	double GetParameter(int index) const;
	void Sample(int num_samples, double* result);
	double SingleSample();

	//void Initialize(ConfigParams &conf_reader);
	void Initialize(int type, double avg = 0.0, const double* params = NULL, bool deep_copy = false);
	void CopyFrom(const Noise &old_noise, bool deep_copy = false);
	int NumParams(bool force_calc = false, const double* params = NULL);

	void Reset();

private:
	double _avg;
	int _type;
	int _num_params;
	double* _params = NULL;
	int64_t _num_sampled;

};


#endif // !NOISE_H
