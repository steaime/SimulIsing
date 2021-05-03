#include "SharedFunctions.h"

int ComplexRoots_P3(double* coeffs, std::complex<double>* res){
	double _D0 = coeffs[1] * coeffs[1] - 3 * coeffs[0] * coeffs[2];
	double _D1 = 2 * coeffs[1] * coeffs[1] * coeffs[1] - 9 * coeffs[0] * coeffs[1] * coeffs[2] + 27 * coeffs[0] * coeffs[0] * coeffs[3];
	std::complex<double> _D2 = _D1 * _D1 - 4 * _D0 * _D0 * _D0;
	std::complex<double> _Cmod = std::pow(0.5 * (_D1 + std::pow(_D2, 0.5)), 1.0 / 3);
	for (int i = 0; i < 3; i++) {
		res[i] = -(1.0 / (3 * coeffs[0])) * (coeffs[1] + _Cmod * ComplexCubicRoot1[i] + _D0 / (_Cmod * ComplexCubicRoot1[i]));
	}
	return 0;
}

int RealRoots_P3(double* coeffs, double* res) {
	return RealRoots_P3(coeffs[0], coeffs[1], coeffs[2], coeffs[3], res);
}

int RealRoots_P3(double C3, double C2, double C1, double C0, double* res)
{
	/// DEBUG THIS: IT CLEARLY DOESN'T WORK!!!!

	int RootNumber = 0;
	double _D0 = C2 * C2 - 3 * C3 * C1;
	double _D1 = 2 * C2 * C2 * C2 - 9 * C3 * C2 * C1 + 27 * C3 * C3 * C0;
	std::complex<double> _D2 = _D1 * _D1 - 4 * _D0 * _D0 * _D0;
	std::complex<double> _Cmod = std::pow(0.5 * (_D1 + std::pow(_D2, 0.5)), 1.0 / 3);
	for (int i = 0; i < 3; i++) {
		std::complex<double> cur_root = -(1.0 / (3 * C3)) * (C2 + _Cmod * ComplexCubicRoot1[i] + _D0 / (_Cmod * ComplexCubicRoot1[i]));
		if (IsReal(cur_root)) {
			res[i] = std::real(cur_root);
			RootNumber++;
		}
	}
	return RootNumber;
}

int RealElements(std::complex<double>* v, int len, double* res) {
	int iCount = 0;
	for (int i = 0; i < len; i++) {
		if (IsReal(v[i])) {
			res[iCount] = std::real(v[i]);
			iCount++;
		}
	}
	return iCount;
}

bool IsReal(std::complex<double> z) {
	return (std::abs(std::imag(z)) < MIN_IMAG_PART);
}

double CalcSum(double* v, int len) {
	if (len >= 0) {
		double res = 0;
		for (int i = 0; i < len; i++) {
			res += v[i];
		}
		return res;
	}
	else {
		return NULL;
	}
}

double CalcSum(double** v, int len1, int* len2) {
	if (len1 >= 0) {
		double res = 0;
		for (int i = 0; i < len1; i++) {
			for (int j = 0; j < len2[i]; j++) {
				res += v[i][j];
			}
		}
		return res;
	}
	else {
		return NULL;
	}
}

double CalcAverage(double v1, double v2, int type)
{
	switch (type) {
	case AvgType::ARITMETIC:
		return 0.5 * (v1 + v2);
		break;
	case AvgType::HARMONIC:
		return 2 / (1 / v1 + 1 / v2);
		break;
	case AvgType::GEOMETRIC:
		return pow(v1 * v2, 0.5);
		break;
	case AvgType::MAXIMUM:
		return std::max<double>(v1, v2);
		break;
	case AvgType::MINIMUM:
		return std::min<double>(v1, v2);
		break;
	default:
		return BAD_VALUE;
	}
}

double CalcAverage(double* v, int len, int type) {
	if (len >= 0) {
		double res = 0;
		if (type == AvgType::ARITMETIC) {
			for (int i = 0; i < len; i++) {
				res += v[i];
			}
			return res / len;
		}
		else if (type == AvgType::HARMONIC) {
			for (int i = 0; i < len; i++) {
				res += 1 / v[i];
			}
			return len / res;
		}
		else if (type == AvgType::GEOMETRIC) {
			// CHECK
			res = 1;
			for (int i = 0; i < len; i++) {
				res *= v[i];
			}
			return pow(res, 1.0 / len);
		}
	}
	else {
		return NULL;
	}
}

double CalcAverage(double** v, int len1, int len2, int type) {
	if (len1 >= 0) {
		double res = 0;
		int count = 0;
		if (type == AvgType::ARITMETIC) {
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2; j++) {
					res += v[i][j];
					count++;
				}
			}
		}
		else if (type == AvgType::HARMONIC) {
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2; j++) {
					res += 1 / v[i][j];
					count++;
				}
			}
		}
		else if (type == AvgType::GEOMETRIC) {
			res = 1;
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2; j++) {
					res *= v[i][j];
					count++;
				}
			}
		}
		if (count > 0) {
			if (type == AvgType::ARITMETIC) {
				return res / count;
			}
			else if (type == AvgType::HARMONIC) {
				return count / res;
			}
			else if (type == AvgType::GEOMETRIC) {
				return pow(res, 1 / count);
			}
		}
		else {
			return -1;
		}
	}
	else {
		return NULL;
	}
}

double CalcAverage(double*** v, int len1, int len2, int len3, int type) {
	if (len1 >= 0) {
		double res = 0;
		int count = 0;
		if (type == AvgType::ARITMETIC) {
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2; j++) {
					for (int k = 0; k < len3; k++) {
						res += v[i][j][k];
						count++;
					}
				}
			}
		}
		else if (type == AvgType::HARMONIC) {
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2; j++) {
					for (int k = 0; k < len3; k++) {
						res += 1 / v[i][j][k];
						count++;
					}
				}
			}
		}
		else if (type == AvgType::GEOMETRIC) {
			res = 1;
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2; j++) {
					for (int k = 0; k < len3; k++) {
						res *= v[i][j][k];
						count++;
					}
				}
			}
		}
		if (count > 0) {
			if (type == AvgType::ARITMETIC) {
				return res / count;
			}
			else if (type == AvgType::HARMONIC) {
				return count / res;
			}
			else if (type == AvgType::GEOMETRIC) {
				return pow(res, 1 / count);
			}
		}
		else {
			return -1;
		}
	}
	else {
		return NULL;
	}
}

double CalcAverage(double** v, int len1, int* len2, int type) {
	if (len1 >= 0) {
		double res = 0;
		int count = 0;
		if (type == AvgType::ARITMETIC) {
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2[i]; j++) {
					res += v[i][j];
					count++;
				}
			}
		}
		else if (type == AvgType::HARMONIC) {
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2[i]; j++) {
					res += 1 / v[i][j];
					count++;
				}
			}
		}
		else if (type == AvgType::GEOMETRIC) {
			res = 1;
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2[i]; j++) {
					res *= v[i][j];
					count++;
				}
			}
		}
		if (count > 0) {
			if (type == AvgType::ARITMETIC) {
				return res / count;
			}
			else if (type == AvgType::HARMONIC) {
				return count / res;
			}
			else if (type == AvgType::GEOMETRIC) {
				return pow(res, 1 / count);
			}
		}
		else {
			return -1;
		}
	}
	else {
		return NULL;
	}
}

double CalcAvgPower(double* v, int len, double power, int type) {
	if (len >= 0) {
		double res = 0;
		if (type == AvgType::ARITMETIC) {
			for (int i = 0; i < len; i++) {
				res += pow(v[i], power);
			}
			return res / len;
		}
		else if (type == AvgType::HARMONIC) {
			for (int i = 0; i < len; i++) {
				res += 1 / pow(v[i], power);
			}
			return len / res;
		}
		else if (type == AvgType::GEOMETRIC) {
			res = 1;
			for (int i = 0; i < len; i++) {
				res *= pow(v[i], power);
			}
			return pow(res, 1 / len);;
		}
	}
	else {
		return NULL;
	}
}


double CalcVariance(double* v, int len) {
	if (len <= 0) {
		return -1;
	}
	else if (len == 1) {
		return 0;
	}
	else {
		double mean = CalcAverage(v, len);
		double res = 0;
		for (int i = 0; i < len; i++) {
			res += (v[i] - mean)*(v[i] - mean);
		}
		return res / len;
	}
}

double CalcVariance(double** v, int len1, int len2)
{
	int num_tot = len1 * len2;
	if (num_tot <= 0) {
		return -1;
	}
	else if (num_tot == 1) {
		return 0;
	}
	else {
		double mean = CalcAverage(v, len1, len2);
		double res = 0;
		for (int i = 0; i < len1; i++) {
			for (int j = 0; j < len2; j++) {
				res += (v[i][j] - mean) * (v[i][j] - mean);
			}
		}
		return res / num_tot;
	}
}

double CalcVariance(double** v, int len1, int* len2)
{
	int num_tot = 0;
	for (int i = 0; i < len1; i++) {
		num_tot += len2[i];
	}
	if (num_tot <= 0) {
		return -1;
	}
	else if (num_tot == 1) {
		return 0;
	}
	else {
		double mean = CalcAverage(v, len1, len2);
		double res = 0;
		for (int i = 0; i < len1; i++) {
			for (int j = 0; j < len2[i]; j++) {
				res += (v[i][j] - mean) * (v[i][j] - mean);
			}
		}
		return res / num_tot;
	}
}

double CalcStd(double* v, int len) {
	if (len <= 0) {
		return -1;
	}
	else {
		return sqrt(CalcVariance(v, len));
	}
}

int CalcMaxMin(double* v, int len, double &max, double &min, double exclude) {
	if (len > 0) {
		max = v[0];
		min = v[0];
		bool valid = (exclude == BAD_VALUE || max != exclude);
		for (int i = 0; i < len; i++) {
			if (exclude == BAD_VALUE || v[i] != exclude) {
				if (valid) {
					if (max < v[i]) {
						max = v[i];
					}
					if (min > v[i]) {
						min = v[i];
					}
				}
				else {
					max = v[i];
					min = v[i];
					valid = (exclude == BAD_VALUE || max != exclude);
				}
			}
		}
		return 0;
	}
	else {
		return -1;
	}
}

int CalcGlobalMaxMin(double** v, int len1, int len2, double &max, double &min, double exclude) {
	if (len1 > 0 && len2 > 0) {
		max = v[0][0];
		min = v[0][0];
		bool valid = (exclude == BAD_VALUE || max != exclude);
		for (int i = 0; i < len1; i++) {
			for (int j = 0; j < len2; j++) {
				if (exclude == BAD_VALUE || v[i][j] != exclude) {
					if (valid) {
						if (max < v[i][j]) {
							max = v[i][j];
						}
						if (min > v[i][j]) {
							min = v[i][j];
						}
					}
					else {
						max = v[i][j];
						min = v[i][j];
						valid = (exclude == BAD_VALUE || max != exclude);
					}
				}
			}
		}
		return 0;
	}
	else {
		return -1;
	}
}

int CalcGlobalMaxMin(double** v, int len1, int* len2, double& max, double& min, double exclude)
{
	if (len1 > 0) {
		max = v[0][0];
		min = v[0][0];
		bool valid = (exclude == BAD_VALUE || max != exclude);
		for (int i = 0; i < len1; i++) {
			for (int j = 0; j < len2[i]; j++) {
				if (exclude == BAD_VALUE || v[i][j] != exclude) {
					if (valid) {
						if (max < v[i][j]) {
							max = v[i][j];
						}
						if (min > v[i][j]) {
							min = v[i][j];
						}
					}
					else {
						max = v[i][j];
						min = v[i][j];
						valid = (exclude == BAD_VALUE || max != exclude);
					}
				}
			}
		}
		return 0;
	}
	else {
		return -1;
	}
}

int CalcGlobalMaxMin(double*** v, int len1, int* len2, int** len3, double &max, double &min, double exclude) {
	max = v[0][0][0];
	min = v[0][0][0];
	bool valid = (exclude == BAD_VALUE || max != exclude);
	for (int i = 0; i < len1; i++) {
		for (int j = 0; j < len2[i]; j++) {
			for (int k = 0; k < len3[i][j]; k++) {
				if (exclude == BAD_VALUE || v[i][j][k] != exclude) {
					if (valid) {
						if (max < v[i][j][k]) {
							max = v[i][j][k];
						}
						if (min > v[i][j][k]) {
							min = v[i][j][k];
						}
					}
					else {
						max = v[i][j][k];
						min = v[i][j][k];
						valid = (exclude == BAD_VALUE || max != exclude);
					}
				}
			}
		}
	}
	return 0;
}

int CalcProduct(int* v, int len) {
	int res = 1;
	for (int i = 0; i < len; i++) {
		res *= v[i];
	}
	return res;
}

int CountGreaterThan(double * v, int len, double dThr)
{
	int res = 0;
	for (int i = 0; i < len; i++) {
		if (v[i] > dThr) {
			res++;
		}
	}
	return res;
}

int CountSmallerThan(double * v, int len, double dThr)
{
	int res = 0;
	for (int i = 0; i < len; i++) {
		if (v[i] < dThr) {
			res++;
		}
	}
	return res;
}

double CalcAvgSquareDiff(double* v1, double* v2, int len) {
	double res = 0;
	for (int i = 0; i < len; i++) {
		res += (v1[i] - v2[i])*(v1[i] - v2[i]);
	}
	return res;
}

void CalcColumnAverage(double ** v, int len1, int len2, double * res)
{
	for (int j = 0; j < len2; j++) {
		res[j] = v[0][j] / len1;
	}
	for (int i = 1; i < len1; i++) {
		for (int j = 0; j < len2; j++) {
			res[j] += v[i][j] / len1;
		}
	}
}

int LogSpace_CalcNum(double start, double end, int ppd) {
	return (int)(abs(log10(end / start) * ppd)) + 1;
}

double LogSpace_CalcPPD(double start, double end, int num) {
#if DEBUG_MODE:
	std::cout << end << " - " << start << " - " << log10(end / start) << " - " << ((num - 1) * 1.0 / abs(log10(end / start))) << std::endl;
#endif
	return ((num - 1) * 1.0 / abs(log10(end / start)));
}

int LogSpace(double start, double end, double ppd, double* res, int maxnum, bool force_num) {
#if DEBUG_MODE:
	int num = LogSpace_CalcNum(start, end, ppd);
#endif
	double cur_v = start;
	bool keep_on = true;
	int count = 0;
	while (keep_on) {
		res[count] = cur_v;
		if (start < end) {
			cur_v *= pow(10, 1.0 / ppd);
			keep_on = (cur_v < end);
		}
		else {
			cur_v /= pow(10, 1.0 / ppd);
			keep_on = (cur_v > end);
		}
		count++;
		if (maxnum > 0) {
			if (force_num) {
				keep_on = (count < maxnum);
			}
			else if (count >= maxnum) {
				keep_on = false;
			}
		}
	}
#if DEBUG_MODE:
	if (count != num) {
		std::cout << "[DEBUG]: warning: number of elements returned by LogSpace() function changed from " << num << " to " << count << std::endl;
	}
#endif
	return count;
}

double LogSpaceNum(double start, double end, int num, double* res) {
	double ppd = LogSpace_CalcPPD(start, end, num);
	LogSpace(start, end, ppd, res, num, true);
	return ppd;
}

int LinSpace(double start, double end, double step, double* res, int maxnum) {
	if (step == 0) {
		return -1;
	}
	else if (step > 0 && end < start) {
		return -2;
	}
	else if (step < 0 && end > start) {
		return -3;
	}
	else {
		double cur_val = start;
		int count = 0;
		if (step > 0) {
			while (cur_val <= end) {
				res[count] = cur_val;
				count++;
				cur_val += step;
				if (maxnum > 0 && count >= maxnum) {
					break;
				}
			}
		}
		else {
			while (cur_val >= end) {
				res[count] = cur_val;
				count++;
				cur_val += step;
				if (maxnum > 0 && count >= maxnum) {
					break;
				}
			}
		}
		return count;
	}
}

double LinSpaceNum(double start, double end, int num, double* res) {
	double step = (end - start) * 1.0 / num;
	LinSpace(start, end, step, res, num);
	return step;
}

void ArrayGreaterThan(double* v, int len, double threshold, bool* res) {
	for (int i = 0; i < len; i++) {
		res[i] = (v[i] > threshold);
	}
}

int Log10Array(double* v, int len, double* res) {
	int count_bad = 0;
	for (int i = 0; i < len; i++) {
		if (v[i] > 0) {
			res[i] = log10(v[i]);
		}
		else {
			res[i] = BAD_VALUE;
			count_bad++;
		}
	}
	return count_bad;
}

int Log10Array(double** v, int len1, int len2, double** res) {
	int count_bad = 0;
	for (int i = 0; i < len1; i++) {
		for (int j = 0; j < len2; j++) {
			if (v[i][j] > 0) {
				res[i][j] = log10(v[i][j]);
			}
			else {
				res[i][j] = BAD_VALUE;
				count_bad++;
			}
		}
	}
	return count_bad;
}

int Log10Array(double** v, int len1, int* len2, double** res)
{
	int count_bad = 0;
	for (int i = 0; i < len1; i++) {
		for (int j = 0; j < len2[i]; j++) {
			if (v[i][j] > 0) {
				res[i][j] = log10(v[i][j]);
			}
			else {
				res[i][j] = BAD_VALUE;
				count_bad++;
			}
		}
	}
	return count_bad;
}

int Log10Array(double*** v, int len1, int len2, int len3, double*** res) {
	int count_bad = 0;
	for (int i = 0; i < len1; i++) {
		for (int j = 0; j < len2; j++) {
			for (int k = 0; k < len3; k++) {
				if (v[i][j][k] > 0) {
					res[i][j][k] = log10(v[i][j][k]);
				}
				else {
					res[i][j][k] = BAD_VALUE;
					count_bad++;
				}
			}
		}
	}
	return count_bad;
}

int Log10Array(double*** v, int len1, int* len2, int** len3, double*** res) {
	int count_bad = 0;
	for (int i = 0; i < len1; i++) {
		for (int j = 0; j < len2[i]; j++) {
			for (int k = 0; k < len3[i][j]; k++) {
				if (v[i][j][k] > 0) {
					res[i][j][k] = log10(v[i][j][k]);
				}
				else {
					res[i][j][k] = BAD_VALUE;
					count_bad++;
				}
			}
		}
	}
	return count_bad;
}

double CalcIQR(double* v, int len, bool sorted) {
	double *vect;
	if (sorted == false) {
		vect = new double[len];
		for (int i = 0; i < len; i++) {
			vect[i] = v[i];
		}
		std::sort(vect, vect + len);
	}
	else {
		vect = v;
	}
	double first_quart = vect[len/4];
	double third_quart = vect[3*len/4];
	if (sorted == false) {
		delete[] vect;
	}
	vect = NULL;
	return third_quart - first_quart;
}

double HistogramAutoBinSize(double* vals, int nvals, bool sorted) {
	return 2 * CalcIQR(vals, nvals, sorted) / pow(nvals, 1.0/3);
}

int HistogramAutoBinNumber(double* vals, int nvals, bool sorted) {
	double max, min;
	CalcMaxMin(vals, nvals, max, min);
	int res = (int)((max - min) / HistogramAutoBinSize(vals, nvals, sorted)) + 1;
	return res;
}

int HistogramAutoBins(double* vals, int nvals, double* res, bool sorted) {
	double max, min;
	CalcMaxMin(vals, nvals, max, min);
	double h = HistogramAutoBinSize(vals, nvals, sorted);
	int num = (int)((max - min) / h) + 1;
	LinSpaceNum(min, max+h, num, res);
	return num;
}

int FindBin(double val, double* bins, int nbins) {
	int idx = 0;
	for (idx = 0; idx < nbins - 1; idx++) {
		if (bins[idx] <= val && bins[idx + 1] > val) {
			break;
		}
	}
	if (idx == nbins - 1 && val == bins[nbins - 2]) {
		idx = nbins - 2;
	}
	return idx;
}

int BuildHistogram(double* vals, int nvals, double* bins, double* cols, int &num_bins, bool norm, bool logx, bool auto_bins, bool sorted) {
#if DEBUG_MODE:
	std::cout << "       BuildHistogram() function called" << std::endl;
#endif
	double* data;
	if (logx) {
		data = new double[nvals];
		Log10Array(vals, nvals, data);
	}
	else {
		data = vals;
	}
	if (auto_bins) {
		num_bins = HistogramAutoBins(data, nvals, bins, sorted);
	}
	for (int i = 0; i < num_bins; i++) {
		cols[i] = 0;
	}
	for (int i = 0; i < nvals; i++) {
		int cur_idx = FindBin(vals[i], bins, num_bins);
		cols[cur_idx] += 1;
	}
	if (norm) {
		for (int i = 0; i < num_bins; i++) {
			cols[i] *= 1.0/nvals;
		}
	}
	if (logx) {
		delete[] data;
	}
	data = NULL;
	return 0;
}

/*
Correlation: corr(x) = <v[i]*v[i+x]>-<v[i]><v[i]>
*/
int AutoCorrelation(double* v, int* shape, double* res, int* x, int num_x, bool pbc, bool normalize, bool log10, int dimensions) {
	double *y;
	if (log10) {
		int totlen = CalcProduct(shape, dimensions);
		y = new double[totlen];
		Log10Array(v, totlen, y);
	}
	else {
		y = v;
	}
	int myres = -1;
	if (dimensions == 1) {
		myres = Correlate1D(y, shape[0], res, x, num_x, pbc, normalize);
	}
	else if (dimensions == 2) {
		myres = Correlate2D(y, shape, res, x, num_x, pbc, normalize);
	}
	else if (dimensions == 3) {
		myres = Correlate3D(y, shape[0], shape[1], shape[2], res, x, num_x, pbc, normalize);
	}
	if (log10) {
		delete[] y;
	}
	y = NULL;
	return myres;
}

int AutoCorrelation(bool* v, int* shape, double* res, int* x, int num_x, bool pbc, bool normalize, int dimensions) {
	int totlen = CalcProduct(shape, dimensions);
	double* data = new double[totlen];
	for (int i = 0; i < totlen; i++) {
		if (v[i]) {
			data[i] = 1;
		}
		else {
			data[i] = 0;
		}
	}
	int result = AutoCorrelation(data, shape, res, x, num_x, pbc, normalize, false, dimensions);
	delete[] data;
	data = NULL;
	return result;
}

int Correlate1D(double* v, int len, double* res, int* x, int num_x, bool pbc, bool normalize)
{
	if (num_x < 0) {
		if (pbc) {
			num_x = len/2;
		}
		else {
			num_x = len;
		}
		for (int i = 0; i < num_x; i++) {
			x[i] = i;
		}
	}

	double avg = CalcAverage(v, len);
	for (int j = 0; j < num_x; j++) {
		res[j] = 0;
		int count_val = 0;
		for (int i = 0; i < len; i++) {
			if ((i + x[j] >= len) && (pbc == false)) {
				break;
			}
			res[j] += v[i] * v[(i + x[j]) % len];
			count_val++;
			if (pbc) {
				res[j] += v[i] * v[(i + len - x[j]) % len];
				count_val++;
			}
		}
		res[j] = res[j] / count_val - avg * avg;
		if (normalize && avg != 0) {
			res[j] = res[j] / (avg * avg);
		}
	}
	return num_x;
}

/*
len1 is the size of first index (row), len2 is the size of second index (column)
*/
int Correlate2D(double* v, int* shape, double* res, int* x, int num_x, bool pbc, bool normalize)
{
	int len1 = shape[0];
	int len2 = shape[1];
	if (num_x < 0) {
		if (pbc) {
			num_x = len1 / 2;
		}
		else {
			num_x = len1;
		}
		for (int i = 0; i < num_x; i++) {
			x[i] = i;
		}
	}

	double avg = CalcAverage(v, len1*len2);
	for (int k = 0; k < num_x; k++) {
		res[k] = 0;
		int count_val = 0;
		for (int i = 0; i < len1; i++) {
			if ((i + x[k] >= len1) && (pbc == false)) {
				break;
			}
			for (int j = 0; j < len2; j++) {
				res[k] += v[i*len2+j] * v[((i + x[k]) % len1)*len2+j];
				count_val++;
			}
			if (pbc) {
				for (int j = 0; j < len2; j++) {
					res[k] += v[i*len2+j] * v[((i + len1 - x[k]) % len1)*len2+j];
					count_val++;
				}
			}
		}
		for (int j = 0; j < len2; j++) {
			if ((j + x[k] >= len2) && (pbc == false)) {
				break;
			}
			for (int i = 0; i < len1; i++) {
				res[k] += v[i*len2+j] * v[i*len2 + (j + x[k]) % len2];
				count_val++;
			}
			if (pbc) {
				for (int i = 0; i < len1; i++) {
					res[k] += v[i*len2+j] * v[i*len2 + (j + len2 - x[k]) % len2];
					count_val++;
				}
			}
		}
		res[k] = res[k] / count_val - avg * avg;
		if (normalize && avg != 0) {
			res[k] = res[k] / (avg * avg);
		}
	}
	return num_x;
}

int Correlate2D(double** v, int len1, int len2, double* res, int* x, int num_x, bool pbc, bool normalize)
{
	if (num_x < 0) {
		if (pbc) {
			num_x = len1 / 2;
		}
		else {
			num_x = len1;
		}
		for (int i = 0; i < num_x; i++) {
			x[i] = i;
		}
	}

	double avg = CalcAverage(v, len1, len2);
	for (int k = 0; k < num_x; k++) {
		res[k] = 0;
		int count_val = 0;
		for (int i = 0; i < len1; i++) {
			if ((i + x[k] >= len1) && (pbc == false)) {
				break;
			}
			for (int j = 0; j < len2; j++) {
				res[k] += v[i][j] * v[(i + x[k]) % len1][j];
				count_val++;
			}
			if (pbc) {
				for (int j = 0; j < len2; j++) {
					res[k] += v[i][j] * v[(i + len1 - x[k]) % len1][j];
					count_val++;
				}
			}
		}
		for (int j = 0; j < len2; j++) {
			if ((j + x[k] >= len2) && (pbc == false)) {
				break;
			}
			for (int i = 0; i < len1; i++) {
				res[k] += v[i][j] * v[i][(j + x[k]) % len2];
				count_val++;
			}
			if (pbc) {
				for (int i = 0; i < len1; i++) {
					res[k] += v[i][j] * v[i][(j + len2 - x[k]) % len2];
					count_val++;
				}
			}
		}
		res[k] = res[k] / count_val - avg * avg;
		if (normalize && avg != 0) {
			res[k] = res[k] / (avg * avg);
		}
	}
	return num_x;
}

int Correlate3D(double* v, int len1, int len2, int len3, double* res, int* x, int num_x, bool pbc, bool normalize) {
	//TODO: to be implemented!
	return -1;
}

int Correlate3D(double*** v, int len1, int len2, int len3, double* res, int* x, int num_x, bool pbc, bool normalize) {
	//TODO: to be implemented!
	return -1;
}

void ExportDataToFile(double* x, double* y, int num_data, std::string out_file) {
	std::ofstream fout(out_file);
	for (int i = 0; i < num_data; i++) {
		fout << x[i] << "\t" << y[i] << "\n";
	}
	fout.close();
}

void ExportDataToRawBinary(bool* x, int* shape, int dimensions, std::string out_file) {
	int32_t n_imgs;
	if (shape[dimensions - 1] % 8 != 0) {
		std::cout << "WARNING: number of pixels in a row (" << shape[dimensions - 1] << ") is not multiple of 8. Raw output will be trimmed." << std::endl;
	}
	int64_t tot_num_px = shape[dimensions - 1];
	int64_t trim_num_px = (long long)(shape[dimensions - 1] / 8) * 8;
	if (dimensions > 1) {
		n_imgs = shape[0];
		for (int i = 1; i < dimensions - 1; i++) {
			trim_num_px *= shape[i];
			tot_num_px *= shape[i];
		}
	}
	else {
		n_imgs = 1;
	}
	int64_t bytes_per_frame = trim_num_px / 8;
	int64_t num_bytes = n_imgs * bytes_per_frame;

#if VERBOSE:
	if (dimensions > 1) {
		std::cout << "Writing " << n_imgs << " images, " << trim_num_px << " (out of " << tot_num_px << ") pixels per image, total " << num_bytes << " bytes to file " << out_file << std::endl;
	}
	else {
		std::cout << "Writing " << trim_num_px << "/" << tot_num_px << " pixels (" << num_bytes << " bytes) to file " << out_file << std::endl;
	}
#endif

	std::fstream myfile;
	myfile = std::fstream(out_file, std::ios::out | std::ios::binary);
	myfile.write(reinterpret_cast<const char *>(&n_imgs), sizeof(n_imgs));
	myfile.write(reinterpret_cast<const char *>(&trim_num_px), sizeof(trim_num_px));
	int64_t cur_el_idx;
	int64_t cur_byte_idx;
	int cur_num;
	char cur_char;
	char* byte_data = new char[num_bytes];
	for (int i = 0; i < num_bytes; i++) {
		byte_data[i] = 0;
	}
	for (int i = 0; i < n_imgs; i++) {
		cur_char = 0;
		cur_num = 0;
		cur_byte_idx = i * bytes_per_frame;
		cur_el_idx = i * tot_num_px;
		for (int j = 0; j < trim_num_px; j++) {
			cur_char += x[cur_el_idx] << cur_num;
			cur_num++;
			cur_el_idx++;
			if (cur_num >= 8) {
				byte_data[cur_byte_idx] = cur_char;
				cur_char = 0;
				cur_num = 0;
				cur_byte_idx++;
			}
		}
	}
	myfile.write(byte_data, num_bytes);
	delete[] byte_data;
	byte_data = NULL;
	myfile.close();

}


int Flatten2DArray(double** v, int len1, int len2, double* res) {
	int shape[2] = { len1, len2 };
	return Flatten2DArray(v, shape, res);
}

int Flatten2DArray(double** v, int* shape, double* res) {
	int* len2 = new int[shape[0]];
	for (int i = 0; i < shape[0]; i++) {
		len2[i] = shape[1];
	}
	int ret = Flatten2DArray(v, shape[0], len2, res);
	delete[] len2;
	return ret;
}

int Flatten2DArray(double** v, int len1, int* len2, double* res){
	int count = 0;
	for (int i = 0; i < len1; i++) {
		for (int j = 0; j < len2[i]; j++) {
			res[count] = v[i][j];
			count++;
		}
	}
	return count;
}

int Flatten3DArray(double*** v, int len1, int len2, int len3, double* res) {
	int shape[3] = { len1, len2, len3 };
	return Flatten3DArray(v, shape, res);
}

int Flatten3DArray(double*** v, int* shape, double* res) {
	int* len2 = new int[shape[0]];
	int** len3 = new int*[shape[0]];
	for (int i = 0; i < shape[0]; i++) {
		len2[i] = shape[1];
		len3[i] = new int[shape[1]];
		for (int j = 0; j < shape[1]; j++) {
			len3[i][j] = shape[2];
		}
	}
	int ret = Flatten3DArray(v, shape[0], len2, len3, res);
	delete[] len2;
	for (int i = 0; i < shape[0]; i++) {
		delete[] len3[i];
	}
	delete[] len3;
	return ret;
}

int Flatten3DArray(double*** v, int len1, int* len2, int** len3, double* res) {
	int count = 0;
	for (int i = 0; i < len1; i++) {
		for (int j = 0; j < len2[i]; j++) {
			for (int k = 0; k < len3[i][j]; k++) {
				res[count] = v[i][j][k];
				count++;
			}
		}
	}
	return count;
}


void DoubleToCharArray(double value, char* bytes) {
	char* b = reinterpret_cast<char*>(&value);
	for (int i = 0; i < 8; i++) {
		bytes[i] = b[i];
	}
	// TODO: test this!!!
}

void ExportDataToRawBinary(double* x, int* shape, int dimensions, std::string out_file) {
	int32_t n_imgs;
	int64_t tot_num_px;
	if (dimensions > 1) {
		n_imgs = shape[0];
		tot_num_px = shape[dimensions - 1];
		for (int i = 1; i < dimensions - 1; i++) {
			tot_num_px *= shape[i];
		}
	}
	else {
		n_imgs = 1;
		tot_num_px = *shape;
	}
	const int bytes_per_pixel = 8;
	int64_t bytes_per_frame = tot_num_px * bytes_per_pixel;
	int64_t num_bytes = n_imgs * bytes_per_frame;

	std::fstream myfile;
	myfile = std::fstream(out_file, std::ios::out | std::ios::binary);
	myfile.write(reinterpret_cast<const char *>(&n_imgs), sizeof(n_imgs));
	myfile.write(reinterpret_cast<const char *>(&tot_num_px), sizeof(tot_num_px));
	int64_t cur_el_idx;
	int64_t cur_byte_idx;
	int cur_num;
	char cur_char;
	char* byte_data = new char[num_bytes];
	for (int i = 0; i < n_imgs; i++) {
		for (int j = 0; j < tot_num_px; j++) {
			DoubleToCharArray(x[i * tot_num_px + j], &byte_data[(i * tot_num_px + j)*bytes_per_pixel]);
			// TODO: test this!!
		}
	}
	myfile.write(byte_data, num_bytes);
	delete[] byte_data;
	byte_data = NULL;
	myfile.close();

}

// to read all values set num_vals to -1
// after the function has finished, num_vals will be equal to the values actually read
void ImportDataFromRawBinary(std::string in_file, double* x, int64_t& num_vals, int64_t max_vals)
{
	int32_t num_imgs = 0;
	int64_t num_px = 0;
	std::fstream fin = OpenRawBinary(in_file, num_imgs, num_px);
	if (num_vals < 0) num_vals = num_imgs * num_px;
	ImportDataFromRawBinary(&fin, x, num_vals, max_vals, true);
}

void ImportDataFromRawBinary(std::fstream *in_file, double* x, int64_t& num_vals, int64_t max_vals, bool close_after, char px_format, bool swap_endian)
{
	if (max_vals > 0 && num_vals > max_vals) num_vals = max_vals;
	int64_t buflen = num_vals * sizeof(double) / sizeof(char);
	char* buffer = new char[buflen];
	in_file->read(buffer, buflen);
	if (close_after) in_file->close();
	if (px_format == 'd') {
		double* double_values = (double*)buffer;//reinterpret as doubles
		for (int i = 0; i < num_vals; i++) {
			if (swap_endian) x[i] = EndianSwap<double>(double_values[i]);
			else x[i] = double_values[i];
		}
	}
	else if (px_format == 'f') {
		float* float_values = (float*)buffer;
		for (int i = 0; i < num_vals; i++) {
			if (swap_endian) x[i] = (double)(EndianSwap<float>(float_values[i]));
			else x[i] = (double)float_values[i];
		}
	}
	else if (px_format == 'I') {
		uint32_t* uint32_values = (uint32_t*)buffer;
		for (int i = 0; i < num_vals; i++) {
			if (swap_endian) x[i] = (double)(EndianSwap<uint32_t>(uint32_values[i]));
			else x[i] = (double)uint32_values[i];
		}
	}
	else if (px_format == 'i') {
		int32_t* int32_values = (int32_t*)buffer;
		for (int i = 0; i < num_vals; i++) {
			if (swap_endian) x[i] = (double)(EndianSwap<int32_t>(int32_values[i]));
			else x[i] = (double)int32_values[i];
		}
	}
	else if (px_format == 'H') {
		uint16_t* uint16_values = (uint16_t*)buffer;
		for (int i = 0; i < num_vals; i++) {
			if (swap_endian) x[i] = (double)(EndianSwap<uint16_t>(uint16_values[i]));
			else x[i] = (double)uint16_values[i];
		}
	}
	else if (px_format == 'h') {
		int16_t* int16_values = (int16_t*)buffer;
		for (int i = 0; i < num_vals; i++) {
			if (swap_endian) x[i] = (double)(EndianSwap<int16_t>(int16_values[i]));
			else x[i] = (double)int16_values[i];
		}
	}
	else if (px_format == 'B') {
		uint8_t* uint8_values = (uint8_t*)buffer;
		for (int i = 0; i < num_vals; i++) {
			x[i] = (double)uint8_values[i];
		}
	}
	else if (px_format == 'b') {
		int8_t* int8_values = (int8_t*)buffer;
		for (int i = 0; i < num_vals; i++) {
			x[i] = (double)int8_values[i];
		}
	}
	else {
		std::cout << "ERROR in ImportDataFromRawBinary(): pixel format '" << px_format << 
			"' not recognized. Allowed formats: d (default), f, I, i, H, h, B, b. Using default format." << std::endl;
		double* def_values = (double*)buffer;
		for (int i = 0; i < num_vals; i++) {
			x[i] = def_values[i];
		}
	}
	delete[] buffer;
	buffer = NULL;
}

std::fstream OpenRawBinary(std::string in_file, char* hdr_bytes, int hdr_len) {
	std::fstream res = std::fstream(in_file, std::ios::in | std::ios::binary);
	if (hdr_len > 0) res.read(hdr_bytes, hdr_len);
	return res;
}


std::fstream OpenRawBinary(std::string in_file, int32_t& num_imgs, int64_t& num_px, int hdr_len) {
	const int hdr_size = sizeof(num_imgs) + sizeof(num_px);
	if (hdr_len < 0) hdr_len = hdr_size;
	char* hdr_bytes = NULL;
	if (hdr_len > 0) hdr_bytes = new char[hdr_len];
	std::fstream res = OpenRawBinary(in_file, hdr_bytes, hdr_len);
	if (hdr_len == 4) {
		num_imgs = 1;
		num_px = (int64_t)Int32FromCharArray(hdr_bytes);
	}
	else if (hdr_len == 8) {
		num_imgs = 1;
		num_px = Int64FromCharArray(hdr_bytes);
	}
	else if (hdr_len == 12) {
		num_imgs = Int32FromCharArray(hdr_bytes);
		num_px = Int64FromCharArray(&hdr_bytes[4]);
	}
	else {
		std::cout << "ERROR in OpenRawBinary(): invalid header size (accepted: 0|4|8|12, given:" << hdr_len << ")" << std::endl;
	}
#if VERBOSE:
	std::cout << "Raw image " << in_file << " contains " << num_imgs << " images with " << num_px << " pixels each." << std::endl;
#endif
	delete[] hdr_bytes;
	hdr_bytes = NULL;
	return res;
}

int32_t Int32FromCharArray(char* bytes) {
	return *(int32_t *)bytes;
}

int64_t Int64FromCharArray(char* bytes) {
	return *(int64_t *)bytes;
}

void ByteToBool(unsigned char c, bool b[8])
{
	for (int i = 0; i < 8; ++i)
		b[i] = (c & (1 << i)) != 0;
}

double DoubleFromCharArray(char* bytes) {
	double res;
	memcpy(&res, bytes, sizeof(double));
	return res;
}

void DoubleFromCharArray(char* bytes, double &val) {
	memcpy(&val, bytes, sizeof(double));
}

cimg_library::CImg<unsigned char> ShowRawBinaryImage(std::string in_file, int* shape, int dimensions) {
	int num_imgs = -1;
	long long num_px = -1;
	const int dim = 2;
	int* cur_shape = new int[dim];
	bool blnByte[8];
	std::fstream fin = OpenRawBinary(in_file, num_imgs, num_px);
	int num_bytes = num_imgs * num_px / 8;
	char* bytes = new char[num_bytes];
	fin.read(bytes, num_bytes);
	fin.close();
	bool* blnData = new bool[8* num_bytes];
	for (int i = 0; i < num_bytes; i++) {
		ByteToBool(bytes[i], blnByte);
		for (int j = 0; j < 8; j++) {
			blnData[i*8+j] = blnByte[j];
		}
	}
	if (shape == NULL) {
		cur_shape[0] = num_imgs;
		cur_shape[1] = num_px;
	}
	else {
		for (int i = 0; i < dimensions; i++) {
			cur_shape[i] = shape[i];
		}
		for (int i = dimensions; i < dim; i++) {
			cur_shape[i] = 1;
		}
	}
	cimg_library::CImg<unsigned char> img(cur_shape[0], cur_shape[1], 1, 1, 0);
	for (int i = 0; i < cur_shape[0]; i++) {
		for (int j = 0; j < cur_shape[1]; j++) {
			if (blnData[i*cur_shape[1]+j]) {
				img(i, j) = 255;
			}
		}
	}
	delete[] blnData;
	delete[] bytes;
	blnData = NULL;
	bytes = NULL;
	return img;
}

cimg_library::CImg<unsigned char> ShowRawBinaryDoubleImage(std::string in_file, double min, double max, int* shape, int dimensions) {
	int num_imgs = -1;
	long long num_px = -1;
	const int dim = 2;
	int* cur_shape = new int[dim];
	bool blnByte[8];
	std::fstream fin = OpenRawBinary(in_file, num_imgs, num_px);
	int num_bytes = num_imgs * num_px * 8;
	char* bytes = new char[num_bytes];
	fin.read(bytes, num_bytes);
	fin.close();
	double* mydata = new double[num_imgs * num_px];
	for (int i = 0; i < num_imgs * num_px; i++) {
		DoubleFromCharArray(&bytes[i*8], mydata[i]);
	}
	if (shape == NULL) {
		cur_shape[0] = num_imgs;
		cur_shape[1] = num_px;
	}
	else {
		for (int i = 0; i < dimensions; i++) {
			cur_shape[i] = shape[i];
		}
		for (int i = dimensions; i < dim; i++) {
			cur_shape[i] = 1;
		}
	}
	cimg_library::CImg<unsigned char> img(cur_shape[0], cur_shape[1], 1, 1, 0);
	for (int i = 0; i < cur_shape[0]; i++) {
		for (int j = 0; j < cur_shape[1]; j++) {
			if (mydata[i*cur_shape[1] + j] > max) {
				img(i, j) = 255;
			}
			else if (mydata[i*cur_shape[1] + j] < min) {
				img(i, j) = 0;
			}
			else {
				img(i, j) = 255*(mydata[i*cur_shape[1] + j] - min)/(max - min);
			}
		}
	}
	delete[] mydata;
	delete[] bytes;
	mydata = NULL;
	bytes = NULL;
	return img;
}

void LoadRawBinaryImages(std::string* file_list, int* shape, int dimensions, cimg_library::CImg<unsigned char>* res) {

	int tot_img_num = 0;
	for (int fidx = 0; fidx < shape[0]; fidx++) {

		// read raw
		int num_imgs = -1;
		long long num_px = -1;
		std::fstream fin = OpenRawBinary(file_list[fidx], num_imgs, num_px);
		long num_bytes = num_imgs * num_px / 8;
		char* bytes = new char[num_bytes];
		fin.read(bytes, num_bytes);
		fin.close();

		// convert to bool array
		bool blnByte[8];
		bool* blnData = new bool[8 * num_bytes];
		for (long i = 0; i < num_bytes; i++) {
			ByteToBool(bytes[i], blnByte);
			for (int j = 0; j < 8; j++) {
				blnData[i * 8 + j] = blnByte[j];
			}
		}

		// fill the array of cimg
		int nimgs_to_create = CalcProduct(shape, dimensions - 2);
		for (int im = 0; im < nimgs_to_create; im++) {
			res[tot_img_num].assign(shape[dimensions - 2], shape[dimensions - 1], 1, 1, 0);
			for (int i = 0; i < shape[dimensions - 2]; i++) {
				for (int j = 0; j < shape[dimensions - 1]; j++) {
					if (blnData[im * shape[dimensions - 1] * shape[dimensions - 2] + i * shape[dimensions - 1] + j]) {
						res[tot_img_num](i, j) = 255;
					}
				}
			}
			tot_img_num++;
		}
		delete[] blnData;
		delete[] bytes;
		blnData = NULL;
		bytes = NULL;
	}
}

int FileCountLines(std::string sFileName) {
	int number_of_lines = 0;
	std::string line;
	std::ifstream myfile(sFileName);
	while (std::getline(myfile, line))
		++number_of_lines;
	return number_of_lines;
	myfile.close();
}

int FileLoadValues(std::string sFileName, double* pdValues) {
	std::ifstream myfile(sFileName);
	if (!myfile.is_open()) {
		return 0; 
	}
	int iCount = 0;
	double dDummy;
	while (myfile >> dDummy)
	{
		pdValues[iCount] = dDummy;
		iCount++;
	}
	myfile.close();
	return iCount;
}

int FileLoadSingleColumn(std::string sFileName, int iNumCol, int iColIdx, double* pdValues)
{
	std::ifstream myfile(sFileName);
	if (!myfile.is_open()) {
		return 0;
	}
	int iCount = 0;
	double dDummy;
	while (myfile >> dDummy)
	{
		if ((iCount % iNumCol) == iColIdx) {
			pdValues[iCount / iNumCol] = dDummy;
		}
		iCount++;
	}
	myfile.close();
	return iCount / iNumCol;
}

bool CopyFileBinary(const char* SRC, const char* DEST)
{
	std::ifstream src(SRC, std::ios::binary);
	std::ofstream dest(DEST, std::ios::binary);
	dest << src.rdbuf();
	return src && dest;
}

bool CopyFileToFolder(std::string source_filepath, std::string dest_folder)
{
	std::string dest_fpath = dest_folder + FilenameFromPath(source_filepath);
	return CopyFileBinary((char*)source_filepath.c_str(), (char*)dest_fpath.c_str());
}

std::string FilenameFromPath(std::string full_path, bool remove_ext)
{
	std::string res(full_path);

	// Remove directory if present.
	// Do this before extension removal incase directory has a period character.
	const size_t last_slash_idx = res.find_last_of("\\/");
	if (std::string::npos != last_slash_idx) res.erase(0, last_slash_idx + 1);
	
	// Remove extension if present.
	if (remove_ext) {
		const size_t period_idx = res.rfind('.');
		if (std::string::npos != period_idx) res.erase(period_idx);
	}
	
	return res;
}

std::string GetExePath()
{
	char path[MAX_PATH];
	GetModuleFileNameA(NULL, path, MAX_PATH);
	return std::string(path);
}

std::string GetExeDir()
{
	return GetParentDirectory(GetExePath(), 1);
}

std::string GetParentDirectory(std::string folder_path, int num_generations)
{
	if (num_generations == 0) return folder_path;
	else {
		std::vector<int> chloc;
		for (int i = 0; i < folder_path.size(); i++)
			if (folder_path[i] == '/' || folder_path[i] == '\\')
				chloc.push_back(i);
		if (chloc.size() >= num_generations) return folder_path.substr(0, chloc[chloc.size() - num_generations]);
		else return "";
	}
}

std::string JoinPath(std::string root, std::string parent, std::string subfolder)
{
	std::string res = root;
	if (parent.size() > 0) {
		if (res.back() != '/' && res.back() != '\\') res = res + "\\";
		res = res + parent;
	}
	if (subfolder.size() > 0) {
		if (res.back() != '/' && res.back() != '\\') res = res + "\\";
		res = res + subfolder;
	}
	return res;
}

bool CheckPathRelative(std::string path_string)
{
	//std::experimental::filesystem::path my_path(path_string); // Construct the path from a string.
	//return my_path.is_relative();
	return PathIsRelativeA(path_string.c_str());
}

bool CheckFileExists(std::string file_path)
{
	std::ifstream ifile;
	ifile.open(file_path);
	if (ifile) return true;
	else return false;
}

std::string NowToString()
{
	time_t     now = std::time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *std::localtime(&now);
	// ref: http://en.cppreference.com/w/cpp/chrono/c/strftime
	std::strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
	return buf;
}

int FileLoadMultiColumn(std::string sFileName, int iNumCol, double** ppdValues)
{
	std::ifstream myfile(sFileName);
	if (!myfile.is_open()) {
		return 0; 
	}
	int iCount = 0;
	double dDummy;
	while (myfile >> dDummy)
	{
		ppdValues[iCount / iNumCol][iCount % iNumCol] = dDummy;
		iCount++;
	}
	myfile.close();
	return iCount / iNumCol;
}