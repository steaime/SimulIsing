#pragma once
#ifndef SHARED_FUNCTIONS_H
#define SHARED_FUNCTIONS_H

#include "stdafx.h"
#include "Constants.h"

int ComplexRoots_P3(double* coeffs, std::complex<double>* res);
int RealRoots_P3(double* coeffs, double* res);
int RealRoots_P3(double C3, double C2, double C1, double C0, double* res);
int RealElements(std::complex<double>* v, int len, double* res);
bool IsReal(std::complex<double> z);

double CalcSum(double* v, int len);
double CalcSum(double** v, int len1, int* len2);
double CalcAverage(double v1, double v2, int type = AvgType::ARITMETIC);
double CalcAverage(double* v, int len, int type = AvgType::ARITMETIC);
double CalcAverage(double** v, int len1, int len2, int type = AvgType::ARITMETIC);
double CalcAverage(double** v, int len1, int* len2, int type = AvgType::ARITMETIC);
double CalcAverage(double*** v, int len1, int len2, int len3, int type = AvgType::ARITMETIC);
double CalcAvgPower(double* v, int len, double power, int type = AvgType::ARITMETIC);
double CalcVariance(double* v, int len);
double CalcVariance(double** v, int len1, int len2);
double CalcVariance(double** v, int len1, int* len2);
double CalcStd(double* v, int len);
int CalcMaxMin(double* v, int len, double &max, double &min, double exclude = BAD_VALUE);
int CalcGlobalMaxMin(double** v, int len1, int len2, double& max, double& min, double exclude = BAD_VALUE);
int CalcGlobalMaxMin(double** v, int len1, int* len2, double& max, double& min, double exclude = BAD_VALUE);
int CalcGlobalMaxMin(double*** v, int len1, int* len2, int** len3, double &max, double &min, double exclude = BAD_VALUE);
int CalcProduct(int* v, int len);
int CountGreaterThan(double* v, int len, double dThr);
int CountSmallerThan(double* v, int len, double dThr);

double CalcAvgSquareDiff(double* v1, double* v2, int len);
void CalcColumnAverage(double** v, int len1, int len2, double* res);

int LogSpace_CalcNum(double start, double end, int ppd);
double LogSpace_CalcPPD(double start, double end, int num);
int LogSpace(double start, double end, double ppd, double* res, int maxnum = -1, bool force_num = false);
double LogSpaceNum(double start, double end, int num, double* res);
int LinSpace(double start, double end, double step, double* res, int maxnum = -1);
double LinSpaceNum(double start, double end, int num, double* res);

void ArrayGreaterThan(double* v, int len, double threshold, bool* res);

int Log10Array(double* v, int len, double* res);
int Log10Array(double** v, int len1, int len2, double** res);
int Log10Array(double** v, int len1, int* len2, double** res);
int Log10Array(double*** v, int len1, int len2, int len3, double*** res);
int Log10Array(double*** v, int len1, int* len2, int** len3, double*** res);

int Flatten2DArray(double** v, int len1, int* len2, double* res);
int Flatten2DArray(double** v, int* shape, double* res);
int Flatten2DArray(double** v, int len1, int len2, double* res);
int Flatten3DArray(double*** v, int len1, int* len2, int** len3, double* res);
int Flatten3DArray(double*** v, int* shape, double* res);
int Flatten3DArray(double*** v, int len1, int len2, int len3, double* res);

double CalcIQR(double* v, int len, bool sorted = false);
double HistogramAutoBinSize(double* vals, int nvals, bool sorted = false);
int HistogramAutoBinNumber(double* vals, int nvals, bool sorted = false);
int HistogramAutoBins(double* vals, int nvals, double* res, bool sorted = false);
int FindBin(double val, double* bins, int nbins);
int BuildHistogram(double* vals, int nvals, double* bins, double* cols, int &num_bins, bool norm = true, bool logx = true, bool auto_bins = true, bool sorted = false);

int AutoCorrelation(bool* v, int* shape, double* res, int* x = NULL, int num_x = -1, bool pbc = true, bool normalize = true, int dimensions = N_DIM);
int AutoCorrelation(double* v, int* shape, double* res, int* x = NULL, int num_x = -1, bool pbc = true, bool normalize = true, bool log10 = false, int dimensions = N_DIM);
int Correlate1D(double* v, int len, double* res, int* x = NULL, int num_x = -1, bool pbc = true, bool normalize = true);
int Correlate2D(double* v, int* shape, double* res, int* x = NULL, int num_x = -1, bool pbc = true, bool normalize = true);
int Correlate2D(double** v, int len1, int len2, double* res, int* x = NULL, int num_x = -1, bool pbc = true, bool normalize = true);
int Correlate3D(double* v, int len1, int len2, int len3, double* res, int* x = NULL, int num_x = -1, bool pbc = true, bool normalize = true);
int Correlate3D(double*** v, int len1, int len2, int len3, double* res, int* x = NULL, int num_x = -1, bool pbc = true, bool normalize = true);

void ExportDataToFile(double* x, double* y, int num_data, std::string out_file);
void ExportDataToRawBinary(bool* x, int* shape, int dimensions, std::string out_file);
void ExportDataToRawBinary(double* x, int* shape, int dimensions, std::string out_file);
void ImportDataFromRawBinary(std::string in_file, double* x, int64_t &num_vals, int64_t max_vals = -1);

std::fstream OpenRawBinary(std::string in_file, char* hdr_bytes, int hdr_len);
std::fstream OpenRawBinary(std::string in_file, int32_t &num_imgs, int64_t &num_px);
int32_t Int32FromCharArray(char* bytes);
int64_t Int64FromCharArray(char* bytes);
void DoubleToCharArray(double value, char* bytes);
double DoubleFromCharArray(char* bytes);
void DoubleFromCharArray(char* bytes, double &val);

void ByteToBool(unsigned char c, bool b[8]);
cimg_library::CImg<unsigned char> ShowRawBinaryImage(std::string in_file, int* shape = NULL, int dimensions = N_DIM);
cimg_library::CImg<unsigned char> ShowRawBinaryDoubleImage(std::string in_file, double min, double max, int* shape = NULL, int dimensions = N_DIM);
void LoadRawBinaryImages(std::string* in_file, int* shape, int dimensions = -1, cimg_library::CImg<unsigned char>* res = NULL);

int FileCountLines(std::string sFileName);
int FileLoadValues(std::string sFileName, double* pdValues);
int FileLoadMultiColumn(std::string sFileName, int iNumCol, double** ppdValues);
int FileLoadSingleColumn(std::string sFileName, int iNumCol, int iColIdx, double* pdValues);

#endif // !SHARED_FUNCTIONS_H
