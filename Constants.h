#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

// Global simulation parameters (cannot be changed at runtime)

#define VERBOSE 0
#define DEBUG_MODE 0
#define FORCE_DEEP_COPY 1	// if 1: When copying/initializing classes such as Noise, IsingParameters or MontecarloParameters,
							//       always deep copy source vectors (don't transfer reference)

#define N_DIM 2			// Number of lattice dimensions.

#define RATE_EPS 1e-20  // Minimum rate accepted (rates cannot be zero or negative).
						// When rates get 0 or negative, they are set to RATE_EPS
#define ALPHA_EPS 0     // Minimum alpha value accepted (alphas cannot be negative).
						// When alphas get 0 or negative, they are set to ALPHA_EPS
#define GLASSY_EPS 1e-7 // Minimum normalized distance from maximum glassiness
#define MAX_RANDOM_ATTEMPTS 100 // Number of attempts to generate a valid random number
								// from distribution before giving up
#define CLIP_ALPHA_POSMIN 0 // eventually, "regularize" histogram by clipping the alphas to their strictly positive global minimum
#define MEANFIELD_USE_LOCALALPHA 1  // When using force_meanfield in a simulation, evaluate local target rate using local alphas
									// (but globally averaged current rate instead of local rate)
									// If set to 0, evaluate local rate by using globally averaged rates AND alpha

#define BAD_VALUE -100
#define BAD_VALUE_2 -101


#define INI_MAXNUMDIM 4
#define INI_MAXPARAMS 201		// maximum number of noise parameters that can be read directly from INI file
#define NOISE_MAXPARAMS 2097161 // maximum number of noise parameters allowed in general (either from INI file or from RAW file)
								// this number is 512*512 (our largest simulations) *8 (number of nearest neighbors in d=4) +1 (to store the number of parameters)
#define INI_TPROT_MAXPARAMS 5



// Correlation parameters

#define CORR_BINARY 1		// Correlate the rates (0) or the solid/fluid binary state (1)
#define CORR_NORM_X0 1		// Normalize (1) / do not normalize (0) all correlation functions by the dx=0 value






#define MAX_SOLID_RATE_RATIO 5.0 // Default threshold separating "slow/solid" sites from "fast/fluid" sites
						// Site i will be slow/solid site if _rates[i] < MAX_SOLID_RATE_RATIO/DEF_TAU0
						// (or the analogous ratio calculated using runtime defined values)






// Simulation parameters

/* Random thermal jump settings: every simulation step
- takes val_current as the current value of the cell,
- computes val_target using the expression based on actual values of the neighbors
- if RANDOM_NOISE_LOG==0:	generates a random number with Gaussian probability distribution
							centered around val_target with variance proportional to effective temperature
							specify proportionality constant. Sets the new rate to this random number
  if RANDOM_NOISE_LOG==1:	generates a random number with Gaussian probability distribution
							centered around log10(val_target) with variance proportional to effective temperature
							specify proportionality constant. Sets the new rate to pow(10, this random number)
*/
#define RANDOM_NOISE_COEFF 1.0
#define RANDOM_NOISE_LOG 1

/* Random simulation step default parameters */

#define DEF_SIM_TEMP 0.0			// Effective temperature: sets the random jump 
									// on top of deterministic dynamics (see RANDOM_NOISE_COEFF)
#define DEF_SIMSTEP_MAX_ITER 10000	// Maximum number of simulation steps
#define DEF_SIM_CONV_PREC 1e-40		// Convergence precision, to decide if simulation has converged
									// Set to 0 not to check for convergence (in this case every step
									// will have the maximum number of iterations)
									// Suggestion: for DEF_DEVIATION_LOG==0, try DEF_SIM_CONV_PREC=1e-6
#define DEF_DEVIATION_LOG 1			// If 0: deviation between current configuration and target configuration
									//       is computed using the difference of relaxation rates
									// If 1: deviation is computed using the difference of log10(rate)
									//       (to give equal weight to slow and fast mode)
//#define DEF_EQ_PROTOCOL 0			// Temperature profile during equilibration
									// check enum TemperatureProtocol in MontecarloParameters.h
#define DEF_EQ_RUN_NUM 1			// Run the whole equilibration protocol (until convergence or for DEF_SIMSTEP_MAX_ITER steps)
									// DEF_EQ_RUN_NUM times. Useful for DEF_EQ_PROTOCOL!=0

/* Parameters involving complex numbers */

#define MIN_IMAG_PART 1e-10			// Largest imaginary part for a complex number to be considered a real

static std::complex<double> ComplexCubicRoot1[3] = { 1.0, std::complex<double>{-0.5, pow(0.75, 0.5)}, std::complex<double>{-0.5, -pow(0.75, 0.5)} };





// Default output parameters
#define DEF_INI_FILE "D:\\steaime\\Documents\\Research\\Projects\\SoftGlasses\\Ising\\config\\simParams.ini"


enum TemperatureProtocol {
	CONSTANT, LINRAMP, EXPRAMP, PWRRAMP
};


enum NoiseTypes
{
	NONOISE, WHITE_FLAT, WHITE_GAUSS, WHITE_LOGNORM, WHITE_CUSTOM, SORTED_LIST
};

/*
NoiseTypes::NONOISE			No parameter required, all values equal to _avg
NoiseTypes::WHITE_FLAT		One parameter required (the noise relative variance)
							NOTE:   if average==0 the parameter will be interpreted as the
									absolute variance, instead of the relative variance.
									Random numbers generated with flat distribution in the
									[_avg-amplitude, _avg+amplitude] interval, with
									amplitude = sqrt(3 * relative_variance * _avg)
NoiseTypes::WHITE_GAUSS		One parameter required (the noise relative variance)
							Random numbers generated following a Gauss distribution
							with the given variance
NoiseTypes::WHITE_LOGNORM	One parameter required (the noise relative variance)
							Random numbers generated following a Lognormal distribution
							with the given variance
NoiseTypes::WHITE_CUSTOM    3N+2 parameters required: the first parameter is the median value
							(i.e. the value for which the cumulative distribution is 0.5)
							The second parameter is N, the number of points
							characterizing the custom cumulative distribution, followed by
							x1, C1, Cerr1, x2, C2, Cerr2, ..., xN, CN, CerrN satisfying : 
							(1) 0 <= Ci <= 1; (2) Ci <= Ci+1; (3) xi < xi+i; (4) Cerri >= 0
							Note: the maximum N is set by the hard-coded constant NOISE_MAXPARAMS.
							      Moreover, maximum number of parameters that can be read from INI file
								  is set by the hard-coded constant INI_MAXPARAMS < NOISE_MAXPARAMS
							Parameters are randomly generated in such a way that their
							cumulative distribution interpolates linearly the points given by {(xi, Ci)}
							with a standard error on each point given by Cerri
							NOTE:   To avoid risky extrapolations, make sure that C0=0 and CN=1
NoiseTypes::SORTED_LIST		N+1 parameter required: the first is the number of values in the list,
							then there is an ordered list of values. 
*/


enum SimStatus {
	NOT_INITIALIZED, INITIALIZED
};

enum EqDetails {
	AVRATE, DEV, RATECHANGE, STDRATE, FASTSITES
};

enum AvgType {
	ARITMETIC, HARMONIC, GEOMETRIC, MAXIMUM, MINIMUM
};

#define NUM_EQ_DETAILS 5
const std::string strEqLabels[NUM_EQ_DETAILS] = { "avrate", "dev", "ratechange", "stdrate", "fastsites" };





#endif
