import sys
import os
import subprocess
import numpy as np
import bisect
from collections.abc import Iterable

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

### SAMPLE SELECTION ###

SAMPLE_SEL = 'LUDOX45'

if (SAMPLE_SEL == 'LUDOX45'):
    logparams = {
        'K':    {'val':np.log(30.74), 'bounds':None, 'sigma':0.06, 'label':r'$\ln(R)$'},
        'n':    {'val':1.31,          'bounds':[1,2],'sigma':0.05, 'label':r'$\ln(n)$'},
        'Tau0': {'val':8.97,          'bounds':None, 'sigma':0.03, 'label':r'$\ln(\tau_0)$'},
        'Gnorm':{'val':-2.58,         'bounds':None, 'sigma':0.35, 'label':r'$\ln(1-g/g_{max})$'},
        }
    exp_omega = 3.14
    FastSitesThr = 7
elif (SAMPLE_SEL == 'LUDOX41_S6'):
    logparams = {}
    exp_omega = 3.14
    FastSitesThr = 2
elif (SAMPLE_SEL == 'LUDOX41_S6_Nopresh'):
    exp_omega = 3.14
    FastSitesThr = 2
elif (SAMPLE_SEL == 'PNIPAM_PRESH_T2'):
    exp_omega = 3.14
    FastSitesThr = 2
elif (SAMPLE_SEL == 'PNIPAM_PRESH_T40'):
    exp_omega = 0.157
    FastSitesThr = 2
elif (SAMPLE_SEL == 'EM65'):
    exp_omega = 6.28
    FastSitesThr = 4
elif (SAMPLE_SEL == 'EM70'):
    exp_omega = 6.28
    FastSitesThr = 4
elif (SAMPLE_SEL == 'EM74'):
    exp_omega = 6.28
    FastSitesThr = 4
elif (SAMPLE_SEL == 'EM82'):
    exp_omega = 6.28
    FastSitesThr = 4
elif (SAMPLE_SEL == 'EM88'):
    exp_omega = 6.28
    FastSitesThr = 4
else:
    print('ERROR: sample not recognized')
    sys.exit()

# INPUT/OUTPUT SETTINGS
data_root       = r'D:\steaime\Documents\Research\Projects\SoftGlasses\Ising\out'
SaveRoot        = os.path.join(data_root, 'v4', SAMPLE_SEL, 'run13_Cust01_NewEq_HarmAvg')
AlphaDistrFile  = os.path.join(SaveRoot, 'alpha_cumul_distr.txt')
ExpDataPath     = os.path.join(data_root, 'v4', SAMPLE_SEL, 'Exp.txt') # None not to compare
ConfigRoot      = r'D:\steaime\Documents\Research\Projects\SoftGlasses\Ising\config'
TemplateFname   = os.path.join(ConfigRoot, 'simParamsTemplate.ini')
SaveFname       = os.path.join(SaveRoot, 'simParams_NNNN.ini')
SimProgPath     = r'C:\Users\steaime\Documents\Visual Studio 2019\SimulIsing_v3\x64\Release\SimulIsing_v3.exe'
sCommFit        = SimProgPath + " " + SaveFname + " -COMPARE " + ExpDataPath + " -SILENT"


def GenParamFile(dict_params, idx=None, key_firstlast='%', key_idx='XXXX', save_name=SaveFname):
    # READ TEMPLATE DATA
    template_data = []
    with open(TemplateFname) as f:
        for line in f:
            template_data.append(line)
    
    # WRITE PARAM FILE
    with open(save_name, 'w') as f:
        for line in template_data:
            line_str = str(line).replace('\n', '')#.replace(' ', '')
            for key in dict_params.keys():
                line_str = line_str.replace(key_firstlast+key+key_firstlast, str(dict_params[key]))
                if (idx != None):
                    line_str = line_str.replace(key_idx, str(idx).zfill(len(key_idx)))
            f.write(line_str)
            f.write('\n')

def SweepAndCompare(params, fOut=None, PrintKeys=[]):
    global iGlobalCount
    iGlobalCount += 1
    GenParamFile(params)
    proc = subprocess.Popen(sCommFit, stdout=subprocess.PIPE)
    res = float(str(proc.stdout.read(), 'utf-8'))
    if (fOut!=None):
        fOut.write('\n' + str(iGlobalCount))
        for key in PrintKeys:
            fOut.write('\t' + str(params[key]))
        fOut.write('\t' + str(res))
        fOut.flush()
    return res

def CalcNormGlassiness(Gamma0, Omega, K, Alpha, Beta=2.0):
    old_psi_c = 4.0*Beta / (Beta**2 - 1)
    old_psi_min = ((Beta+1.0) / Beta)**Beta / (Beta - 1)
    old_psi = ((Beta+1.0)/(Beta-1.0))**Beta * (Gamma0/Omega)**(Beta-1) * (K / Alpha)
    return (old_psi_c - old_psi) / (old_psi_c - old_psi_min)
    
def CalcNormParams(modelParams):
    res = modelParams.copy()
    if res['UseNewEqn']:
        res['Beta'] = 2.0
    res['Gamma0'] = 1.0 / res['Tau0']
    res['xi'] = (res['Beta'] - 1.0) / (res['Beta'] + 1.0)
    res['Gammac'] = res['Gamma0'] * (1.0 + 2.0 / (res['Beta'] - 1.0))
    res['strainc'] = np.power((res['xi'] / res['AAvg']) * np.power(res['xi']*res['Gamma0'] /\
                   res['Omega'] , res['Beta']), 1.0 / res['n'])
    res['psi'] = (res['K'] / (res['AAvg'] * np.power(res['xi'], res['Beta']) * np.power(res['Tau0'] * res['Omega'], res['Beta'] - 1.0))) /\
                 ((4.0 * res['Beta']) / (np.power(res['Beta'], 2.0) - 1.0))
    res['psic'] = 1.0
    res['g'] = 1.0 - res['psi']
    return res

def ModelParam_SetG_strainc(normParams, g_val, strainc_val):
    res = CalcNormParams(normParams)
    res['Tau0'] = 1.0 / res['Gamma0']
    res['strainc'] = strainc_val
    res['g'] = g_val
    res['AAvg'] = (res['xi'] / np.power(res['strainc'], res['n'])) *\
                    np.power(res['xi'] * res['Gamma0'] / res['Omega'], res['Beta'])
    res['K'] = (1.0 - res['g']) * (np.power(res['Tau0'] * res['Omega'], res['Beta'] - 1.0) *\
                                   res['AAvg'] * np.power(res['xi'], res['Beta']))
    return res

def ModelParam_SetG(params, g_val):
    res = params.copy()
    res['psi'] = 1.0 - g_val
    res['AAvg'] = (res['K'] / (res['psi'] * ((4.0 * res['Beta']) / (np.power(res['Beta'], 2.0) - 1.0)) * np.power(res['xi'], res['Beta']))) *\
                    np.power(res['Gamma0'] / res['Omega'], res['Beta'] - 1.0)
    return res

def AlphaFromNormGlass(gnorm, _params):
    psi = 1 - gnorm * (1.0 -  (9.0 / 4.0) / (8.0 / 3.0))
    return _params['K'] * 9.0 / (psi * (8.0 / 3.0) * _params['Omega'] * _params['Tau0'])

def MFstrain(rate, params):
    normrate = np.divide(rate, params['Omega'])
    return np.power(np.subtract(np.divide(params['K'], np.subtract(normrate, 1.0/(params['Omega']*params['Tau0']))),
                                params['Alpha'] / np.square(normrate)), -1.0/params['n'])

def CubicRoots(coeffs):
    # coeffs: largest-to-lowest power
    # reference: https://en.wikipedia.org/wiki/Cubic_equation#General_cubic_formula
    _D0 = coeffs[1]**2 - 3*coeffs[0]*coeffs[2]
    _D1 = 2*coeffs[1]**3 - 9*coeffs[0]*coeffs[1]*coeffs[2] + 27*coeffs[0]**2*coeffs[3]
    _cbr1 = [1, complex(-.5,.75**.5), complex(-.5,-.75**.5)]
    _Cmod = complex(0.5*(_D1 + complex(_D1**2 - 4*_D0**3)**(1/2)))**(1/3)
    return [-(1.0/(3*coeffs[0]))*(coeffs[1]+_Cmod*xi+_D0/(_Cmod*xi)) for xi in _cbr1]

def FilterReal(complex_list, imag_thr=1e-10):
    return [np.real(x) for x in complex_list if np.abs(np.imag(x)) < imag_thr]

def MF_jumpStrain(params):
    _n, _t0, _w, _K, _a, _gn = params['n'], params['Tau0'], params['Omega'], params['K'], params['Alpha'], params['Gnorm']
    coeffs = [_K*_w, -2*_a*_w**2, 4*_a*_w**2/_t0, -2*_a*_w**2/_t0**2]
    ok_sol = [x for x in FilterReal(CubicRoots(coeffs)) if x >= 1.0/_t0]
    if len(ok_sol)>0:
        return  MFstrain(np.min(ok_sol), params)
    else:
        return np.nan

def JumpOrAccelStrain(params, threshold):
    try_js =  MF_jumpStrain(params)
    if np.isnan(try_js):
        return MFstrain(threshold * 1.0 / params['Tau0'], params)
    else:
        return try_js
    
def MFrate(strain, params, min_root=True):
    '''  Calculates mean field rate
    
    Parameters
    ----------
    strain: float or array of floats, strain values (in strain units)
    params: dict with model parameters
    min_root: array of bool (same size as strain) or None. 
                If None or true, the smallest physically acceptable solution is returned
                if False, the largert physically acceptable solution is returned
                physically acceptable solutions must be real and larger than the spontaneous rate
    '''
    _n, _t0, _w, _K, _a, _gn = params['n'], params['Tau0'], params['Omega'], params['K'], params['Alpha'], params['Gnorm']
    if not isinstance(strain, Iterable):
        strain = [strain]
    strain = np.asarray(strain)
    if not isinstance(min_root, Iterable):
        min_root = [min_root] * len(strain)
    res = np.empty(len(strain), dtype=float)
    for i in range(len(strain)):
        _g0 = strain[i]
        coeffs = [_g0**(-_n), -_K*_w-_g0**(-_n)/_t0, _a*_w**2, -_a*_w**2/_t0]
        ok_sol = [x for x in FilterReal(CubicRoots(coeffs)) if x >= 1.0/_t0]
        if len(ok_sol)>0:
            if min_root[i]:
                res[i] = np.min(ok_sol)
            else:
                res[i] = np.max(ok_sol)
        else:
            res[i] = np.nan 
    return res

def MinJumpStrain(params, rate_threshold):
    bkp_alpha = params['Alpha']
    params['Alpha'] = 0
    min_strain = MFstrain(rate_threshold * 1.0 / params['Tau0'], params)
    params['Alpha'] = bkp_alpha
    return min_strain    

def IsJumpStrainPhysical(params, strain, rate_threshold):
    return (MinJumpStrain(params, rate_threshold) < strain)

def GetAlphaFromYieldStrain(params, yieldstrain, rate_threshold, rel_precision=1e-4, step_size=1e-2, debug=False, verbose=0, max_steps=1000, auto_restart=True):
    if not IsJumpStrainPhysical(params, yieldstrain, rate_threshold):
        if debug:
            if verbose:
                print('Min strain: {0:.4f}. Required {1:.4f}. Returning 0'.format(MinJumpStrain(params, rate_threshold), yieldstrain))
            return 0.0, []
        else:
            return 0.0
    # Back-up original alpha:
    alpha_bkp = params['Alpha']
    # Let's check that alpha does not exceet alpha_max. If it's the case, let's restart from an arbitrary glassiness
    alpha_max = AlphaFromNormGlass(1, params)
    if debug and verbose:
        print(alpha_max)
    if params['Alpha'] >= alpha_max:
        params['Alpha'] = AlphaFromNormGlass(0.5, params)
    count_steps = 0
    cur_ys = JumpOrAccelStrain(params, rate_threshold)
    rel_diff = (cur_ys - yieldstrain) / yieldstrain
    if debug:
        trace = [[params['Alpha'], cur_ys, rel_diff, (1 - rel_diff*step_size)]]
    while(np.abs(rel_diff) > rel_precision):
        params['Alpha'] *= (1 - rel_diff*step_size)
        params['Alpha'] = max(1e-15, min(params['Alpha'], alpha_max*0.9999999999))
        cur_ys = JumpOrAccelStrain(params, rate_threshold)
        rel_diff = (cur_ys - yieldstrain) / yieldstrain
        count_steps += 1
        if debug:
            if verbose:
                print(trace[-1])
            trace.append([params['Alpha'], cur_ys, rel_diff, (1 - rel_diff*step_size)])
        if count_steps > max_steps:
            if auto_restart:
                params['Alpha'] = alpha_bkp
                if debug:
                    print('retrying with step size reduced from {0} to {1}'.format(step_size, step_size*0.1))
                    params['Alpha'], trace = GetAlphaFromYieldStrain(params, yieldstrain, rate_threshold, rel_precision, step_size=step_size*0.1, 
                                                                     debug=debug, verbose=verbose, max_steps=max_steps)
                else:
                    params['Alpha'] = GetAlphaFromYieldStrain(params, yieldstrain, rate_threshold, rel_precision, step_size=step_size*0.1, 
                                                                     debug=debug, verbose=verbose, max_steps=max_steps)
            else:
                params['Alpha'] = np.nan
            break
    res_alpha = params['Alpha']
    params['Alpha'] = alpha_bkp
    if debug:
        return res_alpha, np.asarray(trace)
    else:
        return res_alpha
    



######################################################

if __name__ == '__main__':

    paramVarKeys = list(logparams.keys())

    ### INITIALIZING PARAMETERS ###
    
    NSites = 32
    MaxAlphaParam = 201
    NumSim = 4
    
    NumRandomDraws = 20
    
    params = {
                'NSim'              : NumSim,
                'NSites'            : NSites,
                'Omega'             : exp_omega,
                'ANoiseType'        : 4, # 2: Gauss, 3: LogNorm, 4: Custom
                'RateAvg'           : 0, # 0: Aritmetic, 1: Harmonic
                'UseNewEqn'         : True,
                'NoiseOnG'          : False,
                'SBNoise'           : True,
                'SBNoiseAvg'        : 2,
                'ForceMF'           : False,
                'gFromFile'         : True,
                'gFile'             : os.path.join(data_root, 'v4', SAMPLE_SEL, 'GammaList.txt'),
                'outFolder'         : SaveRoot,
                'outPre'            : 'outXXXX_',
                'statPre'           : 'statXXXX_',
                'parName'           : 'paramsXXXX.txt',
                'alphaHistName'     : 'alphaHistXXXX',
                'chiThr'            : FastSitesThr,
                'saveRateHist'      : True,
                'saveAlphaHist'     : True,
            }
    
    for i in range(1, MaxAlphaParam):
        params['AlphaParam{0}'.format(i)] = 0.0
    if params['ANoiseType'] == 4:
        cust_alpha, cust_CD = np.loadtxt(AlphaDistrFile, unpack=True)
        cust_median = cust_alpha[cust_CD >= 0.5][0]
        params['AlphaParam0'] = cust_median
        params['AlphaParam1'] = len(cust_alpha)
        for i in range(len(cust_alpha)):
            params['AlphaParam{0}'.format(3*i+2)] = cust_alpha[i]
            params['AlphaParam{0}'.format(3*i+3)] = cust_CD[i]
            params['AlphaParam{0}'.format(3*i+4)] = 0.0
            
    MFpars = []
    Sample_pars = []
    prodList = np.empty((NumRandomDraws, len(paramVarKeys)), dtype=float)
    for i in range(len(paramVarKeys)):
        prodList[:,i] = np.random.normal(logparams[paramVarKeys[i]]['val'], logparams[paramVarKeys[i]]['sigma'], NumRandomDraws)
    prodList = np.exp(prodList)
    
    if ExpDataPath is not None:
        exp_gamma, exp_chi, exp_chierr, exp_rate, exp_rateerr = np.loadtxt(ExpDataPath, skiprows=1, unpack=True)

    CALC = True
    
    n_parall_proc = 4
    if CALC:
        for iCount in range(0, len(prodList), n_parall_proc):
            proc_list = []
            for iProc in range(iCount, min(iCount+n_parall_proc, len(prodList))):
                for iKey in range(len(prodList[iProc])):
                    params[paramVarKeys[iKey]] = prodList[iProc][iKey]
                params['alphaHistName'] = str('alphaHistXXXX').replace('XXXX', str(iProc).zfill(4))
                cur_pars = {'Omega':params['Omega'], 'n':params['n'], 'K':params['K'], 
                               'Tau0':params['Tau0'], 'Gnorm':params['Gnorm']}
                params['AAvg'] = AlphaFromNormGlass(params['Gnorm'], cur_pars)
                cur_pars['Alpha'] = params['AAvg']
                MFpars.append(cur_pars)
                print(cur_pars)
                #print(iProc, params['alphaHistName'])
        
                if False:
                    test_alpha, test_trace = GetAlphaFromYieldStrain(cur_pars, 0.285, FastSitesThr, step_size=1.0, debug=True, verbose=0, auto_restart=True)
                    print(test_alpha)
                    cur_pars['Alpha'] = test_alpha
                    print(MF_jumpStrain(cur_pars))
                    print(JumpOrAccelStrain(cur_pars, FastSitesThr))
                    print('Conveged in ' + str(len(test_trace)) + ' steps')
                    sys.exit()
                cur_pars['Alpha'] = params['AAvg']
                
                if False:
                    # Generate cumulative distribution of alpha_ij from chi(gamma)
                    cumdistr_x, cumdistr_y, cumdistr_yerr = [], [], []
                    bln_useless = True
                    for i in range(len(exp_gamma)):
                        # leading 0s and trailing 1s are just useless...
                        if bln_useless:
                            if (exp_chi[i] > 0.0001 and exp_chi[i] < 0.9999):
                                bln_useless = False
                            if (i > 0):
                                if (exp_chi[i-1] > 0.0001 and exp_chi[i-1] < 0.9999):
                                    bln_useless = False
                            if (i < len(exp_gamma)-1):
                                if (exp_chi[i+1] > 0.0001 and exp_chi[i+1] < 0.9999):
                                    bln_useless = False
                        #print([i, exp_chi[i], bln_useless])
                        if not bln_useless:
                            cumdistr_x.append(GetAlphaFromYieldStrain(cur_pars, exp_gamma[i], FastSitesThr, step_size=1.0, debug=False, verbose=0, auto_restart=True))
                            if i==0:
                                cumdistr_y.append(0)
                            elif i== len(exp_gamma)-1:
                                cumdistr_y.append(1)
                            elif len(cumdistr_y)==0:
                                cumdistr_y.append(1 - exp_chi[i])
                            else:
                                cumdistr_y.append(max(cumdistr_y[-1], 1 - exp_chi[i]))
                            cumdistr_yerr.append(exp_chierr[i])
                            if exp_chi[i]<=0:
                                bln_useless = True
                                break
                    #print(bisect.bisect(cumdistr_y, 0.5))
                    median_idx = bisect.bisect(cumdistr_y, 0.5)
                    if median_idx<=0:
                        median_idx = 1
                    elif median_idx==len(cumdistr_y):
                        median_idx = len(cumdistr_y) - 1
                    #print(median_idx)
                    if cumdistr_y[median_idx-1]==cumdistr_y[median_idx]:
                        if median_idx < len(cumdistr_y)-1:
                            median_idx+=1
                        else:
                            median_idx-=1
                    #print(median_idx)
                    median_val = cumdistr_x[median_idx-1] + (cumdistr_x[median_idx] - cumdistr_x[median_idx-1]) * (0.5 - cumdistr_y[median_idx-1]) / (cumdistr_y[median_idx] - cumdistr_y[median_idx-1])
                    params['AlphaParam0'] = median_val
                    params['AlphaParam1'] = len(cumdistr_x)
                    for i in range(len(cumdistr_x)):
                        params['AlphaParam{0}'.format(3*i+2)] = cumdistr_x[i]
                        params['AlphaParam{0}'.format(3*i+3)] = cumdistr_y[i]
                        params['AlphaParam{0}'.format(3*i+4)] = cumdistr_yerr[i]
                        
                    if True:
                        test_ys, test_jumporacc = [], []
                        test_alpha, test_chi = [], []
                        for tmp in range(len(exp_gamma)):
                            if tmp < len(cumdistr_x):
                                cur_pars['Alpha'] = cumdistr_x[tmp]
                                test_jumporacc.append(JumpOrAccelStrain(cur_pars, FastSitesThr))
                                test_ys.append(MF_jumpStrain(cur_pars))
                            if (exp_chi[tmp] < 1.0):
                                cur_pars['Alpha'] = params['AAvg']
                                test_alpha.append(GetAlphaFromYieldStrain(cur_pars, exp_gamma[tmp], FastSitesThr, step_size=1.0, debug=False, verbose=0, auto_restart=True))
                                test_chi.append(1-exp_chi[tmp])
                        plt.plot(test_jumporacc, cumdistr_y, 'kx:')
                        plt.plot(test_ys, cumdistr_y, 'y--')
                        plt.plot(exp_gamma, 1-exp_chi, 'r.')
                        plt.figure()
                        plt.plot(test_alpha, test_chi, 'g^')
                        plt.plot(cumdistr_x, cumdistr_y, 'b--')
                        plt.show()
                        sys.exit()
                    
                cur_savename = str(SaveFname).replace('NNNN', str(iProc))
                GenParamFile(params, iProc, save_name=cur_savename)
                if False:
                    sys.exit()
                proc = subprocess.Popen(SimProgPath + " " + cur_savename, stdout=subprocess.PIPE)
                proc_list.append(proc)
            for cur_p in proc_list:
                cur_p.communicate()
    col_max_num = 20000
    fpar = open(os.path.join(params['outFolder'], '__par.txt'), 'w')
    fpar.write('#')
    for key in paramVarKeys:
        fpar.write('\t' + str(key))
    if ExpDataPath is not None:
        fpar.write('\tloglike_rate\tloglike_chi\tloglike_tot')
    data_first = []
    with PdfPages(os.path.join(params['outFolder'], 'figres.pdf')) as fig_pdf:
        for iFile in range(0, len(prodList), col_max_num):
            fgamma = open(os.path.join(params['outFolder'], '_all_gamma_' + str(iFile//col_max_num).zfill(4) + '.txt'), 'w')
            fchi = open(os.path.join(params['outFolder'], '_all_chi_' + str(iFile//col_max_num).zfill(4) + '.txt'), 'w')
            fgamma.write('g')
            fchi.write('g')
            data_g = []
            data_c = []
            for iCol in range(iFile, min(iFile+col_max_num, len(prodList))):
                cur_col_g = []
                cur_col_c = []
                with open(os.path.join(params['outFolder'], str(params['outPre']).replace('XXXX', str(iCol).zfill(4))+'ps_all.txt'), 'r') as fin:
                    nhead = 0
                    for line in fin:
                        if (nhead >= 1):                      
                            words = line.split()
                            cur_col_c.append(float(words[2]))
                            if (NumSim == 1):
                                cur_col_g.append(float(words[1]))
                            if (iCol == 0):
                                data_first.append(float(words[0]))
                        nhead += 1
                if (NumSim > 1):
                    with open(os.path.join(params['outFolder'], str(params['statPre']).replace('XXXX', str(iCol).zfill(4))+'ps_all.txt'), 'r') as fin:
                        nhead = 0
                        for line in fin:
                            if (nhead >= 1):                      
                                words = line.split()
                                cur_col_g.append(float(words[1]))
                            nhead += 1
                if (True or ((cur_col_c[0] == 0) and (cur_col_c[-1] > 0) and (cur_col_g[0] < 0.003))):# and (cur_col_g[-1] > 0.3) and (cur_col_g[-1] < 30))):
                    fgamma.write('\t' + str(iCol))
                    fchi.write('\t' + str(iCol))
                    data_g.append(cur_col_g)
                    data_c.append(cur_col_c)
                    fpar.write('\n' + str(iCol))
                    for iKey in range(len(prodList[iCol])):
                        fpar.write('\t' + str(prodList[iCol][iKey]))
                if ExpDataPath is not None:
                    cur_MF_idx = iFile*col_max_num + iCol
                    cur_MF = MFpars[cur_MF_idx]
                    fig, ax = plt.subplots(figsize=(12,8))
                    ax.errorbar(exp_gamma, exp_rate, yerr=exp_rateerr, fmt='ko')
                    rate_MFlist = np.geomspace(np.min(exp_rate), np.max(exp_rate), 300)
                    ax.plot(MFstrain(rate_MFlist, cur_MF), rate_MFlist, 'k--')
                    ax.axhline(y=1.0/cur_MF['Tau0'], c='k', linestyle='-.')
                    ax.plot(data_first, cur_col_g, 'kx-')
                    ax.set_xscale('log')
                    ax.set_yscale('log')
                    ax.set_xlabel(r'$\gamma_0$')
                    ax.set_ylabel(r'$\Gamma [s^{-1}]$', color='k')
                    ax.tick_params(axis='y', labelcolor='k')
                
                    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
                    ax2.set_ylabel(r'$\chi$', color='g')  # we already handled the x-label with ax1
                    ax2.errorbar(exp_gamma, exp_chi, yerr=exp_chierr, fmt='gs')
                    ax2.plot(data_first, cur_col_c, 'g+-')
                    ax2.tick_params(axis='y', labelcolor='g')
                    str_title = r'[#{6}]:: $K={0:.1f} - n={1:.1f} - \tau_0={2:.2e} - \bar\alpha={3:.4e} - \sigma_\alpha^2={4:.3f} - g_n={5:.3f}$'.format(cur_MF['K'],
                                            cur_MF['n'], cur_MF['Tau0'], cur_MF['Alpha'], cur_MF['AVar'], cur_MF['GNorm'], cur_MF_idx)
                
                    if (np.max(np.abs(exp_gamma - data_first)) < 1e-3):
                        chi_chisq = np.nansum(np.true_divide(np.square(np.subtract(exp_chi, cur_col_c)), np.square(exp_chierr)))
                        rate_chisq = np.nansum(np.true_divide(np.square(np.subtract(exp_rate, cur_col_g)), np.square(exp_rateerr)))
                        chi_loglike = np.sum(np.log(np.true_divide(1.0, np.sqrt(np.multiply(2.0*np.pi,np.square(exp_chierr)))))) - 0.5*chi_chisq
                        rate_loglike = np.sum(np.log(np.true_divide(1.0, np.sqrt(np.multiply(2.0*np.pi,np.square(exp_rateerr)))))) - 0.5*rate_chisq
                        fpar.write('\t' + str(rate_loglike) + '\t' + str(chi_loglike) + '\t' + str(rate_loglike + chi_loglike))
                        str_title += r' - llike={0:.3e}'.format(rate_loglike + chi_loglike)
                    else:
                        fpar.write('\t-\t-\t-')

                    fig.suptitle(str_title)
                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    fig_pdf.savefig(fig)
                    plt.close('all')
                    
            if (len(data_g) > 0):
                for iRow in range(len(data_g[0])):
                    fgamma.write('\n' + str(data_first[iRow]))
                    fchi.write('\n' + str(data_first[iRow]))
                    for iCol in range(len(data_g)):
                        fgamma.write('\t' + str(data_g[iCol][iRow]))
                        if (data_c[iCol][-1] != 0):
                            div = data_c[iCol][-1]
                        else:
                            div = 1.0
                        fchi.write('\t' + str(data_c[iCol][iRow] / div))
            fgamma.close()
            fchi.close()
        fpar.close()
        print('... all done!')
                        