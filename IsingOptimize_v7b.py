import sys
import os
import subprocess
import numpy as np
import itertools

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

### SAMPLE SELECTION ###

SAMPLE_SEL = 'LUDOX45'

plot_range = None
defNoiseG = False
defgFromFile = False
defGammaPPD = 20
defPrintRateSTD = True
defPreshear = True
defAscDesc = False
if (SAMPLE_SEL == 'LUDOX45'):
    defK = [30.74]
    defn = [3.702]
    defTau0 = [7891.356]
    defAAvg = [None]
    defGnorm = [0.7]#np.linspace(0.9, 0.99, 10)
    defAVar = [0.0001, 0.22, 0.25, 0.3]#np.geomspace(0.17, 0.23, 4)
    defNoiseG = False
    if defNoiseG:
        defAVar = defAVar * 1000
        defAVar = defAVar * 10000
    defOmega = 3.14
    defgFromFile = False
    defGammaMin = 0.005
    defGammaMax = 0.5
    defGammaPPD = 500
    FastSitesThr = 7
    defPreshear = False
    defAscDesc = True
elif (SAMPLE_SEL == 'LUDOX41_S6'):
    defK = [0.14]
    defn = [1.8]
    defTau0 = [8000]
    defAAvg = [None]
    defGnorm = [0.05]#[0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
    defAVar = [0.01]#[0.01, 0.03, 0.1, 0.3]
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 1.0
    FastSitesThr = 2
    plot_range = [[0.01, 0.2], [1e-4, 1e-2]]
elif (SAMPLE_SEL == 'LUDOX41_S6_Nopresh'):
    defK = [61.6]
    defn = [3.0]
    defTau0 = [10000]
    defAAvg = [None]
    defGnorm = [0.985]
    defAVar = [0.005]
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 0.2
    FastSitesThr = 2
    plot_range = [[0.02, 0.04], [1e-4, 1e-2]]
elif (SAMPLE_SEL == 'PNIPAM_PRESH_T2'):
    defK = [0.02]
    defn = [3.2]
    defTau0 = [70569]
    defAAvg = [None]
    defGnorm = [0.988]#[0.984, 0.988]
    defAVar = [0.15]#[0.15, 0.17]
    defOmega = 3.14
    defgFromFile = False
    defGammaPPD = 100
    defGammaMin = 0.001
    defGammaMax = 1.0
    FastSitesThr = 2
    #plot_range = [[0.1, 0.7], [1e-5, 1e-1]]
elif (SAMPLE_SEL == 'PNIPAM_PRESH_T40'):
    defK = [0.012]
    defn = [2.0]
    defTau0 = [55000]
    defAAvg = [None]
    defGnorm = [0.4]#[0.2, 0.4, 0.6, 0.8]
    defAVar = [0.2]#[0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
    defOmega = 0.157
    defGammaMin = 0.01
    defGammaMax = 1.0
    FastSitesThr = 2
elif (SAMPLE_SEL == 'EM65'):
    defK = [5e11]#[1e18]#[1.5e14]#[1e9]#[2e24]
    defn = [10.0]#[15.0]#[12.0]#[8.0]#[20.0]
    defTau0 = [617.0]#[650.0]
    defAAvg = [None]
    defGnorm = [0, 0.9, 0.95, 0.99]
    defAVar = [0.001, 0.01, 0.1, 0.2]
    defOmega = 6.28
    defgFromFile = False
    defGammaMin = 0.01
    defGammaMax = 0.08
    defGammaPPD = 200
    FastSitesThr = 4
elif (SAMPLE_SEL == 'EM70'):
    defK = [3.1e23]#[2.41e8]
    defn = [20.0]#[8.0]
    defTau0 = [8000]
    defAAvg = [None]
    defGnorm = [0.998]
    defAVar = [0.16]
    defOmega = 6.28
    defgFromFile = False
    defGammaMin = 0.01
    defGammaMax = 0.1
    defGammaPPD = 400
    FastSitesThr = 4
elif (SAMPLE_SEL == 'EM74'):
    defK = [1.8e7]
    defn = [8.0]
    defTau0 = [1630]
    defAAvg = [7000]
    defAVar = [0.006]
    defOmega = 6.28
    defGammaMin = 0.03
    defGammaMax = 0.2
    FastSitesThr = 4
elif (SAMPLE_SEL == 'EM82'):
    defK = [7.08e5]
    defn = [8.0]
    defTau0 = [1000]
    defAAvg = [441]
    defAVar = np.logspace(-4, -3, 11)
    defOmega = 6.28
    defGammaMin = 0.03
    defGammaMax = 0.2
    FastSitesThr = 4
elif (SAMPLE_SEL == 'EM88'):
    defK = [1.15e6]
    defn = [8.0]
    defTau0 = [1500]
    defAAvg = [488.2]
    defAVar = [0.0002]
    defOmega = 6.28
    defGammaMin = 0.03
    defGammaMax = 0.3
    FastSitesThr = 4
else:
    print('ERROR: sample not recognized')
    sys.exit()

# INPUT/OUTPUT SETTINGS
data_root       = r'D:\steaime\Documents\Research\Projects\SoftGlasses\Ising\out'
SaveRoot        = os.path.join(data_root, 'v4', SAMPLE_SEL, 'run22_hysteresis_v3')
ExpDataPath     = os.path.join(data_root, 'v4', SAMPLE_SEL, 'Exp.txt') # None not to compare
ConfigRoot      = r'D:\steaime\Documents\Research\Projects\SoftGlasses\Ising\config'
TemplateFname   = os.path.join(ConfigRoot, 'simParamsTemplate.ini')
SaveFname       = os.path.join(SaveRoot, 'simParams_NNNN.ini')
SimProgPath     = r'C:\Users\steaime\Documents\Visual Studio 2019\SimulIsing_v3\x64\Release\SimulIsing_v3.exe'
sCommFit        = SimProgPath + " " + SaveFname + " -COMPARE " + str(ExpDataPath) + " -SILENT"


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

if __name__ == '__main__':

    paramVarKeys = ['K', 'n', 'Tau0', 'AAvg', 'AVar', 'Gnorm']

    
    ### INITIALIZING PARAMETERS ###
    
    NSites = 64
    NumSim = 1
    
    params = {
                'NSim'              : NumSim,
                'NSites'            : NSites,
                'K_list'            : defK,
                'n_list'            : defn,
                'Tau0_list'         : defTau0,
                'AAvg_list'         : defAAvg,
                'Gnorm_list'        : defGnorm,
                'AVar_list'         : defAVar,
                'K'                 : defK[0],
                'n'                 : defn[0],
                'Tau0'              : defTau0[0],
                'AAvg'              : defAAvg[0],
                'NormGlass'         : None,
                'ANoiseType'        : 2, # 2: Gauss, 3: LogNorm
                'AVar'              : defAVar[0],
                'Omega'             : defOmega,
                'Beta'              : 2.0,
                'RateAvg'           : 1, # 0: Aritmetic, 1: Harmonic
                'UseNewEqn'         : True,
                'NoiseOnG'          : defNoiseG,
                'ForbidUnphys'      : False,
                'ForceMF'           : False,
                'SBNoise'           : False,
                'SBNoiseAvg'        : 2,
                'gFromFile'         : defgFromFile,
                'gFile'             : os.path.join(data_root, 'v4', SAMPLE_SEL, 'GammaList.txt'),
                'gMin'              : defGammaMin,
                'gMax'              : defGammaMax,
                'gPPD'              : defGammaPPD,
                'Preshear'          : defPreshear,
                'AscDesc'           : defAscDesc,
                'outFolder'         : SaveRoot,
                'outPre'            : 'outXXXX_',
                'statPre'           : 'statXXXX_',
                'parName'           : 'paramsXXXX.txt',
                'alphaHistName'     : 'alphaHistXXXX',
                'rateHistName'      : 'rateHistXXXX',
                'chiThr'            : FastSitesThr,
                'wRate'             : 1.0,
                'wChi'              : 100.0,
                'dComb'             : True,
                'saveRateHist'      : True,
                'saveAlphaHist'     : True,
                'printRateSTD'      : defPrintRateSTD,
                'saveRawSnapshots'  : True,
            }
    
    for i in range(201):
        params['AlphaParam{0}'.format(i)] = 0.0

    
    if False:
        npar = CalcNormParams(params)
        params = ModelParam_SetG(npar, 0.8)
        print(params)
        sys.exit()
    
    if False:
        npar = CalcNormParams(params)
        print(npar)
        print('\n\n\n')
        npar = ModelParam_SetG_strainc(params, 0.93, 0.002)
        print(npar)
        #int_test = input('exit?')
        sys.exit()
        
    if False:
        npar = CalcNormParams(params)
        print(npar)
        sys.exit()
        
    varParamList = []
    for key in paramVarKeys:
        varParamList.append(params[key+'_list'])
    prodList = list(itertools.product(*varParamList))
    MFpars = []
    
    print(len(prodList))
    
    CALC = True
    
    n_parall_proc = 4
    if CALC:
        for iCount in range(0, len(prodList), n_parall_proc):
            proc_list = []
            for iProc in range(iCount, min(iCount+n_parall_proc, len(prodList))):
                for iKey in range(len(prodList[iProc])):
                    params[paramVarKeys[iKey]] = prodList[iProc][iKey]
                params['alphaHistName'] = str('alphaHistXXXX').replace('XXXX', str(iProc).zfill(4))
                params['rateHistName'] = str('rateHistXXXX').replace('XXXX', str(iProc).zfill(4))
                cur_pars = {'Omega':params['Omega'], 'n':params['n'], 'K':params['K'], 
                               'Tau0':params['Tau0'], 'Alpha':params['AAvg'], 'Gnorm':params['Gnorm']}
                if params['AAvg'] is None:
                    params['AAvg'] = AlphaFromNormGlass(params['Gnorm'], cur_pars)
                    cur_pars['Alpha'] = params['AAvg']
                else:
                    params['Gnorm'] = CalcNormGlassiness(1.0/params['Tau0'], params['Omega'], params['K'], params['AAvg'], params['Beta'])
                    cur_pars['Gnorm'] = params['Gnorm']
                cur_pars['AVar'] = params['AVar']
                if (params['ANoiseType'] == 4):
                    print('ANoiseType==4 not supported. Use program version 8 or later')
                else:
                    params['AlphaParam0'] = params['AVar']
                MFpars.append(cur_pars)
                print(cur_pars)
                #print(iProc, params['alphaHistName'])
                cur_savename = str(SaveFname).replace('NNNN', str(iProc))
                GenParamFile(params, iProc, save_name=cur_savename)
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
            data_g_desc = []
            data_c_desc = []
            for iCol in range(iFile, min(iFile+col_max_num, len(prodList))):
                if defPreshear:
                    res_suffix = 'ps_all.txt'
                else:
                    res_suffix = 'ad_all.txt'
                cur_col_g = []
                cur_col_c = []
                cur_col_g_desc = []
                cur_col_c_desc = []
                with open(os.path.join(params['outFolder'], str(params['outPre']).replace('XXXX', str(iCol).zfill(4))+res_suffix), 'r') as fin:
                    nhead = 0
                    for line in fin:
                        if (nhead >= 1):                      
                            words = line.split()
                            if defPrintRateSTD:
                                colc_idx = 3
                            else:
                                colc_idx = 2
                            cur_col_c.append(float(words[colc_idx]))
                            if (NumSim == 1):
                                cur_col_g.append(float(words[1]))
                            if (iCol == 0):
                                data_first.append(float(words[0]))
                            if defAscDesc:
                                if defPrintRateSTD:
                                    cur_col_c_desc.append(float(words[colc_idx+3]))
                                else:
                                    cur_col_c_desc.append(float(words[colc_idx+2]))
                                if (NumSim == 1):
                                    cur_col_g_desc.append(float(words[colc_idx+1]))
                        nhead += 1
                if (NumSim > 1):
                    with open(os.path.join(params['outFolder'], str(params['statPre']).replace('XXXX', str(iCol).zfill(4))+res_suffix), 'r') as fin:
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
                    data_g_desc.append(cur_col_g_desc)
                    data_c_desc.append(cur_col_c_desc)
                    fpar.write('\n' + str(iCol))
                    for iKey in range(len(prodList[iCol])):
                        fpar.write('\t' + str(prodList[iCol][iKey]))
                
                fig, ax = plt.subplots(figsize=(12,8))
                if ExpDataPath is not None:
                    exp_gamma, exp_chi, exp_chierr, exp_rate, exp_rateerr = np.loadtxt(ExpDataPath, skiprows=1, unpack=True)
                    cur_MF_idx = iFile*col_max_num + iCol
                    cur_MF = MFpars[cur_MF_idx]
                    ax.errorbar(exp_gamma, exp_rate, yerr=exp_rateerr, fmt='ko')
                    rate_MFlist = np.geomspace(np.min(exp_rate), np.max(exp_rate), 300)
                    ax.plot(MFstrain(rate_MFlist, cur_MF), rate_MFlist, 'k--')
                    ax.axhline(y=1.0/cur_MF['Tau0'], c='k', linestyle='-.')
                    ax.plot(data_first, cur_col_g, 'kx-')
                    if defAscDesc:
                        ax.plot(data_first, cur_col_g_desc, 'k+-')
                    ax.set_xscale('log')
                    ax.set_yscale('log')
                    ax.set_xlabel(r'$\gamma_0$')
                    ax.set_ylabel(r'$\Gamma [s^{-1}]$', color='k')
                    ax.tick_params(axis='y', labelcolor='k')
                
                    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
                    ax2.set_ylabel(r'$\chi$', color='g')  # we already handled the x-label with ax1
                    ax2.errorbar(exp_gamma, exp_chi, yerr=exp_chierr, fmt='gs')
                    ax2.plot(data_first, np.subtract(1, cur_col_c), 'g+-')
                    ax2.tick_params(axis='y', labelcolor='g')
                    str_title = r'[#{6}]:: $K={0:.1f} - n={1:.1f} - \tau_0={2:.2e} - \bar\alpha={3:.4e} - \sigma_\alpha^2={4:.3f} - g_n={5:.3f}$'.format(cur_MF['K'],
                                            cur_MF['n'], cur_MF['Tau0'], cur_MF['Alpha'], cur_MF['AVar'], cur_MF['Gnorm'], cur_MF_idx)
                
                    if plot_range is not None:
                        ax.set_xlim(plot_range[0])
                        ax.set_ylim(plot_range[1])
                    ax2.set_ylim([-0.1, 1.1])
                
                    if len(exp_gamma) == len(data_first):
                        if (np.max(np.abs(exp_gamma - data_first)) < 1e-3):
                            chi_chisq = np.nansum(np.true_divide(np.square(np.subtract(exp_chi, cur_col_c)), np.square(exp_chierr)))
                            rate_chisq = np.nansum(np.true_divide(np.square(np.subtract(exp_rate, cur_col_g)), np.square(exp_rateerr)))
                            chi_loglike = np.sum(np.log(np.true_divide(1.0, np.sqrt(np.multiply(2.0*np.pi,np.square(exp_chierr)))))) - 0.5*chi_chisq
                            rate_loglike = np.sum(np.log(np.true_divide(1.0, np.sqrt(np.multiply(2.0*np.pi,np.square(exp_rateerr)))))) - 0.5*rate_chisq
                            fpar.write('\t' + str(rate_loglike) + '\t' + str(chi_loglike) + '\t' + str(rate_loglike + chi_loglike))
                            str_title += r' - llike={0:.3e}'.format(rate_loglike + chi_loglike)
                        else:
                            fpar.write('\t-\t-\t-')
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
                        