import sys
import subprocess
import numpy as np
import itertools

### SAMPLE SELECTION ###

SAMPLE_SEL = 'LUDOX45_n3'

if (SAMPLE_SEL == 'LUDOX45_n4'):
    defK = [110.0]
    defn = [4.0]
    defTau0 = [2400.0]
    defAAvg = [0.05]
    defAVar = [0.5]
    defOmega = 3.14
    defGammaMin = 0.005
    defGammaMax = 0.5
    FastSitesThr = 5
    defAVar_bounds = [0.0, 1.0]
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run1_Ludox45\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'LUDOX45_n3'):
    defK = [8.0]
    defn = [3.0]
    defTau0 = [9000.0]
    defAAvg = [1.112E-3]
    defAVar = [0.08]
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 0.3
    FastSitesThr = 7
    defAVar_bounds = [0.0, 1.0]
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run1_Ludox45\\run10\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'LUDOX41_S6'):
    defK = [1.0]#[0.6]#[1.9]
    defn = [3.0]
    defTau0 = [2400]
    defAAvg = [4.21e-4]
    defAVar = [0.8]
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 1.0
    FastSitesThr = 2
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run8_Ludox41\\out5\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'LUDOX41_S6_bis'):
    defK = np.linspace(3.7, 4.4, 8)
    defn = [3.0]
    defTau0 = [4800]
    defAAvg = [1e-3]
    defAVar = [0.8]
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 1.0
    FastSitesThr = 2
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run8_Ludox41\\out7\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'LUDOX41_S6_Nopresh'):
    defK = [61.6]
    defn = [3.0]
    defTau0 = [10000]
    defAAvg = [7.79e-3]
    defAVar = [0.05]
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 0.2
    FastSitesThr = 2
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run10_Ludox41Presh\\out4\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'LUDOX41_S6_n4'):
    defK = [9.1]#[1.8]
    defn = [4.0]
    defTau0 = [2400]
    defAAvg = [0.004]#[8e-4]
    defAVar = np.logspace(-2, 0, 20)
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 1.0
    FastSitesThr = 2
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run8_Ludox41\\out2\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'PNIPAM_PRESH_T2'):
    defK = [0.00125]
    defn = [2.0]
    defTau0 = [90000]
    defAAvg = [1.67e-8]
    defAVar = [0.11]
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 1.0
    FastSitesThr = 2
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run2_PnipamPresh\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'PNIPAM_PRESH_T40'):
    defK = [0.012]
    defn = [2.0]
    defTau0 = [54000]
    defAAvg = np.linspace(3e-6, 7e-6, 16)
    defAVar = [0.5]
    defOmega = 0.157
    defGammaMin = 0.01
    defGammaMax = 1.0
    FastSitesThr = 2
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run9_PnipamPreshT40\\out6\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'EM65_n20'):
    defK = [3e23]
    defn = [20.0]
    defTau0 = [617.0]
    defAAvg = [2e20]
    defAVar = [0.5]
    defOmega = 6.28
    defGammaMin = 0.001
    defGammaMax = 0.08
    FastSitesThr = 4
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run3_Em65\\'
    ExpDataPath = SaveRoot + "Exp_Em65.txt"
elif (SAMPLE_SEL == 'EM65'):
    defK = [1e9]
    defn = [8.0]
    defTau0 = [617.0]
    defAAvg = [8.8e5]
    defAVar = [0.35]
    defOmega = 6.28
    defGammaMin = 0.001
    defGammaMax = 0.08
    FastSitesThr = 4
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run3_Em65\\'
    ExpDataPath = SaveRoot + "Exp_Em65.txt"
elif (SAMPLE_SEL == 'EM70'):
    defK = np.logspace(9, 10, 5)
    defn = [8.0]
    defTau0 = np.logspace(np.log10(7000), np.log10(20000), 5)
    defAAvg = np.logspace(5, 6.5, 21)
    defAVar = np.logspace(-1.2, 0.5, 11)
    defOmega = 6.28
    defGammaMin = 0.01
    defGammaMax = 0.1
    FastSitesThr = 4
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run4_Em70\\out3\\'
    ExpDataPath = SaveRoot + "Exp_Em70.txt"
elif (SAMPLE_SEL == 'EM70_2'):
    defK = [1e9]
    defn = [8.0]
    defTau0 = [12000]
    defAAvg = [53830]
    defAVar = [0.009]
    defOmega = 6.28
    defGammaMin = 0.005
    defGammaMax = 0.2
    FastSitesThr = 4
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run4_Em70\\out5\\'
    ExpDataPath = SaveRoot + "Exp_Em70.txt"
elif (SAMPLE_SEL == 'EM70_3'):
    defK = [2.41e8]
    defn = [8.0]
    defTau0 = [7500]
    defAAvg = [20550]
    defAVar = [0.002]
    defOmega = 6.28
    defGammaMin = 0.01
    defGammaMax = 0.2
    FastSitesThr = 4
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run4_Em70\\out8\\'
    ExpDataPath = SaveRoot + "Exp_Em70.txt"
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
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run5_Em74\\out3\\'
    ExpDataPath = SaveRoot + "Exp_Em74.txt"
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
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run6_Em82\\out2\\'
    ExpDataPath = SaveRoot + "Exp_Em82.txt"
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
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run7_Em88\\out2\\'
    ExpDataPath = SaveRoot + "Exp_Em88.txt"
else:
    print('ERROR: sample not recognized')
    sys.exit()

# INPUT/OUTPUT SETTINGS
ConfigRoot      = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\config\\'
TemplateFname   = ConfigRoot + 'simParamsTemplate.ini'
SaveFname       = SaveRoot + 'simParams_NNNN.ini'
SimProgPath     = "C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\SimulIsing_v3\\x64\\Release\\SimulIsing_v3_190524.exe"
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

def CalcNormParams(modelParams):
    res = modelParams.copy()
    res['Gamma0'] = 1.0 / res['Tau0']
    res['xi'] = (res['Beta'] - 1.0) / (res['Beta'] + 1.0)
    res['Gammac'] = res['Gamma0'] * (1.0 + 2.0 / (res['Beta'] - 1.0))
    res['strainc'] = np.power((res['xi'] / res['AAvg']) * np.power(res['xi']*res['Gamma0'] /\
                   res['Omega'] , res['Beta']), 1.0 / res['n'])
    res['psi'] = res['K'] / (res['AAvg'] * np.power(res['xi'], res['Beta']) * np.power(res['Tau0'] * res['Omega'], res['Beta'] - 1.0))
    res['psic'] = (4.0 * res['Beta']) / (np.power(res['Beta'], 2.0) - 1.0)
    res['g'] = res['psi'] / res['psic']
    return res

def ModelParam_SetG_strainc(normParams, g_val, strainc_val):
    res = CalcNormParams(normParams)
    res['Tau0'] = 1.0 / res['Gamma0']
    res['strainc'] = strainc_val
    res['g'] = g_val
    res['AAvg'] = (res['xi'] / np.power(res['strainc'], res['n'])) *\
                    np.power(res['xi'] * res['Gamma0'] / res['Omega'], res['Beta'])
    res['K'] = res['g'] * (4.0 * np.power(res['Tau0'] * res['Omega'], res['Beta'] - 1.0) *\
                                   res['AAvg'] * res['Beta'] * np.power(res['xi'], res['Beta'])) /\
                                   (np.power(res['Beta'], 2.0) - 1.0)
    return res

def ModelParam_SetG(params, g_val):
    res = params.copy()
    res['psi'] = g_val * res['psic']
    res['AAvg'] = (res['K'] / (res['psi'] * np.power(res['xi'], res['Beta']))) *\
                    np.power(res['Gamma0'] / res['Omega'], res['Beta'] - 1.0)
    return res

if __name__ == '__main__':

    paramVarKeys = ['K', 'n', 'Tau0', 'AAvg', 'AVar']

    
    ### INITIALIZING PARAMETERS ###
    
    NSites = 32
    NumSim = 10
    
    params = {
                'NSim'              : NumSim,
                'NSites'            : NSites,
                'K_list'            : defK,
                'n_list'            : defn,
                'Tau0_list'         : defTau0,
                'AAvg_list'         : defAAvg,
                'AVar_list'         : defAVar,
                'K'                 : defK[0],
                'n'                 : defn[0],
                'Tau0'              : defTau0[0],
                'AAvg'              : defAAvg[0],
                'AVar'              : defAVar[0],
                'Omega'             : defOmega,
                'Beta'              : 2.0,
                'gFromFile'         : False,
                'gFile'             : 'XXX',
                'gMin'              : defGammaMin,
                'gMax'              : defGammaMax,
                'gPPD'              : 100,
                'outFolder'         : SaveRoot,
                'outPre'            : 'outXXXX_',
                'statPre'           : 'statXXXX_',
                'parName'           : 'paramsXXXX.txt',
                'chiThr'            : FastSitesThr,
                'wRate'             : 1.0,
                'wChi'              : 100.0,
                'dComb'             : True
            }
    
    
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
    
    print(len(prodList))
    
    CALC = True
    
    n_parall_proc = 4
    if CALC:
        for iCount in range(0, len(prodList), n_parall_proc):
            proc_list = []
            for iProc in range(iCount, min(iCount+n_parall_proc, len(prodList))):
                print(prodList[iProc])
                for iKey in range(len(prodList[iProc])):
                    params[paramVarKeys[iKey]] = prodList[iProc][iKey]
                cur_savename = str(SaveFname).replace('NNNN', str(iProc))
                GenParamFile(params, iProc, save_name=cur_savename)
                proc = subprocess.Popen(SimProgPath + " " + cur_savename, stdout=subprocess.PIPE)
                proc_list.append(proc)
            for cur_p in proc_list:
                cur_p.communicate()
    col_max_num = 20000
    fpar = open(params['outFolder']+'__par.txt', 'w')
    fpar.write('#')
    for key in paramVarKeys:
        fpar.write('\t' + str(key))
    data_first = []
    for iFile in range(0, len(prodList), col_max_num):
        fgamma = open(params['outFolder']+'_all_gamma_' + str(iFile//col_max_num).zfill(4) + '.txt', 'w')
        fchi = open(params['outFolder']+'_all_chi_' + str(iFile//col_max_num).zfill(4) + '.txt', 'w')
        fgamma.write('g')
        fchi.write('g')
        data_g = []
        data_c = []
        for iCol in range(iFile, min(iFile+col_max_num, len(prodList))):
            cur_col_g = []
            cur_col_c = []
            strname = str(params['outFolder']+params['outPre']+'ps_all.txt').replace('XXXX', str(iCol).zfill(4))
            with open(strname) as fin:
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
                strname = str(params['outFolder']+params['statPre']+'ps_all.txt').replace('XXXX', str(iCol).zfill(4))
                with open(strname) as fin:
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
                    