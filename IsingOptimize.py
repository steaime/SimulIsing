import subprocess
import sys
import os
import numpy as np

# INPUT/OUTPUT SETTINGS
froot           = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\config\\'
TemplateFname   = froot + 'simParamsTemplate.ini'
SaveFname       = froot + 'simParams_0.ini'
SimProgPath     = "C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\SimulIsing_v3\\x64\\Release\\SimulIsing_v3.exe"
ExpDataPath     = froot + "ExpLudox.txt"
GlobalOut       = froot + "ChiOut.txt"
sComm           = SimProgPath + " " + SaveFname + " -COMPARE " + ExpDataPath + " -SILENT"

# PARAMETERS
NumSites = 128
K = 8
n = 3
Tau0 = 9000
Omega = 3.14
Beta = 2
AlphaAvg = 0.0011
AlphaVar = 0.08
NumSim = 3
gamma_file = "XXX"
out_folder = "XXX"
FastSitesThr = 5
DiffWRate = 1.0
DiffWChi = 1.0


def GenParamFile(dict_params, template_fname=TemplateFname, save_fname=SaveFname):
    
    # READ TEMPLATE DATA
    template_data = []
    with open(template_fname) as f:
        for line in f:
            template_data.append(line)
    
    # WRITE PARAM FILE
    with open(save_fname, 'w') as f:
        for line in template_data:
            line_str = str(line).replace('\n', '')#.replace(' ', '')
            for key in dict_params.keys():
                line_str = line_str.replace(key, str(dict_params[key]))
            f.write(line_str)
            f.write('\n')

def SweepAndCompare(params):
    GenParamFile(params)
    proc = subprocess.Popen(sComm, stdout=subprocess.PIPE)
    #print("Result: {0}".format(str(proc.stdout.read(), 'utf-8')))
    return float(str(proc.stdout.read(), 'utf-8'))


if __name__ == '__main__':
    
    alphaavg_list = np.logspace(np.log10(AlphaAvg)-1, np.log10(AlphaAvg)+1, 21)
    alphavar_list = np.logspace(np.log10(AlphaVar)-2, np.log10(AlphaVar)+1, 31)
    Tau0_list = np.logspace(np.log10(Tau0)-0.2, np.log10(Tau0)+0.2, 7)
    K_list = np.logspace(np.log10(K)-0.4, np.log10(K)+0.4, 13)
    
    with open(GlobalOut, 'w') as fOut:
        
        fOut.write('K\tTau0\tAlphaAvg\tAlphaVar\tDiffRate\tDiffChi')

        kw_dict_base = {
                '%NUMSITES%'        : str(NumSites),
                '%K%'               : str(K),
                '%n%'               : str(n),
                '%Tau0%'            : str(Tau0),
                '%OMEGA%'           : str(Omega),
                '%BETA%'            : str(Beta),
                '%ALPHA_AVG%'       : str(AlphaAvg),
                '%ALPHA_VAR%'       : str(AlphaVar),
                '%NUM_SIM%'          : str(NumSim),
                '%GAMMA_FILE%'      : str(gamma_file),
                '%OUT_FOLDER%'      : str(out_folder),
                '%CHITHR%'          : str(FastSitesThr),
                '%DIFFWRATE%'       : str(DiffWRate),
                '%DIFFWCHI%'        : str(DiffWChi)
                }
        
        total_dict_list = []
        for cur_alpha in alphaavg_list:
            total_dict_list.append(kw_dict_base.copy())
            total_dict_list[-1]['%ALPHA_AVG%'] = cur_alpha
        for cur_alphavar in alphavar_list:
            total_dict_list.append(kw_dict_base.copy())
            total_dict_list[-1]['%ALPHA_VAR%'] = cur_alphavar
        for cur_K in K_list:
            total_dict_list.append(kw_dict_base.copy())
            total_dict_list[-1]['%K%'] = cur_K
        for cur_tau0 in Tau0_list:
            total_dict_list.append(kw_dict_base.copy())
            total_dict_list[-1]['%Tau0%'] = cur_tau0
        
        print("Total number of simulations: {0}".format(2*len(total_dict_list)))
        
        for cur_dict in total_dict_list:
            cur_dict['%DIFFWRATE%'] = 1.0
            cur_dict['%DIFFWCHI%'] = 0
            cur_wRate = SweepAndCompare(cur_dict)
            cur_dict['%DIFFWRATE%'] = 0
            cur_dict['%DIFFWCHI%'] = 1.0
            cur_wChi = SweepAndCompare(cur_dict)
            fOut.write('\n' + str(cur_dict['%K%']) + '\t' + str(cur_dict['%Tau0%']) + '\t' +\
                       str(cur_dict['%ALPHA_AVG%']) + '\t' + str(cur_dict['%ALPHA_VAR%']) + '\t' +\
                       str(cur_wRate) + '\t' + str(cur_wChi))
            fOut.flush()
            print([cur_wRate, cur_wChi])