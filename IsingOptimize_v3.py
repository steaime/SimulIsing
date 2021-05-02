import sys
import subprocess
import numpy as np
import datetime

### SAMPLE SELECTION ###

SAMPLE_SEL = 'EMULSION65'

if (SAMPLE_SEL == 'LUDOX45'):
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run1_Ludox45\\'
    ExpDataPath     = SaveRoot + "ExpLudox.txt"
    if False: #run 1
        defK = 8.0
        defn = 3.0
        defTau0 = 9000.0
        defAAvg = 0.0011
        defAVar = 0.08
        varK = 0.2
        varn = 0.1        
        varTau0 = 0.1
        varAAvg = 0.2
        varAVar = 0.5
    else: #run 2
        defK = 7.93718519
        defn = 2.767168107
        defTau0 = 8217.835594
        defAAvg = 0.001224168
        defAVar = 0.065487347
        varK = 0.01
        varn = 0.005       
        varTau0 = 0.005
        varAAvg = 0.01
        varAVar = 0.02
    defOmega = 3.14
    FastSitesThr = 7
elif (SAMPLE_SEL == 'PNIPRESH'):
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run2_PnipamPresh\\'
    ExpDataPath     = SaveRoot + "Exp_Pnipam.txt"
    defK = 0.00135
    defn = 1.9
    defTau0 = 70000.0
    defAAvg = 1.8e-8
    defAVar = 0.11
    varK = 0.2
    varn = 0.1        
    varTau0 = 0.1
    varAAvg = 0.2
    varAVar = 0.5
    defOmega = 3.14
    FastSitesThr = 7
elif (SAMPLE_SEL == 'EMULSION65'):
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run3_Em65\\'
    ExpDataPath     = SaveRoot + "Exp_Em65.txt"
    defK = 1e9
    defn = 8.0
    defTau0 = 617.0
    defAAvg = 8.8e5
    defAVar = 0.35
    varK = 0.5
    varn = 0.5
    varTau0 = 0.2
    varAAvg = 0.5
    varAVar = 0.5
    defOmega = 6.28
    FastSitesThr = 5
else:
    print('ERROR: sample not recognized')
    sys.exit()

# INPUT/OUTPUT SETTINGS
ConfigRoot      = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\config\\'
TemplateFname   = ConfigRoot + 'simParamsTemplate.ini'
SaveFname       = ConfigRoot + 'simParams_0.ini'
SimProgPath     = "C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\SimulIsing_v3\\x64\\Release\\SimulIsing_v3_190524.exe"
GlobalOut       = SaveRoot + "OutFit.txt"
LogOut          = SaveRoot + "LogFit.txt"
sCommFit        = SimProgPath + " " + SaveFname + " -COMPARE " + ExpDataPath + " -SILENT"
sCommSave        = SimProgPath + " " + SaveFname
iGlobalCount    = 0



def GenParamFile(dict_params, template_fname=TemplateFname, save_fname=SaveFname, key_firstlast='%'):
    
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
                line_str = line_str.replace(key_firstlast+key+key_firstlast, str(dict_params[key]))
            f.write(line_str)
            f.write('\n')

def SweepAndCompare(params, fOut=None, PrintKeys=[]):
    global iGlobalCount
    iGlobalCount += 1
    GenParamFile(params)
    proc = subprocess.Popen(sCommFit, stdout=subprocess.PIPE)
    strRes = str(proc.stdout.read(), 'utf-8')
    resList = strRes.split()
    for i in range(len(resList)):
        resList[i] = float(resList[i])
    if (fOut!=None):
        fOut.write('\n' + str(iGlobalCount))
        for key in PrintKeys:
            fOut.write('\t' + str(params[key]))
        for i in range(len(resList)):
            fOut.write('\t' + str(resList[i]))
        fOut.flush()
    return resList





if __name__ == '__main__':
    
    ### INITIALIZING PARAMETERS ###
    
    LatticeSize     = 128
    NumSim          = 3
    
    LatticeSize     = 64
    NumSim          = 1
    
    constPar = {
                'Omega'             : defOmega,
                'Beta'              : 2.0,
                'gFromFile'         : False,
                'gFile'             : 'XXX',
                'gMin'              : 0.005,
                'gMax'              : 0.5,
                'gPPD'              : 100,
                'outFolder'         : SaveRoot,
                'chiThr'            : FastSitesThr,
                'wRate'             : 1.0,
                'wChi'              : 10.0,
                'dComb'             : False,
                'NSim'              : NumSim,
                'NSites'            : LatticeSize
            }
    initPar = {
                'K'                 : defK,
                'n'                 : defn,
                'Tau0'              : defTau0,
                'AAvg'              : defAAvg,
                'AVar'              : defAVar
            }
    varPar = {
                'K'                 : varK,
                'n'                 : varn,
                'Tau0'              : varTau0,
                'AAvg'              : varAAvg,
                'AVar'              : varAVar
            }

    ### INITIALIZING OUTPUT FILES ###

    fLog = open(LogOut, 'w')
    fOut = open(GlobalOut, 'w')

    strLog  = 'Program started on: ' + str(datetime.datetime.now())
    strLog += '\n\nLatticeSize\tNumSimulations : '
    strLog += '\n' + str(LatticeSize) + '\t' + str(NumSim)
    strLog += '\n\nProgram parameters:'
    strLog += '\n - Detailed output saved to file : ' + str(GlobalOut)
    strLog += '\n - Config root folder : ' + str(ConfigRoot)
    strLog += '\n - Template *.ini filename : ' + str(TemplateFname)
    strLog += '\n - Save root folder : ' + str(SaveRoot)
    strLog += '\n - Work *.ini filename : ' + str(SaveFname)
    strLog += '\n - Simulation exe : ' + str(SimProgPath)
    strLog += '\n - Experimental data file : ' + str(ExpDataPath)
    strLog += '\n - Simulation fit command line : ' + str(sCommFit)
    strLog += '\n - Simulation save command line : ' + str(sCommFit)
    strLog += '\n\nConstant simulation parameters:'
    for key in constPar.keys():
        strLog += '\n - ' + str(key) + '\t: ' + str(constPar[key])
    strLog += '\n\nInitial variable simulation parameters:'
    for key in initPar.keys():
        strLog += '\n - ' + str(key) + '\t: ' + str(initPar[key])
    fLog.write(strLog)
    fLog.flush()
    print(strLog)

    PrintKeys = []
    for key in initPar.keys():
        PrintKeys.append(str(key))
    fOut.write('#')
    for key in PrintKeys:
        fOut.write('\t' + str(key))
    fOut.write('\tDiffRate\tDiffChi')
    
        
    for i in range(10000):
        
        rndPar = initPar.copy()
        for key in initPar.keys():
            cur_rnd = np.random.normal(0, varPar[key])
            rndPar[key] = np.exp(np.log(rndPar[key])+cur_rnd)
    
        curDiff = SweepAndCompare({**constPar, **rndPar}, fOut, PrintKeys)
        print(curDiff)

    ### TERMINATE PROGRAM

    fLog.close()
    fOut.close()