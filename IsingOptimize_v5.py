import sys
import subprocess
import numpy as np
import datetime

### SAMPLE SELECTION ###

SAMPLE_SEL = 'LUDOX45_n3'

if (SAMPLE_SEL == 'LUDOX45_n4'):
    defK = 110.0
    defn = 4.0
    defTau0 = 2400.0
    defAAvg = 0.05
    defAVar = 0.5
    defOmega = 3.14
    defGammaMin = 0.005
    defGammaMax = 0.5
    FastSitesThr = 5
    defAVar_bounds = [0.0, 1.0]
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run1_Ludox45\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"
elif (SAMPLE_SEL == 'LUDOX45_n3'):
    defK = 8.0
    defn = 3.0
    defTau0 = 9000.0
    defAAvg = 1.112E-3
    defAVar = 0.08
    defOmega = 3.14
    defGammaMin = 0.001
    defGammaMax = 1.0
    FastSitesThr = 7
    defAVar_bounds = [0.0, 1.0]
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run1_Ludox45\\'
    ExpDataPath = SaveRoot + "ExpLudox.txt"

elif (SAMPLE_SEL == 'EM65_n20'):
    defK = 3e23
    defn = 20.0
    defTau0 = 617.0
    defAAvg = 2e20
    defAVar = 0.5
    defOmega = 6.28
    defGammaMin = 0.001
    defGammaMax = 0.08
    FastSitesThr = 4
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run3_Em65\\'
    ExpDataPath = SaveRoot + "Exp_Em65.txt"
elif (SAMPLE_SEL == 'EM65'):
    defK = 1e9
    defn = 8.0
    defTau0 = 617.0
    defAAvg = 8.8e5
    defAVar = 0.35
    defOmega = 6.28
    defGammaMin = 0.001
    defGammaMax = 0.08
    FastSitesThr = 4
    SaveRoot = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run3_Em65\\'
    ExpDataPath = SaveRoot + "Exp_Em65.txt"
else:
    print('ERROR: sample not recognized')
    sys.exit()

# INPUT/OUTPUT SETTINGS
ConfigRoot      = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\config\\'
TemplateFname   = ConfigRoot + 'simParamsTemplate.ini'
SaveFname       = SaveRoot + 'simParams_0.ini'
SimProgPath     = "C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\SimulIsing_v3\\x64\\Release\\SimulIsing_v3_190524.exe"
GlobalOut       = SaveRoot + "OutFit.txt"
LogOut          = SaveRoot + "LogFit.txt"
sCommFit        = SimProgPath + " " + SaveFname + " -COMPARE " + ExpDataPath + " -SILENT"
sCommSave        = SimProgPath + " " + SaveFname
iGlobalCount    = 0

"""
WORKFLOW:
- Program input:
    - set of initial parameters
    - set of initial increment steps
    - set of initial increment direction (+/-)
    - increment step decrease factor
    - set of ultimate increment steps
    - set of boolean flags "parameter enabled"
    - list of lattice sizes
- For each lattice size, small to large: 
    - Run simulations on starting set of parameters, get reference deviation
    - Print to output current set of parameters, current increment step and direction, 
      if parameter is enabled, current reference deviation
    - Starting from current set of parameters, for each enabled parameter...
        - try to increment parameter in previous direction
        - run simulations, get deviation, print to output
        - if deviation is decreased from reference deviation:
            - set next parameter to incremented one
        - else:
            - try to increment parameter in opposite direction
            - run simulations, get deviation, print to output
            - if deviation is decreased from reference deviation:
                - set next parameter to incremented one
                - change increment direction
            - else:
                - keep next parameter to current one
                - compare the two deviations obtained by incrementing parameter
                - set increment direction to the most promising one
                - decrease increment step by defined decrease factor
                - if increment step falls below ultimate increment step:
                    - fix parameter to current value (= disable parameter)
    - If parameters are all disabled: break cycle, get to larger lattice size
"""

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
    
    #paramFitKeys = ['AAvg', 'AVar']
    paramFitKeys = ['AVar']
   

    
    ### INITIALIZING PARAMETERS ###
    
    LatticeSizes    = [16]
    NumSim          = [1]
    MaxLoopNum      = 100
    
    
    ParamIncrement = 0.3
    ParamIncDecrease = 0.7
    ParamIncThr = 0.01
    
    initPar = {
                'K'                 : defK,
                'n'                 : defn,
                'Tau0'              : defTau0,
                'AAvg'              : defAAvg,
                'AVar'              : defAVar,
                'Omega'             : defOmega,
                'Beta'              : 2.0,
                'gFromFile'         : False,
                'gFile'             : 'XXX',
                'gMin'              : defGammaMin,
                'gMax'              : defGammaMax,
                'gPPD'              : 100,
                'outFolder'         : SaveRoot,
                'outPre'            : 'outXXXX_',
                'chiThr'            : FastSitesThr,
                'wRate'             : 1.0,
                'wChi'              : 100.0,
                'dComb'             : True,
                'AVar_bounds'       : defAVar_bounds
            }
    
    
    if False:
        npar = CalcNormParams(initPar)
        initPar = ModelParam_SetG(npar, 0.8)
        print(initPar)
        sys.exit()
    
    if True:
        npar = CalcNormParams(initPar)
        print(npar)
        print('\n\n\n')
        npar = ModelParam_SetG_strainc(initPar, 0.98, 0.031)
        print(npar)
        sys.exit()
        
    if True:
        npar = CalcNormParams(initPar)
        print(npar)
        sys.exit()
    
    fitPar = initPar.copy()

    initParIncr = {}
    parIncrDir = {}
    parIncrDecr = {}
    parIncrThr = {}
    parEnabled = {}
    parIncrGrow = {}
    for key in paramFitKeys:
        initParIncr[key]    = ParamIncrement
        parIncrDir[key]     = 1.0
        parIncrDecr[key]    = ParamIncDecrease
        parIncrThr[key]     = ParamIncThr
        parEnabled[key]     = True
        parIncrGrow[key]    = True
    parIncr = initParIncr.copy()

    ### INITIALIZING OUTPUT FILES ###

    fLog = open(LogOut, 'w')
    fOut = open(GlobalOut, 'w')

    strLog  = 'Program started on: ' + str(datetime.datetime.now())
    strLog += '\n\nLatticeSize\tNumSimulations : '
    for i in range(len(LatticeSizes)):
        strLog += '\n' + str(LatticeSizes[i]) + '\t' + str(NumSim[i])
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
    strLog += '\n - Parameter initial relative increment  : ' + str(ParamIncrement)
    strLog += '\n - Increment relative reduction  : ' + str(ParamIncDecrease)
    strLog += '\n - Parameter increment threshold : ' + str(ParamIncThr)
    strLog += '\n\nSimulation parameters:'
    for key in fitPar.keys():
        strLog += '\n - ' + str(key) + '\t: ' + str(fitPar[key])
    strLog += '\n\nFit parameters:'
    strLog += '\nParam\tValue\tMin\tMax\tIncrement\tDirection\tEnabled\tStepDecr\tThreshold'
    for key in paramFitKeys:
        strLog += '\n' + str(key) + '\t' + str(fitPar[key]) + '\t' 
        if (key+'_bounds' in fitPar):
            strLog += str(fitPar[key+'_bounds'][0]) + '\t' + str(fitPar[key+'_bounds'][1]) 
        else:
            strLog += '-\t-'
        strLog += '\t' + str(initParIncr[key]) +\
                    '\t' + str(parIncrDir[key]) + '\t' + str(parEnabled[key]) +\
                    '\t' + str(parIncrDecr[key]) + '\t' + str(parIncrThr[key])
    fLog.write(strLog)
    fLog.flush()
    print(strLog)

    PrintKeys = ['ID', 'NSites', 'nLoop', 'K', 'AAvg', 'Tau0']
    for key in paramFitKeys:
        PrintKeys.append(str(key))
    fOut.write('#')
    for key in PrintKeys:
        fOut.write('\t' + str(key))
    fOut.write('\tDiff')
    
    ### LOOP ON LATTICE SIZES ###
        
    for iSizeIdx in range(len(LatticeSizes)):
    
        ### INITIALIZE NEW LATTICE SIZE ###
        
        fitPar['nLoop'] = 0
        fitPar['ID'] = 'FIT'
        bkpPar = fitPar.copy()
        bkpPar['ID'] = 'BKP'
        if (iSizeIdx > 0):
            for key in parIncr.keys():
                parEnabled[key]  = True
                parIncrGrow[key] = True
                if False:
                    # Do this if you want to start again from a small increment
                    # It saves a bit of time, but it is more likely trapped
                    # in secondary minima
                    parIncr[key] = parIncrThr[key] * np.power(parIncrDecr[key], -3)
                else:
                    parIncr[key] = ParamIncrement
        
        fitPar['NSites'] = LatticeSizes[iSizeIdx]
        fitPar['NSim']   = NumSim[iSizeIdx]        
        refDiff = SweepAndCompare(fitPar, fOut, PrintKeys)
        
        strLog  = '\n\n[' + str(datetime.datetime.now()) + ']'
        strLog += '\nStarting optimization with NSites=' + str(fitPar['NSites'])
        strLog += '\nNumber of parallel simulations : ' + str(fitPar['NSim'])
        strLog += '\nParam\tValue\tIncrement\tDirection'
        for key in paramFitKeys:
            strLog += '\n' + str(key) + '\t' + str(fitPar[key]) + '\t' + str(parIncr[key])
            strLog += '\t' + str(parIncrDir[key])
        strLog += '\nInitial (reference) difference: ' + str(refDiff)
        strLog += '\n\n#\tDiff\tnActive'
        for key in paramFitKeys:
            strLog += '\t' + str(key) + '\t' + str(key) + '_inc'
        fLog.write(strLog)
        fLog.flush()
        print(strLog)
        
        # curDiff: difference at the beninning of optimization loop
        curDiff = refDiff
        # nLoopCount: count number of optimization loops
        nLoopCount = 0
        
        while (sum(1 for x in parEnabled.values() if x==True) > 0):
            
            for key in paramFitKeys:
                if (parEnabled[key]):
                    workPar = fitPar.copy()
                    workPar['ID'] = 'WRK_' + str(key)
                    # Try to increment parameter in previous direction
                    cur_par = np.exp(np.log(fitPar[key]) + parIncrDir[key]*parIncr[key])
                    if (key+'_bounds' in fitPar) : cur_par = np.clip(cur_par, fitPar[key+'_bounds'][0], fitPar[key+'_bounds'][1])
                    # Only do that if clipping didn't keep parameter the same
                    if (workPar[key] != cur_par):
                        workPar[key] = cur_par
                        # Run simulations, get deviation, print to output
                        workDiff = SweepAndCompare(workPar, fOut, PrintKeys)
                    else:
                        # Just skip this
                        workDiff = curDiff + 1
                    # if deviation is decreased from reference deviation:
                    if (workDiff < curDiff):
                        # This is our new reference difference
                        curDiff = workDiff
                        # Set next parameter to incremented one
                        fitPar[key] = workPar[key]
                        # Check if increment has to be amplified
                        if (parIncrGrow[key]) :  parIncr[key] /= parIncrDecr[key]
                    else:
                        # Try to increment parameter in opposite direction
                        cur_par = np.exp(np.log(fitPar[key]) - parIncrDir[key]*parIncr[key])
                        # Eventually clip with parameter bounds
                        if (key+'_bounds' in fitPar) : cur_par = np.clip(cur_par, fitPar[key+'_bounds'][0], fitPar[key+'_bounds'][1])
                        if (workPar[key] != cur_par):
                            workPar[key] = cur_par
                            # Run simulations, get deviation, print to output
                            workDiffRev = SweepAndCompare(workPar, fOut, PrintKeys)
                        else:
                            workDiffRev = curDiff + 1
                        # if deviation is decreased from reference deviation:
                        if (workDiffRev < curDiff):
                            # This is our new reference difference
                            curDiff = workDiffRev
                            # Set next parameter to incremented one
                            fitPar[key] = workPar[key]
                            # Change increment direction
                            parIncrDir[key] *= -1.0
                            # If this is not the first loop, disable eventual increment growth
                            if (nLoopCount > 0):
                                parIncrGrow[key] = False
                            # Check if increment has to be amplified
                            elif (parIncrGrow[key]):
                                parIncr[key] /= parIncrDecr[key]
                        else:
                            # Set increment direction to the most promising one
                            if (workDiffRev < workDiff) : parIncrDir[key] *= -1.0
                            # Decrease increment step by defined decrease factor
                            # Don't go below parIncrThr[key]
                            parIncr[key] = np.clip(parIncr[key]*parIncrDecr[key], parIncrThr[key], parIncr[key]*parIncrDecr[key]+1)
                            # Uncomment this to only enable parameter growth once
                            # Probably faster, but potentially less reliable
                            if False:
                                # Disable increment growth
                                parIncrGrow[key] = False
                            # Set to True to disable parameters once they reached the threschold
                            # Probably faster, but potentially less reliable
                            if True:
                                # if increment step falls below ultimate increment step, disable parameter
                                if (parIncr[key] <= parIncrThr[key]) : parEnabled[key] = False
                            

            nLoopCount += 1
            fitPar['nLoop'] = nLoopCount
            
            # Check if all parameters are good
            nActive = 0
            for key in paramFitKeys:
                if (parIncr[key] > parIncrThr[key]) : nActive += 1
            if (nActive == 0):
                strLog = '\n\nLoop #' + str(nLoopCount) + ': fit converged! Fit parameters:'
                strLog += '\nParam\tInitial\tFinal'
                for key in paramFitKeys:
                    strLog += '\n' + str(key) + '\t' + str(bkpPar[key]) + '\t' + str(fitPar[key])
                strLog += '\nInitial difference: ' + str(refDiff) + ' | Final : ' + str(curDiff)
            else :
                strLog  = '\nLoop #' + str(nLoopCount) + ': difference lowered to ' + str(curDiff)
                strLog +=' (' + str(nActive) + '/' + str(len(paramFitKeys)) + ' active parameters)'
                strLog = '\n' + str(nLoopCount) + '\t' + str(curDiff) + '\t' + str(nActive)
                for key in paramFitKeys:
                    strLog += '\t' + str(fitPar[key]) + '\t' + str(parIncr[key])
            fLog.write(strLog)
            fLog.flush()
            print(strLog)
            
            if (nActive == 0):
                break
            elif (nLoopCount >= MaxLoopNum):
                strLog = '\nMaximum number of loops reached. Optimization loop ended.'
                fLog.write(strLog)
                print(strLog)
                break
            
    strLog  = '\n\n[' + str(datetime.datetime.now()) + ']'
    strLog += '\nProcedure terminated!\n\nInitial and final model parameters:'
    strLog += '\nParam\tInitial\tFinal'
    for key in paramFitKeys:
        strLog += '\n' + str(key) + '\t' + str(initPar[key]) + '\t' + str(fitPar[key])
    strLog += '\nFinal difference: ' + str(refDiff)
    strLog += '\n\n\nNow running full simulation and saving output...'
    fLog.write(strLog)
    fLog.flush()
    print(strLog)
    
    ### RUN FINAL SIMULATION
    if (fitPar['NSim'] > 3):
        fitPar['NSim'] = 3
    GenParamFile(fitPar)
    proc = subprocess.Popen(sCommSave, stdout=subprocess.PIPE)
    
    strLog  = '\n\n[' + str(datetime.datetime.now()) + ']'
    strLog += '\n\nResult:\n'
    strLog += str(proc.stdout.read(), 'utf-8')
    strLog += '\n\n...done!'
    fLog.write(strLog)
    fLog.flush()
    print(strLog)

    ### TERMINATE PROGRAM

    fLog.close()
    fOut.close()
    
    
    print(CalcNormParams(fitPar))