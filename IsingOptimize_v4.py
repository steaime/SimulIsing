import sys
import subprocess
import numpy as np
import datetime

# INPUT/OUTPUT SETTINGS
ConfigRoot      = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\config\\'
TemplateFname   = ConfigRoot + 'simParamsTemplate.ini'
SaveRoot        = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v4\\run1_Ludox45\\'
SaveFname       = SaveRoot + 'simParams_0.ini'
SimProgPath     = "C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\SimulIsing_v3\\x64\\Release\\SimulIsing_v3_190524.exe"
ExpDataPath     = SaveRoot + "ExpLudox.txt"
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

def CalcIndependentParam(modelParams):
    res = modelParams.copy()
    res['Gamma0'] = 1.0 / res['Tau0']
    res['csi'] = (res['Beta'] - 1.0) / (res['Beta'] + 1.0)
    res['Gammac'] = res['Gamma0'] * (1.0 + 2.0 / (res['Beta'] - 1.0))
    res['strainc'] = np.power((res['csi'] / res['AAvg']) * np.power(res['csi'] /\
                   (res['Omega'] * res['Tau0']), res['Beta']), 1.0 / res['n'])
    res['g'] = (res['K'] * (np.power(res['Beta'], 2.0) - 1.0)) /\
            (4.0 * np.power(res['Tau0'] * res['Omega'], res['Beta'] - 1.0) *\
             res['AAvg'] * res['Beta'] * np.power(res['csi'], res['Beta']))
    return res

def CalcModelParam(independentParams):
    res = independentParams.copy()
    res['Tau0'] = 1.0 / res['Gamma0']
    res['AAvg'] = (res['csi'] / np.power(res['strainc'], res['n'])) *\
                    np.power(res['csi'] * res['Gamma0'] / res['Omega'], res['Beta'])
    res['K'] = res['g'] * (4.0 * np.power(res['Tau0'] * res['Omega'], res['Beta'] - 1.0) *\
                                   res['AAvg'] * res['Beta'] * np.power(res['csi'], res['Beta'])) /\
                                   (np.power(res['Beta'], 2.0) - 1.0)
    return res

if __name__ == '__main__':
    
    #paramFitKeys = ['n', 'g', 'Gamma0', 'strainc', 'AVar']
    paramFitKeys = ['K', 'Gamma0', 'strainc']
   
    ### SAMPLE SELECTION ###
    
    SAMPLE_SEL = 'LUDOX45'
    
    if (SAMPLE_SEL == 'LUDOX45'):
        defK = 110.0
        defn = 4.0
        defTau0 = 2400.0
        defAAvg = 0.019
        defAVar = 0
        defOmega = 3.14
        defGammaMin = 0.005
        defGammaMax = 0.5
        FastSitesThr = 5
        defn_bounds = [1.0, 10.0]
        defg_bounds = [0.5, 1.3]
        defGamma0_bounds = [1.0/30000.0, 1.0/300.0]
        defstrainc_bounds = [0.0, 0.08]
        defAVar_bounds = [0.0, 1.0]
    else:
        print('ERROR: sample not recognized')
        sys.exit()
    
    ### INITIALIZING PARAMETERS ###
    
    #LatticeSizes    = [ 4,  8, 16, 32, 64, 128]
    #NumSim          = [32, 16, 8,  4,  2,   1]
    LatticeSizes    = [1]
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
                'chiThr'            : FastSitesThr,
                'wRate'             : 1.0,
                'wChi'              : 20.0,
                'dComb'             : True,
                'n_bounds'          : defn_bounds,
                'g_bounds'          : defg_bounds,
                'Gamma0_bounds'     : defGamma0_bounds,
                'strainc_bounds'    : defstrainc_bounds,
                'AVar_bounds'       : defAVar_bounds
            }
    
    # Add to initPar the derived parameters on which to fit
    initPar = CalcModelParam(CalcIndependentParam(initPar))
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
        strLog += '\n' + str(key) + '\t' + str(fitPar[key]) + '\t' + str(fitPar[key+'_bounds'][0]) +\
                    '\t' + str(fitPar[key+'_bounds'][1]) + '\t' + str(initParIncr[key]) +\
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
                    cur_par = np.clip(np.exp(np.log(fitPar[key]) + parIncrDir[key]*parIncr[key]), fitPar[key+'_bounds'][0], fitPar[key+'_bounds'][1])
                    # Only do that if clipping didn't keep parameter the same
                    if (workPar[key] != cur_par):
                        workPar[key] = cur_par
                        # Recalculate model parameters
                        workPar = CalcModelParam(workPar)
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
                        cur_par = np.clip(np.exp(np.log(fitPar[key]) - parIncrDir[key]*parIncr[key]), fitPar[key+'_bounds'][0], fitPar[key+'_bounds'][1])
                        if (workPar[key] != cur_par):
                            workPar[key] = cur_par
                            # Recalculate model parameters
                            workPar = CalcModelParam(workPar)
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
            
                # Next loop we will start from the best parameters found in currend loop
                fitPar = CalcModelParam(fitPar)
                

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