import sys
import subprocess
import numpy as np
import datetime

# INPUT/OUTPUT SETTINGS
ConfigRoot      = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\config\\'
TemplateFname   = ConfigRoot + 'simParamsTemplate.ini'
SaveRoot        = 'C:\\Users\\STEAIME\\Documents\\Montpellier\\Ludox\\Ising\\SimulIsing\\out\\v3\\run1_Ludox45\\'
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





if __name__ == '__main__':
    
    ### SAMPLE SELECTION ###
    
    SAMPLE_SEL = 'LUDOX45'
    
    if (SAMPLE_SEL == 'LUDOX45'):
        defK = 7.0
        defn = 3.0
        defTau0 = 8500.0
        defAAvg = 0.001
        defAVar = 0.08
        defOmega = 3.14
        defGammaMin = 0.005
        defGammaMax = 0.5
        FastSitesThr = 5
    else:
        print('ERROR: sample not recognized')
        sys.exit()
    
    ### INITIALIZING PARAMETERS ###
    
    LatticeSizes    = [ 4,  8, 16, 32, 64, 128]
    NumSim          = [64, 32, 16,  4,  1,   1]
    MaxLoopNum      = 20
    
    ParamIncrement = 0.1
    ParamIncDecrease = 0.5
    ParamIncThr = 0.001
    
    constPar = {
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
                'wChi'              : 10.0,
                'dComb'             : True
            }
    sitesPar = {
                'NSim'              : NumSim[0],
                'NSites'            : LatticeSizes[0]
            }
    initPar = {
                'K'                 : defK,
                'n'                 : defn,
                'Tau0'              : defTau0,
                'AAvg'              : defAAvg,
                'AVar'              : defAVar,
            }
    fitPar = initPar.copy()
    initParIncr = {}
    parIncrDir = {}
    parIncrDecr = {}
    parIncrThr = {}
    parEnabled = {}
    parIncrGrow = {}
    for key in initPar.keys():
        initParIncr[key]    = ParamIncrement
        parIncrDir[key]     = 1.0
        parIncrDecr[key]    = ParamIncDecrease
        parIncrThr[key]     = ParamIncThr
        parEnabled[key]     = True
        parIncrGrow[key]    = False
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
    strLog += '\n\nConstant simulation parameters:'
    for key in constPar.keys():
        strLog += '\n - ' + str(key) + '\t: ' + str(constPar[key])
    strLog += '\n\nInitial variable simulation parameters:'
    for key in initPar.keys():
        strLog += '\n - ' + str(key) + '\t: ' + str(initPar[key])
    fLog.write(strLog)
    print(strLog)

    PrintKeys = ['NSites']
    for key in initPar.keys():
        PrintKeys.append(str(key))
    fOut.write('#')
    for key in PrintKeys:
        fOut.write('\t' + str(key))
    fOut.write('\tDiff')
    
    ### LOOP ON LATTICE SIZES ###
        
    for iSizeIdx in range(len(LatticeSizes)):
    
        ### INITIALIZE NEW LATTICE SIZE ###
        
        bkpPar = fitPar.copy()
        if (iSizeIdx > 0):
            for key in parIncr.keys():
                parEnabled[key]  = True
                parIncrGrow[key] = True
        
        sitesPar['NSites'] = LatticeSizes[iSizeIdx]
        sitesPar['NSim']   = NumSim[iSizeIdx]        
        refDiff = SweepAndCompare({**constPar, **sitesPar, **fitPar}, fOut, PrintKeys)
        
        strLog  = '\n\n[' + str(datetime.datetime.now()) + ']'
        strLog += '\nStarting optimization with NSites=' + str(sitesPar['NSites'])
        strLog += '\nNumber of parallel simulations : ' + str(sitesPar['NSim'])
        strLog += '\nParam\tValue\tIncrement\tDirection\tEnabled'
        for key in fitPar.keys():
            strLog += '\n' + str(key) + '\t' + str(fitPar[key]) + '\t' + str(parIncr[key])
            strLog += '\t' + str(parIncrDir[key]) + '\t' + str(parEnabled[key])
        strLog += '\nInitial (reference) difference: ' + str(refDiff)
        fLog.write(strLog)
        fLog.flush()
        print(strLog)
        
        # curDiff: difference at the beninning of optimization loop
        curDiff = refDiff
        # minDiff: minimum difference encountered in single parameter variation
        minDiff = curDiff
        # newPar: each parameter minimizing difference after optimization loop
        newPar = fitPar.copy()
        # minPar: same as fitPar, with only one parameter change, the one that produces minDiff        
        minPar = fitPar.copy()
        # nLoopCount: count number of optimization loops
        nLoopCount = 0
        
        while (sum(1 for x in parEnabled.values() if x==True) > 0):
            
            if (nLoopCount > 0):
                curDiff = SweepAndCompare({**constPar, **sitesPar, **fitPar}, fOut, PrintKeys)
                if (curDiff > minDiff):
                    fitPar = minPar.copy()
                    curDiff = minDiff
                    
            minDiff = curDiff
            newPar = fitPar.copy()
            minPar = fitPar.copy()
            
            for key in initPar.keys():
                if (parEnabled[key]):
                    workPar = fitPar.copy()
                    # Try to increment parameter in previous direction
                    workPar[key] = np.exp(np.log(fitPar[key]) + parIncrDir[key]*parIncr[key])
                    # Run simulations, get deviation, print to output
                    workDiff = SweepAndCompare({**constPar, **sitesPar, **workPar}, fOut, PrintKeys)
                    # if deviation is decreased from reference deviation:
                    if (workDiff < curDiff):
                        # Set next parameter to incremented one
                        newPar[key] = workPar[key]
                        # Check if increment has to be amplified
                        if (parIncrGrow[key]):
                            parIncr[key] /= parIncrDecr[key]
                        # Keep track of minimum difference found
                        if (workDiff < minDiff): 
                            minDiff = workDiff
                            minPar = workPar.copy()
                    else:
                        # Try to increment parameter in opposite direction
                        workPar[key] = np.exp(np.log(fitPar[key]) - parIncrDir[key]*parIncr[key])
                        # Run simulations, get deviation, print to output
                        workDiffRev = SweepAndCompare({**constPar, **sitesPar, **workPar}, fOut, PrintKeys)
                        # if deviation is decreased from reference deviation:
                        if (workDiffRev < curDiff):
                            # Set next parameter to incremented one
                            newPar[key] = workPar[key]
                            # Change increment direction
                            parIncrDir[key] *= -1.0
                            # If this is not the first loop, disable eventual increment growth
                            if (nLoopCount > 0):
                                parIncrGrow[key] = False
                            # Check if increment has to be amplified
                            elif (parIncrGrow[key]):
                                parIncr[key] /= parIncrDecr[key]
                            # Keep track of minimum difference found
                            if (workDiffRev < minDiff): 
                                minDiff = workDiffRev
                                minPar = workPar.copy()
                        else:
                            # Set increment direction to the most promising one
                            if (workDiffRev < workDiff) : parIncrDir[key] *= -1.0
                            # Decrease increment step by defined decrease factor
                            parIncr[key] *= parIncrDecr[key]
                            # Disable increment growth
                            parIncrGrow[key] = False
                            # if increment step falls below ultimate increment step, disable parameter
                            if (parIncr[key] < parIncrThr[key]) : parEnabled[key] = False

            # Next loop we will start from the best parameters found in currend loop
            fitPar = newPar.copy()
            nLoopCount += 1
            
            if (sum(1 for x in parEnabled.values() if x==True) > 0):
                strLog = '\nLoop #' + str(nLoopCount) + ': difference lowered to ' + str(minDiff)
            else:
                strLog = '\nLoop #' + str(nLoopCount) + ': fit converged! Fit parameters:'
                strLog += '\nParam\tInitial\tFinal'
                for key in fitPar.keys():
                    strLog += '\n' + str(key) + '\t' + str(bkpPar[key]) + '\t' + str(fitPar[key])
                strLog += '\nInitial difference: ' + str(refDiff) + ' | Final : ' + str(minDiff)
            fLog.write(strLog)
            fLog.flush()
            print(strLog)
            
            if (nLoopCount >= MaxLoopNum):
                strLog = '\nMaximum number of loops reached. Optimization loop ended.'
                fLog.write(strLog)
                print(strLog)
                break
            
    strLog  = '\n\n[' + str(datetime.datetime.now()) + ']'
    strLog += '\nProcedure terminated!\n\nInitial and final model parameters:'
    strLog += '\nParam\tInitial\tFinal'
    for key in fitPar.keys():
        strLog += '\n' + str(key) + '\t' + str(initPar[key]) + '\t' + str(fitPar[key])
    strLog += '\nFinal difference: ' + str(refDiff)
    strLog += '\n\n\nNow running full simulation and saving output...'
    fLog.write(strLog)
    fLog.flush()
    print(strLog)
    
    ### RUN FINAL SIMULATION
    sitesPar['NSim'] = 1
    GenParamFile({**constPar, **sitesPar, **fitPar})
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