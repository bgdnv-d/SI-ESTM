function results = RunSelfConsTransition(setup,rootDir,projName,systemParams,activeElements,modFiles,Reg,ProduserType,costYear,prosShare,el_cost,feedin_cost,endHour)
%
% FUNCTION Run_Self_Cons_Transition(setup, rootDir, projName, systemParams, activeElements, modFiles, Reg, ProduserType, costYear, prosShare, el_cost, feedin_cost, endHour)
%
% Runs a self consumption transition scenario for prosumers.
%
%
% INPUT:
%            setup:          Structure that contains all necessary settings and data for processing.
%            rootDir:        Main directory where project files are located.
%            projName:       ADD!!!!!!!!!!!!!!
%            systemParams:   Structure with all system parameters.
%            activeElements: Structure with active elements in the system
%            modFiles:       Files used for modifying the scenario setup.
%            Reg:            Region for which the scenario is run.
%            ProduserType:   Type of producer (residential, commercial).
%            costYear:       Year to which all cost values are adjusted.
%            prosShare:      Share of energy covered by prosumers.
%            el_cost:        Electricity cost.
%            feedin_cost:    Feedin tariff for exported electricity.
%            endHour:        Last hour of the simulation period.
%
% OUTPUT:
%            results:        Structure containing the scenario results.
%
%Dmitrii Bogdanov
%last change 24.07.2025


setup.ProsumersRun = 1;
try 
    setup.units = setup.prosumerUnits;
catch
    setup.units = 'kW';
end

datFiles.cost = ['cost_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.dat'];
datFiles.demand = ['demand_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.dat'];
datFiles.feedIn = ['feedin_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.dat'];
datFiles.sets = ['sets_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.dat'];
datFiles.physical = ['physical_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.dat'];
datFiles.limits = ['limits_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.dat'];
datFiles.misc = ['misc_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.dat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write .dat-Files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transmissionCapacity = '';
try
    tempFilename=[rootDir filesep 'projects' filesep projName filesep 'output' filesep 'matlab' filesep 'results_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.mat'];
    f=dir(tempFilename);
    if f.bytes<10^5
        delete(tempFilename);
    end
    load(tempFilename);
catch

    WriteCostFile([rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.cost],setup,systemParams,activeElements,costYear);
    WriteSetsProsumers([rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.sets],setup,systemParams,activeElements)
    WritePhysicalData([rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.physical],setup,systemParams,activeElements,costYear)
    WriteFeedIn([rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.feedIn],setup,systemParams,activeElements,endHour)
    WriteDemand([rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.demand],setup,systemParams,activeElements,costYear,endHour)
    WriteMiscProsumers([rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.misc],setup,systemParams,activeElements,transmissionCapacity,0,0,el_cost,feedin_cost,setup.shareOfSector,setup.shareIndHeatPros,costYear)
    limits = WriteLimits([rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep datFiles.limits],setup,systemParams,activeElements);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Generate Model file %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % paths and files
    modelFile = ['GTM' num2str(Reg) '_' num2str(costYear) '.mod'];
    modelFileLP = strrep(modelFile,'.mod','.lp');

    % concatenation of data files string
    TMPdataFilesStr = fieldnames(datFiles);
    dataFilesString = '';

    for index=1:length(TMPdataFilesStr)
        dataFilesString = [dataFilesString ' -d ' rootDir filesep 'projects' filesep projName filesep 'dat-files' filesep getfield(datFiles,char(TMPdataFilesStr(index)))];
    end

    addString = [];

    if (costYear <= setup.startYear) & ~setup.OvernightFlag

        addString = [addString sprintf('s.t. SBAT_LowerLimitation2{nBat in SBAT_LowerLimNodes}: storageCapacity[nBat,''SBAT'']/10^3 = SBAT_LowerLim[nBat]/10^3;\n')];
        addString = [addString sprintf('s.t. IBAT_LowerLimitation2{nIBat in IBAT_LowerLimNodes}: storageInterface[nIBat,''IBAT'']/10^3 = IBAT_LowerLim[nIBat]/10^3;\n')];
        addString = [addString sprintf('s.t. RPVO_LowerLimitation2{nPV in RPVO_LowerLimNodes}: instCapacity_electFeedIn[nPV,''RPVO'']/10^3 = RPVO_LowerLim[nPV]/10^3;\n')];

        addString = [addString sprintf('s.t. RRSH_LowerLimitation2{nRSH in RRSH_LowerLimNodes}: instCapacity_heatFeedIn[nRSH,''RRSH'']/10^3 = RRSH_LowerLim[nRSH]/10^3;\n')];

        addString = [addString sprintf('s.t. THHR_LowerLimitation2{nHHR in THHR_LowerLimNodes}: instCapacity_heatTransf[nHHR,''THHR'']/10^3 = THHR_LowerLim[nHHR]/10^3;\n')];
        addString = [addString sprintf('s.t. THHP_LowerLimitation2{nHHP in THHP_LowerLimNodes}: instCapacity_heatTransf[nHHP,''THHP'']/10^3 = THHP_LowerLim[nHHP]/10^3;\n')];
        addString = [addString sprintf('s.t. THNG_LowerLimitation2{nHNG in THNG_LowerLimNodes}: instCapacity_heatTransf[nHNG,''THNG'']/10^3 = THNG_LowerLim[nHNG]/10^3;\n')];

        addString = [addString sprintf('s.t. THBP_LowerLimitation2{nHBP in THBP_LowerLimNodes}: instCapacity_heatTransf[nHBP,''THBP'']/10^3 = THBP_LowerLim[nHBP]/10^3;\n')];
        addString = [addString sprintf('s.t. THBG_LowerLimitation2{nHBG in THBG_LowerLimNodes}: instCapacity_heatTransf[nHBG,''THBG'']/10^3 = THBG_LowerLim[nHBG]/10^3;\n')];

        addString = [addString sprintf('s.t. shareEl{t in T, n in N}: gen_heatTransf[t,n,''THHR'']+gen_heatTransf[t,n,''THHP''] + gen_heatFeedIn[t,n,''RRSH''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n]) * (shareIndHeatPros[n])) * shareHeatfromEl[n]+ excess_HeLo_El[t,n];\n')];
        addString = [addString sprintf('s.t. shareGa{t in T, n in N}: gen_heatTransf[t,n,''THNG''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n])) * shareHeatfromGa[n] + excess_HeLo_Ga[t,n];\n')];
        addString = [addString sprintf('s.t. shareOi{t in T, n in N}: gen_heatTransf[t,n,''THOI''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n])) * shareHeatfromOi[n]+ excess_HeLo_Oi[t,n];\n')];
        addString = [addString sprintf('s.t. shareBm{t in T, n in N}: gen_heatTransf[t,n,''THBP''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n])) * shareHeatfromBm[n]+ excess_HeLo_Bm[t,n];\n')];
        addString = [addString sprintf('s.t. shareBg{t in T, n in N}: gen_heatTransf[t,n,''THBG''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n])) * shareHeatfromBg[n]+ excess_HeLo_Bg[t,n];\n')];


    else

        addString = [addString sprintf('s.t. shareEl{t in T, n in N}: gen_heatTransf[t,n,''THHR'']+gen_heatTransf[t,n,''THHP''] + gen_heatFeedIn[t,n,''RRSH''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n]) * (shareIndHeatPros[n])) * shareHeatfromEl[n]+ excess_HeLo_El[t,n];\n')];
        addString = [addString sprintf('s.t. shareGa{t in T, n in N}: gen_heatTransf[t,n,''THNG''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n])) * shareHeatfromGa[n] + excess_HeLo_Ga[t,n];\n')];
        addString = [addString sprintf('s.t. shareOi{t in T, n in N}: gen_heatTransf[t,n,''THOI''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n])) * shareHeatfromOi[n]+ excess_HeLo_Oi[t,n];\n')];
        addString = [addString sprintf('s.t. shareBm{t in T, n in N}: gen_heatTransf[t,n,''THBP''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n])) * shareHeatfromBm[n]+ excess_HeLo_Bm[t,n];\n')];
        addString = [addString sprintf('s.t. shareBg{t in T, n in N}: gen_heatTransf[t,n,''THBG''] = ((demand_LHSP[t,n]+demand_LHDW[t,n])*(1-shareOfDistrHeat[n]) * (shareOfSector[n])) * shareHeatfromBg[n]+ excess_HeLo_Bg[t,n];\n')];


    end



    if setup.BlockFossil.Flag
        if costYear >= setup.BlockFossil.Year
            addString = [addString sprintf('s.t. fossilLim1BF{t in T,n in N}: inputFuel[t,n,''RHAR''] = 0;\n')];
            addString = [addString sprintf('s.t. fossilLim2BF{t in T,n in N}: inputFuel[t,n,''RPET''] = 0;\n')];
        end
    end


    if costYear == setup.endYear

        addString = [addString sprintf('s.t. RRSH_LowerLimitation2{nRSH in RRSH_LowerLimNodes}: instCapacity_heatFeedIn[nRSH,''RRSH'']/10^3 <= 1.01*RRSH_LowerLim[nRSH]/10^3;\n')];

        addString = [addString sprintf('s.t. THBG_LowerLimitation2{nHBG in THBG_LowerLimNodes}: instCapacity_heatTransf[nHBG,''THBG'']/10^3 <= 1.01*THBG_LowerLim[nHBG]/10^3;\n')];
        addString = [addString sprintf('s.t. THBP_LowerLimitation2{nHBP in THBP_LowerLimNodes}: instCapacity_heatTransf[nHBP,''THBP'']/10^3 <= 1.01*THBP_LowerLim[nHBP]/10^3;\n')];
        addString = [addString sprintf('s.t. THNG_LowerLimitation2{nHNG in THNG_LowerLimNodes}: instCapacity_heatTransf[nHNG,''THNG'']/10^3 <= 1.01*THNG_LowerLim[nHNG]/10^3;\n')];
        addString = [addString sprintf('s.t. THOI_LowerLimitation2{nHOI in THOI_LowerLimNodes}: instCapacity_heatTransf[nHOI,''THOI'']/10^3 <= 1.01*THOI_LowerLim[nHOI]/10^3;\n')];

    end

    try
        setup.ProsumersRun = 1;
        eval(setup.additionalConstraints);
    catch
        warning('nNo additional constraints added');
    end
    setup.ProsumersRun = 0;
    % concatenation of model string
    modelString = [sprintf('param hourEnd := ') num2str(endHour) sprintf(';\n') ...
        fileread([rootDir  filesep 'projects' filesep projName filesep 'models' filesep modFiles.paramsAndSets]) ...
        fileread([rootDir filesep 'projects' filesep projName filesep 'models' filesep 'varsMain.mod']) ...
        fileread([rootDir filesep 'projects' filesep projName filesep 'models' filesep modFiles.objAndCons]) ...
        addString, ...
        fileread([rootDir filesep 'projects' filesep projName filesep 'models' filesep modFiles.end])];


    % write concatenated string to model file
    fModel = fopen([rootDir filesep 'projects' filesep projName filesep 'input-data' filesep modelFile],'w');
    fprintf(fModel,'%s',modelString);
    fclose(fModel);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Process model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    glpsolArgs = [' --check --math ' rootDir filesep 'projects' filesep projName filesep 'input-data' filesep modelFile ' --wlp ' rootDir filesep 'projects' filesep projName filesep modelFileLP ' '];
	GLPKlog = [];
    if ismac
		% Code to run on Mac platform
		[~,GLPKlog] = system(['"' setup.GLPKPath '"' glpsolArgs dataFilesString]); % used to run a zsh command
	elseif isunix
		% Code to run on Linux platform
		disp('Linux platform is not supported yet')
		[~,GLPKlog] = system(['"' setup.GLPKPath '"' glpsolArgs dataFilesString]); % used to run a zsh command
	elseif ispc
		% Code to run on Windows platform
		[~,GLPKlog] = system(['"' setup.GLPKPath '"' glpsolArgs dataFilesString]); % ! = used for run a bash command
	else
		disp('Platform not supported')
    end

    if isempty(GLPKlog) || contains(GLPKlog(max(0,end-5):end),'error')
	    display(GLPKlog(max(0,end-500):end));
        error('Problem with compiling the model and writing to LP file')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Solve problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % default is usage of MATLAB MOSEK toolbox

    % read problem to pStruct.prob
    [rDummy,pStruct] = mosekopt(['read(' rootDir filesep 'projects' filesep projName filesep modelFileLP ')']);

    % setting parameters for MOSEK solver
    param.MSK_IPAR_NUM_THREADS = setup.SolverNumThreads; % specifies number of threads used by interior-point solver

    param.MSK_IPAR_INTPNT_BASIS  = 'MSK_BI_NEVER'; % solely using interior-point solver

    param.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON'; % switch infeasibility report to ON
    param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 6000;
    param.MSK_DPAR_INTPNT_TOL_PFEAS = 1.0e-8;
    param.MSK_DPAR_INTPNT_TOL_DFEAS = 1.0e-8;
    param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0e-9;
    param.MSK_IPAR_INFEAS_REPORT_LEVEL = 1000;

    % settings for writing log to text file
    fid                = fopen([rootDir filesep 'projects' filesep projName filesep 'output' filesep 'log_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.log'],'w');
    callback.log       = 'LogPrint';
    callback.loghandle = fid;

    %  solve problem via mosek solver call
    [r, res] = mosekopt('minimize info',pStruct.prob,param,callback);

    fclose(fid);

    display('importing results to MATLAB...')
    results=PrepareResults(setup,pStruct.prob.names.var,res.sol.itr.xx,res.sol.itr.pobjval);

    results.SC_share = prosShare;
    results.SC_demand = systemParams.ValueLoad(find(ismember(systemParams.IndexIDL,'LELE')),:,1,1)';
    results.shareIndHeatElPros = setup.shareIndHeatPros;
    %results.sectorShare =
    display('[done] result import.')

    save([rootDir filesep 'projects' filesep projName filesep 'output' filesep 'results_' num2str(Reg) '_' ProduserType '_' num2str(costYear) '.mat'],'results')
	delete([rootDir filesep 'projects' filesep projName filesep modelFileLP])
end