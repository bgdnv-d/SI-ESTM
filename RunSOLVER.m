function results = RunSOLVER(solver,pName,setup,rootDir,costYear,SupportFossilsFlag,name,varargin)
%
% FUNCTION runSOLVER_TR_regNew(solver, pName, setup, rootDir, costYear, SupportFossilsFlag, name, varargin)
%
% Runs the transition region solver with given parameters and options.
%
%
% INPUT:
%            solver:                Solver function or identifier to be used.
%            pName:                 Project name.
%            setup:                 Structure that contains all necessary settings and data for processing.
%            rootDir:               Main directory path for project files.
%            costYear:              Year to which costs are adjusted.
%            SupportFossilsFlag:    Flag indicating whether fossil fuel support is enabled (true/false).
%            name:                  Name identifier for the run or configuration.
%            varargin:              
%               endHour - hourly profiles length (default 8760)
%               solveTool - solver call approach, matlab interface call ('mtlb') or system call ('sysc'), (default 'mtlb')
%               numThreads - number of threads used by solver in barrier or interior point methods (default 2)
%
% OUTPUT:
%            results:               Structure containing the solver output results.
%
%Dmitrii Bogdanov
%last change 24.07.2025


switch nargin
    case 7
        endHour = 8760;
        solveTool = 'mtlb';
        numThreads = 2;

    case 8
        endHour = varargin{1};
        solveTool = 'mtlb';
        numThreads = 2;

    case 9
        endHour = varargin{1};
        solveTool = varargin{2};
        numThreads = 2;

    case 10
        endHour = varargin{1};
        solveTool = varargin{2};
        numThreads = varargin{3};

end


try
    setup.SolverDemLim;
catch
    setup.SolverDemLim=0;
end

try
    setup.BlockFossil.Limit;
catch
    setup.BlockFossil.Limit = 0;
end

try 
    setup.flexStart;
catch
    try
        setup.flexStart = setup.flex2015;
    catch
        setup.flexStart = 0;
    end
end

if isempty(name)

    resFilename=[rootDir filesep 'projects' filesep pName filesep 'output' filesep 'results.mat'];

else
    resFilename=[rootDir filesep 'projects' filesep pName filesep 'output' filesep 'results_' name '.mat'];

end

s= (dir(resFilename));
if isfile(resFilename) & s.bytes>1000000
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load(resFilename);
    disp('Solver is not run, existing results are used');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(name)
        load([rootDir filesep 'projects' filesep pName filesep 'input-data' filesep 'simulation-input_' pName]);
    else
        load([rootDir filesep 'projects' filesep pName filesep 'input-data' filesep 'simulation-input_' pName '_' name ]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    addString = [];
    addParam = [];



    if costYear==setup.startYear

        if (setup.SupportFossilsFlag)&~(setup.OvernightFlag)



            addString = [addString sprintf('\ns.t. RPVO_UpperLimitation2{nPV in RPVO_LowerLimNodes}: instCapacity_electFeedIn[nPV,''RPVO'']/10^3 = (RPVO_LowerLim[nPV])/10^3;\n')];
            addString = [addString sprintf('s.t. RPVA_UpperLimitation2{nPVA in RPVA_LowerLimNodes}: instCapacity_electFeedIn[nPVA,''RPVA'']/10^3 = (RPVA_LowerLim[nPVA])/10^3;\n')];
            addString = [addString sprintf('s.t. RPVF_UpperLimitation2{nPV in RPVF_LowerLimNodes}: instCapacity_electFeedIn[nPV,''RPVF'']/10^3 = (RPVF_LowerLim[nPV])/10^3;\n')];
            addString = [addString sprintf('s.t. RPBV_UpperLimitation2{nPV in RPBV_LowerLimNodes}: instCapacity_electFeedIn[nPV,''RPBV'']/10^3 = (RPBV_LowerLim[nPV])/10^3;\n')];
            addString = [addString sprintf('s.t. RPBO_UpperLimitation2{nPV in RPBO_LowerLimNodes}: instCapacity_electFeedIn[nPV,''RPBO'']/10^3 = (RPBO_LowerLim[nPV])/10^3;\n')];
            addString = [addString sprintf('s.t. RPBA_UpperLimitation2{nPV in RPBA_LowerLimNodes}: instCapacity_electFeedIn[nPV,''RPBA'']/10^3 = (RPBA_LowerLim[nPV])/10^3;\n')];

            addString = [addString sprintf('s.t. RRRI_UpperLimitation2{nRRI in RRRI_LowerLimNodes}: instCapacity_electFeedIn[nRRI,''RRRI'']/10^3 = RRRI_LowerLim[nRRI]/10^3;\n')];
            addString = [addString sprintf('s.t. HDAM_UpperLimitation2{nDR in HDAM_LowerLimNodes}: instCapacity_Hydro[nDR,''HDAM'']/10^3 = HDAM_LowerLim[nDR]/10^3;\n')];
            addString = [addString sprintf('s.t. RWIN_UpperLimitation2{nRWIN in RWIN_LowerLimNodes}: instCapacity_electFeedIn[nRWIN,''RWIN'']/10^3 = RWIN_LowerLim[nRWIN]/10^3;\n')];
            addString = [addString sprintf('s.t. RWIO_UpperLimitation2{nRWIO in RWIO_LowerLimNodes}: instCapacity_electFeedIn[nRWIO,''RWIO'']/10^3 = RWIO_LowerLim[nRWIO]/10^3;\n')];
            addString = [addString sprintf('s.t. ROWI_UpperLimitation2{nROWI in ROWI_LowerLimNodes}: instCapacity_electFeedIn[nROWI,''ROWI'']/10^3 = ROWI_LowerLim[nROWI]/10^3;\n')];
            addString = [addString sprintf('s.t. RCSP_UpperLimitation2{nRCSP in RCSP_LowerLimNodes}: instCapacity_heatFeedIn[nRCSP,''RCSP'']/10^3 = RCSP_LowerLim[nRCSP]/10^3;\n')];
            addString = [addString sprintf('s.t. RDSH_UpperLimitation2{nRDSH in RDSH_LowerLimNodes}: instCapacity_heatFeedIn[nRDSH,''RDSH'']/10^3 = RDSH_LowerLim[nRDSH]/10^3;\n')];

            addString = [addString sprintf('s.t. IPHS_LowerLimitation2{nIPS in IPHS_LowerLimNodes}: storageInterface[nIPS,''IPHS''] /10^3 = IPHS_LowerLim[nIPS]/10^3;\n')];


            if ~setup.Mobility
                addString = [addString sprintf('\ns.t. TWEL_UpperLimitation2{n in N}: instCapacity_ptgTransf[n,''TWEL'']/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. TMET_UpperLimitation2{n in N}: instCapacity_ptgTransf[n,''TMET'']/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. TCOS_UpperLimitation2{n in N}: instCapacity_ptgTransf[n,''TCOS'']/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. TFTU_UpperLimitation2{n in N}: instCapacity_ptgTransf[n,''TFTU'']/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. TBGU_UpperLimitation2{n in N}: instCapacity_ptgTransf[n,''TBGU'']/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. FTCon0{t in T, n in N}: gen_ptgTransf[t,n,''TFTU''] = 0;\n')];
                addString = [addString sprintf('s.t. FTCon1{t in T, n in N}: gen_ptgTransf[t,n,''TWEL''] = 0;\n')];
                addString = [addString sprintf('s.t. FTCon3{t in T, n in N}: gen_ptgTransf[t,n,''TCOS''] = 0;\n')];
            end

            if setup.flexGT
                addString = [addString sprintf('s.t. TICM_UpperLimitation2{n in N}: instCapacity_electTransf[n,''TICM'']/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. TGCS_UpperLimitation2{n in N}: instCapacity_electTransf[n,''TGCS'']/10^3 = 0/10^3;\n')];
            end

            %% allows extra OCG and ICG capacity to be installed


            if ~setup.flexStart
                addString = [addString sprintf('s.t. TOCG_UpperLimitation2{nOCG in TOCG_LowerLimNodes}: instCapacity_electTransf[nOCG,''TOCG'']/10^3 = TOCG_LowerLim[nOCG]/10^3;\n')];
                addString = [addString sprintf('s.t. TICG_UpperLimitation2{nICG in TICG_LowerLimNodes}: instCapacity_electTransf[nICG,''TICG'']/10^3 = TICG_LowerLim[nICG]/10^3;\n')];
            end


            addString = [addString sprintf('s.t. existingTransmissionCapacityDC2{l in L}: instCapacitiyTransmLines_DC[l]/10^3 = TransmissionLines_DC_LowerLim[l]/10^3;\n')];
            addString = [addString sprintf('s.t. existingTransmissionCapacityAC2{l in L}: instCapacitiyTransmLines_AC[l]/10^3 = TransmissionLines_AC_LowerLim[l]/10^3;\n')];


        end

    else



    end



    if setup.SolverDemLim == 1
        if setup.Heat.Flag
            addString = [addString sprintf('\ns.t. elProdLim{n in N}: 3 * sum{t in T} (demand_LELE[t,n] + demand_LHSP[t,n] + demand_LHDW[t,n] + 2*demand_LIGA[t,n] + 2*demandProsGas[t,n])/10^6 >= sum{t in T} (sum{eT in electTransf} gen_electTransf[t,n,eT] + sum{cT in chpTransf} gen_chpElTransf[t,n,cT] + sum{eF in electFeedIn} gen_electFeedIn[t,n,eF] + sum{H in Hydro} gen_Hydro[t,n,H])/10^6;\n')];
        else
            addString = [addString sprintf('\ns.t. elProdLim{n in N}: 3 * sum{t in T} (demand_LELE[t,n])/10^6 >= sum{t in T} (sum{eT in electTransf} gen_electTransf[t,n,eT] + sum{cT in chpTransf} gen_chpElTransf[t,n,cT] + sum{eF in electFeedIn} gen_electFeedIn[t,n,eF] + sum{H in Hydro} gen_Hydro[t,n,H])/10^6;\n')];
        end

    end

    if costYear >= setup.BlockFossil.Year
        try
            if setup.BlockFossil.Flag

                addString = [addString sprintf(['\ns.t. fossilLim1{t in T,n in N}: inputFuel[t,n,''RNGA'']/1000 <= ' num2str(setup.BlockFossil.Limit) '/1000;\n'])];
                addString = [addString sprintf(['\ns.t. fossilLim2{t in T,n in N}: inputFuel[t,n,''RHAR'']/1000 <= ' num2str(setup.BlockFossil.Limit) '/1000;\n'])];
                addString = [addString sprintf(['\ns.t. fossilLim3{t in T,n in N}: inputFuel[t,n,''RPET'']/1000 <= ' num2str(setup.BlockFossil.Limit) '/1000;\n'])];
            end
        catch
            warning('No BlockFossil Data - add setup.BlockFossil.Flag and setup.BlockFossil.Value')
        end

    end

    if systemParams.GridMax==0
        addString = [addString sprintf('\ns.t. gridMaxConstrainUp{t in T, n in N}: transPower_AC[t,n]+transPower_DC[t,n] <= 0;\n')];
        addString = [addString sprintf('s.t. gridMaxConstrainDown{t in T, n in N}: transPower_AC[t,n]+transPower_DC[t,n] >= 0;\n')];
    end

    if setup.Heat.Flag

        if (costYear~=2015)&(costYear <= 2040)
            if ~setup.Mobility
                try
                    setup.SolverPeakLoadLim;
                catch
                    setup.SolverPeakLoadLim=0;
                end
                if setup.SolverPeakLoadLim
                    addString = [addString sprintf('\ns.t. peakLoadLim{n in N}: (instCapacity_chpTransf[n,''TMSW'']+instCapacity_chpTransf[n,''TCHP'']+instCapacity_chpTransf[n,''TCBP'']+instCapacity_chpTransf[n,''TCNG'']+instCapacity_chpTransf[n,''TCOI'']+instCapacity_chpTransf[n,''TCCO'']+instCapacity_electTransf[n,''TGEO'']+instCapacity_electTransf[n,''TBPP'']+instCapacity_electTransf[n,''THPP'']+instCapacity_electTransf[n,''TCCG'']+instCapacity_electTransf[n,''TOCG'']+instCapacity_electTransf[n,''TICG'']+instCapacity_electTransf[n,''TNUC'']+instCapacity_electTransf[n,''TSTU'']) <= 5*sum{t in T}(demand_LELE[t,n])/hourEnd;\n')];
                end
            end
        end

        addString = [addString sprintf('\n# biomass cooking\n')];
        addString = [addString sprintf('\ns.t. BCHDemandSatis{t in T, n in N}: (BGAToIndCook[t,n] + BioToCook[t,n] + MBBioToCook[t,n])/10^3 = (demand_LHBC[t,n])/10^3;\n')];

        addString = [addString sprintf('\n# centralised system\n')];
        addString = [addString sprintf('\ns.t. CentrHeatSatis{t in T, n in N}: (demand_LHIN[t,n] + (demand_LHSP[t,n]+demand_LHDW[t,n])*shareOfDistrHeat[n] + MidHetoInd[t,n] + HighHetoInd[t,n] + co2Scrubbing_He[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n] + desalination_Heat_demand[t,n] + heatToPower[t,n])/10^3\n')];
        addString = [addString sprintf('\n= DistrHeatEff[n] * (sum{chT in heatDTransf} gen_heatTransf[t,n,chT] + sum{chF in heatDFeedIn} gen_heatFeedIn[t,n,chF] + sum{cT in chpTransf} gen_chpHeTransf[t,n,cT] + dischargePower[t,n,''SSHS''] + dischargePower[t,n,''SDHS''] + dischargePower[t,n,''SHOT''] + LossToHeat[t,n] - chargePower[t,n,''SSHS''] - chargePower[t,n,''SDHS''] - chargePower[t,n,''SHOT''] - excess_He[t,n])/10^3+\n')];
        addString = [addString sprintf('(0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n]))/10^3;\n')];

        addString = [addString sprintf('\n# centralised system hight heat\n')];

        %High heat directly with fuels
        addString = [addString sprintf('\ns.t. CentrHighHeatSatis{t in T, n in N}: (demand_LHIN[t,n]*shareOfHighHeatInd[n]+ HighHetoInd[t,n])/10^3 <= 0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n] )/10^3;\n')];



        %end


    else

        addString = [addString sprintf('\n# biomass cooking\n')];
        addString = [addString sprintf('\ns.t. BCHDemandSatis{t in T, n in N}: (BGAToIndCook[t,n] + BioToCook[t,n] + MBBioToCook[t,n])/10^3 = 0;\n')];


        addString = [addString sprintf('\n# heat cycle \n')];

        addString = [addString sprintf('\ns.t. elHeatProd1{t in T, n in N}: gen_heatTransf[t,n,''TDHR'']  = HighHeat[t,n];\n')];
        addString = [addString sprintf('\ns.t. elHeatProd2{t in T, n in N}: 0 = LowHeat[t,n];\n')];
        addString = [addString sprintf('\ns.t. CentrLowHeatSatis{t in T, n in N}: LossToHeat[t,n] + LowHeat[t,n] + gen_heatTransf[t,n,''TDHP''] + dischargePower[t,n,''SSHS''] + dischargePower[t,n,''SDHS''] - chargePower[t,n,''SSHS''] - chargePower[t,n,''SDHS''] = co2Scrubbing_He[t,n] + excess_He[t,n];\n')];
        addString = [addString sprintf('\ns.t. CentrHighHeatSatis{t in T, n in N}: HighHeat[t,n] + gen_heatFeedIn[t,n,''RCSP''] + dischargePower[t,n,''SHOT''] - chargePower[t,n,''SHOT''] = heatToPower[t,n];\n')];
        addString = [addString sprintf('\ns.t. ExHeatLim{t in T, n in N}: LossToHeat[t,n] >= excess_He[t,n];\n')];
        addString = [addString sprintf('\ns.t. ExHeatLim2{t in T, n in N}: 0 = excess_HeLo[t,n];\n')];


    end
    %% Scenario specific
    try
        setup.Standard;
    catch
        setup.Standard=1;
    end

    if setup.Standard == 0
        if isfield(setup,'additionalConstraints')
            eval(setup.additionalConstraints);
        else
            warning('Please add additional constraints configuration, setup switches to standard');
            setup.Standard = 1;
        end
    end

    if setup.Standard

        if costYear == setup.startYear
            addString = [addString sprintf('\ns.t. electricalTransformersFLH_NUC{n in N}: (sum{t in T} gen_electTransf[t,n,''TNUC''])/10^3= hourEnd*0.7*instCapacity_electTransf[n,''TNUC'']/10^3;\n')];
        else
            addString = [addString sprintf('\ns.t. electricalTransformersFLH_NUC{n in N}: (sum{t in T} gen_electTransf[t,n,''TNUC''])/10^3= hourEnd*0.7*instCapacity_electTransf[n,''TNUC'']/10^3;\n')];
        end

        if setup.Mobility

            if costYear > 2015


                if costYear == 2030
                    addString = [addString sprintf('s.t. DieselForMobility3{n in N}: sum{t in T} (diesToMobility[t,n])*0.026/10^6 <= sum{t in T}(FTdiesel[t,n] + BIOdiesel[t,n] + importFT_diesel[t,n])/10^6;\n')];
                    addString = [addString sprintf('s.t. KeroseneForMobility3{n in N}: sum{t in T} (kersToMobility[t,n])*0.026/10^6 <= sum{t in T}(FTkerosene[t,n] + BIOkerosene[t,n] + importFT_kerosene[t,n])/10^6;\n')];
                end
                if costYear == 2035
                    addString = [addString sprintf('s.t. DieselForMobility3{n in N}: sum{t in T} (diesToMobility[t,n])*0.1/10^6 <= sum{t in T}(FTdiesel[t,n] + BIOdiesel[t,n]+ importFT_diesel[t,n])/10^6;\n')];
                    addString = [addString sprintf('s.t. KeroseneForMobility3{n in N}: sum{t in T} (kersToMobility[t,n])*0.1/10^6 <= sum{t in T}(FTkerosene[t,n] + BIOkerosene[t,n] + importFT_kerosene[t,n])/10^6;\n')];
                end
                if costYear == 2040
                    addString = [addString sprintf('s.t. DieselForMobility3{n in N}: sum{t in T} (diesToMobility[t,n])*0.43/10^6 <= sum{t in T}(FTdiesel[t,n] + BIOdiesel[t,n]+ importFT_diesel[t,n])/10^6;\n')];
                    addString = [addString sprintf('s.t. KeroseneForMobility3{n in N}: sum{t in T} (kersToMobility[t,n])*0.43/10^6 <= sum{t in T}(FTkerosene[t,n] + BIOkerosene[t,n] + importFT_kerosene[t,n])/10^6;\n')];
                end
                if costYear == 2045
                    addString = [addString sprintf('s.t. DieselForMobility3{n in N}: sum{t in T} (diesToMobility[t,n])*0.76/10^6 <= sum{t in T}(FTdiesel[t,n] + BIOdiesel[t,n]+ importFT_diesel[t,n])/10^6;\n')];
                    addString = [addString sprintf('s.t. KeroseneForMobility3{n in N}: sum{t in T} (kersToMobility[t,n])*0.76/10^6 <= sum{t in T}(FTkerosene[t,n] + BIOkerosene[t,n] + importFT_kerosene[t,n])/10^6;\n')];
                end
                if costYear == 2050
                end
            end
            if costYear <=2050
                if costYear == setup.startYear
                    addString = [addString sprintf('s.t. FuelFlowRBDS_STR{n in N}: ((sum{t in T} inputFuel[t,n,''RBDS'']))/10^6 <= (hourEnd/8760*(totalFuelLimit[n,''RBDS'']))/10^6;\n')];
                else
                    addString = [addString sprintf('s.t. FuelFlowRBDS_STR{n in N}: ((sum{t in T} inputFuel[t,n,''RBDS'']))/10^6 = (hourEnd/8760*(totalFuelLimit[n,''RBDS'']))/10^6;\n')];

                end
            end
        end
        %% Hydrogen use for power and heat

        if setup.flexGT
            if costYear <2030
                addString = [addString sprintf('s.t. h2Use1{t in T,n in N}: h2toCCGT[t,n]/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. h2Use2{t in T,n in N}: h2toOCGT[t,n]/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. h2Use3{t in T,n in N}: h2toTGCS[t,n]/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. h2Use4{t in T,n in N}: h2toTCNG[t,n]/10^3 = 0/10^3;\n')];
                addString = [addString sprintf('s.t. h2Use5{t in T,n in N}: h2toTICM[t,n]/10^3 = 0/10^3;\n')];

                addString = [addString sprintf('s.t. TGCS_UpperLimitationPre2030{n in N}: instCapacity_electTransf[n,''TGCS'']/10^3 = 0/10^3;\n')];

            elseif costYear ==2030
                addString = [addString sprintf('s.t. h2Use1{t in T,n in N}: h2toCCGT[t,n]/10 <= 0.3*(importGasToCCGT[t,n]+SNGGasToCCGT[t,n]+OilToCCGT[t,n]+h2toCCGT[t,n])/10;\n')];
                addString = [addString sprintf('s.t. h2Use2{t in T,n in N}: h2toOCGT[t,n]/10 <= 0.3*(importGasToOCGT[t,n]+SNGGasToOCGT[t,n]+OilToOCGT[t,n]+h2toOCGT[t,n])/10;\n')];
                addString = [addString sprintf('s.t. h2Use3{t in T,n in N}: h2toTGCS[t,n]/10 <= 0.3*(importGasToTGCS[t,n]+SNGGasToTGCS[t,n]+OilToTGCS[t,n]+h2toTGCS[t,n])/10;\n')];
                addString = [addString sprintf('s.t. h2Use4{t in T,n in N}: h2toTCNG[t,n]/10 <= 0.3*(importGasToCHP[t,n]+SNGGasToCHP[t,n]+OilToTCNG[t,n]+h2toTCNG[t,n])/10;\n')];
                addString = [addString sprintf('s.t. h2Use5{t in T,n in N}: h2toTICM[t,n]/10 <= 0.3*(importGasToTICM[t,n]+SNGGasToTICM[t,n]+OilToTICM[t,n]+h2toTICM[t,n])/10;\n')];
            elseif costYear ==2035
                addString = [addString sprintf('s.t. h2Use1{t in T,n in N}: h2toCCGT[t,n]/10 <= 0.7*(importGasToCCGT[t,n]+SNGGasToCCGT[t,n]+OilToCCGT[t,n]+h2toCCGT[t,n])/10;\n')];
                addString = [addString sprintf('s.t. h2Use2{t in T,n in N}: h2toOCGT[t,n]/10 <= 0.7*(importGasToOCGT[t,n]+SNGGasToOCGT[t,n]+OilToOCGT[t,n]+h2toOCGT[t,n])/10;\n')];
                addString = [addString sprintf('s.t. h2Use3{t in T,n in N}: h2toTGCS[t,n]/10 <= 0.7*(importGasToTGCS[t,n]+SNGGasToTGCS[t,n]+OilToTGCS[t,n]+h2toTGCS[t,n])/10;\n')];
                addString = [addString sprintf('s.t. h2Use4{t in T,n in N}: h2toTCNG[t,n]/10 <= 0.7*(importGasToCHP[t,n]+SNGGasToCHP[t,n]+OilToTCNG[t,n]+h2toTCNG[t,n])/10;\n')];
                addString = [addString sprintf('s.t. h2Use5{t in T,n in N}: h2toTICM[t,n]/10 <= 0.7*(importGasToTICM[t,n]+SNGGasToTICM[t,n]+OilToTICM[t,n]+h2toTICM[t,n])/10;\n')];
            end
        else
            addString = [addString sprintf('s.t. h2Use1{t in T,n in N}: h2toCCGT[t,n]/10^3 = 0/10^3;\n')];
            addString = [addString sprintf('s.t. h2Use2{t in T,n in N}: h2toOCGT[t,n]/10^3 = 0/10^3;\n')];
            addString = [addString sprintf('s.t. h2Use3{t in T,n in N}: h2toTGCS[t,n]/10^3 = 0/10^3;\n')];
            addString = [addString sprintf('s.t. h2Use4{t in T,n in N}: h2toTCNG[t,n]/10^3 = 0/10^3;\n')];
            addString = [addString sprintf('s.t. h2Use5{t in T,n in N}: h2toTICM[t,n]/10^3 = 0/10^3;\n')];
        end

        if setup.Heat.Flag
            addString = [addString sprintf('\n# centralised system hight and medium heat\n')];
            %High heat directly with fuels with spilover to medium heat

            switch costYear
                case 2020

                case 2025
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n] + HighHetoInd[t,n] + heatToPower[t,n])/10^3 <= DistrHeatEff[n] * (gen_heatFeedIn[t,n,''RCSP''] + gen_heatTransf[t,n,''TDHR''] + dischargePower[t,n,''SHOT''] - chargePower[t,n,''SHOT''] )/10^3+0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis2{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]+ HighHetoInd[t,n])/10^3 >= 0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];

                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis3{n in N}: 0.20*(sum{t in T} (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n]-shareOfHighHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]))/10^6 <= DistrHeatEff[n] * (sum{t in T} gen_heatTransf[t,n,''TDHR''])/10^6;\n')];
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis3u{n in N}: 0.25*(sum{t in T} (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n]-shareOfHighHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]))/10^6 >= DistrHeatEff[n] * (sum{t in T} gen_heatTransf[t,n,''TDHR''])/10^6;\n')];
                case 2030
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n] + HighHetoInd[t,n] + heatToPower[t,n])/10^3 <= DistrHeatEff[n] * (gen_heatFeedIn[t,n,''RCSP''] + gen_heatTransf[t,n,''TDHR''] + dischargePower[t,n,''SHOT''] - chargePower[t,n,''SHOT''] )/10^3+0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis2{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]+ HighHetoInd[t,n])/10^3 >= 0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];

                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis3{n in N}: 0.35*(sum{t in T} (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n]-shareOfHighHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]))/10^6 <= DistrHeatEff[n] * (sum{t in T} gen_heatTransf[t,n,''TDHR''])/10^6;\n')];
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis3u{n in N}: 0.40*(sum{t in T} (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n]-shareOfHighHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]))/10^6 >= DistrHeatEff[n] * (sum{t in T} gen_heatTransf[t,n,''TDHR''])/10^6;\n')];
                case 2035
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n] + HighHetoInd[t,n] + heatToPower[t,n])/10^3 <= DistrHeatEff[n] * (gen_heatFeedIn[t,n,''RCSP''] + gen_heatTransf[t,n,''TDHR''] + dischargePower[t,n,''SHOT''] - chargePower[t,n,''SHOT''] )/10^3+0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis2{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]+ HighHetoInd[t,n])/10^3 >= 0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];

                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis3{n in N}: 0.50*(sum{t in T} (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n]-shareOfHighHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]))/10^6 <= DistrHeatEff[n] * (sum{t in T} gen_heatTransf[t,n,''TDHR''])/10^6;\n')];
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis3u{n in N}: 0.60*(sum{t in T} (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n]-shareOfHighHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]))/10^6 >= DistrHeatEff[n] * (sum{t in T} gen_heatTransf[t,n,''TDHR''])/10^6;\n')];
                case 2040
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n] + HighHetoInd[t,n] + heatToPower[t,n])/10^3 <= DistrHeatEff[n] * (gen_heatFeedIn[t,n,''RCSP''] + gen_heatTransf[t,n,''TDHR''] + dischargePower[t,n,''SHOT''] - chargePower[t,n,''SHOT''] )/10^3+0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis2{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]+ HighHetoInd[t,n])/10^3 >= 0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];

                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis3{n in N}: 0.5*(sum{t in T} (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n]-shareOfHighHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]))/10^6 <= DistrHeatEff[n] * (sum{t in T} gen_heatTransf[t,n,''TDHR''])/10^6;\n')];
                otherwise
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n] + HighHetoInd[t,n] + heatToPower[t,n])/10^3 <= DistrHeatEff[n] * (gen_heatFeedIn[t,n,''RCSP''] + gen_heatTransf[t,n,''TDHR''] + dischargePower[t,n,''SHOT''] - chargePower[t,n,''SHOT''] )/10^3+0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];
                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis2{t in T, n in N}: (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]+ HighHetoInd[t,n])/10^3 >= 0.9*(importGasToIndHeat[t,n] + SNGGasToIndHeat[t,n] + OilToIndHeat[t,n] + CoalToIndHeat[t,n] + h2toIndHeat[t,n] +BioToIndHeat[t,n]+MBBioToIndHeat[t,n]+ BGAToIndHeat[t,n])/10^3;\n')];

                    addString = [addString sprintf('\ns.t. CentrMeanHighHeatSatis3{n in N}: 0.5*(sum{t in T} (demand_LHIN[t,n]*(1-shareOfLowHeatInd[n]-shareOfHighHeatInd[n])+ MidHetoInd[t,n] + He_to_TPSC[t,n] + He_to_TPSP[t,n] + He_to_TRME[t,n] + He_to_TRSi[t,n] + He_to_TRGE[t,n]))/10^6 <= DistrHeatEff[n] * (sum{t in T} gen_heatTransf[t,n,''TDHR''])/10^6;\n')];

            end
        end
    end

    %% Curtailment limitations
    try 
        setup.ExcessLim;
    catch
        setup.ExcessLim = 0;
    end

    if setup.ExcessLim;
        addString = [addString sprintf('s.t. elelExcessLim{n in N}: 0.055 * sum{t in T} (demand_LELE[t,n] + TotalDesalinationElDemand[t,n] + powerToHeat[t,n] +electrolysis[t,n] + co2Scrubbing_El[t,n]+EltoLH2[t,n] + EltoLNG[t,n] + EltoMobility[t,n]+EltoInd[t,n]+ EltoMET[t,n]+ EltoMEO[t,n]+ EltoNH3[t,n] + EltoHyS[t,n])/10^6 >= sum{t in T} (excess_El[t,n])/10^6;\n')];
    end

    %% Peak load
    try 
        setup.SolverPeakLoadLim;
    catch
        setup.SolverPeakLoadLim = 0;
    end
    if setup.SolverPeakLoadLim
        addString = [addString sprintf('\ns.t. peakLoadLim2{n in N}: (instCapacity_chpTransf[n,''TMSW'']+instCapacity_chpTransf[n,''TCHP'']+instCapacity_chpTransf[n,''TCBP'']+instCapacity_chpTransf[n,''TCNG'']+instCapacity_chpTransf[n,''TCOI'']+instCapacity_chpTransf[n,''TCCO'']+instCapacity_electTransf[n,''TGEO'']+instCapacity_electTransf[n,''TBPP'']+instCapacity_electTransf[n,''THPP'']+instCapacity_electTransf[n,''TCCG'']+instCapacity_electTransf[n,''TOCG'']+instCapacity_electTransf[n,''TICG'']+instCapacity_electTransf[n,''TICM'']+instCapacity_electTransf[n,''TNUC'']+instCapacity_electTransf[n,''TSTU'']) <= 2.5*sum{t in T}(demand_LELE[t,n])/hourEnd;\n')];
    end

    %% Support technology, reinvestments
    try 
        setup.fixTech.Flag;
    catch
        setup.fixTech.Flag = 0;
    end

    if (setup.fixTech.Flag)&&(costYear>setup.startYear)

        for tt = 1:length(setup.fixTech.Tech)

            if sum(ismember(activeElements.labels.elFeedIn,setup.fixTech.Tech{tt}))

                if strcmp(setup.fixTech.Tech{tt},'RWIN')
                    addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: (instCapacity_electFeedIn[n,''' setup.fixTech.Tech{tt} ''']+instCapacity_electFeedIn[n,''RWIO''])/10^3 >= (minimum_' setup.fixTech.Tech{tt} '[n]+minimum_RWIO[n])/10^3;\n'])];
                    addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
                    addParam = [addParam sprintf('param minimum_RWIO{n in N};\n')];
                else
                    addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: instCapacity_electFeedIn[n,''' setup.fixTech.Tech{tt} ''']/10^3 >= minimum_' setup.fixTech.Tech{tt} '[n]/10^3;\n'])];
                    addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
                end
            end

            if sum(ismember(activeElements.labels.elTransformer,setup.fixTech.Tech{tt}))
                addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: instCapacity_electTransf[n,''' setup.fixTech.Tech{tt} ''']/10^3 >= minimum_' setup.fixTech.Tech{tt} '[n]/10^3;\n'])];
                addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
            end

            if sum(ismember(activeElements.labels.heatFeedIn,setup.fixTech.Tech{tt}))
                addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: instCapacity_heatFeedIn[n,''' setup.fixTech.Tech{tt} ''']/10^3 >= minimum_' setup.fixTech.Tech{tt} '[n]/10^3;\n'])];
                addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
            end

            if sum(ismember(activeElements.labels.distrHeatTransformer,setup.fixTech.Tech{tt}))
                addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: instCapacity_heatTransf[n,''' setup.fixTech.Tech{tt} ''']/10^3 >= minimum_' setup.fixTech.Tech{tt} '[n]/10^3;\n'])];
                addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
            end

            if sum(ismember(activeElements.labels.chpTransformer,setup.fixTech.Tech{tt}))
                addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: instCapacity_chpTransf[n,''' setup.fixTech.Tech{tt} ''']/10^3 >= minimum_' setup.fixTech.Tech{tt} '[n]/10^3;\n'])];
                addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
            end

            if sum(ismember(activeElements.labels.gasTransformer,setup.fixTech.Tech{tt}))
                addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: instCapacity_ptgTransf[n,''' setup.fixTech.Tech{tt} ''']/10^3 >= minimum_' setup.fixTech.Tech{tt} '[n]/10^3;\n'])];
                addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
            end

            if sum(ismember(activeElements.labels.storage,setup.fixTech.Tech{tt}))
                addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: storageCapacity[n,''' setup.fixTech.Tech{tt} ''']/10^3 >= minimum_' setup.fixTech.Tech{tt} '[n]/10^3;\n'])];
                addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
            end

            if sum(ismember(activeElements.labels.storageInterface,setup.fixTech.Tech{tt}))
                addString = [addString sprintf(['s.t. ' (setup.fixTech.Tech{tt}) '_minimum{n in N}: storageInterface[n,''' setup.fixTech.Tech{tt} ''']/10^3 >= minimum_' setup.fixTech.Tech{tt} '[n]/10^3;\n'])];
                addParam = [addParam sprintf(['param minimum_' (setup.fixTech.Tech{tt}) '{n in N};\n'])];
            end
        end
    end

    %% Model compilation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Generate Model file %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % paths and files
    modelFile = ['GTM' name '.mod'];
    modelFileLP = strrep(modelFile,'.mod','.lp');
    modelFileSol = strrep(modelFile,'.mod','.sol');
    modelFileBas = strrep(modelFile,'.mod','.bas');
    modelFileMPS = strrep(modelFile,'.mod','.mps');

    % concatenation of data files string
    TMPdataFilesStr = fieldnames(datFiles);
    dataFilesString = '';

    for i=1:length(TMPdataFilesStr)
        dataFilesString = [dataFilesString ' -d "' rootDir filesep 'projects' filesep pName filesep 'dat-files' filesep getfield(datFiles,char(TMPdataFilesStr(i))) '"'];
        %dataFilesString = [dataFilesString ' -d ' rootDir filesep 'projects' filesep pName filesep 'dat-files' filesep getfield(datFiles,char(TMPdataFilesStr(i))) ''];
    end

    if setup.Heat.Flag
        heat = [];
    else
        heat = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'withoutHeat.mod']);
    end

    if setup.Mobility
        mobility = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'withMobility.mod']);
        mobilityVars = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'varsMobility.mod']);
    else
        mobility = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'withoutMobility.mod']);
        mobilityVars = [];
    end

    if setup.DesalinationFlag
        desalination = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'withDesalination.mod']);
    else
        desalination = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'withoutDesalination.mod']);
    end

    if setup.IndustryFlag
        industry = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'withIndustry.mod']);
        industryVars = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'varsIndustry.mod']);
    else
        industry = fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'withoutIndustry.mod']);
        industryVars = [];
    end

    % concatenation of model string
    modelString = [sprintf('param hourEnd := ') num2str(endHour) sprintf(';\n') ...
        fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep modFiles.paramsAndSets]) ...
        addParam,...
        fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'varsMain.mod']) ...
        mobilityVars,...
        industryVars,...
        fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep 'objAndConsNewNDMI.mod']) ...
        heat,...
        desalination,...
        mobility,...
        industry,...
        addString, ...
        fileread([rootDir filesep 'projects' filesep pName filesep 'dat-files' filesep modFiles.grid]) ...
        fileread([rootDir filesep 'projects' filesep pName filesep 'models' filesep modFiles.end])];

    % write concatenated string to model file
    fModel = fopen([rootDir filesep 'projects' filesep pName filesep 'input-data' filesep modelFile],'w');
    fprintf(fModel,'%s',modelString);
    fclose(fModel);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Process model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Building .LP file for a solver...\n');
    glpsolArgs = [' --check --math "' rootDir filesep 'projects' filesep pName filesep 'input-data' filesep modelFile '" --wlp "' rootDir filesep 'projects' filesep pName filesep modelFileLP '"'];

    %if isunix
    %    eval(['!' 'glpsol ' glpsolArgs dataFilesString])	% ! = used for run a bash command
    %else
    %    system(['"' setup.GLPKPath '"' glpsolArgs dataFilesString '"'],'-echo'); % ! = used for run a bash command
    %end
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
    else
        fprintf('[done] Model .LP file is compiled.\n\n');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Solve problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    try
        setup.SolverFeasibilityTol;
    catch
        setup.SolverFeasibilityTol = 1.0e-6;
    end

    try
        setup.SolverOptimalityTol;
    catch
        setup.SolverOptimalityTol = 1.0e-8;
    end

    % default is usage of MATLAB MOSEK toolbox
    if strcmp(solver,'Mosek')

        
        if (~exist('solveTool','var') || strcmp(solveTool,'mtlb'))
            
            disp('Solve the model using Mosek matlab function')
            
            % read problem to pStruct.prob
            [~,pStruct] = mosekopt(['read(' rootDir filesep 'projects' filesep pName filesep modelFileLP ')']);

            % setting parameters for MOSEK solver
            param.MSK_IPAR_NUM_THREADS = numThreads; % specifies number of threads used by interior-point solver
            try setup.SolverPreSolve
            catch
                setup.SolverPreSolve = 0;
            end
            if ~setup.SolverPreSolve
                param.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';
            end
            param.MSK_IPAR_INTPNT_BASIS  = 'MSK_BI_NEVER'; % solely using interior-point solver
            param.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON'; % switch infeasibility report to ON
            param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 600000;
            param.MSK_DPAR_INTPNT_TOL_PFEAS = setup.SolverFeasibilityTol;
            param.MSK_DPAR_INTPNT_TOL_DFEAS = setup.SolverFeasibilityTol;
            param.MSK_DPAR_INTPNT_TOL_REL_GAP = setup.SolverOptimalityTol;
            param.MSK_IPAR_INFEAS_REPORT_LEVEL = 10;



            % settings for writing log to text file
            if setup.endHour==8760
                fid                = fopen([rootDir filesep 'projects' filesep pName filesep 'output' filesep 'log_' name '.log'],'w');
            else
                fid                = fopen([rootDir filesep 'projects' filesep pName filesep 'output' filesep 'log_' name 'sm.log'],'w');
            end
            callback.log       = 'LogPrint';
            callback.loghandle = fid;

            %  solve problem via mosek solver call

            [r, res] = mosekopt('minimize info',pStruct.prob,param,callback);

            fclose(fid);
            
            if r==0
                disp('Solution is feasible');
            else
                disp('Some solution provided');
            end

            disp('importing results to MATLAB...');
            try
                if isMATLABReleaseOlderThan('R2022b')
                    [results]=PrepareResults(setup,pStruct.prob.names.var,res.sol.itr.xx,res.sol.itr.pobjval);
                else
                    [results] = PrepareResultsD(setup,pStruct.prob.names.var,res.sol.itr.xx,res.sol.itr.pobjval, activeElements,systemParams.IndexID,length(systemParams.IndexNumNodes),length(systemParams.gridMat),endHour);
                end
            catch
                [results]=PrepareResults(setup,pStruct.prob.names.var,res.sol.itr.xx,res.sol.itr.pobjval);
            end

            results.SolvParams = param;
            results.code = res.rcode;
            if ~setup.Mobility
                results = results_WO_Mobility(results,length(results.OPT_SIZE_RWIN),endHour);
            end

            if ~setup.IndustryFlag
                results = results_WO_Industry(results,length(results.OPT_SIZE_RWIN),endHour);
            end
            disp('[done] result import.');

        elseif strcmp(solveTool,'sysc')
            %  (for large optimisation-problems): system call

            disp('Solve the model using Mosek system call')

			% Define solution + log file names
            solFile = [rootDir filesep 'projects' filesep pName filesep modelFileSol];
            logFile = [rootDir filesep 'projects' filesep pName filesep 'output' filesep strrep(modelFileLP,'lp','log')];
			idcs   = strfind(setup.SolverPath,filesep);
            if ismac
				% Code to run on Mac platform
                system(['"' setup.SolverPath(1:idcs(end-1)-1) filesep 'tools' filesep 'platform' filesep 'osxaarch64' filesep 'bin' filesep 'mosek"' ...
                    ' -d MSK_IPAR_NUM_THREADS ' num2str(numThreads) ...
                    ' -d MSK_IPAR_INTPNT_BASIS MSK_BI_NEVER' ...
                    ' -d MSK_IPAR_INFEAS_REPORT_AUTO MSK_ON' ...
                    ' -d MSK_IPAR_INTPNT_MAX_ITERATIONS ' num2str(600000) ...
                    ' -d MSK_DPAR_INTPNT_TOL_PFEAS ' num2str(setup.SolverFeasibilityTol) ...
                    ' -d MSK_DPAR_INTPNT_TOL_DFEAS ' num2str(setup.SolverFeasibilityTol) ...
                    ' -d MSK_DPAR_INTPNT_TOL_REL_GAP ' num2str(setup.SolverOptimalityTol) ...
                    ' -d MSK_IPAR_INFEAS_REPORT_LEVEL ' num2str(1000) ...
                    ' "' rootDir filesep 'projects' filesep pName filesep modelFileLP '"'...
                    ' -itro "' solFile '"' ...
                    ' -q "' logFile '"' ' -silent']);       
			elseif isunix
				% Code to run on Linux platform
				disp('Linux platform is not supported yet')
				%eval(['!mosek -d MSK_IPAR_INTPNT_NUM_THREADS ' num2str(numThreads) ' '  rootDir filesep pName filesep modelFileLP ' -q ' rootDir filesep pName filesep strrep(modelFileLP,'lp','log')])
			elseif ispc
				% Code to run on Windows platform
				%system(['"' erase(setup.SolverPath,[filesep 'toolbox' filesep 'r2017aom']) filesep 'tools' filesep 'platform' filesep 'win64x86' filesep 'bin' filesep 'mosek.exe"' ...
                system(['"' setup.SolverPath(1:idcs(end-1)-1) filesep 'tools' filesep 'platform' filesep 'win64x86' filesep 'bin' filesep 'mosek.exe"' ...
                    ' -d MSK_IPAR_NUM_THREADS ' num2str(numThreads) ...
                    ' -d MSK_IPAR_INTPNT_BASIS MSK_BI_NEVER' ...
                    ' -d MSK_IPAR_INFEAS_REPORT_AUTO MSK_ON' ...
                    ' -d MSK_IPAR_INTPNT_MAX_ITERATIONS ' num2str(600000) ...
                    ' -d MSK_DPAR_INTPNT_TOL_PFEAS ' num2str(setup.SolverFeasibilityTol) ...
                    ' -d MSK_DPAR_INTPNT_TOL_DFEAS ' num2str(setup.SolverFeasibilityTol) ...
                    ' -d MSK_DPAR_INTPNT_TOL_REL_GAP ' num2str(setup.SolverOptimalityTol) ...
                    ' -d MSK_IPAR_INFEAS_REPORT_LEVEL ' num2str(1000) ...
                    ' "' rootDir filesep 'projects' filesep pName filesep modelFileLP '"'...
                    ' -itro "' solFile '"' ...
                    ' -q "' logFile '"' ' -silent']);
			else
				disp('Platform not supported')
			end
            
            disp('importing results to MATLAB...');
            
            disp('reading values from solution file');
            [varNames, solVec, objVal] = PrepareValuesFromSolD(solFile);
            disp('values read');
            
            
                       
            results = PrepareResultsD(setup, varNames, solVec, objVal, activeElements,systemParams.IndexID,length(systemParams.IndexNumNodes),length(systemParams.gridMat),endHour);

            if ~setup.Mobility
                results = results_WO_Mobility(results, length(results.OPT_SIZE_RWIN), endHour);
            end
            if ~setup.IndustryFlag
                results = results_WO_Industry(results, length(results.OPT_SIZE_RWIN), endHour);
            end
            
            disp('[done] result import.');

        else
            error('Wrong solve option chosen or MOSEK Toolbox not installed (or path to this toolbox is not added)')
        end
    elseif strcmp(solver,'Gurobi')

        model = gurobi_read([rootDir filesep 'projects' filesep pName filesep modelFileLP]);

        params.method = 2;
        params.threads = 30;
        params.ScaleFlag = 0;
        params.FeasibilityTol = setup.SolverFeasibilityTol;
        params.OptimalityTol = setup.SolverOptimalityTol;
        params.Crossover = 0;
        params.BarHomogeneous = 1;

        params.method = 2;
        params.threads = numThreads;
        params.ScaleFlag = 0;
        params.BarConvTol = 1.e-5
        params.FeasibilityTol = 1.e-4;
        params.OptimalityTol = 1.e-4;
        params.ObjScale = -0.5
        params.Crossover = 0;
        params.BarHomogeneous = 1;

        params.LogFile = [rootDir filesep 'projects' filesep pName filesep 'output' filesep 'log_' name 'gur.log']
        res=gurobi(model,params)
        if strcmp(res.status,'OPTIMAL')
            resultsStat = 1; %Optimal
        elseif strcmp(res.status,'INF_OR_UNBD')
            resultsStat = -1; %Optimal
        elseif strcmp(res.status,'NUMERIC')
            resultsStat = 11; %NUMERIC
        elseif strcmp(res.status,'SUBOPTIMAL')
            resultsStat = 12; %SUBOPTIMAL
        elseif strcmp(res.status,'INTERRUPTED')
            resultsStat = 10; %Interrupted
        else
            resultsStat = 0; %Unknown
        end
        disp('importing results to MATLAB...');

        try
            if isMATLABReleaseOlderThan('R2022b')
                [results]=PrepareResults(setup,model.varnames,res.x,res.objval);
            else
                [results] = PrepareResultsD(setup,model.varnames,res.x,res.objval, activeElements,systemParams.IndexID,length(systemParams.IndexNumNodes),length(systemParams.gridMat),endHour);
            end
        catch
            [results]=PrepareResults(setup,model.varnames,res.x,res.objval);
        end

        %results=prepareResults(setup,model.varnames,res.x,res.objval);
        if ~setup.Mobility
            results = ResultsNoMobility(results,length(results.OPT_SIZE_RWIN),endHour);
        end

        if ~setup.IndustryFlag
            results = ResultsNoIndustry(results,length(results.OPT_SIZE_RWIN),endHour);
        end
        disp('[done] result import.');

        results.resultsStat =resultsStat;
        results.StatusGur =res.status;
    else
        disp('unknown solver option');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Delete lp-model file %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist([rootDir filesep 'projects' filesep pName filesep modelFileLP])
        delete([rootDir filesep 'projects' filesep pName filesep modelFileLP])
    else warning('No model file (.lp) found.')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save or return results %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(name)
        save([rootDir filesep 'projects' filesep pName filesep 'output' filesep 'results.mat'],'results')
    else
        save([rootDir filesep 'projects' filesep pName filesep 'output' filesep 'results_' name '.mat'],'results')
    end
    disp('Results saved')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end