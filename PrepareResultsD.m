function [results] = prepareResultsAcc(setup,namesAndInd,values,totalCost, activeElements,IndexID,NumNodes,NumLines,endHour)

disp('read raw results')

try setup.ProsumersRun;
catch
    setup.ProsumersRun = 0;
end
%tic 
A = dictionary(namesAndInd,values);
%toc
% input data in kW, kWh, kg for CO2 and tonnes of product for industry, m3
% per day for desalination
try setup.units;
catch
    setup.units = 'kW';
end
switch setup.units
    case 'kW'
        m_factor = 1;
    case 'MW'
        m_factor = 0.001;
    case 'GW'
        m_factor = 0.000001;
end
warning('off', 'all')
     
%% EL Feed-in

tech = [IndexID(activeElements.activeElFeedIn)];
for t = 1:length(tech)
    name ='gen_electFeedIn';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.([tech{t} '_EL']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='instCapacity_electFeedIn';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end

   


%% El Transformer        
if ~setup.ProsumersRun
    tech = [IndexID(activeElements.activeElTransformer)];
    for t = 1:length(tech)
        name = 'gen_electTransf';
        key = {};
        for i=1:endHour
            for r = 1:NumNodes   
                    key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
            end    
        end 
        try
            results.([tech{t} '_EL']) = A(key)/m_factor;
            %A = remove(A,key);
        catch
            warning([key{1,1} ' not found'])
        end
    
        name ='instCapacity_electTransf';
        key = {};
        for r = 1:NumNodes   
                key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
        end    
        try
            results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
            %A = remove(A,key);
        catch
            warning([key{1,1} ' not found'])
        end
    end
end




%% Heat
if ~setup.ProsumersRun
    tech = [IndexID(activeElements.activeHeatTransformer)];
    for t = 1:length(tech)
        name = 'gen_heatTransf';
        key = {};
        for i=1:endHour
            for r = 1:NumNodes   
                    key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
            end    
        end 
        try 
            results.([tech{t} '_HE']) = A(key)/m_factor;
            %A = remove(A,key);
        catch
            warning([key{1,1} ' not found'])
        end
    
        name ='instCapacity_heatTransf';
        key = {};
        for r = 1:NumNodes   
                key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
        end    
        try
            results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
            %A = remove(A,key);
        catch
            warning([key{1,1} ' not found'])
        end
    end

else
    tech = [IndexID(activeElements.activeLHeatTransformer)];
    for t = 1:length(tech)
        name = 'gen_heatTransf';
        key = {};
        for i=1:endHour
            for r = 1:NumNodes   
                    key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
            end    
        end 
        try 
            results.([tech{t} '_HE']) = A(key)/m_factor;
            %A = remove(A,key);
        catch
            warning([key{1,1} ' not found'])
        end
    
        name ='instCapacity_heatTransf';
        key = {};
        for r = 1:NumNodes   
                key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
        end    
        try
            results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
            %A = remove(A,key);
        catch
            warning([key{1,1} ' not found'])
        end
    end
end

%% Heat Feed-In

tech = [IndexID(activeElements.activeHeatFeedIn)];
for t = 1:length(tech)
    name = 'gen_heatFeedIn';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try 
        results.([tech{t} '_HE']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='instCapacity_heatFeedIn';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end
tech = [IndexID(activeElements.activeLHeatFeedIn)];
for t = 1:length(tech)
    name = 'gen_heatFeedIn';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try 
        results.([tech{t} '_HE']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='instCapacity_heatFeedIn';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end


%% Storage

tech = [IndexID(activeElements.activeStorage)];
for t = 1:length(tech)
    name ='dischargePower';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.([tech{t} '_EL']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='chargePower';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.(['EL_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='storageState';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.(['SoC_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='storageCapacity';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end

tech = [IndexID(activeElements.activeStorageInterface)];
for t = 1:length(tech)
    
    name ='storageInterface';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end


%% Power to X        

tech = [IndexID(activeElements.activeGasTransformer)];
for t = 1:length(tech)
    name = 'gen_ptgTransf';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.([tech{t} '_GA']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='instCapacity_ptgTransf';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end
 %% Fuels inflow        

tech = [IndexID(activeElements.activeResource)];
for t = 1:length(tech)
    name = 'inputFuel';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.([tech{t} '_FU']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='instCapacity_Resource';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end             





%% CO2 exhaust
key = {};  name = 'fossilCO2_RNGA';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
	results.RNGA_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'fossilCO2_RPET';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
	results.RPET_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'fossilCO2_RHAR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
	results.RHAR_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'fossilCO2';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
	results.Fossil_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'geoheat_out';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RGEO_HE = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'excessGeoHeat';

    for r = 1:NumNodes   
        key(r) = (cellstr([name '(' num2str(r) ')']));       
    end    

try
		results.GHE_EXCESS = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

  


	 % NG flows
key = {};  name = 'importGasToOCGT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TOCG_GAS_FOS = A(key)/m_factor; 
    %A = remove(A,key); 
catch     
    warning([key{1,1} ' not found']) 
end
key = {};  name = 'importGasToCCGT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TCCG_GAS_FOS = A(key)/m_factor; 
catch     
    warning([key{1,1} ' not found']) 
end
key = {};  name = 'importGasToTGCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TGCS_GAS_FOS = A(key)/m_factor; 
    %A = remove(A,key); 
catch     
    warning([key{1,1} ' not found']) 
end
key = {};  name = 'importGasToTICM';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TICM_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToCHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TCNG_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TDNG_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToIndHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.LHIN_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToIndivHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Indiv_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToInd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.LIGA_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToDes';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WDES_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end       
key = {};  name = 'demandProsGasFos';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Pros_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToLNG';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TLNG_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToStRef';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TSMR_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToTSMR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TSMR_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'importGasToTSMC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TSMC_GAS_FOS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

        %SNG flows   
key = {};  name = 'SNGGasToOCGT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TOCG_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'SNGGasToCCGT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TCCG_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'SNGGasToTGCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TGCS_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'SNGGasToTICM';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TICM_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'SNGGasToCHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TCNG_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'SNGGasToHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TDNG_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'SNGGasToIndHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.LHIN_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'SNGGasToInd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.LIGA_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'SNGGasToDes';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WDES_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'demandProsGasRen';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Pros_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'SNGGasToLNG';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TLNG_GAS_REN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

        %Coal flows   
key = {};  name = 'CoalToPower';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RHAR_THPP = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'CoalToTHCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RHAR_THCS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'CoalToCHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RHAR_TCCO = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'CoalToHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RHAR_TDCO = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'CoalToIndHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RHAR_LHIN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

     %Oil flows  
key = {};  name = 'OilToPower';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_TICG = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'OilToTICM';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_TICM = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'OilToCCGT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_TCCG = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'OilToOCGT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_TOCG = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'OilToTGCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_TGCS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'OilToTCNG';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_TCNG = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'OilToCHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_TCOI = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'OilToHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_TDOI = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'OilToIndHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_LHIN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'FOSdiesel';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_DI = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'FOSkerosene';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RPET_KE = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
     %RWOO flows  
key = {};  name = 'BioToPower';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_TBPP = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BioToPowerCCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_TBCS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BioToCHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_TCBP = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BioToCHPCCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_TCBC = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BioToHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_TDBP = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BioToIndHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_LHIN = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BioToCook';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_COOK = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BioToBDS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_RBDS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BioToBiocharCDR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWOO_TRBC = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
     %RWWO flows  
key = {};  name = 'MBBioToPower';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWWO_TBPP = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MBBioToPowerCCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWWO_TBCS = A(key)/m_factor; 
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MBBioToCHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWWO_TCBP = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MBBioToCHPCCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWWO_TCBC = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MBBioToHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWWO_TDBP = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MBBioToIndHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWWO_LHIN = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MBBioToCook';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWWO_COOK = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MBBioToBDS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RWWO_RBDS = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
     %RBMW
key = {};  name = 'WastesToCHPCCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBMW_TCWC = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'WastesToCHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBMW_TMSW = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

     %RBDS   
key = {};  name = 'BIOdiesel';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBDS_DI = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'BIOkerosene';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBDS_KE = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
     %RBGA   
key = {};  name = 'BGAToCHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBGA_TCHP = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
% key = {};  name = 'BGAToHeat'
% 		results.RBGA_TDBP = A(key)/m_factor;
key = {};  name = 'BGAToIndHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBGA_LHIN = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

 key = {};  name = 'BGAToIndCook';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBGA_COOK = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'BGAToBME';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBGA_TBGU = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'demandProsBGA';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBGA_THBG = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'unusedBiomass';

for r = 1:NumNodes   
    key(r) = cellstr([name '(' num2str(r) ')']);       
end    

try
    results.UnusedBiomass = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end





     %Power to heat

key = {};  name = 'powerToTHHR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_THHR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'powerToTHHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_THHP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end







key = {};  name = 'excess_El';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_EXCESS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'excess_He';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_EXCESS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'excess_HeLo';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_EXCESS_Local = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'excess_HeLo_EP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_EXCESS_Local_EP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'excess_HeLo_HP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_EXCESS_Local_HP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
%%
if ~setup.ProsumersRun
%% CHP

tech = [IndexID(activeElements.activeCHPTransformer)];
for t = 1:length(tech)
    name = 'gen_chpElTransf';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try 
        results.([tech{t} '_EL']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
    
    name = 'gen_chpHeTransf';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try 
        results.([tech{t} '_HE']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='instCapacity_chpTransf';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end

%% Hydro

tech = [IndexID(activeElements.activeHydro)];
for t = 1:length(tech)
    name ='gen_Hydro';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.([tech{t} '_EL']) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='storageStateDam';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.(['SoC_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='chargePowerDam';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.(['Ch_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='dischargePowerDam';
    key = {};
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
        results.(['DCh_' tech{t} ]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='instCapacity_Hydro';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end

    name ='storageDam';
    key = {};
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_ST_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end 

%% Transport
tech = [{'MRLI'};{'MRWI'};{'MRBI'};{'MRMI'};{'MRHI'}];
for t = 1:length(tech)
    key = {};  name = 'iceFu_Dem';
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
         results.(['FU_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end

tech = [{'MRLB'};{'MRWB'};{'MRBB'};{'MRMB'};{'MRHB'}];
for t = 1:length(tech)
    key = {};  name = 'bevEl_Dem';
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
         results.(['EL_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end

tech = [{'MRLF'};{'MRWF'};{'MRBF'};{'MRMF'};{'MRHF'}];
for t = 1:length(tech)
    key = {};  name = 'fceFu_Dem';
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
         results.(['HY_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end

tech = [{'MRLP'};{'MRWP'};{'MRBP'};{'MRMP'};{'MRHP'}];
for t = 1:length(tech)
    key = {};  name = 'phvFu_Dem';
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
         results.(['FU_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end


    key = {};  name = 'phvEl_Dem';
    for i=1:endHour
        for r = 1:NumNodes   
                key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' tech{t} ')']));       
        end    
    end 
    try
         results.(['EL_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end

%% Grids

key = {}; 
name = 'fromgrid';
for i = 1:endHour   
        key(i) = (cellstr([name '(' num2str(i) ')']));       
end
try
    results.EL_GRID= A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
     
key = {}; 
name = 'fromgridToEl';
for i = 1:endHour   
        key(i) = (cellstr([name '(' num2str(i) ')']));       
end  
try
results.EL_GRIDforEl= A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

     
key = {}; 
name = 'fromgridToHe';
for i = 1:endHour   
        key(i) = (cellstr([name '(' num2str(i) ')']));       
end  
try
results.EL_GRIDforHe = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
     
key = {}; 
name = 'instCapacitiyTransmLines_DC';
for r = 1:NumLines   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end
try
results.OPT_SIZE_TRTL = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
     
key = {}; 
name = 'instCapacitiyTransmLines_AC';
for r = 1:NumLines   
    key(r) = (cellstr([name '(' num2str(r)  ')']));       
end
try
results.OPT_SIZE_THAO = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
	 
key = {}; 
name = 'transPower_DC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
results.GRID_DC = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

     
key = {}; 
name = 'transPower_AC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
results.GRID_AC = A(key)/m_factor;
%A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

	 
key1 = {};
key2 = {};
name = 'line_DC';
for i=1:endHour
    for r = 1:NumLines   
        key1(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' 'pos' ')']));
        key2(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' 'neg' ')']));
    end    
end
try
results.LINEpos_DC= A(key1)/m_factor;
%A = remove(A,key1);
catch
    warning([key1{1,1} ' not found'])
end
try
results.LINEneg_DC= A(key2)/m_factor;
%A = remove(A,key2);
catch
    warning([key2{1,1} ' not found'])
end

	 
key = {}; 
name = 'line_AC';
for i=1:endHour
    for r = 1:NumLines   
        key1(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' 'pos' ')']));
        key2(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ',' 'neg' ')']));
    end    
end
try
results.LINEpos_AC= A(key1)/m_factor;
%A = remove(A,key1);
catch
    warning([key1{1,1} ' not found'])
end
try
results.LINEneg_AC= A(key2)/m_factor;
%A = remove(A,key2);
catch
    warning([key2{1,1} ' not found'])
end

    % e-fuels import

key = {};  name = 'importH2';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importH2 = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importH2_regas';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importH2_regas = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importH2_liq';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importH2_liq = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importLNG';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importLNG = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importLNG_regas';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importLNG_regas = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importLNG_liq';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importLNG_liq = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importFT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importFT = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importFT_kerosene';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importFT_kerosene = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importFT_diesel';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importFT_diesel = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importEL';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importEL = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importMeOH';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importMeOH = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'importNH3';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.importNH3 = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

% Power-to-Heat
key = {};  name = 'powerToTDHR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TDHR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'powerToTDHP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TDHP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'powerToTDGE';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TDGE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'electrolysis';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TWEL = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'HefromFT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTU_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HefromMET';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TMET_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HefromWEL';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TWEL_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HefromNH3';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TNH3_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HefromMEO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TMEO_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'co2Scrubbing_El';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TCOS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'co2Scrubbing_He';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TCOS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'NH_LINH';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.NH_LINH = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MeOtoMobility';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MeOtoMobility = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'MeOtoMethanolFT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.ME_TFTM = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
   
key = {};  name = 'ME_LIME';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.ME_LIME = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'NH3toMobility';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.NH3toMobility = A(key)/m_factor; 
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
  

key = {};  name = 'h2toMET';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TMET = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toFT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TFTU = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toLH2';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TLH2 = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toMobility';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TRSP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toMEO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TMeO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toNH3';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TNH3 = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toTISH';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TISH = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toCCGT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TCCG = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toOCGT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TOCG = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toTGCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TGCS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toTCNG';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TCNG = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toIndHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_LHIN = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toTICM';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TICM = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'h2toTRSi';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HY_TRSi = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'CO2toMET';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.CO_TMET = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CO2toFT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.CO_TFTU = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CO2toMEO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.CO_TMeO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CO2toCCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.CO_LCCS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'CO2fromSMC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TSMC_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CO2fromGCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TGCS_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CO2fromHCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.THCS_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CO2fromCBC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TCBC_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CO2fromBCS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TBCS_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CO2fromCWC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TCWC_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end



key = {};  name = 'FTkeroseneProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTU_KE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'FTnaphtaProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTU_NA = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'FTdieselProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTU_DI = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
  
key = {};  name = 'BioFTkeroseneProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTB_KE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'BioFTnaphtaProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTB_NA = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'BioFTdieselProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTB_DI = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
  
key = {};  name = 'methanolFTkeroseneProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTM_KE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'methanolFTnaphtaProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTM_NA = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'methanolFTdieselProd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TFTM_DI = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'BIOdiesel';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.BIOdiesel = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'BIOkerosene';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.BIOkerosene = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'FTdiesel';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.FTdiesel = A(key)/m_factor; 
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
 
key = {};  name = 'FTkerosene';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.FTkerosene = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'BioFTdiesel';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.BioFTdiesel = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'BioFTkerosene';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.BioFTkerosene = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
  
key = {};  name = 'methanolFTdiesel';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.methanolFTdiesel = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'methanolFTkerosene';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.methanolFTkerosene = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'GasToLNG';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.GA_TLNG = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'EltoFT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TFTU = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EltoLNG';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TLNG = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EltoLH2';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TLH2 = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
   
key = {};  name = 'EltoMET';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TMET = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EltoMEO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TMeO = A(key)/m_factor; 
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
  
key = {};  name = 'EltoNH3';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TNH3 = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EltoHyS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_THyS = A(key)/m_factor; 
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
  
key = {};  name = 'EltoBioFT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TFTB = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EltoMethanolFT';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TFTM = A(key)/m_factor; 
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'heatToPower';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TSTU = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'powerToHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'biomethane';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.RBME = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'geoheat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.SoC_RGEO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

% key = {};  name = 'LowHeat'
%          results.LowT_HE(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo);
% key = {};  name = 'HighHeat'
%          results.HigT_HE(str2num(varIndicesSplit{1}), str2num(varIndicesSplit{2})) = values(varNo);
key = {};  name = 'LossToHeat';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Loss_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


%% Transport

key = {};  name = 'elToLDV';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToLDV = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'LDV_chRate';
for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.LDV_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToW23';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToW23 = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'W23_chRate';
for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.W23_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToBus';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToBus = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'BUS_chRate';
for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.BUS_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMDV';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMDV = A(key)/m_factor;
key = {};  name = 'MDV_chRate';
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.MDV_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToHDV';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToHDV = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HDV_chRate';
for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.HDV_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'elToLDV_V2G';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToLDV_V2G = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToW23_V2G';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToW23_V2G = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToBus_V2G';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToBus_V2G = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMDV_V2G';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMDV_V2G = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToHDV_V2G';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToHDV_V2G = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


 key = {};  name = 'elToLDV_G2V';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToLDV_G2V = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToW23_G2V';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToW23_G2V = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToBus_G2V';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToBus_G2V = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMDV_G2V';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMDV_G2V = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToHDV_G2V';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToHDV_G2V = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'elToLDV_DumpCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToLDV_DumpCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToW23_DumpCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToW23_DumpCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToBus_DumpCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToBus_DumpCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMDV_DumpCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMDV_DumpCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToHDV_DumpCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToHDV_DumpCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'elToLDV_SmartCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToLDV_SmartCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToW23_SmartCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToW23_SmartCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToBus_SmartCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToBus_SmartCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMDV_SmartCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMDV_SmartCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToHDV_SmartCh';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToHDV_SmartCh = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'elToLDV_DumpDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToLDV_DumpDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToW23_DumpDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToW23_DumpDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToBus_DumpDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToBus_DumpDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMDV_DumpDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMDV_DumpDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToHDV_DumpDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToHDV_DumpDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'elToLDV_SmartDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToLDV_SmartDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToW23_SmartDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToW23_SmartDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToBus_SmartDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToBus_SmartDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMDV_SmartDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMDV_SmartDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToHDV_SmartDisch';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToHDV_SmartDisch = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'MRPF_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MRPF_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MRPE_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MRPE_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MRFF_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MRFF_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MRFE_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MRFE_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'MMPF_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMPF_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMPE_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMPE_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMPH_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMPH_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMPG_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMPG_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMPA_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMPA_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMPM_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMPM_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMFF_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMFF_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMFE_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMFE_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMFH_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMFH_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMFG_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMFG_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMFA_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMFA_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MMFM_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MMFM_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'elToMF';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToMF = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MP_chRate';
for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.MP_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MF_chRate';
for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.MP_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end



key = {};  name = 'MAPF_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MAPF_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MAPE_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MAPE_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MAPH_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MAPH_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MAFF_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MAFF_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MAFE_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MAFE_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'MAFH_Dem';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MAFH_Dem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'elToAP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToAP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'elToAF';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.elToAF = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'AP_chRate';
for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.AP_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'AF_chRate';
for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.AF_chRate = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'EltoMobility';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EltoMobility = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
% key = {};  name = 'h2toMobility'
% 		results.h2toMobility = A(key)/m_factor;
key = {};  name = 'diesToMobility';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.diesToMobility = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'kersToMobility';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.kersToMobility = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
key = {};  name = 'LNGtoMobility';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.LNGtoMobility = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'LH2toMobility';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.LH2toMobility = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


        %%% Desalination

        % Water production
key = {};  name = 'waterDesalination';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Desalination = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'waterRODesalination';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WROD_Desalination = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'waterMSSDesalination';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WMSS_Desalination = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'waterMSCDesalination';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WMSC_Desalination = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'waterMDSDesalination';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WMDS_Desalination = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'waterMDCDesalination';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WMDC_Desalination = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'desalination_Gas_demand';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.GA_WDES = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_Gas_demand_MSC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.GA_WMSC = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_Gas_demand_MDC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.GA_WMDC = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'desalination_Heat_demand';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_WDES = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_Heat_demand_MSS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_WMSS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_Heat_demand_MDS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_WMDS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'desalination_El_demand';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_WDES = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_El_demand_RO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_WROD = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_El_demand_MSS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_WMSS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_El_demand_MSC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_WMSC = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_El_demand_MDS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_WMDS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_El_demand_MDC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_WMDC = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'Vpump_demand';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_WVPU = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'Hpump_demand';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_WHPU = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'TotalDesalinationElDemand';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TotDesDem = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_El_production';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WDES_EL = A(key)/m_factor; 
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
  
key = {};  name = 'desalination_El_production_MSC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WMSC_EL = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'desalination_El_production_MDC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.WMDC_EL = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'waterStorageState';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.SoC_SWAT = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

tech = [IndexID(activeElements.activeDesalination)];
name = 'instCapacitiy_Desalination';
for t = 1:length(tech)
    key = {};  
    for r = 1:NumNodes   
            key(r) = (cellstr([name '(' num2str(r) ',' tech{t} ')']));       
    end    
    try
        results.(['OPT_SIZE_' tech{t}]) = A(key)/m_factor;
        %A = remove(A,key);
    catch
        warning([key{1,1} ' not found'])
    end
end
  %% industry      
key = {};  name = 'TICB_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TICB_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TICI_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TICI_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TISB_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TISB_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TISH_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TISH_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TISR_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TISR_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TISE_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TISE_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIAA_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIAA_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIAM_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIAM_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIAR_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIAR_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIPP_PR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIPP_PR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'EL_TICB';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TICB = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EL_TICI';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TICI = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EL_TISB';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TISB = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EL_TISH';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TISH = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EL_TISR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TISR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EL_TISE';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TISE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EL_TIAM';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TIAM = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EL_TIAR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TIAR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'EL_TIPP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TIPP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'HE_TICB';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TICB = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HE_TICI';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TICI = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HE_TISB';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TISB = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HE_TISH';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TISH = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HE_TISR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TISR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HE_TISE';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TISE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HE_TIAA';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TIAA = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HE_TIAR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TIAR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HE_TIPP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TIPP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'HC_TISB';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HC_TISB = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CC_TISH';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.CC_TISH = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CC_TISR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.CC_TISR = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'CC_TISE';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.CC_TISE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'TIAM_HE';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIAM_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIAR_HE';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIAR_HE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'El_to_TPSC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TPSC = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TPSP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TPSP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'He_to_TPSC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TPSC = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'He_to_TPSP';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TPSP = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'Wa_to_TRAD';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Wa_to_TRAD = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'Wa_to_TRMI';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Wa_to_TRMI = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'Wa_to_TRMO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Wa_to_TRMO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'Wa_to_TRBC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.Wa_to_TRBC = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRSS';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRSS = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRSL';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRSL = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRMI';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRMI = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRMO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRMO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRME';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRME = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRSi';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRSi = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRAD';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRAD = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TREW';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TREW = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRBC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRBC = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'El_to_TRGE';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EL_TRGE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'He_to_TRME';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TRME = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'He_to_TRSi';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TRSi = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'He_to_TRGE';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HE_TRGE = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TICB_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TICB_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TICI_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TICI_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TISB_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TISB_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TISR_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TISR_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TISH_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TISH_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TISE_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TISE_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIAA_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIAA_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIAR_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIAR_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIAM_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIAM_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'TIPP_CO';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.TIPP_CO = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'CO2EmissionsfromInd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.CO2EmissionsfromInd = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'EltoInd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.EltoInd = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'HighHetoInd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.HighHetoInd = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MidHetoInd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MidHetoInd = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'MidHefromInd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.MidHefromInd = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end

key = {};  name = 'LowHefromInd';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.LowHefromInd = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end


key = {};  name = 'TIPP_WWtot';

for r = 1:NumNodes   
    key(r) = (cellstr([name '(' num2str(r) ')']));       
end    
try
    results.TIPP_WWtot = A(key)/m_factor;
    %A = remove(A,key);
catch
    warning([key{1,1} ' not found'])
end
 

key = {};  name = 'addElPlus';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.addElPlus = A(key)/m_factor;
    %A = remove(A,key);
catch
    %warning([key{1,1} ' not found'])
end
 
key = {};  name = 'addElMinus';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.addElMinus = A(key)/m_factor;
    %A = remove(A,key);
catch
    %warning([key{1,1} ' not found'])
end

key = {};  name = 'addElSoC';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.addElSoC = A(key)/m_factor;
    %A = remove(A,key);
catch
    %warning([key{1,1} ' not found'])
end
 

key = {};  name = 'extraHydroRoR';
for i=1:endHour
    for r = 1:NumNodes   
        key(i,r) = (cellstr([name '(' num2str(i) ',' num2str(r) ')']));       
    end    
end
try
    results.extraHydroRoR = A(key)/m_factor;
    %A = remove(A,key);
catch
    %warning([key{1,1} ' not found'])
end
   
end



results.totalCost = totalCost;
warning('on', 'all')
end

%  [pathstr, name, ext] = fileparts(file);
%  save([dir name '.mat'])
