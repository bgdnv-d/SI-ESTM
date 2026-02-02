function writeDemand(fileName,setup,systemParams,activeElements,costYear,endHour)
%
% FUNCTION writeDemand(fileName, setup, systemParams, activeElements, costYear, endHour)
%
% Writes demand data to a file for the given scenario.
%
% INPUT:
%            fileName: Desired filename without extension
%            setup: Structure that contains all necessary settings and data for processing
%            systemParams: Structure with system parameters for the scenario
%            activeElements: Structure with active elements in the system
%            costYear: Year to which costs are adjusted
%            endHour: Last hour of the simulation period
%
% OUTPUT:
%
%
%Dmitrii Bogdanov
%last change 24.07.2025


% writes demand time series in gmpl conform data format
%endYear = 8760;
%endYear = 168;

try setup.ProsumersRun;
catch
    setup.ProsumersRun = 0;
end

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

% reshaping of data structure
for h=1:length(activeElements.labels.load)
    for i=1:size(systemParams.ValueLoad,4)

      	 demand(:,i,h) = (systemParams.ValueLoad(ismember(systemParams.IndexIDL,activeElements.labels.load(h)),1:endHour,2,i)');
    end
end


WriteParamFileBegin(fileName);


for i = 1:size(activeElements.labels.load,1)

    WriteParams(fileName,['demand_' activeElements.labels.load{i}],round(demand(:,:,i)*m_factor,6),[1:size(demand,1)],[1:size(demand,2)]);

end

for i = 1:size(activeElements.labels.load,1)

    WriteParams(fileName,['PeakDemand_' activeElements.labels.load{i}],round(max(demand(:,:,i)*m_factor),6)',[1:length(systemParams.IndexNodes)]);

end



try setup.ModelType;
catch
    setup.ModelType = 'Main';
end
switch setup.ModelType
    case 'Main'
        if setup.Mobility
            WriteParams(fileName,'LDV_connProf' ,round(systemParams.Mobility.LDV_connProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'W23_connProf' ,round(systemParams.Mobility.W23_connProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'BUS_connProf' ,round(systemParams.Mobility.BUS_connProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'MDV_connProf' ,round(systemParams.Mobility.MDV_connProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'HDV_connProf' ,round(systemParams.Mobility.HDV_connProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);


            WriteParams(fileName,'LDV_chProf' ,round(systemParams.Mobility.LDV_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'W23_chProf' ,round(systemParams.Mobility.W23_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'BUS_chProf' ,round(systemParams.Mobility.BUS_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'MDV_chProf' ,round(systemParams.Mobility.MDV_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'HDV_chProf' ,round(systemParams.Mobility.HDV_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'MP_chProf' ,round(systemParams.Mobility.MP_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'MF_chProf' ,round(systemParams.Mobility.MF_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'AP_chProf' ,round(systemParams.Mobility.AP_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
            WriteParams(fileName,'AF_chProf' ,round(systemParams.Mobility.AF_chProf(:,1:endHour),3)',[1:endHour],[1:length(systemParams.IndexNodes)]);
        end
        try
            WriteParams(fileName,'demandProsGas' ,round(systemParams.ProsDemand_Gas(:,1:endHour)*m_factor,6)',[1:endHour],[1:size(systemParams.ProsDemand_Gas,1)]);
        catch
            'No Prosumers Gas consumption data'
            WriteParams(fileName,'demandProsGas' ,round(zeros(length(systemParams.IndexNodes),endHour)*m_factor,6)',[1:endHour],[1:length(systemParams.IndexNodes)]);
        end
        try
            WriteParams(fileName,'demandProsOil' ,round(systemParams.ProsDemand_Oil(:,1:endHour)*m_factor,6)',[1:endHour],[1:size(systemParams.ProsDemand_Oil,1)]);
        catch
            'No Prosumers Oil consumption data'
            WriteParams(fileName,'demandProsOil' ,round(zeros(length(systemParams.IndexNodes),endHour)*m_factor,6)',[1:endHour],[1:length(systemParams.IndexNodes)]);
        end
        try
            WriteParams(fileName,'demandProsBGA' ,round(systemParams.ProsDemand_BGA(:,1:endHour)*m_factor,6)',[1:endHour],[1:size(systemParams.ProsDemand_BGA,1)]);
        catch
            'No Prosumers BGA consumption data'
            WriteParams(fileName,'demandProsBGA' ,round(zeros(length(systemParams.IndexNodes),endHour)*m_factor,6)',[1:endHour],[1:length(systemParams.IndexNodes)]);
        end
        try
            WriteParams(fileName,'demandProsWOO' ,round(systemParams.ProsDemand_WOO(:,1:endHour)*m_factor,6)',[1:endHour],[1:size(systemParams.ProsDemand_WOO,1)]);
        catch
            'No Prosumers Biomass consumption data'
            WriteParams(fileName,'demandProsWOO' ,round(zeros(length(systemParams.IndexNodes),endHour)*m_factor,6)',[1:endHour],[1:length(systemParams.IndexNodes)]);
        end
end

WriteParams(fileName,'annualDemand',round(shiftdim(round(systemParams.Instalations(:,find(systemParams.IndexYears==costYear),find(ismember(systemParams.IndexID,activeElements.labels.load)))),2)'*m_factor,6),[1:length(systemParams.IndexNodes)],activeElements.labels.load);

if ~setup.ProsumersRun
    %% for Macro simulations
    if setup.MacroReg
        try
            for i =1:length(systemParams.IndexNodes)
                try
                    temp = load([setup.rootDir '\projects\Base\input-data\Macro\RegProf\SP\Profile_AC_' num2str(systemParams.IndexNumNodes(i)) '_tot_' num2str(costYear) '.mat']);
                catch
                    temp.profAC = zeros(8760,1);
                end
                totProfilesAC(1:endHour,i) = temp.profAC;
                try
                    temp = load([setup.rootDir '\projects\Base\input-data\Macro\RegProf\SP\Profile_DC_' num2str(systemParams.IndexNumNodes(i)) '_tot_' num2str(costYear) '.mat']);
                catch
                    temp.profDC = zeros(8760,1);
                end

                totProfilesDC(1:endHour,i) = temp.profDC;

            end
            WriteParams(fileName,'addDemand' ,round((systemParams.totProfilesAC(1:endHour,:) + systemParams.totProfilesDC(1:endHour,:))*m_factor,6),[1:endHour],[1:length(systemParams.IndexNodes)]);
        catch
            WriteParams(fileName,'addDemand' ,(zeros(endHour,length(systemParams.IndexNodes))),[1:endHour],[1:length(systemParams.IndexNodes)]);

            display('No data for totProfilesAC or struct.totProfilesDC;')
        end

        try
            load([setup.rootDir '\projects\Base\input-data\Macro\RegProf\SP\EfuelsFlows_' num2str(systemParams.IndexNumNodes(1)) '.mat']);

            WriteParams(fileName,'FTkeroseneInflow',EfuelsFlows.TFTU_KE(EfuelsFlows.years==costYear)*m_factor);
            WriteParams(fileName,'FTdieselInflow',EfuelsFlows.TFTU_DI(EfuelsFlows.years==costYear)*m_factor);

        catch
            WriteParams(fileName,'FTkeroseneInflow' ,0);
            WriteParams(fileName,'FTdieselInflow' ,0);

            display('No data for EfuelsFlows')

        end



        try
            MeOHtemp=zeros(length(systemParams.IndexNumNodes),1);
            load([setup.rootDir '\projects\Base\input-data\Macro\RegProf\SP\EfuelsFlows_' num2str(systemParams.IndexNumNodes(1)) '.mat']);

            WriteParams(fileName,'MeOHInFlow',0.99*EfuelsFlows.TMeO_GA(EfuelsFlows.years==costYear)*m_factor);
        catch
            WriteParams(fileName,'MeOHInFlow' ,0);
            display('No data for e-Methanol Flows')
        end



        try
            NH3temp=zeros(endHour,length(systemParams.IndexNumNodes));
            load([setup.rootDir '\projects\Base\input-data\Macro\RegProf\SP\EfuelsFlows_' num2str(systemParams.IndexNumNodes(1)) '.mat']);

            for i=1:length(systemParams.IndexNumNodes)
                try
                    temp.NH3=getfield(EfuelsFlows,['TNH3_' num2str(systemParams.IndexNumNodes(i)) '_prof']);
                catch
                    temp.NH3=zeros(7,endHour);
                end
                NH3temp(:,i)=temp.NH3(EfuelsFlows.years==costYear,:)';
            end

            WriteParams(fileName,'NH3InFlow',NH3temp(1:endHour,:)*m_factor,[1:endHour],[1:length(systemParams.IndexNodes)]);

        catch
            WriteParams(fileName,'NH3InFlow',(zeros(endHour,length(systemParams.IndexNodes))),[1:endHour],[1:length(systemParams.IndexNodes)]);
            display('No data for e-Ammonia Flows')
        end


        try
            LNGtemp=zeros(endHour,length(systemParams.IndexNumNodes));
            load([setup.rootDir '\projects\Base\input-data\Macro\RegProf\SP\EfuelsFlows_' num2str(systemParams.IndexNumNodes(1)) '.mat']);

            for i=1:length(systemParams.IndexNumNodes)
                try
                    temp.LNG=getfield(EfuelsFlows,['TLNG_' num2str(systemParams.IndexNumNodes(i)) '_prof']);
                catch
                    temp.LNG=zeros(7,endHour);
                end
                LNGtemp(:,i)=temp.LNG(EfuelsFlows.years==costYear,:)';
            end

            WriteParams(fileName,'LNGInFlow',LNGtemp(1:endHour,:)*m_factor,[1:endHour],[1:length(systemParams.IndexNodes)]);

        catch
            WriteParams(fileName,'LNGInFlow',(zeros(endHour,length(systemParams.IndexNodes))),[1:endHour],[1:length(systemParams.IndexNodes)]);
            display('No data for e-LNG Flows')
        end
    else
        WriteParams(fileName,'addDemand' ,(zeros(endHour,length(systemParams.IndexNodes))),[1:endHour],[1:length(systemParams.IndexNodes)]);
        WriteParams(fileName,'FTkeroseneInflow' ,0);
        WriteParams(fileName,'FTdieselInflow' ,0);
        WriteParams(fileName,'MeOHInFlow' ,0);
        WriteParams(fileName,'NH3InFlow',(zeros(endHour,length(systemParams.IndexNodes))),[1:endHour],[1:length(systemParams.IndexNodes)]);
        WriteParams(fileName,'LNGInFlow',(zeros(endHour,length(systemParams.IndexNodes))),[1:endHour],[1:length(systemParams.IndexNodes)]);
    end
end

WriteParamFileEnd(fileName);
