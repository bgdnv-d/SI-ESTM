function writeFeedIn(fileName,setup,struct,activeElements,endYear)
%
% FUNCTION writeFeedIn(fileName, setup, struct, activeElements, endYear)
%
% Writes feed-in data to a file for the specified scenario.
%
% INPUT:
%            fileName: Desired filename without extension
%            setup: Structure that contains all necessary settings and data for processing
%            struct: Parameter structure obtained from Excel sheets
%            activeElements: Structure with active elements in the system
%            endYear: Last year for which data is written
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 24.07.2025


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

try setup.ProsumersRun;
catch
    setup.ProsumersRun = 0;
end

% reshaping of data structure

for h=1:length(activeElements.labels.elFeedIn)
    for i=1:size(struct.ValueResource,4)
      	 elFeedIn(:,i,h) = struct.ValueResource(ismember(struct.IndexIDR,activeElements.labels.elFeedIn(h)),1:endYear,2,i)';
    end
end

for h=1:length(activeElements.labels.heatFeedIn)
    for i=1:size(struct.ValueResource,4)
      	 heatFeedIn(:,i,h) = round(struct.ValueResource(ismember(struct.IndexIDR,activeElements.labels.heatFeedIn(h)),1:endYear,2,i)'); %% in kWh/m2 (RCSP) or kWh/kW (RRSH)
    end
end

for h=1:length(activeElements.labels.resource)
    for i=1:size(struct.ValueResource,4)
      	 resource(:,i,h) = struct.ValueResource(ismember(struct.IndexIDR,activeElements.labels.resource(h)),1:endYear,1,i)';
    end
    resource(isinf(resource)) = 10^10;
    resource(isnan(resource)) = 0;
end



WriteParamFileBegin(fileName);

WriteParams(fileName,'maxPowerNorm',round(elFeedIn,3),[1:size(elFeedIn,1)],[1:size(elFeedIn,2)],activeElements.labels.elFeedIn)
try
    WriteParams(fileName,'maxHeatNorm',round(heatFeedIn,3),[1:size(heatFeedIn,1)],[1:size(heatFeedIn,2)],activeElements.labels.heatFeedIn)
catch
    warning('no Heat Feedin CF data');
end

try
    WriteParams(fileName,'feedInFuelInput',round(resource.*m_factor,0),[1:size(resource,1)],[1:size(resource,2)],activeElements.labels.resource)
catch
    warning('no Fuels data');
end

WriteParamFileEnd(fileName);
