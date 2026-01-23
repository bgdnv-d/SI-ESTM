function writeHydro(fileName,setup,struct,activeElements,endYear)
%
% FUNCTION writeHydro(fileName, setup, struct, activeElements, endYear)
%
% Writes hydro power data to a file for the specified scenario.
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


% feed-in time series for hydro resources
% endYear = 8760;
% endYear = 168;
% reshaping of data structure
for h=1:length(activeElements.labels.hydro)
    for i=1:size(struct.ValueHydroDam,4)
   	 hydro(:,i,h) = struct.ValueHydroDam(ismember(struct.IndexIDH,activeElements.labels.hydro(h)),1:endYear,2,i)';
     dischargePowerDamMin(i) = 0.5*min(struct.ValueHydroDam(ismember(struct.IndexIDH,activeElements.labels.hydro(h)),1:endYear,2,i)');
    end
end



WriteParamFileBegin(fileName);

WriteParams(fileName,'precipitations',round(hydro,3),[1:size(hydro,1)],[1:size(hydro,2)])
WriteParams(fileName,'dischargePowerDamMin',round(dischargePowerDamMin',3),[1:length(struct.IndexNodes)]')

WriteParamFileEnd(fileName);