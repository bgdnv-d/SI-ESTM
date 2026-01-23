function CollectAndWriteResultsAll(setup,fileName,varargin)
%
% FUNCTION collectAndWriteResults_D(setup, fileName, varargin)
%
% Collects and writes result data to a file for all regions.
%
%
% INPUT:
%            setup:     Structure that contains all necessary settings and data for processing.
%            fileName:  Desired filename without extension
%            varargin:  ADD!!!!!!!!!!!!
%
% OUTPUT:
%            -
%
%Dmitrii Bogdanov
%last change 22.07.2025


if nargin ==3
    numReg = varargin{1};
end

projName = [setup.projType];

if setup.MacroMc
    projName = [projName '_Mc'];
end

if setup.MacroReg
    projName = [projName '_R'];
end

if setup.OvernightFlag
    projName = [projName '_Overnight'];
end

if setup.SC.Flag
    projName = [projName '_EL_SC'];
end

if setup.Heat.Flag
    projName = [projName '_HE'];
end

if setup.Mobility
    projName = [projName '_TR'];
end

if setup.IndustryFlag
    projName = [projName '_IND'];
end

if setup.GasFlag
    projName = [projName '_GAS'];
end

if setup.DesalinationFlag
    projName = [projName '_DES'];
end

if setup.Only.Flag
    projName = [projName '_Only'];
end



yearToRead = [setup.startYear];


baseData = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'simulation-input_' 'Base' '.mat'])


yearInd = 2050*ones(1,numReg);
Reg = 1
Reg1 = baseData.systemParams.IndexNumNodes(Reg);
Reg2 = baseData.systemParams.IndexNumNodes(numReg);
setup.Regions = baseData.systemParams.IndexNumNodes(Reg:numReg)';


costYear(Reg)=yearToRead(Reg);
while costYear(Reg)<=setup.endYear
    try
        calc = PrepareScenarioResultsNewAll(setup,costYear(Reg),baseData,Reg1,Reg2,'Temp','all')


        PrepareScenarioResultsNewAll(setup,costYear(Reg),baseData,Reg1,Reg2,'Hist','all')

        try setupPrint = setup.setupPrint;
        end

        yearInd(Reg) = 1+(costYear(Reg)-setup.startYear)/setup.stepYear;

        formatOut = 'yymmdd';
        if strcmp(setup.projType,'Countries')|strcmp(setup.projType,'Regions')
            
            if strcmp(setup.projType,'Countries')
                CountryName = strsplit(baseData.systemParams.IndexNodes{Reg1},'-');
            else
                CountryName = baseData.systemParams.IndexNodes(Reg1);
            end
                        
            exelFile = [setup.rootDir filesep projName '_' fileName '_' CountryName '_' datestr(now,formatOut) '.xlsx' ];
            exelFileSm = [setup.rootDir filesep projName '_' fileName '_' CountryName '_Sm_' datestr(now,formatOut) '.xlsx' ];
            
        end
        if strcmp(setup.projType,'Area')
            exelFile = [setup.rootDir filesep projName '_' fileName '_' datestr(now,formatOut)  '.xlsx' ];
        end
        try
            exelFile = regexprep(exelFile,'''','');
        end

        WriteResultsToExelAll(setup,setup.rootDir,setup.projType,exelFile,costYear(Reg),0,0,'all',calc)

        costYear(Reg) = costYear(Reg) + setup.stepYear;
        yearToRead(Reg) = costYear(Reg);
    catch
        display(['Region ' num2str(Reg1) ' ' num2str(costYear(Reg)) ' is not ready'])
        yearToRead(Reg) = costYear(Reg);
        break
    end

end

end

