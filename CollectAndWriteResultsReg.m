function collectAndWriteResultsAll(setup,fileName,varargin)
%
% FUNCTION collectAndWriteResultsAll(setup, fileName, varargin)
%
% Collects and writes result data to a file according to the regions specified.
%
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




if strcmp(setup.projType,'Regions')
    setup.Countries = [[1:numReg]',[1:numReg]'];
end
if strcmp(setup.projType,'Area')
    setup.Countries = [1,numReg];
end

numReg = size(setup.Countries,1);

yearToRead = [setup.startYear*ones(numReg)];


baseData = load([setup.rootDir filesep 'projects' filesep 'Base' filesep 'input-data' filesep 'simulation-input_' 'Base' '.mat'])

yearInd = ones(1,numReg);

yearInd = 2050*ones(1,numReg);
startReg = 1

for Reg = startReg:numReg

    if iscell(setup.Countries)
        Reg1 = setup.Countries{Reg}(1);
        Reg2 = setup.Countries{Reg}(end);

        setup.Regions = setup.Countries{Reg};

    else
        Reg1 = setup.Countries(Reg,1);
        Reg2 = setup.Countries(Reg,2);

        setup.Regions = [setup.Countries(Reg,1):Reg2];

    end



    costYear(Reg)=yearToRead(Reg);
    while costYear(Reg)<=setup.endYear
        try
            calc = PrepareScenarioResultsNewAll(setup,costYear(Reg),baseData,Reg1,Reg2,'Temp')
            PrepareScenarioResultsNewAll(setup,costYear(Reg),baseData,Reg1,Reg2,'Hist')

            yearInd(Reg) = 1+(costYear(Reg)-setup.startYear)/setup.stepYear;

            formatOut = 'yymmdd';
            if strcmp(setup.projType,'Countries')|strcmp(setup.projType,'Regions')
                if strcmp(setup.projType,'Countries')
                    CountryName = strsplit(baseData.systemParams.IndexNodes{Reg1},'-');
                else
                    CountryName = baseData.systemParams.IndexNodes(Reg1);
                end

                if strcmp(CountryName{1},'''BKN')
                    CountryName{1} = baseData.systemParams.IndexNodes{Reg1};
                end

                exelFile = [setup.rootDir filesep projName '_' CountryName{1} '_' fileName '_' datestr(now,formatOut) '.xlsx' ];
                exelFileSm = [setup.rootDir filesep projName '_' CountryName{1} '_Sm_' fileName '_' datestr(now,formatOut) '.xlsx' ];
            end
            if strcmp(setup.projType,'Area')
                exelFile = [setup.rootDir filesep projName '_' datestr(now,formatOut) '_' fileName '_Area.xlsx' ];
            end
            try
                exelFile = regexprep(exelFile,'''','');
            end


            WriteResultsToExelAll(setup,setup.rootDir,setup.projType,exelFile,costYear(Reg),0,0,num2str(Reg1),calc)
			%WriteResultsToExelAll(setup,setup.rootDir,setup.projType,exelFile,costYear(Reg),0,0,'all',calc)
            costYear(Reg) = costYear(Reg) + setup.stepYear;
            yearToRead(Reg) = costYear(Reg);
        catch
            display(['Region ' num2str(Reg1) ' ' num2str(costYear(Reg)) ' is not ready'])
            yearToRead(Reg) = costYear(Reg);
            break
        end

    end

end

end










