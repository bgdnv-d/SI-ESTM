function writeRampingCostData(fileName,setup,systemParams,activeElements)
%
% FUNCTION writeRampingCostData(fileName, setup, systemParams, activeElements)
%
% Writes ramping cost data for technologies to a file.
%
% INPUT:
%            fileName: Desired filename without extension
%            setup: Structure that contains all necessary settings and data for processing
%            systemParams: Structure with system parameters
%            activeElements: Structure with active elements in the system
%
% OUTPUT:
%            - 
%
%Dmitrii Bogdanov
%last change 24.07.2025


WriteParamFileBegin(fileName);

WriteParams(fileName,'rampingCost',round(systemParams.rampingCost(activeElements.activeAllTransformer),3),systemParams.IndexID(activeElements.activeAllTransformer))

WriteParamFileEnd(fileName);