function writeParamFileBegin(fileName)
%
% FUNCTION writeParamFileBegin(fileName)
%
% Starts a new parameter file by writing initial formatting or structure.
%
% INPUT:
%            fileName: Desired filename without extension
%
% OUTPUT:
%            - (writes beginning of parameter file, no output variables)
%
%Dmitrii Bogdanov
%last change 24.07.2025


fid = fopen(fileName,'w');
fprintf(fid,'data;\n\n');
fclose(fid);