function writeParamFileEnd(fileName)
%
% FUNCTION writeParamFileEnd(fileName)
%
% Finalizes a parameter file by writing any required ending statements.
%
% INPUT:
%            fileName: Desired filename without extension
%
% OUTPUT:
%            - (writes end of parameter file, no output variables)
%
%Dmitrii Bogdanov
%last change 24.07.2025


fid = fopen(fileName,'a');
fprintf(fid,'end;\n');
fclose(fid);