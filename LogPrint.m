function logPrint(handle,str)
%
% LOGPRINT write log file information to specified log file
% handle: Is user defined data structure
% str   : Is a log string.
%
%Dmitrii Bogdanov
%last change 23.07.2025


fprintf(handle,'%s',str);