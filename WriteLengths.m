function writeLengths(fileName,setup,struct)

% writes lengths of transmission lines to data file

%last change 28th March 2018 by Dmitrii Bogdanov

WriteParamFileBegin(fileName);

WriteParams(fileName,'length',round(struct.TLlength),[1:length(struct.TLlength)])

WriteParamFileEnd(fileName);