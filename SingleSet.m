function singleSet(fileName,setName, elements)
% writes a line (to an existing file) which specifies the elements of a set
% setName: name of set
% elements: elements of a set
%
%Dmitrii Bogdanov
%last change 24.07.2025


fid = fopen(fileName,'a');
fprintf(fid,'set %s :=',setName);
for k=1:length(elements)
    fprintf(fid,' %s',elements{k});
end
fprintf(fid,';\n');
fclose(fid);
