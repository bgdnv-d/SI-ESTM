function WriteParams(fileName,paramName,values,varargin)
%
% FUNCTION WRITEPARAMS(FILENAME,PARAMNAME,VALUES,VARARGIN)
%
% Writes Parameter files for arbitrary numbers of indices given via varargin
% (indices represent the number of dimensions)
% case 3: scalar (no index)
% case 4: 1 dimension (1 index)
% case 5: 2 dimensions (2 indices)
% case 6: 3 dimensions (3 indices), eg. dimension 1: time, dimension 2: nodes, dimension 3: technology
%
% INPUT
%	FILENAME: parameter file name (incl. ext.)
%	PARAMNAME: name of parameter
%	VALUES: value(s) which will be assigned to the specified parameter
%	VARARGIN: variable number of indices (zero to three)
%
%last change June 18th 2012 by Guido Plessmann


% consistency checks on input data:
% if more than one index is given appropriate values are needed

switch nargin
    case 3
    case 4
        index1 = varargin{1};
        if size(values,1) ~= length(index1)
            error('Values can''t assigned to indices.')
        end
    case 5
        index1 = varargin{1};
        index2 = varargin{2};
        if size(values,1) ~= length(index1) || size(values,2) ~= length(index2)
            error('Values can''t assigned to indices.')
        end
    case 6
        % USE try,catch,error OR warning,error
        index1 = varargin{1};
        index2 = varargin{2};
        index3 = varargin{3};
        if size(values,1) ~= length(index1) || size(values,2) ~= length(index2) || size(values,3) ~= length(index3)
            error('Values can''t assigned to indices.')
        end
end

fid = fopen(fileName,'a');
% line ONE
fprintf(fid,'param %s := \n',paramName);
% line TWO
switch nargin
    case 3
        fprintf(fid,'%.4f \n',values);
        fprintf(fid,';\n');

    case 4
        % line FOUR
        for m = 1:length(index1)
            if isnumeric(index1)
                fprintf(fid,'%d %.4f\n', index1(m), values(m));
            else
                fprintf(fid,'''%s'' %.4f\n', index1{m}, values(m));
            end
        end
        fprintf(fid,';\n');

    case 5
        fprintf(fid,'[*,*]: \n');

        % line THREE
        if isnumeric(index2)
            fprintf(fid, '%d ', index2);
        else
            if ~isempty(index2)
                fprintf(fid, '''%s'' ', index2{:});
            end
        end
        fprintf(fid, ':= \n');

        % line FOUR
        for m = 1:length(index1)
            if isnumeric(index1)
                IndexString = sprintf('%d ', index1(m));
            else
                if ~isempty(index2)
                    fprintf(fid, '''%s'' ', index2{:});
                end
            end
            ValueString = sprintf('%.4f ', values(m, :));
            fprintf(fid, '%s%s\n', IndexString, ValueString);
        end
        fprintf(fid, ';\n');

    case 6
        for k = 1:length(index3)
            if isnumeric(index3)
                fprintf(fid,'[*,*,%d]: \n', index3(k));
            else
                fprintf(fid,'[*,*,''%s'']: \n', index3{k});
            end

            % line THREE
            if isnumeric(index2)
                fprintf(fid, '%d ', index2);
            else
                if ~isempty(index2)
                    fprintf(fid, '''%s'' ', index2{:});
                end
            end
            fprintf(fid, ':= \n');

            % line FOUR
            if all(isnumeric(index1))
                paramValues = [index1', values(:,:,k)];
                fprintfString = repmat([repmat('%.4f ',1,size(values,2)+1) '\n'],1,size(values,1));
                fprintf(fid, fprintfString, paramValues');
                clear paramValues;
            else
                error('currently not available!')
            end
        end
        fprintf(fid,';\n');
end

fclose(fid);
