function [varNames, outVect, objVal] = PrepareValuesFromSolD(file)
    % read variables from  MOSEK .sol file 
    % INPUTS:
    % 
    % OUTPUTS:
    %   variable names
    %   variables values, 
    %   objective value
    %   
    % Dmitrii Bogdanov, Omar Elyamany, 2025

    % Open file
    fileID = fopen(file,'r');

    % --- Extract objective value ---
    objVal = [];
    while isempty(objVal)
        tempLine = fgetl(fileID);
        objLine = regexp(tempLine, 'PRIMAL OBJECTIVE\s*:\s*([\-0-9\.Ee+]+)', 'tokens', 'once');
        if ~isempty(objLine)
            objVal = str2double(objLine{1});
        end
    end

    % --- Locate VARIABLES section ---
    variablesSection = 0;
    while ~variablesSection
        tempLine = fgetl(fileID);
        variablesSection = startsWith(tempLine, 'VARIABLES');
    end

    outVect   = [];
    varNames  = {};
    
    fgetl(fileID);
    tempLine = fgetl(fileID);

    while tempLine~=-1

        % Stop when VARIABLES section ends
        if isempty(tempLine) || startsWith(tempLine, 'CONSTRAINTS') || startsWith(tempLine, 'PROBLEM')
            break;
        end

        % Split line into tokens
        tokens = regexp(tempLine, '\s+', 'split');

        % Expected format:
        % {Index, Name, Status, Value, Lower, Upper, DualLower, DualUpper}
        if numel(tokens) >= 4
            % Variable name
            varNames{end+1} = tokens{2}; %#ok<AGROW>

            % Numeric value
            valStr = tokens{4};
            val = str2double(valStr);
            if ~isnan(val)
                outVect(end+1) = val; %#ok<AGROW>
            else
                outVect(end+1) = NaN; %#ok<AGROW>
            end
        end
        tempLine = fgetl(fileID);
    end
    fclose(fileID);
    delete(file);
end


