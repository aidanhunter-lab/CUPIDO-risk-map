function [] = extractVarargin(z, varargin)
extraArguments = ~isempty(z); % any extra arguments to extract?
if extraArguments
    % extra arguments must be formatted as name-value pairs
    nv = length(z);
    isEven = mod(nv, 2) == 0;
    allNamed = all(cellfun(@(x) ischar(x), z(1:2:nv)));
    if ~isEven || ~allNamed
        error('Error in extractVarargin.m -- varargin must be formatted as a (1,2*n) cell array of name-value pairs')
    end
    % assign variables to workspace
    workspace = 'caller';
    if ~isempty(varargin)
        j = strcmp(varargin, 'workspace');
        if any(j), workspace = varargin{find(j) + 1}; end
        if ~ismember(workspace, {'base', 'caller'})
            warning('The "workspace" argument in extractVarargin.m should be either "base" or the default "caller". It has been reset to the default value.')
            workspace = 'caller';
        end
    end
    for j = 1:2:nv
       assignin(workspace, z{j}, z{j+1})
    end
end