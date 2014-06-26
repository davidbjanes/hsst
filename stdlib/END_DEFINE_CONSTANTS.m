function varargout = END_DEFINE_CONSTANTS(throwError,ignoreCase)
%END_DEFINE_CONSTANTS  Overrides default constant values
%
%   extras = END_DEFINE_CONSTANTS(*throwError,*ignoreCase)
%
%   END_DEFINE_CONSTANTS(*throwError,*ignoreCase) Will not print "extras"
%   output (optional output)
%
%   INPUTS
%   =======================================================================
%   throwError : (default true), if false will return a non-empty error message
%   ignoreCase : (default true), if false, case sensitive differences are respected
%
%   OUTPUT
%   =======================================================================
%   extras : (struture)
%   	.errMsg - an error message if any occurred
%       .initDefaultStruct  - structure of constant names and values before overriding
%       .finalDefaultStruct - structure of constant names and values at end
%       .newValuesStruct    - structure that only includes new values
%
%   EXAMPLE:
%   ==================================================
%   1) The Call:
%      myFunc(2,3,'C1',2,'C2',-3)
%               or
%      params.C1 = 2
%      params.C2 = -3
%      myFunc(2,3,params)
%      %NOTE: One might be preferred over another
%
%   2) The function definition
%       function myFunc(a,b,varargin)
%       DEFINE_CONSTANTS
%           C1 = 5;
%           C2 = 10;
%           NAME = 'temp';
%       END_DEFINE_CONSTANTS
%
%
%       END RESULT:
%       C1 = 2
%       C2 = -3
%       NAME = 'temp' %still, no change
%
%   ==================================================
%
%   See Also:
%       DEFINE_CONSTANTS
%       variablesToStruct
%       toCellArrayString

%Moved this up, goal is to quit as early as possible ...
[initialVariables,params] = gsr_DEFINE_CONSTANTS('get');

%PARAMS FORMAT
%==========================================================================
%either structure or prop/value cell

%Quit if there is nothing to do ...
if isempty(params) && nargout == 0;
   return
end

if ~exist('throwError','var')
   throwError = true; 
end

if ~exist('ignoreCase','var')
   ignoreCase = true; 
end

%Grabbing new variables that have been defined since DEFINE_CONSTANTS
%--------------------------------------------------------------------------
temp  = evalin('caller','whos');
vars2 = {temp.name};

%Need to remove '(unassigned)' -> nested functions
classValue = {temp.class};
vars2(strcmpi(classValue,'(unassigned)')) = [];

%OLD CODE: used setdiff
%constantNames = setdiff(vars2,initialVariables);
%I think this tends to be faster since we tend to have few variables
keepMask = true(1,length(vars2));
for iVar = 1:length(initialVariables)
   keepMask(strcmp(vars2,initialVariables{iVar})) = false; 
end
constantNames = vars2(keepMask);

if isempty(constantNames)
    %This might occur in a script
    formattedWarning('No constants detected, please debug')
    %myCaller = evalin('caller','mfilename'); -> if empty
    %indicates access to the base workspace
end

if nargout
    extras = struct;
    fString = ['variablesToStruct(' toCellArrayString(constantNames) ')'];
    extras.initDefaultStruct = evalin('caller',fString);
end

%CONVERTING THE PROPERTY VALUE PAIRS TO A STRUCTURE
%==========================================================================
if iscell(params)
    params = params(:)';
    params = cell2struct(params(2:2:end),params(1:2:end),2);
end

if isempty(params) && nargout
    extras.newValuesStruct    = struct([]);
    extras.finalDefaultStruct = extras.initDefaultStruct;
    extras.errMsg = [];
    varargout{1} = extras;
    return
end

%DETERMING BAD VARIABLES
%===========================================================
fn = fieldnames(params);
if ignoreCase
    [isPresent,loc] = ismember_str(upper(fn),upper(constantNames));
else
    [isPresent,loc] = ismember_str(fn,constantNames);
end

errMsg = [];
if ~all(isPresent)
    badVariables = fn(~isPresent);
    errMsg = ['Bad variable names given in input structure: ' cellArrayToString(badVariables,',')];
    if throwError
        error(errMsg)
    else
        for iBad = 1:length(badVariables)
            params = rmfield(params,badVariables{iBad});
        end
        loc(~isPresent) = [];
        fn(~isPresent)  = [];
    end
end

%ASSIGNMENT IN CALLER
%==========================================================================
for i = 1:length(fn)
    %NOTE: By using constantNames we ensure case matching
    assignin('caller',constantNames{loc(i)},params.(fn{i}))
end

if nargout
    extras.newValuesStruct = params;
    extras.errMsg = errMsg;
    extras.finalDefaultStruct = evalin('caller',fString);
end


if nargout
    varargout{1} = extras;
end

end




