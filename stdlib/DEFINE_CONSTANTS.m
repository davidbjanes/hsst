function DEFINE_CONSTANTS(params)
%DEFINE_CONSTANTS Defines variables in memory (use all caps)
%
%   DEFINE_CONSTANTS(*params) 
%
%   Defines all variables in memory in the caller
%   so that a call to END_DEFINE_CONSTANTS will be able to tell
%   what constants have been defined
%
%   OPTIONAL INPUTS
%   =======================================================================
%   params: (default, grabs "varargin" from caller)
%           Formats can be as follows:
%           1) a structure (field names are properties, values are values)
%           2) the more typical property/value pairs
%           
%   See Also: 
%       END_DEFINE_CONSTANTS
%
%   CAVEATS: 
%   1) This function only works inside functions or classes where there
%   is not access to defined variables on repeat calls (like in a script)
%   2) This function does not support threaded calls or nested calls. In
%       general this means that the constants can not call functions which
%       themselves call this function (could be changed). In addition
%       parallel computation could screw this up as well.
%

%NOTE: Internal storage of variables done using gsr_DEFINE_CONSTANTS

%Grab variables in caller's workspace
%--------------------------------------------------------------------------
temp = evalin('caller','whos');
VARS_CONSTANTS = {temp.name};

%Need to remove '(unassigned)' -> nested functions
classValue = {temp.class};
VARS_CONSTANTS(strcmpi(classValue,'(unassigned)')) = [];

if exist('params','var')
    v = params;
elseif any(strcmp(VARS_CONSTANTS,'varargin'))
    v = evalin('caller','varargin');
else
    v = [];
end

if isempty(v)
    params = [];
elseif isstruct(v)
    params = v;
elseif isstruct(v{1}) && length(v) == 1
    params = v{1};
elseif iscell(v) && length(v) == 1 && isempty(v{1})
    params = [];
else
    params = v;
    isStr  = cellfun('isclass',v,'char');
    if ~all(isStr(1:2:end))
        error('Unexpected format for varargin, not all properties are strings')
    end
    
    if mod(length(v),2) ~= 0
        error('Property/value pairs are not balanced, length of input: %d',length(v))
    end
end

VALUES_REPLACE_CONSTANTS = params;

%NOTE: Could add on caller support, store caller's name
%and then check on get that last setter was the same
gsr_DEFINE_CONSTANTS('set',VARS_CONSTANTS,VALUES_REPLACE_CONSTANTS)



