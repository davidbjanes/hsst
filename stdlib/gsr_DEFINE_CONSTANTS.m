function varargout = gsr_DEFINE_CONSTANTS(option,variablesAtStart,params)
%gsr_DEFINE_CONSTANTS Holds variables for END_DEFINE_CONSTANTS
%
%   USAGE FORMS
%   ================================================================
%   called by DEFINE_CONSTANTS:
%       gsr_DEFINE_CONSTANTS('set',variablesAtStart,params)
%
%   called by END_DEFINE_CONSTANTS:
%       [variablesAtStart,params] = gsr_DEFINE_CONSTANTS('get')
%       
%       gsr_DEFINE_CONSTANTS() - clears variables in memory
%
%   INPUTS/OUTPUTS
%   ============================================================
%   variablesAtStart : (cellstr) list of variables in place during call
%                       to DEFINE_CONSTANTS, for later exclusion when
%                       figuring out what variable names are constants
%   params : (param/value pairs or structure), cleaned by DEFINE_CONSTANTS
%             but parsed by END_DEFINE_CONSTANTS
%
%   NOTE: This approach was chosen instead of creating a global variable as
%   it was thought that it made the sharing between variables more obvious
%   and was something that could easily be checked for conflicts. That
%   being said it might not be the best way of handling this ...

%NOTE: Could make this more complex to support file name filtering
%i.e. store these variables based on the name of the function calling them

persistent VARS_AT_START
persistent PARAMS_REPLACE

%NOTE: Since this is only called by two functions we should
%be able to do little checking on this
switch option(1)
    case 's'
        VARS_AT_START  = variablesAtStart;
        PARAMS_REPLACE = params;
    case 'g'
        varargout{1} = VARS_AT_START;
        varargout{2} = PARAMS_REPLACE;
    otherwise
        VARS_AT_START = {};
        PARAMS_REPLACE = {};
end