function args = sanitizeVarargin(args,argFilter)
% SANITIZEVARARGIN Format varargins just the way you like them!
%
%   args = sanitizeVarargin(args,*argFilter)
%
%   - sets all varargin names lower case
%   - replaces all spaces with the underscore character '_'
%   - allows the user to filter out varargins
%
% INPUTS
% =========================================================================
%   args - (cell) varargin as passed to a function
%   argFilter - (cell) list of varargins to be removed.
%

args(1:2:end) = regexprep(lower(args(1:2:end)),'\s+','_');

if nargin > 1 && ~isempty(argFilter)
    argFilter  = regexprep(argFilter,'\s+','_');
    [mask,idx] = ismember(argFilter,args(1:2:end));
    if any(mask)
        idx = [idx(mask) idx(mask)+1];
        args(idx) = [];
    end
end

end