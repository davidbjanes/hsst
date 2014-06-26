function [tf,loc] = ismember_str(input_1,input_2,isUnique)
%ismember_str  Generally faster version of ismember for strings
%
%   [tf,loc] = ismember_str(input_1, input_2, *isUnique)
%
%   INPUTS
%   =================================================
%   input_1  : (cellstr or char)
%   input_2  : (cellstr or char)
%   isUnique : (default false), if true indicates that both inputs are
%              unique, if this is true it can speed up the code a bit
%
%   OUTPUTS
%   =================================================
%   tf  : whether or not each element of input_1 is in input_2
%   loc : location of each element of input_1 in input_2, is equal to zero
%         for no match and gives the location of the last match when a 
%         match is present

%POSSIBLE IMPROVEMENTS
%------------------------------------------------------------------
%1) Speedup sort of cell array of strings, not sure why it is so slow 
%       compared to sort rows
%2) Write function that on sort also tests neighbor equality, this has to
%   be done already, just not returned :/, then we would essentially be
%   done
%3) Write a sort function with two inputs, instead of concatenating
%4) Sort function with 2 inputs that checks if both are already sorted
%   and then allows for a merge sort

%I didn't make these optional to try and speed up general usage

SEARCH_LIMIT_RATIO_MAX = 0.1; %If one of the inputs is significantly shorter
%than the other, then we search for the smaller set and only strcmp on that
%set, ignoring a lot of string comparisons between members of the same
%group
SEARCH_DIRECT_MAX = 5; %If one of the inputs is this size, then a direct search
%for each element is employed

nIn = nargin;

if nIn == 2
    isUnique = false;
elseif nIn == 3
    if ischar(isUnique)
        isUnique = strcmp(isUnique,'true');
    end
else
    error('Incorrect # of inputs')
end

if ischar(input_1)
    input_1 = cellstr(input_1); 
end
if ischar(input_2)
    input_2 = cellstr(input_2); 
end

%We'll skip this for now
%This is if we wanted to submit this to Matlab as a substitute for
%cell.ismember
% if ~(iscellstr(input_1) && iscellstr(input_2))
%     error('Incorrect type of inputs to function')
% end

%Some input processing
%------------------------------------
L1 = length(input_1);
L2 = length(input_2);

if L1 ~= numel(input_1) || L2 ~= numel(input_2)
    error('Function only works on 1d arrays')
end

if isempty(input_1) && isempty(input_2)
    tf  = logical([]); % indefinite result.
    loc = []; % indefinite result.
    return
elseif L1 <= SEARCH_DIRECT_MAX  %Any restrictions on L2????
    %COULD IMPROVE THIS AND NEXT IF WE KNEW INPUTS WERE SORTED
    loc = zeros(1,L1);
    tf  = false(1,L1);
    for iTF = 1:L1
        temp_I = find(strcmp(input_1{iTF},input_2),1,'last');
        if ~isempty(temp_I)
            loc(iTF) = temp_I;
            tf(iTF)  = true;
        end
    end
    return
elseif L2 <= SEARCH_DIRECT_MAX
    loc = zeros(1,L1);
    tf  = false(1,L1);
    for iTF = 1:L2
        temp_I = find(strcmp(input_2{iTF},input_1));
        if ~isempty(temp_I)
            loc(temp_I) = iTF;
            tf(temp_I)  = true;
        end
    end
    return
end

%We need to sort, make rows
if size(input_1,1) > 1, input_1 = input_1'; end
if size(input_2,1) > 1, input_2 = input_2'; end


%CASE 1: We are not sure if the inputs are unique
%NOTE: I think checking this is and using the other technique isn't as
%good as
if ~isUnique %Default case
    
    L1 = length(input_1);
    [s_inputs,I]   = sort([input_1 input_2]);
    
    matchesNextGuy = strcmp(s_inputs(1:end-1),s_inputs(2:end));
    transZones     = find(diff([false matchesNextGuy false] ~= 0));
    
    Istart = transZones(1:2:end);
    Iend   = transZones(2:2:end);
    
    
    %Potetial room for improvement
% %     if all(Iend - Istart == 2)
% %         %keyboard
% %     end
    
    
    if nargout == 1
        badRegions = find(I(Istart) > L1 | I(Iend) <= L1);
    else
        %This construct drives me nuts, why can't we split
        badMask     = I(Istart) > L1 | I(Iend) <= L1;
        badRegions  = find(badMask);
        goodRegions = find(~badMask);
    end
    
    %Remove all cases in which matches were only from the same input
    for iBad = badRegions(:)'
        matchesNextGuy(Istart(iBad):Iend(iBad)) = false;
    end
    
    %We won't grab matches from input_2
    matchesNextGuy(I > L1) = false;
    
    tf = false(1,length(input_1)); %Initialize output
    if nargout == 1
        tf(I(matchesNextGuy)) = true;
    else
        %might also be able to speedup with search on mask
        %Should put check in place here ..., but really, more than
        %4294967295 elements ...
        locTemp = zeros(1,length(s_inputs),'uint32'); %trying to save on memory
        for iGood = goodRegions(:)'
            locTemp(Istart(iGood):Iend(iGood)) = I(Iend(iGood))-L1;
        end
        matchesNextGuy_I = find(matchesNextGuy);
        loc = zeros(1,L1);
        loc(I(matchesNextGuy_I)) = locTemp(matchesNextGuy_I);
        tf(I(matchesNextGuy_I))  = true;
    end
    
%UNIQUE SORT
%==========================================================================
else
    [s_inputs,I] = sort([input_1 input_2]);
    tf = false(1,L1);
    
    %EXAMINE DOING LESS STRING COMPARISONS
    if L1/L2 <= SEARCH_LIMIT_RATIO_MAX || L2/L1 <= SEARCH_LIMIT_RATIO_MAX
        if L1 < L2
            I_1 = find(I(1:end-1) <= L1);
        else
            I_1 = find(I(2:end) > L1);
            %This is a bit tricky, only look for I_1 where the guy 
            %to the right is from input_2
        end
        if nargout == 1
            tf(I(I_1(strcmp(s_inputs(I_1),s_inputs(I_1+1))))) = true;
        else
            isMatch = I_1(strcmp(s_inputs(I_1),s_inputs(I_1+1)));
            tf(I(isMatch)) = true;
            loc = zeros(1,L1);
            loc(tf) = I(isMatch+1) - L1;
        end
    else
        if nargout == 1
            tf(I(strcmp(s_inputs(1:end-1),s_inputs(2:end)))) = true;
        else
            isMatch = strcmp(s_inputs(1:end-1),s_inputs(2:end)); %again we are short
            tf(I(isMatch)) = true;
            loc = zeros(1,length(input_1));
            loc(tf) = I([false isMatch]) - length(input_1); %the [false tf] shifts everything
            %to the right to get the index from the next guy
        end
    end
end

%If we want location we'll need to hang onto the strcmp ...