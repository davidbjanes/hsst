%% Duke Sorting Algorithm 

function [y,aligned_data] = dukeSpikeSorting(wf, varargin)

%% Input Variables
    
varargin = sanitizeVarargin(varargin);
DEFINE_CONSTANTS
    burnin = 400;                   % (1000) variable in the Bayesian statistic model
    num = 200;                      % variable in the Bayesian statistic model
    space = 1;                      % variable in the Bayesian statistic model
    Ncentres = 10;                  % initial assumed number of clusters
    K = 32;                         % Dimensionality of Waveform
    NegativePeak = 11;              % Default Alignment function Option

    ANP = 10^4;                     % Default Additive Noise Parameter 

    verbose = false;                % Outputs Raw Updates to fprintf 
    showplots = false;              % plot figures from iterator
END_DEFINE_CONSTANTS


%% Initialize Varibles

data = wf;      % Snip Waveforms

%% Alignment

% NegetivePeak varable
% when the samples of spike is 22, the most waveforms with negetive peak are around 6,  
% aligning waveforms with negetive peak at this point will keep the most information,
% the original function of "alignment_function" align the waveforms which have the negetive peak at the
% point 11 based on the previous data with length 32 by default. 
%%%--------------------------------------
%-------------------------------------------
% Maxposition = 11; % the number of samples =32;
[aligned_data,Bmatrix] = alignment_function(data, NegativePeak);% alignment of spike waveform from plexon software 

% Normalize by the mean of the maximum values across each snip
datass_train = aligned_data./mean(max(abs(aligned_data)));
datass_test{1} = data;

    
%% Tunable Params ---------------------------------------------------------
% This parameter which is very important controls noise variance,In
% general, when the spike number is less than 10000, we set the value around 10^-6
% when the spike number is more than 10000 but less than 20000,always (10^-5
% ~10^-4); when spike number is huge (> 20000),we always set it
% around(10^-2 or 10^-4);
% these are emprical value, we will analyze the relationship between noise
% variance and this parameter.
    
    TuningParameter1 = ANP;       %1*10^6;

% This parameter controls the cluster variance the range always be set
% around(1~1000);in most  case the emprical parameter is 100.
    
    %TuningParameter2 = KnobB;       %100;

    
%% Main Code
 
if verbose, fprintf('Sorting '); end
                                      
spl = Testing_DPBPSVD_Time_depedent_model_v3(datass_train, ...
                                             K,...
                                             Ncentres, ...
                                             TuningParameter1, ...
                                             burnin,...
                                             num,...
                                             space, ... 
                                             showplots);
                                      

%------------------------------------------------------
% find the most probablity sample
I = find(spl.Likelihood == max(spl.Likelihood));
% I=num;    % taking the last sample 

if verbose, fprintf(' Done. \r'); end


%% SORTING RESULTS

% Saves Data for importing back to the database
y = spl.H_z{I};

end % dukeSpikeSorting function