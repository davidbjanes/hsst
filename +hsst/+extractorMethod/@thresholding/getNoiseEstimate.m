function  [noise,info] = getNoiseEstimate(data,fs,varargin)
% iterativeNoiseEstimation - remove spikes iteratively to estimate the noise.
% Written by Matt Bauman (mbauman@gmail.com)
% 
% Algorithm:
%  0. Dropouts are removed
%  1. A normal gaussian is fit to the waveform data, estimating mu and sigma
%  2. Spikes are removed at DETECT_SIGMA multiples of sigma (default: 4)
%  3. Compute a new normal gaussian fit, re-estimating mu and sigma
%  4. If the new sigma is outside the previous confidence interval, GOTO 1
%  5. Otherwise, NOISE is returned at NOISE_SIGMA multiples of the final sigma
% 
% Optional Key-value Argument pairs
%   DETECT_SIGMA   = 4      the threshold to detect spikes, a multiple of sigma
%   SIGMA_CONVERGE = 0.05	alpha level for the confidence interval about sigma
%   NOISE_SIGMA    = 3      multiple of sigma that determines the noise level
%   TIME_RANGE     = []     seconds to use for estimate, all data when empty
%   MAX_ITERATIONS = 50     maximum number of iterations to use before bailing
%   THROW_ERROR    = false  throw errors when things go screwy
% 
% Outputs
%   NOISE - NOISE_SIGMA multiples of the final sigma estimate
%   DATA  - struct of interal data about the computation
% 
% See also: RawWaveform.estimateNoise PDecWaveform.estimateNoise

% Argument handling
% varargin = sanitizeVarargin(varargin); %#ok<NASGU>
% DEFINE_CONSTANTS
max_iterations  = 50;
sigma_converge  = 0.05;
detect_sigma    = 4;
noise_sigma     = 3;
debug_plot      = false;
time_range      = [];
throw_errors    = false;
max_remove_pct  = 66;
dropout_epsilon = 0; % microvolts
dropout_min_duration = 0.001; % seconds
% args = END_DEFINE_CONSTANTS;
% notice = @formattedWarning;
% if throw_errors
%     notice = @error; %#ok<*UNRCH>
% end

% Restrict data by time_range if requested
if isempty(time_range)
    % Use all the data
    time_range = [0 (length(data)-1)/fs];
else
    if time_range(1) < 0
        time_range(1) = 0;
    end
    if time_range(2) > (length(data)-1)/fs
        time_range(2) = (length(data)-1)/fs;
    end
    
    sample_range = floor(time_range*fs)+1;
    data = data(sample_range(1):sample_range(2));
end
original_length = length(data);

time = []; % only track time and the original in debug mode
if debug_plot
    original_data   = data;
    time = time_range(1):1/fs:time_range(2);
end

% Setup
iiIter = 1;
prev_sigmaci = [inf inf];
spikes_removed = 0;

% Algorithm
[data,dropouts] = removeDropouts(data,dropout_epsilon,round(dropout_min_duration*fs));
dropout_fraction = 1-length(data)/original_length;
while iiIter <= max_iterations
    [muhat,sigmahat,muci,sigmaci] = normfit(data, sigma_converge);
    
    if sigmaci(2) < prev_sigmaci(1)
        threshold     = detect_sigma * sigmahat;
        rawSpikeIdxs = detectSpikes(data, threshold);
        if ~isempty(rawSpikeIdxs)
            spikes_removed = spikes_removed + numel(rawSpikeIdxs);
            [data,time] = removeSpikesfromRaw(data, fs, rawSpikeIdxs, time);
        end
        iiIter = iiIter + 1;
        prev_sigmaci = sigmaci;
    else
        break
    end
end

% Setup output
noise            = noise_sigma*sigmahat;
info.noise       = noise;

info.method      = mfilename;
% info.args        = args.finalDefaultStruct;
% info.args.time_range = time_range; % This gets updated programmatically
info.date        = now;

info.iterations  = iiIter;
info.muhat       = muhat;
info.sigmahat    = sigmahat;
info.muci        = muci';
info.sigmaci     = sigmaci';
info.dropouts    = dropouts;
info.percentDropped = dropout_fraction*100;
info.spikesRemoved = spikes_removed;
info.percentRemaining = (length(data)/original_length)*100;

info.noise_wf = data;

if iiIter > max_iterations || info.percentRemaining < (100-max_remove_pct)
    warning('Removed %.1f%% of the waveform after %d iterations',...
        100-info.percentRemaining,iiIter)
    noise = nan;
    info.noise = nan;
end

if debug_plot
    figure %#ok<UNRCH>
    plot(time_range(1):1/fs:time_range(2),original_data);
    line(time,data,'Color',[.5 .5 .5]);
    hold on
    line(xlim,  noise_sigma*sigmahat(end)*[1 1],'Color','red','LineWidth',2);
    line(xlim, -noise_sigma*sigmahat(end)*[1 1],'Color','red','LineWidth',2);
    line(xlim,  noise_sigma*sigmahat(1)*[1 1],  'Color','g','LineWidth',2);
    line(xlim, -noise_sigma*sigmahat(1)*[1 1],  'Color','g','LineWidth',2);
end
end

function [data, count] = removeDropouts(data,eps,dur)
% Remove chunks of data where it is near zero for extended periods of time
if nargin == 1
    eps = 1.5; % microvolts, 9 possible values in Plexon between -1.5 and 1.5 uV
    dur = 40;  % samples, 1 ms at 40kHz
else
    assert(nargin == 3);
end

nearzero = abs(data(:)) < eps;

diffdrop = [nearzero(1); diff(nearzero)];
% Diffdrop is a vector of 0,1,-1.
% I want segments where 1 and -1 are separated by dur (or more) idxs.
% Prefix the diff by the first element so that the first nonzero element
% demarcates the beginning of a run of nearzero values
idxs = [find(diffdrop); length(diffdrop)+1];
runlengths = diff(idxs);
% runlengths contains both lengths in and out of nearzero, starting with the
% first run of nearzero values. All even runlengths are not nearzero
runlengths(2:2:end) = 0; % set the even runlengths to zero

remove_idxs = idxs(runlengths >= dur);
remove_lens = runlengths(runlengths >= dur);

mask = false(size(data));
for i = 1:length(remove_idxs)
    idx = remove_idxs(i);
    len = remove_lens(i);
    mask(idx:idx+len-1) = true;
end

data(mask) = [];
count = length(remove_idxs);
end

function [rawData,time] = removeSpikesfromRaw(rawData, rawFs, spikeIdxs, time)

remove_window = [-0.0003 0.0007]; % The amount of data to remove in seconds
remove_samples = floor(remove_window(1)*rawFs):floor(remove_window(2)*rawFs);

remove = bsxfun(@plus,remove_samples,spikeIdxs(:));

% Ensure that the indices to remove are within the bounds of rawNoise
remove(remove < 1 | remove > length(rawData)) = [];
rawData(remove) = [];

if ~isempty(time)
    time(remove) = [];
end

end

function spike_mask = detectSpikes(data, threshold)
%Look for positive and negative going spikes.  The spikeTimes vector that
%gets returned will be the one that with the most spikes in it.

posData = data > threshold;
negData = data < -threshold;

%Find positive spikes
posThresholds = find(diff(posData) == 1);
posThresholds = posThresholds + 1;

%Find negative spikes
negThresholds = find(diff(negData) == 1);
negThresholds = negThresholds + 1;

if length(posThresholds) >= length(negThresholds)
    spike_mask = posThresholds;
else
    spike_mask = negThresholds;
end
end
