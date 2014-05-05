close all

cat_name  = 'Erol';
block     = 144;
channel   = 65;
tr        = getDBObj('trial',cat_name,block);
if isprop(tr,PDecData.listAs) && isempty(tr.pdecData.wf(1).noiseEstimate)
    estimateNoise(tr.pdecData);
end
spike_data = tr.spikeData(1);
operators = {SNROperator(1.2) AlignmentOperator(3) VPPOperator(0,5e2)};
setOperators(spike_data,operators);
% plotRecordingQualitySummary(spike_data,'use operators',true)

Fs       = spike_data.fs;

spikes = ss_default_params(Fs);
% HANDJAM =================================================================
if ~isprop(tr,RawData.listAs)
    [ts,spike_mask] = getSpikesUsingOperators(spike_data.st(channel));
    %     wf       = spike_data.st(channel).snips.wf(:,spike_mask);
    wf       = alignment_function(spike_data.st(channel).snips.wf(:,spike_mask),11);
    nSnips   = size(wf,2);
    nSamples = size(wf,1);

    spikes.waveforms  = wf';
    spikes.spiketimes = ts;
    spikes.trials     = ones(1,nSnips);
    
    spikes.params.max_jitter   = 3/Fs*1e3;
    spikes.params.window_size  = (nSamples)/Fs*1e3;   % ms, width of a spike
    spikes.params.shadow       = 2/Fs*1e3;          % ms, enforced dead region after each spike
    spikes.params.cross_time   = 0/Fs*1e3;         % ms, alignment point for peak of waveform
    
    % spikes.info.detect.stds          = 28.8141
    spikes.info.detect.thresh        = spike_data.st(channel).threshold;
    spikes.info.detect.event_channel = ones(nSnips,1);
    spikes.info.detect.dur           = tr.duration;
    spikes.info.detect.align_sample  = spike_data.st(channel).snips.nSamplesPreThreshold+1;
    spikes.unwrapped_times           = ts;
%     spikes = ss_align(spikes);
else
    spikes = ss_detect({tr.rawData.wf(channel).data},spikes);
    spikes = ss_align(spikes);
end
% SORT =================================================================

[pca.u,pca.s,pca.v] = svd(detrend(spikes.waveforms(:,:),'constant'), 0);             % SVD the data matrix
spikes.info.pca = pca;
spikes = ss_kmeans(spikes);
spikes = ss_energy(spikes);

dim_red       = DimensionalityReduction;
[~, coeff, ~] = PCA(dim_red, spikes.waveforms');
    
agg_values = spikes.params.agg_cutoff*[7.5e-4 1e-3 .01 .1 .25 .5 1 1.1 2 5 10 50 100];
agg_values = spikes.params.agg_cutoff*[1:10 50 100]*1e-4;

[hfig,haxis] = PLOT_layoutGridByN(length(agg_values), ...
    'axis style',   'outer', ...
    'labels',       arrayfun(@num2str,agg_values,'un',0));
haxis = haxis{1};
for ii = 1:length(agg_values)
    parent = haxis(ii);
    spikes_copy = spikes;
    spikes_copy.params.agg_cutoff = agg_values(ii);
    spikes_copy = ss_aggregate(spikes_copy);

    [uId,~,iiB] = unique_2011b(spikes_copy.assigns);
    nUnits = length(uId)
    
    if nUnits < 7
        colors = lines(length(uId));
    else
        colors = hsv(length(uId));
    end
    
    hold(parent,'all')
    for ii = 1:length(uId)
        mask = ii == iiB;
        plot(coeff(mask,1),coeff(mask,2),'.','color',colors(ii,:),'parent',parent)
    end
end
    
adjustAllAxisLimits

spikes.params.agg_cutoff = agg_values(1);
spikes= ss_aggregate(spikes);
% splitmerge_tool(spikes)

[uId,~,iiB] = unique_2011b(spikes.assigns);
nUnits = length(uId);

if nUnits < 7
    colors = lines(length(uId));
else
    colors = hsv(length(uId));
end

[hfig,haxis] = PLOT_layoutGridByN(nUnits,'axis style','outer');
haxis = haxis{1};
for ii = 1:length(uId)
    parent = haxis(ii);
    hold(parent,'on')
    mask = ii == iiB;
    plot(spikes.waveforms(mask,:)', ...
        'parent',parent, ...
        'color',colors(ii,:))
end
adjustAllAxisLimits

dim_red       = DimensionalityReduction;
[~, coeff, ~] = PCA(dim_red, spikes.waveforms');

figure
hold all
for ii = 1:length(uId)
    mask = ii == iiB;
    plot(coeff(mask,1),coeff(mask,2),'.','color',colors(ii,:))
end
adjustAllAxisLimits

sort_ids = zeros(size(spike_mask),'int8');
sort_ids(spike_mask) = iiB;
saveSortCodes(tr.spikeData(1), ...
    'UltraMegaSort2000', ...
    {sort_ids}, ...
    channel, ...
    'force_replace', true);
plot(tr.spikeData(1).st(channel).snips, ...
    'SortCode','UltraMegaSort2000',...
    'align method','global_min')