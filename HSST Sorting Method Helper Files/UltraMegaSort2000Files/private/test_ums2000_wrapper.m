sortParameters          = 5e-6*[.01 .1 1:10 100 1e3 1e4 1e5 1e6]*1e-4;
cat_name  = 'Erol';
sort_name = 'UltraMegaSort2000_rethreshold';
operators = {VPPOperator(50,700) SNROperator(1.2) AlignmentOperator(3) BlankingOperator};
exp = getDBObj('exp','Erol');
blocks = exp.trial.id;
blocks = [blocks{:}];
blocks = 290;
for iiBlock = 1:length(blocks)
    try
        this_block  = blocks(iiBlock);
        fprintf('%d: %d/%d\n',this_block,iiBlock,length(blocks));
        tr        = getDBObj('trial',cat_name,this_block);
        spike_data = tr.spikeData(1);
        channels = spike_data.st.channel;
        channels = [channels{:}];
        if isprop(tr,PDecData.listAs) && isempty(tr.pdecData.wf(1).noiseEstimate)
            estimateNoise(tr.pdecData);
        end
        setOperators(spike_data,operators);
        
        % [~,spike_mask] = getSpikesUsingOperators(spike_data.st(channels));
        % plotRecordingQualitySummary(spike_data.st(channels), ...
        %     'snip_align_method','')
        
        tic
        sorts = UltraMegaSort2000_wrapper(cat_name, this_block, channels, ...
            'force_score',          false,  ...
            'force_replace',        false,   ...
            'save_to_qnap',         true,   ...
            'sort name',            sort_name,      ...
            'sortParameters',       sortParameters, ...
            'use_operators',        true);
        toc
        sort_mask  = arrayfun(@(x)any(strcmp(sort_name,x.sortCodes)),tr.spikeData);
        spike_data = tr.spikeData(sort_mask);

        if length(channels) == 1
            nUnits = cellfun(@(x)length(unique(x)),sorts.sortCodes{1}(1,:));
            
            [nUnits; [sorts.sortScores{1}.score]]
            plot(spike_data.st(channels).snips,...
                'align method', 'global_min',...
                'sortcode',     sort_name);
            [~,idx]         = max(nUnits);
            % [~,spike_mask]  = getSpikesUsingOperators(spike_data.st(channels));
            spike_mask = sorts.spike_mask{1};
            sort_ids  = zeros(length(spike_mask),1,'uint8');
            [~,~,iiB] = unique_2011b(sorts.sortCodes{1}{1,idx});
            sort_ids(spike_mask) = iiB;
            plot(spike_data.st(channels).snips, ...
                'ClusterIDX',   sort_ids,...
                'Align Method', 'global_min', ...
                'Spike Mask',   spike_mask);
            
            wf       = spike_data.st(channels).snips.wf;
            wf       = wf(:,spike_mask);
            align_wf = alignment_function(wf,10);
            
            %     figure;
            %     subplot(1,2,1);
            %     plot(wf);
            %     subplot(1,2,2);
            %     plot(align_wf);
        else
            plotRecordingQualitySummary(spike_data, ...
                'use operators',        true,...
                'recalculate_noise',    true,...
                'sort_name',            sort_name)
        end
        close(tr)
    catch ME
        simpleExceptionDisplay(ME)
    end
end