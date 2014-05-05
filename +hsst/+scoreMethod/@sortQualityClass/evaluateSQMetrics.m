
% Metric Evaluation
function [this] = evaluateSQMetrics(this)

    bool_output = zeros(this.num_units, this.num_of_metrics);
    raw_output  = {};
    for i = 1:this.num_of_metrics,
        [bool_output(:,i), raw_output{i}] = ...
            this.METRICS{i}.scoreAllUnits(this.property_struct.wf, ...
                                          this.property_struct.ts, ...
                                          this.noiseEstimate, ...
                                          this.sort_IDs, ...
                                          this.property_struct);
    end       
    this.final_scores = boolean(bool_output);
    this.raw_values = raw_output;
    
    % DPrime Dissimilar Peaks
    dissi_pks = this.final_scores(:,5);
    for cur_pk_ind = find(dissi_pks == 0)',
        similar_pk_ind = raw_output{5}{cur_pk_ind};
        all_ind = [similar_pk_ind, cur_pk_ind];
        [~, ind] = max(this.num_wf(all_ind));
        if ind == length(all_ind),
            dissi_pks(cur_pk_ind) = 1;
        end               
    end
    this.final_scores(:,5) = dissi_pks;
    
    
    % DPrime Dissimilar Waveforms
    dprime_dissim = this.final_scores(:,6);
    for cur_pk_ind = find(dprime_dissim == 0)',
        similar_pk_ind = raw_output{6}{cur_pk_ind};
        all_ind = [similar_pk_ind, cur_pk_ind];
        [~, ind] = max(this.num_wf(all_ind));
        if ind == length(all_ind),
            dprime_dissim(cur_pk_ind) = 1;
        end               
    end
    this.final_scores(:,6) = dprime_dissim;
    
end