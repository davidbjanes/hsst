function [ this ] = calculateFinalScore( this )

    if sum(this.num_wf > 100) > this.max_num_units
       this.score = 0;
       return
    end

    % Find Noise Units
    snr_val   = this.final_scores(:,1)';
    sing_pk   = this.final_scores(:,2)';
    isi_viol  = this.final_scores(:,3)';
    exp_Fit   = this.final_scores(:,4)';
    dissi_pks = this.final_scores(:,5)';
    dprime    = this.final_scores(:,6)';
    xCoor_val = this.final_scores(:,7)';
    temp_mat  = this.final_scores(:,8)';
    bimodal_pk= this.final_scores(:,9)';
    th_slope  = this.final_scores(:,10)';

    this.oversorted_units  = ~dissi_pks | ~dprime;
    this.undersorted_units = ~xCoor_val | ~temp_mat | ~bimodal_pk | ... 
                             ~th_slope  | ~isi_viol;

    this.noise_units = (~snr_val | ~exp_Fit | ~sing_pk | ~isi_viol) ...
                       & ~(~xCoor_val | ~temp_mat | ~bimodal_pk);


    % Unit Scores
    unit_scores = mean([~this.noise_units; ...
                        ~this.undersorted_units; ...
                        ~this.oversorted_units]);              

                    
    % If only one unit exists, the final score is that unit's score
    if length(this.unitIDs) == 1,
        score_percentage = mean(unit_scores);

    else
        if sum(this.noise_units) > 0,
            [~, max_wf_noise_ind] = max(this.num_wf(this.noise_units));           
            noise_unit_ind = find(this.noise_units > 0);
            unit_scores(noise_unit_ind(max_wf_noise_ind)) = 1;
        end

        num_waveforms = this.num_wf;
        total_num_waveforms = sum(num_waveforms);

        if this.normByWf,
            score_percentage = (unit_scores * num_waveforms');
            score_percentage = score_percentage / total_num_waveforms;
        else
            score_percentage = mean(unit_scores);
        end
    end

    
    this.score = score_percentage;
    if this.property_struct.verbose, fprintf('Sort Score: %2.1f%%\r', mean(this.score)*100); end


    %% Old Version
    % Calcuate score of each unit
%             total_sum = [];
%             for unit = this.unitIDs,
%                 total_sum = [total_sum, this.getUnitScore(unit)];
%             end
%             
%             % If only one unit exists, the final score is that unit's score
%             if length(this.unitIDs) == 1,
%                 score_percentage = total_sum; 
%                 
%             % if more than one unit exist
%             else
%                  
%                 % Find lowest scoring "Noise Unit"
%                 [~, ind] = min(total_sum(this.noise_units));
%                 total_sum_ind = find(this.noise_units == 1);
%                 
%                 noise_unit_mask = boolean(zeros([1 length(this.unitIDs)]));
%                 noise_unit_mask(total_sum_ind(ind)) = 1;
%                                
%                 
%                 % Remove it from the total_sum
%                 total_sum(noise_unit_mask) = [];
% 
%                 if this.normByWf,
%                     weighting_vector = this.num_wf(~noise_unit_mask);
%                     total_sum = total_sum * weighting_vector';
%                     score_percentage = total_sum / sum(this.num_wf(~noise_unit_mask));
% 
%                 else
%                     score_percentage = mean(total_sum); 
%                 end    
%             end
% 
%     this.score = score_percentage;
%     if this.property_struct.verbose, fprintf('Sort Score: %2.1f%%\r', mean(this.score)*100); end

end

