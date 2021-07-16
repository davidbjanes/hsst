classdef KMeans < hsst.sorter
    
    properties (Constant)
        sortMethodLabel = 'K-Means';
        defaultParameters = [1:6];
    end
    
    methods (Static)
        
        function [sort_ids, wf, ts, property] = sort(wf, ts, ~, thresh, ...
                                                 ~, ~, ~, ...
                                                 sortparameter)
            
            % Handle Inputs
            if size(wf, 1) ~= length(ts),
                wf = wf';
                if size(wf, 1) ~= length(ts),
                    error('Time Stamp and Number of Waveform Mismatch');
                end
            end
            
            K = sortparameter;
            debug = false;
            showplots = false;
            V_thresh = thresh;
                        
            num_snips = size(wf, 1);
            num_samples = size(wf, 2);


            %% Method 1 - determine starting cluster location
%             u_k(:,1) = wf(randi(num_snips), :);
%             u_k(:,1) = wf(10, :);
%             for n = 1:num_snips,
%                  dist(n) = norm(wf(n,:) - u_k(:,1)');        
%             end
% 
%             [val, ind] = sort(dist);
%         %     TODO: initalize u_k according to distance distrb.
%         %     dist_distribution = smooth(diff(dist(ind)),50);
%             for k = 2:K,
%                 u_k(:,k) = wf( ind(floor(num_snips/K)*(k-1)) ,:);
%             end
            
            %% Method 2 - initalize starting clusters
            % Initalize U_K
            dist_matrix = squareform(pdist(wf));

            val = max(max(dist_matrix));
            [~, init_indices] = find(dist_matrix == val);
            chosen_indices = 1:num_snips == init_indices(1);

            for k = 2:K,
                total_dist_sum = sum(dist_matrix(chosen_indices, :), 1);
                total_dist_sum(chosen_indices) = 0;

                [~, next_index] = max(total_dist_sum);
                chosen_indices(next_index) = 1;
            end

            [COEFF, score] = pca(wf);
            mean_waveform = mean(wf);
    
            low_d = score(chosen_indices, :);
            low_d(:,3:31) = 0;
            projected = bsxfun(@plus, low_d/COEFF, mean_waveform);
            u_k = projected';
%             u_k = wf(chosen_indices, :)';

            
            %% Sort Algorithm
            % Parameter Initalization
            J = [];
            eps = 1e-20;
            r = [];

            
            % Sort!
            while length(J) < 2 || abs(J(length(J)) - J(length(J)-1)) > eps,

                % E-Step
                for i = 1:num_snips,
                    dist = [];
                    for k = 1:K,
                        dist(k) = norm(wf(i,:) - u_k(:,k)');        
                    end

                    [~,ind] = min(dist);

                    r(i) = ind;
                end

                % M-Step
                for k = 1:K,
                    u_k(:,k) = mean(wf(r == k, :));
                end

                % Calculate J
                nxt_step_J = 0;
                for i = 1:num_snips,
                    nxt_step_J = nxt_step_J + norm(wf(i,:) - u_k(:,r(i))')^2;        
                end

                J(length(J)+1) = nxt_step_J;


                % Plot
                if showplots,
                    figure
                        hold on
                        plot([1:size(wf,2)], ones([size(wf,2) 1])*V_thresh, 'r:') 
                        format = {'r', 'b', 'c', 'g', 'k', 'm'};
                        for k = 1:K,
                            plot([1:size(wf,2)], wf(r == k, :)', format{k})
                        end

                        for k = 1:K,
                            plot([1:size(wf,2)], u_k(:,k), '*', 'Color', format{k+3}, 'LineWidth', 3)
                        end
                        xlabel('Sample Number')
                        ylabel('Amplitude [uV]')
                        title('Voltage vs Time of Snips')

                    figure
                        plot(J, '*-')
                        xlabel('Time Step Number')
                        ylabel('Value of J')
                        title('Inital U_k = InitThreeClusters_2');
                end

                if debug, 
                    plotWaveforms(snips, r, V_thresh, u_k)
                end

            end
                
            property = J(end);
            sort_ids = r;
            
        end
        
    end
    
end
