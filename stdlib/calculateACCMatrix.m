function [total_error, err, mapping, acc, tpr, spc] = calculateACCMatrix(sortCode_known, sortCode_est)
    
    %% Process Inputs
    if ~isrow(sortCode_known),
        sortCode_known = sortCode_known(:)';
    end
    
    if size(sortCode_est,2) ~= length(sortCode_known),
        sortCode_est = sortCode_est';
        if size(sortCode_est,2) ~= length(sortCode_known),
            warning('Data Dimension Mismatch!')
            percent_incorrect_identified = [];
            mapping = [];
            acc = [];
            return;
        end
    end
    
    debug = false;
    
    %% Code  
    
    UnitIDs_known = unique(sortCode_known);
    UnitIDs_est = unique(sortCode_est);
    
    acc = zeros([length(UnitIDs_known), length(UnitIDs_est)]); 
    tpr = zeros([length(UnitIDs_known), length(UnitIDs_est)]); 
    spc = zeros([length(UnitIDs_known), length(UnitIDs_est)]); 
    err = zeros([length(UnitIDs_known), length(UnitIDs_est)]); 
    
    for unit_id_known = UnitIDs_known,
        known_unit_ind = find(UnitIDs_known == unit_id_known);
        
        for unit_id_est = UnitIDs_est, 
            est_unit_ind = find(UnitIDs_est == unit_id_est);

            est_unit_mask = sortCode_est == unit_id_est;
            known_unit_mask = sortCode_known == unit_id_known;
            
%             relvant_unit_mask = est_unit_mask | known_unit_mask;
%             est_unit_mask(~relvant_unit_mask) = [];
%             known_unit_mask(~relvant_unit_mask) = [];
            
            TP = sum(known_unit_mask & est_unit_mask);
            FP = sum(~known_unit_mask & est_unit_mask);
            TN = sum(~known_unit_mask & ~est_unit_mask);
            FN = sum(known_unit_mask & ~est_unit_mask);
            
            ppv = TP / (TP + FP);   % positive predictive value
            npv = TN / (TN + FN);   % negative predictive value  
            
            accuracy = (TP + TN)/(TP + TN + FP + FN);
            sensitivity = TP / (TP + FN);   % True Positive Rate
            specificity = TN / (TN + FP);   % SPeCificity
            prevalence = (TP + FN)/length(sortCode_est);
            
            est_acc = sensitivity*prevalence + specificity*(1-prevalence);
            est_acc = round(est_acc*1000)./1000;
            mea_acc = round(accuracy*1000)./1000;
            if est_acc ~= mea_acc,
                ('WARNING: Measured Accuracy and Estimated Acc are not Equal!');
            end
            
            acc(known_unit_ind, est_unit_ind) = accuracy;
            tpr(known_unit_ind, est_unit_ind) = sensitivity;
            spc(known_unit_ind, est_unit_ind) = specificity;
            err(known_unit_ind, est_unit_ind) = 1 - (TP / sum(known_unit_mask));
        end
    end
    
    % K_known -> K_est(mapping)
    mapping = travelingSalesmanUniqueMapping(1 - acc');
    
    TP = zeros([1 length(UnitIDs_known)]);
    for i = 1:length(UnitIDs_known),
        unit_id_known = UnitIDs_known(i);
        map_index = mapping(i);
        
        if map_index ~= 0,
            unit_id_est = UnitIDs_est(map_index);
            temp_indices_est = sortCode_est == unit_id_est;
            temp_indices_known = sortCode_known == unit_id_known;
            TP(i) = sum( temp_indices_est & temp_indices_known );
        end        
    end
    
    total_error = 1 - (sum(TP) / length(sortCode_known));
    
    if debug == true,
        num_known_ids = length(UnitIDs_known);
        num_est_ids = length(UnitIDs_est);
        
        debug_cell = cell(num_known_ids+1, num_est_ids+1);
        debug_cell{1,1} = 'x';
        debug_cell(1,2:num_est_ids+1) = num2cell(UnitIDs_est);
        debug_cell(2:num_known_ids+1,1) = num2cell(UnitIDs_known);
        debug_cell(2:num_known_ids+1,2:num_est_ids+1) = num2cell(acc);
        debug_cell
    end

end

% Traveling Salesman Problem (only unique combinations)
function [mapping] = travelingSalesmanUniqueMapping(distance_matrix)
    % Distance Matrix is KxN (K - outputs, N - inputs, where K >= N)
    % Mapping is 1xN - mapping each input to an output
    
    size_dist = size(distance_matrix);
    if size_dist(1) < size_dist(2),
        distance_matrix = distance_matrix';
        FLIP_FLAG = true;
    else
        FLIP_FLAG = false;
    end
    
    K = size(distance_matrix,1);
    N = size(distance_matrix,2);
    
    combos = nchoosek([K:-1:1],N);
    N_fact = factorial(N);
    permuation_index = zeros(N_fact*nchoosek(K,N), N);
    for i = 1:size(combos,1),
        permuation_index(((i-1)*N_fact + 1):N_fact*i,:) = perms(combos(i,:));
    end

    index_vector = linspace(0,K*(N-1),N);
    indexing_matrix = repmat(index_vector, [size(permuation_index,1), 1]);

    permuation_matrix = bsxfun(@plus, permuation_index, indexing_matrix);

    total_perm_sum = sum(distance_matrix(permuation_matrix), 2);
    [~, ind] = min(total_perm_sum);

    mapping = permuation_index(ind, :);
        
    if FLIP_FLAG,
%         distance_matrix = distance_matrix';  % For Debugging
        flipped_map = zeros(1, K);
        flipped_map(mapping) = [1:N];
        
        mapping = flipped_map;
    end       
end

