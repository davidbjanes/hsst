function  Merged_final_detection_flag=...
    MergeMinReso_modified(final_detection_flag, final_detection_level, final_detection_value,W)

n=length (final_detection_flag);

Merged_final_detection_flag   = zeros(1,n);

    t = find(final_detection_flag==1); %% change 3
        level_temp = final_detection_level(t);
        level_temp=W(level_temp);
        new_t = [0, t] + [0, level_temp];
        a = (t - new_t(1:end-1)>0);
        b = strfind(a, [0 0 0]);
        d = b;
        if ~isempty(b)
            for i = 1 : length(b)
                if ~isempty(d)
                    c = d(1);
                    a(c+1) = 1;
                    d = strfind(a,[0 0 0]);
                end
            end
        end
        t = t(a==1);

        Merged_final_detection_flag(t)    =1;
        
