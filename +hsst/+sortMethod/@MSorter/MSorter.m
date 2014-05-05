classdef MSorter < hsst.sorter
    
    properties (Constant)
        sortMethodLabel = 'M-Sorter';
        defaultParameters = [1:6];
    end
    
    methods (Static)
        
        function [sort_ids, wf, ts, property] = sort(wf, ts, ~, thresh, ...
                                                 ~, ~, ~, ...
                                                 sortparameter)
                    
            sort_ids = MSorter_DavidMod(wf, thresh, sortparameter);
            property = [];
            
        end
        
    end
    
end