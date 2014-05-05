
classdef DummySorterExample < hsst.sorter
    
    properties (Constant)
        sortMethodLabel = 'Dummy Sorter Example';
        defaultParameters = [1];
    end
    
    methods (Static)      
        
        function [sort_ids, wf, ts, property] = sort(wf, ts, ~, ~, ...
                                                 ~, ~, ~, ~)
            
            %% Insert New Sorting Algorithm Here
            
            % Do Something Here
            %
            %
            %
            %
            
            %% Output
            
            property = [];
            sort_ids = ones([1 length(ts)]);
            
        end
        
    end
    
end
