
classdef optimizer < handle
        
    properties (Access = public)
        verbose = false;
    end
    
    properties (Access = private)
        dataObj                 % dataObject instance with data vars
        SoMObj                  % Sort Method Object
        ScMObj                  % Score Method Object
        optimalParameter_ind    % Optimal Parameter Selection Index
    end
    
    
    %% Public Methods
    methods (Access = public)
        
        % Constructor
        function obj = optimizer(dataObj, sortMethodObj, scoreMethodObj, sortParameters)
            obj.dataObj         = dataObj;
            obj.SoMObj          = sortMethodObj;
            obj.ScMObj          = scoreMethodObj;
                        
            if isempty(sortParameters)
                sortParameters = eval(sprintf('%s.defaultParameters',class(sortMethodObj)));
            end
            
            generateParameterScores(obj, sortParameters);
           
            scoreObj_list = obj.dataObj.sortScoreObjList;
            obj.optimalParameter_ind = obj.chooseOptimalParameter( scoreObj_list );
            
        end
               
        % Getter: Optimal Parameter
        function [optimalParam_val, optimalParam_ind] = returnOptimalParameter(this)
            
            optimalParam_ind = this.optimalParameter_ind;
            optimalParam_val = this.dataObj.sortParameterList{optimalParam_ind};

        end
        
        % Getter: Optimal SortCode (determined by optimal parameter)
        function sortCode = returnOptimalSortCode(this)
                    
            [~,ind] = this.returnOptimalParameter();
            sortCode = this.dataObj.sortCodeList{ind};
            
        end
        
        % Getter: Optimal Score Object (determined by optimal parameter)
        function scoreObj = returnOptimalScoreObj(this)
                    
            [~,ind] = this.returnOptimalParameter();
            scoreObj = this.resultObj(ind).scoreObj;
            
        end
        
        % Getter: All Results
        function resultObj = returnResultObj(this)
                    
            resultObj = this.resultObj;
            
        end
            
    end
    
    
    %% Viewable to all
    methods (Static)
        
        % Pick Optimal Sort from List of ScoreObjects
        function [ index ] = chooseOptimalParameter( scoreObj_list )

            N = length(scoreObj_list);

            score_list = zeros([1 N]);
            for n = 1:length(scoreObj_list)
                score_list(n) = scoreObj_list{n}.score;
            end

            index = find(score_list == max(score_list),1,'last');

        end
        
        function [ sortCode ] = cleanUpSortCode( sortCode )

            MIN_NUM_OF_WAVEFORMS = round(length(sortCode)*0.05);
            
            table = tabulate(sortCode);
            num_occur = table(:,2);
            unitIDs = table(:,1);
            
            if length(unitIDs) > 6,
                unitIDs_to_remove = unitIDs(num_occur < MIN_NUM_OF_WAVEFORMS);

                sortCode_ind_to_remove = ismember(sortCode,unitIDs_to_remove);

                sortCode(sortCode_ind_to_remove) = 0;
            end
        end
        
    end
    
    
    %% Private Helper Methods
    methods (Hidden, Access = private)
        
        generateParameterScores(this, sortParameters)
               
    end
    
end

