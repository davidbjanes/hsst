function [ boolean_output, message ] = isObjectScoreMethod( object )

    %% Check ScoreMethodObj is correct class type (if it was an input)    
    scoreMethodList = hsst.getScoreMethods();
    superClassScorerNames = superclasses(class(scoreMethodList{1}));
    
    inputScorer_superClassNames = superclasses(class(object));
    if ~all(ismember(superClassScorerNames, inputScorer_superClassNames)),
        message = 'Invalid Scoring Class Object';
        boolean_output = false;
    else
        message = 'Valid Scoring Class Object';
        boolean_output = true;
    end  
    
end

