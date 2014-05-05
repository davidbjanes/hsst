function [ boolean_output, message ] = isObjectExtractorMethod( object )

    %% Check ScoreMethodObj is correct class type (if it was an input)    
    extractorMethodList = hsst.getExtractorMethods();
    superClassExtractorNames = superclasses(class(extractorMethodList{1}));
    
    inputExtractor_superClassNames = superclasses(class(object));
    if ~all(ismember(superClassExtractorNames, inputExtractor_superClassNames)),
        message = 'Invalid Extractor Class Object';
        boolean_output = false;
    else
        message = 'Valid Extractor Class Object';
        boolean_output = true;
    end  
    
end