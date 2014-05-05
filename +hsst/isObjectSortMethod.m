function [ boolean_output, message ] = isObjectSortMethod( object )

    %% Check ScoreMethodObj is correct class type (if it was an input)    
    sortMethodList = hsst.getSortMethods();
    superClassSorterNames = superclasses(class(sortMethodList{1}));
    
    inputSorter_superClassNames = superclasses(class(object));
    if ~all(ismember(superClassSorterNames, inputSorter_superClassNames)),
        message = 'Invalid Sorting Class Object';
        boolean_output = false;
    else
        message = 'Valid Sorting Class Object';
        boolean_output = true;
    end  
    
end

