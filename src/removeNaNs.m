function outputMatrix = removeNaNs(inputMatrix)
%***********************************************************************
% function removes the NaN values from the input Matrix and sets to 
% zero values
%
%***********************************************************************
    bNaNMask = isnan(inputMatrix);
    
    outputMatrix = inputMatrix;
    outputMatrix(bNaNMask) = 0;
end