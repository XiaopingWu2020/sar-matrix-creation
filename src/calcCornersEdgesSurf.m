function [corners, edges, surfaces, lNoOfInd] = calcCornersEdgesSurf(lRadius)
%%**********************************************************************
%
%   calculates the corners of cubes with growing diameter. E.g. for a 3x3x3
%   cube the corners are [1 1 1] [1 1 -1] ....
%   for a 5x5x5 cube the corners are [2 2 2] [2 2 -2]... 
%
%
%   INPUT:                                                  [unit]
%   ----------------------------------------------------------------
%   corners:  3D Matrix
%   edges:
%   surfaces
%  
%
%   OUTPUT:
%   ----------------------------------------------------------------
%   corners:  3D Matrix with coordinates of corners. 
%             1st dimension: radius of the shell, e.g. for 3x3x3 matrix 1st 
%                            dim can be 0 and 1; 
%             2nd dimension: counts all possible corners. for corners this number is
%                            always 8 (except for l=1, here corners == 1). 
%             3rd dimension: gives the coordinate. e.g. [ 1 1 1]
%   edges:    same handling as corners. 2nd dim is i.g. larger
%   surfaces  same handling as corners. 2nd dim is i.g. larger
%   lNoOfInd  2D Matrix with size [lRadius x 3] Gives for each radius the
%             number of corner points, edge points and surface points
%
%%**********************************************************************

    [traj, trajval] = buildSARAveTrajChessb(lRadius*2+1);
    trajvaluniq = unique(trajval);
    
    shell = zeros(1,3);
    
    %loop over the radius
    for lL=0:max(trajvaluniq(:))
        shellmask   = (trajval == lL);
        shell       = traj(shellmask,:);
        shabs       = abs(shell);
        %find corners:
        cornind     = (shabs(:,1) == lL) & (shabs(:,2) == lL) & (shabs(:,3) == lL);
        
        edgeind     = ((shabs(:,1) == lL) & (shabs(:,2) == lL)) | ...
                      ((shabs(:,2) == lL) & (shabs(:,3) == lL)) | ...
                      ((shabs(:,1) == lL) & (shabs(:,3) == lL));
        
        corntmp     = shell(cornind,:);
        cs          = size(corntmp);
        edgetmp     = shell(edgeind & ~cornind,:);
        es          = size(edgetmp);
        surftmp     = shell(~edgeind & ~cornind,:);
        ss          = size(surftmp);
        corners(lL+1,1:cs(1),1:cs(2))   = corntmp;
        edges(lL+1,1:es(1),1:es(2))     = edgetmp;
        surfaces(lL+1,1:ss(1),1:ss(2))  = surftmp;
        
        lNoOfInd(lL+1,:) = [cs(1) es(1) ss(1)];
    end

end