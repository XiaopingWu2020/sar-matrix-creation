function [traj, trajval] = buildSARAveTrajChessb(lKernelsize)
%*************************************
% 
% creates a Trajectory, along which the SAR is averaged.
%
% INPUT: lKernelsize (e.g. 20) gives the size diameter of the averaging 3
% dimensional cube
%           lMode:  chessboard (fixed in this function!)
%                  
%
% Output: is a Nx4 vector. 1.-3. element gives karthesian koordinates of
% the repective trajectory element and 4 the distance to the center
%
%*****************************************************


lKern = zeros(lKernelsize, lKernelsize, lKernelsize);
lCenter =   [floor(lKernelsize/2)+mod(lKernelsize,2)...
            floor(lKernelsize/2)+mod(lKernelsize,2)...
            floor(lKernelsize/2)+mod(lKernelsize,2)];
lKern(lCenter(1),lCenter(2),lCenter(3)) = 1;

lMode = 2;

switch lMode
    case 1 
        sMode = 'euclidean';
    case 2
        sMode = 'chessboard';
    case 3 
        sMode = 'cityblock';
    case 4
        sMode = quasi-euclidean'
    otherwise
        sMode = ''
        error('ERROR: wrong mode');
end

%create distance Map with mode specified as defined by argument lMode of the
%function
dDistMap = bwdist(lKern,sMode);

%Nx1 Array with sorted distance:
dSortDist = unique(sort(dDistMap(:)));

traj    = [];
trajval = [];

for lL = 1:length(dSortDist)
    atmpInd = find(dDistMap == dSortDist(lL));
    
    [i1,i2,i3] = ind2sub(size(lKern), atmpInd);
    ival = repmat(dSortDist(lL),[size(i1),1]);
    
    %koordinates of all the trajectory points sorted by distance from
    %center
    %trajtmp = [i1-floor(lKernelsize/2) i2-floor(lKernelsize/2) i3-floor(lKernelsize/2)];
    trajtmp = [i1-lCenter(1) i2-lCenter(2) i3-lCenter(3)];
    
    
    %sort traj to use radial symmetry
    trajlen = length(trajtmp);
    for lM = 2:2:trajlen
        trajtmp = [trajtmp(1:(lM-1),:); circshift(trajtmp(lM:end,:),[1,0])];
    end
    
    %use euclidean distance to sort within one shell (does not make sense
    %if mode == 1 ...)
    dEucDist = trajtmp(:,1).^2+trajtmp(:,2).^2+trajtmp(:,3).^2;
    [dEucDistSort, dEucDistInd] = sort(dEucDist);
    trajtmp = trajtmp(dEucDistInd,:);
    
    %just take one surface of the cube:
    lRadius = max(trajtmp(:,1));
    trajtmpSurf1 = find(trajtmp(:,1) == lRadius);
    
    trajXYZ = trajtmp(trajtmpSurf1,:);
    trajX = circshift(trajXYZ,[0 0]);
    trajY = circshift(trajXYZ,[0 1]);
    trajZ = circshift(trajXYZ,[0 2]);
    
    trajAllNoCut = [trajX, -trajX, trajY, -trajY, trajZ, -trajZ];
    tsize = size(trajAllNoCut);
    tlength = tsize(1)*tsize(2);
    trajAllNoCutResh = reshape(trajAllNoCut',3,tlength/3)';
    [b, m, n] = unique(trajAllNoCutResh,'rows','first');
    trajtmp = trajAllNoCutResh(sort(m),:);
    

    traj = [traj; trajtmp];

    %values of each traj.-point defined by traj
    trajval = [trajval;ival]; 
end

end %function