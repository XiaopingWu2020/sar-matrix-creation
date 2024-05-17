function [Sv_vox, Sv_VOX, sarmask] = construct_Sv_voxsar (myE, mycond, myden, mask, isFDTD)

% CONSTRUCT_VOX_SV Construct characteristic space matrix Sv for voxel SAR
% calculations given precalculated E field and electric properties of a tissue
% model.
%
% Usage: [Sv_vox, Sv_VOX, sarmask] = construct_Sv_voxsar (myE, mycond, myden, mask, isFDTD)
%
% Returns
% -------
% Sv_vox: nchs x nchs x nValidVoxels array, nValidVoxels is determined by mask.
%
% Sv_VOX: (nchs x nchs) x nValidVoxels matrix.
%
% sarmask: a mask for SAR. 
%
% Expects
% -------
% myE: [x y z 3 nchs] array for E field
% mycond: [x y z 3] array for conductivity
% myden: [x y z 3] array for mass density
%
% mask: a mask determined by remcom vox sar calc, representing valid voxels.
% When isFDTD is false, this mask is not used.
%
% isFDTD: a flag indicating if the EM simulation has been done using
% FDTD. If so, given the fact that E fields and tissue properties are defined on edges of Yee
% cells, the voxel SAR calculation takes two steps. The first is to
% calculate edge SAR and the second step to average/filter edge SAR values to
% calculate voxel SAR values for individual Yee cells. If not, only the
% first step is done assuming E fields and tissue properties are defined
% for each voxel. Defaults to true meaning EM simulation is done using
% FDTD.
%
% See also: construct_Sv_avesar construct_S_glosar find_VOPs compute_voxSAR
%
%
% Copyright (C) 2012 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu>
% Created: Thu Jun 14 16:00:18 2012
%

if nargin< 5
    isFDTD= true;
end

cod= removeNaNs(0.5*mycond./myden);
codX= cod(:,:,:,1);
codY= cod(:,:,:,2);
codZ= cod(:,:,:,3);

if isFDTD
    mymask= codX>0| codY>0| codZ>0;
else
    mymask= mask;
end

codx= codX(mymask);
cody= codY(mymask);
codz= codZ(mymask);

siz= size(myE);
nchs= siz(end);
EX= squeeze(myE(:,:,:,1,:));
EY= squeeze(myE(:,:,:,2,:));
EZ= squeeze(myE(:,:,:,3,:));
myMask=repmat(mymask,[1 1 1 nchs]);
Ex = reshape(EX(myMask),[],nchs);
Ey = reshape(EY(myMask),[],nchs);
Ez = reshape(EZ(myMask),[],nchs);

% construct Sv for edge sar calculations
myidx= find(mymask);
disp('-> Edge SAR Sv under construction ...')
Sv_edgeX= complex(zeros([nchs nchs length(myidx)]));
Sv_edgeY= Sv_edgeX;
Sv_edgeZ= Sv_edgeX;
parfor idx=1:length(myidx) % loop over space
    Sv_edgeX(:,:,idx)= codx(idx).* Ex(idx,:)' * Ex(idx,:);
    Sv_edgeY(:,:,idx)= cody(idx).* Ey(idx,:)' * Ey(idx,:);
    Sv_edgeZ(:,:,idx)= codz(idx).* Ez(idx,:)' * Ez(idx,:);
end
clear E* myE
disp('-> Edge SAR Sv construction DONE ...')

% construct Sv for voxel sar calculations
if isFDTD
    Sv_edgeXX= reshape(Sv_edgeX,[],length(myidx));
    Sv_edgeYY= reshape(Sv_edgeY,[],length(myidx));
    Sv_edgeZZ= reshape(Sv_edgeZ,[],length(myidx));
    
    Sv_edgeXXf= complex(zeros(size(Sv_edgeXX,1),length(find(mask))));
    Sv_edgeYYf= Sv_edgeXXf;
    Sv_edgeZZf= Sv_edgeYYf;
    
    % filter edge sar along each of 3 dimensions.
    disp('-> FDTD: Filtering edge SAR Sv started ...')
    
    kernel= designKernel();
    % x
    parfor ind= 1:size(Sv_edgeXXf,1) % loop over elements of Sv
        mymesh= complex(zeros(size(mymask)));
        mymesh(mymask)= Sv_edgeXX(ind,:);
        mytmp= imfilter(mymesh,kernel(:,:,:,1));
        mytmp= circshift(mytmp,[-1 -1]);
        Sv_edgeXXf(ind,:) = mytmp(mask);
    end
    disp('-> Filtering along x done...')
    % y
    parfor ind= 1:size(Sv_edgeYYf,1)
        mymesh= complex(zeros(size(mymask)));
        mymesh(mymask)= Sv_edgeYY(ind,:);
        mytmp= imfilter(mymesh,kernel(:,:,:,2));
        mytmp= circshift(mytmp,[-1 -1]);
        Sv_edgeYYf(ind,:) = mytmp(mask);
    end
    disp('-> Filtering along y done...')
    % z
    parfor ind= 1:size(Sv_edgeZZf,1)
        mymesh= complex(zeros(size(mymask)));
        mymesh(mymask)= Sv_edgeZZ(ind,:);
        mytmp= imfilter(mymesh,kernel(:,:,:,3));
        mytmp= circshift(mytmp,[-1 -1]);
        Sv_edgeZZf(ind,:) = mytmp(mask);
    end
    disp('->Filtering along z done...')
    
    %
    Sv_VOX= Sv_edgeXXf+ Sv_edgeYYf+ Sv_edgeZZf;
    Sv_vox= reshape(Sv_VOX,nchs,nchs,[]);
    sarmask= mask;
else
    disp('-> NOT FDTD: skipping the second step...')
    Sv_vox= Sv_edgeX+ Sv_edgeY+ Sv_edgeZ;
    Sv_VOX= reshape(Sv_vox, nchs*nchs,[]);
    sarmask= mymask;
end

disp('-> Voxel SAR Sv construction DONE...')


%% ====================
%%
%%  Local functions
%%
%% ====================

%
function kernel = designKernel()
% construct the kernel for filtering the 3d sar distribution.

% this first kernel is designed based on the mesh structure of remcom tissue model.

kernel = zeros(2,2,2,3);                             % 3 for x, y, z

% x
kernel(1,1,1,1)=1;
kernel(1,2,1,1)=1;
kernel(1,1,2,1)=1;
kernel(1,2,2,1)=1;

% y
kernel(1,1,1,2)=1;
kernel(2,1,1,2)=1;
kernel(1,1,2,2)=1;
kernel(2,1,2,2)=1;

%z
kernel(1,1,1,3)=1;
kernel(2,1,1,3)=1;
kernel(1,2,1,3)=1;
kernel(2,2,1,3)=1;

kernel = 0.25* kernel;
