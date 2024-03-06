function [Sv_ave, mymask, Sv_global] = construct_Sv_avesar (Sv_vox, mydenm, ...
    res, Ng, threshold)
% CONSTRUCT_SV_AVESAR Construct Sv for local Ng SAR calculations.
%
% Usage: [Sv_ave, mymask, Sv_global] = construct_Sv_avesar (Sv_vox, mydenm, res, Ng)
%
% Returns
% -------
%
% Sv_ave: nchs x nchs x nValidVoxels array of local sar matrices
%
% mymask: a mask defining nValidVoxels.
% Sv_global: nchs x nchs global sar matrix.
%
% Expects
% -------
%
% Sv_vox: nchs x nchs x nValidVoxels array, or (nchsxnchs) x nValidVoxels matrix
% for voxel sar Sv.
%
% mydenm: masked 3D density map.
%
% res: a 1x3 vector specifying resolution in m. defaults to 1e-3*[5 5 5]. if
% anisotropic, interp for isotropic will be performed before local sar is
% calculated.
%
% Ng: target mass in kg. defaults to 0.01 kg (10g)
%
% threshold: specifies the threshold used to define edge tissue. defaults to
% {'THRESHOLD=1'}, meaning no edge. a reasonal threshold for edge definition 
% would be {'THRESHOLD=0.12'}.
%
%
% See also: clean_Qmat construct_Sv_voxsar find_VOPs construct_S_glosar calcNgAveSARAdv3
%
%
% Copyright (C) 2012 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu>
% Created: Fri Jun 15 16:05:44 2012
%

if nargin<3|| isempty(res)
    res= 1e-3*[5 5 5];end

if nargin<4|| isempty(Ng)
    Ng= 0.01;end

if nargin<5
    threshold={'THRESHOLD=1'};
end

%
siz= size(Sv_vox);
Sv_VOX = reshape(Sv_vox,[],siz(end));
mask=~~mydenm;

res_min=min(res);
if any(res-res_min) % anisotropic, thus interp
    disp('-> Converting anisotropic voxels to isotropic via interp...')
    elemsiz= res_min./res;
    mydenm= affine(mydenm,eye(4),elemsiz);
    
    parfor ind=1:size(Sv_VOX,1)
        tmpsarv = complex(zeros(size(mask)));
        tmpsarv(mask) = Sv_VOX(ind,:);
        tmpsarv1= affine(tmpsarv,eye(4),elemsiz);
        Sv_VOX1(ind,:)= tmpsarv1(~~mydenm);
    end
    
    Sv_VOX= Sv_VOX1;
    mask=~~mydenm;
end

vol= res_min.^3;

disp('-> Local Ng SAR Sv construction started...')
% run once mainly to get mymask
ind=1;
tmpsarv = complex(zeros(size(mask)));
tmpsarv(mask) = Sv_VOX(ind,:);

tmpsara = calcNgAveSARAdv3(tmpsarv,mydenm,vol,Ng,threshold);
mymask= ~~abs(tmpsara);
Sv_AVE(ind,:) = tmpsara(mymask);
Sv_GLB(ind)= sum(tmpsarv(mask).*mydenm(mask))./sum(mydenm(mask));
%
parfor ind= 2: size(Sv_VOX,1) % loop over elements of Sv, each treated as if it
    % were a voxel sar map.
    tmpsarv = complex(zeros(size(mask)));
    tmpsarv(mask) = Sv_VOX(ind,:);
    % threshold=1, no edge
    tmpsara = calcNgAveSARAdv3(tmpsarv,mydenm,vol,Ng,threshold);
    Sv_AVE(ind,:) = tmpsara(mymask);
    Sv_GLB(ind)= sum(tmpsarv(mask).*mydenm(mask))./sum(mydenm(mask));
end

if any(res-res_min) % if interp, inv interp
    disp('-> Converting isotropic voxels back to anisotropic via interp...')
    elemsiz= res./res_min;
    % run once to get mymask;
    ind=1;
    tmpsara=complex(zeros(size(mymask)));
    tmpsara(mymask) = Sv_AVE(ind,:);
    tmpsara1= affine(tmpsara,eye(4),elemsiz);
    mymask1= ~~abs(tmpsara1);
    Sv_AVE1(ind,:)= tmpsara1(mymask1);
    
    parfor ind=2:size(Sv_AVE,1)
        tmpsara=complex(zeros(size(mymask)));
        tmpsara(mymask) = Sv_AVE(ind,:);
        tmpsara1= affine(tmpsara,eye(4),elemsiz);
        Sv_AVE1(ind,:)= tmpsara1(mymask1);
    end
    
    Sv_AVE= Sv_AVE1;
    mymask= mymask1;
end

nchs= siz(1);
Sv_ave= reshape(Sv_AVE,nchs,nchs,[]);
Sv_global= reshape(Sv_GLB,nchs,nchs);


disp('-> Local Ng SAR Sv construction DONE...')

