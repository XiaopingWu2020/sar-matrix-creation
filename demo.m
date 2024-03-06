%%% A demo showing how to create voxel and local SAR matrices given the EM simulation of multichannel transmission. 
%%% This workflow is updated to work with EM simulation run by Alireza using CST.
%%%
%%% Created by Xiaoping Wu, 3/5/2024
%%%

clear

% %% data prep
% load('Duke_data_Ch_by_Ch_E_Field_for_Xiaoping.mat')
% load('Sigma_3Components.mat'); load('Rho_v2.mat');
% 
% Efield= complex(zeros([size(Mask) 3 N]));
% iEfield= Efield(:,:,:,:,1);
% iEf= iEfield(:,:,:,1);
% for ich=1:N
%     eval(['iE= E_Ch',num2str(ich),';']);
%     Efield(:,:,:,:,ich)= reshape(iE,[size(Mask) 3]);
% end
% 
% condMap(:,:,:,1)= Sigma1;
% condMap(:,:,:,2)= Sigma2;
% condMap(:,:,:,3)= Sigma3; 
% mdenMap= cat(4,Rho_3D_v2, Rho_3D_v2, Rho_3D_v2); % Rho= sum_i(J_i* E_i)/ pointSAR
% 
% % truncate the data so that it can be shared on github. 
% soi= 81:120;
% Efield= Efield(:,:,soi,:,:);
% save Efield Efield
% condMap= condMap(:,:,soi,:);
% mdenMap= mdenMap(:,:,soi,:);
% Rho= Rho(:,:,soi); % Rho= powerDensity/ pointSAR
% save tissueProperties condMap mdenMap Rho

%
load('Efield.mat')
load('tissueProperties.mat')
res=1e-3*[2 2 2]; % m3
Ng= 0.01; %kg, 10g

%% voxel SAR matrices calc
[Sv_vox,~,mask]= construct_Sv_voxsar(Efield,condMap,mdenMap, logical(Mask), false);
voxSARmatrix.Smatrix= Sv_vox;
voxSARmatrix.mask= mask;
save('voxSARmatrix', 'voxSARmatrix','-v7.3');

%% local SAR matrices calc
sTh={'THRESHOLD=1'}; % no edge
[Sv_ave,mask]= construct_Sv_avesar(voxSARmatrix.Smatrix,Rho,res,Ng,sTh);
NgSARmatrix.Smatrix= Sv_ave;
NgSARmatrix.mask= mask;
save('NgSARmatrix', 'NgSARmatrix','-v7.3');
