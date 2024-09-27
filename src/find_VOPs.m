function [AA, myeps, SS, u, clusterLabel] = find_VOPs (Sv, ubp, bIsatest)

% FINDVOPS find virtual observation points given the characteristic matrices of
% EM model. Algorithm is implemented following G. Eichfelder and M. Gebhardt,
% MRM 2011.
%
% Usage: [AA, myeps, SS, u, clusterLabel] = find_VOPs (Sv, ubp, bIsatest)
%
% Returns
% -------
% AA: VOPs for individual clusters created.
% myeps: overshoots of individual clusters.
% SS: a cell array containing individual clusters.
% u: upper bound for overestimation.
% clusterLabel: labels indicating to which cluster each matrix belongs.
%
%
% Expects
% -------
% Sv: characteristic matrices of size nchs x nchs x nValidVoxels.
% 
% ubp: upper bound percentage, which will be used to determine the
% overestimation of SAR. defaults to 0.05, meaning overestimation will be 5% of
% max eigenvalue of all Sv matrices of the model.
%
% bIsatest: when its true, only info about upper bound will be displayed. defaults to false,
% meaning vops will be created. 
%
% See also: construct_Sv_avesar construct_Sv_voxsar
%
%
% Copyright (C) 2012 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Tue Jun 19 10:53:57 2012
%
if nargin<2
  ubp= 0.05;
end
if nargin<3
    bIsatest= false;
end

siz= size(Sv);
npts = siz(end);

%% STEP0: 
% set u based on largest eigenvalue of the model
eigMAX= zeros(1,npts);
parfor ind= 1:npts,
  d = real(eig(Sv(:,:,ind))); % remove numeric error by taking real part
  eigMAX(ind)= max(d);
end
eigmax = max(eigMAX);
u= ubp* eigmax;

if bIsatest
    fprintf('-> Upper bound is %f for %2.1f%% of worst case max SAR.\n',...
        u,100*ubp);
    return
end

AA= complex(zeros(siz(1),siz(2),1000));
myeps= zeros(1,1000);
SS= cell(1,1000);

k_clusterCounter=0;% set k
BB= Sv;

clusterLabel= zeros(1,npts,'uint16');
ix0= 1:npts;

while (~isempty(BB))
    
    k_clusterCounter= k_clusterCounter+ 1;
    
    %% STEP1:
    % choose B_kstar= argmax(||B||), B in OMEGA.
    npts= size(BB,3);
    mynorm= zeros(1,npts);
    parfor ind= 1:npts,
        mynorm(ind)= norm(BB(:,:,ind),2);
    end
    [~,idxnormmax]= max(mynorm(:));
    B_kstar= BB(:,:,idxnormmax);
    
    %% STEP2: clustering
    % sort all B wrt eigval_min(B_kstar - B) in decreasing order
    dmin = zeros(1,npts);
    parfor ind=1:npts,
        d = real(eig(B_kstar- BB(:,:,ind)));
        dmin(ind) = min(d(:));
    end
    SS_k= BB(:,:,dmin>=0);

    clusterLabel(ix0(dmin>=0))= k_clusterCounter;
    
    [dminsorted,ix] = sort(-dmin);
    dminsorted= -dminsorted;
    dminsorted1= dminsorted(dminsorted<0);
    ix1= ix(dminsorted<0);

    %ix0= ix0(ix);
    
    %% STEP3: solve for Z
    Z_k= 0;
    eps_k= 0;
    
    if ~isempty(ix1)
        for idx=1: (length(ix1))
            len= idx;
            
            Q= B_kstar+ Z_k- BB(:,:,ix1(len));
            % spectral decomposition of Q, i.e., Q = Q+ - Q-.
            [V,E]= eig(Q);
            E= diag(real(diag(E)));
            diagE= diag(E);
            diagE(diagE<0)= 0;
            E_plus= diag(diagE);
            E_minus= E_plus- E;
            % Q_plus= V* E_plus* V';
            Q_minus= V* E_minus* V';
            
            if norm(Z_k+Q_minus,2)> u
                len= len- 1;
                break;
            end
            
            eps_k= -dminsorted1(len);
            Z_k= Z_k+ Q_minus;
            
        end
    else
        len=0;
    end

    %% STEP4:
    % set cluster
    if len> 0
        SS_k= cat(3,SS_k,BB(:,:,ix1(1:len)));
        clusterLabel(ix0(ix1(1:len)))= k_clusterCounter;
        ix1= ix1(len+1:end);
    end
    
    % save VOP, overshoot and cluster
    AA(:,:,k_clusterCounter)= B_kstar+ Z_k;
    myeps(k_clusterCounter) = eps_k;
    SS{k_clusterCounter}= SS_k;
    clear SS_k;
    
    % update set
    BB= BB(:,:,ix1);    
    ix0= ix0(ix1);
    
    fprintf(1,'-> VOP created for cluster %d...\n',k_clusterCounter);
    
end

AA= AA(:,:,1:k_clusterCounter);
myeps= myeps(1:k_clusterCounter);
SS= SS(1:k_clusterCounter);

nvops=size(AA,3);
parfor idx=1:nvops
    iAA= AA(:,:,idx);
    AA(:,:,idx)= iAA- diag(diag(iAA))+ diag(real(diag(iAA)));
end

disp('-> VOP construction DONE!')
