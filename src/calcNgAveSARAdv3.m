function [sarFIN, massAct,outprop] = calcNgAveSARAdv3(sarmaporig,densmaporig,voxvol,dNGrams, varargin)
%% ************************************************************************
%
% Calculates N grams average SAR of a given 3D local SARmap. Uses the
% advanced calculation algorithm described in Caputa et al. IEEE Antennas and Propagation
% Magazine Vol 41 No 4, 1999
% Additionnaly the IEEE recommendation for edge handling is used!
%
%   modified by xwu, Mar 2013:
%   skip the skin treatment if there is no skin specified (ie a large
%   threshold is used). in addition, calculation is made compatible with
%   VOP generation.
%
%   modified March 2012:
%   add optional parameters: theshold for boundary voxels can be set now
%
%   INPUT:                                                  [unit]
%   ----------------------------------------------------------------
%   sarmap:     3D matrix with voxel SAR data               W/kg
%   densmap:    density map, corresponding to sarmap        kg/m^3
%   voxvol:     voxelvolume                                 m^3
%   dNGrams:    Grams to average                            kg
%
%
%   OUTPUT:
%   ----------------------------------------------------------------
%   sarFIN      averaged Ng SAR                             W/kg
%   masssum   	exact averaged mass of each voxel           kg
%
%% ************************************************************************

%ssm 220312
propdefault = struct();
propdefault.THRESHOLD = 0.12;

prop  = catstruct(propdefault,parseVariableInputs(varargin));

%add padding to avoid problems with circshift-function!
nz=10;
sarmap  = addPadding(sarmaporig,nz);
densmap = addPadding(densmaporig,nz);

sarsize = size(sarmap);
denssize = size(densmap);
if (sarsize ~= denssize)
    error('ERROR: size of SAR map and Density must be equal');
end
if (size(sarsize) ~= 3)
    error('ERROR: dimension of SAR map / Density must be 3');
end

%replace NaNvalues in sarmap by Zeros
bNaNVal = isnan(sarmap);
sarmap(bNaNVal) = 0;

densmapmask = ~~(densmap);

%convert some units: calculate everything in SI units:
dNGramsSI = dNGrams;                %in kg
%(%changed in adv. Version: *1E-3;); input is already kg!
voxvolSI    = voxvol;               %in m3

%calculate the growing trajectory for averaging the SAR
[traj, trajval] = buildSARAveTrajChessb(20);
trajvaluniq = unique(trajval);

lCPt        = 1;
aveSAR      = zeros(sarsize);
avemask     = densmapmask; % ones(sarsize);
denssum     = zeros(sarsize);
lNoOfElem   = length(trajval);
lCountMap   = zeros(sarsize);
edgemask    = false(sarsize);
skinmask    = false(sarsize);
lAirVox     = uint16(zeros(sarsize));

%Main loop. It loops over the distance from Isocenter
fprintf('Calculate kernel. Radius: ');
for lL = 1:size(trajvaluniq)
    
    %Inner loop which loops over the individual voxel of a constant
    %distance from Isoventer. By doing this, the average SAR grows in
    %shells
    while(trajvaluniq(lL) == trajval(lCPt))
        %helper matrix, saves time:
        densmaptmp  = circshift(densmap,traj(lCPt,:)).*avemask;
        %this is the average SAR, but in W/m^3
        aveSAR      = aveSAR+circshift(sarmap,traj(lCPt,:)).*densmaptmp;
        %sum of the density in kg/m^3
        denssum     = denssum + densmaptmp;
        
        %mask for the edges; those voxels are treated seperately later.
        edgemask = edgemask | (logical(avemask) & ~logical(densmaptmp));
        %~logical(densmaptmp) is the shifted density map multiplied by
        %avemask and then negated. this is equeivalent to
        % edgemask |  (logical(avemask) & (~densmap | ~avemask) )
        % => edgemask |  (logical(avemask) & (~densmap )) because
        % a&(b|c) = a&b | a&c
        
        %counter for the number of AirVoxels for the average:
        lAirVox     = lAirVox+ uint16((logical(avemask) & ~logical(densmaptmp)));
        
        if(lCPt == lNoOfElem)
            break;
        end
        
        %this has to be commented out. Thus only full cubes are
        %calculated
        %avemask = ((denssum*voxvolSI <  dNGramsSI) & densmapmask);
        
        %lCountMap(avemask) = lCPt;
        
        %running variable:
        lCPt = lCPt+1;
    end
    
    %gives the 'radius' of the cube for averaging:
    lCountMap(avemask) = trajval(lCPt-1);
    
    lThrld   = ceil((trajvaluniq(lL)*2+1)^3*prop.THRESHOLD);
    %determine mask for next round, voxels are not considered anymore, when averaged
    %mass exceeds the Ng mass threshold defined as input parameter, or
    %when voxel is marked as edge:
    skinmask = skinmask | (lAirVox>lThrld);
    avemask = ((denssum*voxvolSI < dNGramsSI) & densmapmask) & ~skinmask;% ~edgemask;
    
    %loops breakes, if all values average SAR values are averaged over
    %mor than N grams:
    if(~avemask)
        break;
    end
    
    fprintf('%d...',lL);
    
end

fprintf('Done\n');
fprintf('Calculate volume fraction of outer shell... ');
%avemask = ((denssum*voxvolSI < dNGramsSI) & densmapmask) & ~(lAirVox>5);%
%mod hier:
edgemask = skinmask;% (lAirVox>lThrld);

%calcuclate final mass and average SAR in correct units
masssum = denssum*voxvolSI;         %mass in kg
%aveSAR  = aveSAR./denssum;          %SAR in W/kg

massCorn    = zeros(size(aveSAR));
massEdge    = zeros(size(aveSAR));
massSurf    = zeros(size(aveSAR));

sarCorn     = zeros(size(aveSAR));
sarEdge     = zeros(size(aveSAR));
sarSurf     = zeros(size(aveSAR));

lMaxRad     = max(lCountMap(:));
[corners, edges, surfaces, lNoOfInd] = calcCornersEdgesSurf(lMaxRad);
%calculate for each voxel the contribution of the corners, edges and
%surfaces of the outermost shell:
for lR = 1:lMaxRad+1
    masktmp = (lCountMap == lR-1) &~edgemask;
    
    for lC=1:lNoOfInd(lR,1)
        %calculate mass of corners
        masstmp     = circshift(densmap, reshape(corners(lR,lC,:),[1 3])).*masktmp;
        massCorn    = massCorn + masstmp;
        
        sarCorn     = sarCorn + ...
            circshift(sarmap, reshape(corners(lR,lC,:),[1 3])).*masstmp;
    end
    %edges:
    for lE=1:lNoOfInd(lR,2)
        %calculate mass of edges
        masstmp     = circshift(densmap, reshape(edges(lR,lE,:),[1 3])).*masktmp;
        massEdge    = massEdge + masstmp;
        
        sarEdge     = sarEdge + ...
            circshift(sarmap, reshape(edges(lR,lE,:),[1 3])).*masstmp;
    end
    %surfaces:
    for lS=1:lNoOfInd(lR,3)
        %calculate mass of surfaces
        masstmp     = circshift(densmap, reshape(surfaces(lR,lS,:),[1 3])).*masktmp;
        massSurf = massSurf + masstmp;
        
        sarSurf     = sarSurf + ...
            circshift(sarmap, reshape(surfaces(lR,lS,:),[1 3])).*masstmp;
    end
end
massCorn    = massCorn*voxvol;
massEdge    = massEdge*voxvol;
massSurf    = massSurf*voxvol;

massShell   = massCorn+massEdge+massSurf;
sarShell    = sarCorn+sarEdge+sarSurf;

%calculate fraction of mass needed from outer shell:
massFrac = (dNGrams*ones(size(masssum))-(masssum-massShell))./ ...
    (massShell);

msize = size(massCorn);
c = massCorn(:);
e = massEdge(:);
s = massSurf(:);
k = massFrac(:).*massShell(:);

fprintf('Done\n');
fprintf('Solve cubic equation... ');

xx = (((k./(2*c) - e.^3./(27*c.^3) + (e.*s)./(6*c.^2)).^2 + (s./(3*c) - ...
    e.^2./(9*c.^2)).^3).^(1/2) + k./(2*c) - e.^3./(27*c.^3) + ...
    (e.*s)./(6*c.^2)).^(1/3) - e./(3*c) - (s./(3*c) - e.^2./(9*c.^2))./(((k./(2*c) ...
    - e.^3./(27*c.^3) + (e.*s)./(6*c.^2)).^2 + (s./(3*c) - e.^2./(9*c.^2)).^3).^(1/2) ...
    + k./(2*c) - e.^3./(27*c.^3) + (e.*s)./(6*c.^2)).^(1/3);

f = reshape(real(xx),msize);
massAct = masssum - massShell + f.^3.*massCorn +f.^2.*massEdge + f.*massSurf;

sarTMP = (aveSAR - sarShell + f.^3.*sarCorn +f.^2.*sarEdge + f.*sarSurf)./massAct*voxvol;

sarAll= sarTMP;
sarNonSkin = removeNaNs(sarTMP).*double(~edgemask);
dCurSarSk   = zeros(size(sarNonSkin));

if any(edgemask(:))
    %handle the edges. calculate SAR for those voxels that fall out
    %of the 'normal' calculation process (above) because they are to close
    %to air. However the local SAR of those voxels are used for the average
    %SAR cacluation of other (neighboring) voxels.
    
    edgeSAR = zeros(size(sarNonSkin));
    %tmp = ((double(~~lCountMap)+1)*2+1).^3;
    %new:
    tmp = (double(lCountMap)*2+1).^3;
    
    %**********************************************************************
    % NEW IMPLEMENTATION:
    % Handle the skin voxels here.
    %**********************************************************************
    
    [trajsk, trajval] = buildSARAveTrajChessb(20);
    %add radius
    trajsk(:,4) = sqrt(trajsk(:,1).^2+trajsk(:,2).^2+trajsk(:,3).^2);
    %sort by radius
    [trajsks,tsind] = sort(trajsk,4);
    %last column will be just the index (is probably not needed..)
    trajsks(:,5) = linspace(1,size(trajsks,1),size(trajsks,1));
    trajsk(:,5) = linspace(1,size(trajsks,1),size(trajsks,1));
    
    [dRUniq,dRUniqInd] = unique(trajsk(:,4));
    
    bIsDone = false;
    lRIndex = 1;
    dCurMass    = zeros(size(sarNonSkin));
    
    lNPxLeft    = sum(edgemask);
    skinmask    = (~~abs(sarmap)-~~abs(sarNonSkin)).*~~abs(sarmap);
    dXGrammMask = skinmask;
    densmapNonSkin = densmap.*~~abs(sarNonSkin);
    
    skinmaskInd   = find(skinmask > 0);
    dCurMassR   = zeros(size(skinmaskInd));
    dXGrammMaskR = zeros(size(skinmaskInd));
    dCurSarSkR   = zeros(size(skinmaskInd));
    lRunMap     = zeros(size(skinmaskInd));
    
    fprintf('\n\n\n\n\n\\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');
    fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');
    
    %loop over the neiboring voxels:
    while ~bIsDone
        
        %shifts the full density map (non reduced)
        
        tmpmass     = circshift(densmapNonSkin,-trajsk(lRIndex,1:3)).*voxvol;
        
        tmpsar      = circshift(sarNonSkin,-trajsk(lRIndex,1:3));
        dCurMassR   = dCurMassR+dXGrammMaskR.*tmpmass(skinmaskInd);
        %sar of reduced matrix:
        dCurSarSkR  = max(dCurSarSkR,dXGrammMaskR.*tmpsar(skinmaskInd));
        
        if(sum(lRIndex == dRUniqInd))
            %dXGrammMask = skinmask.*(dCurMass < dNGramsSI);
            dXGrammMaskR = (dCurMassR < dNGramsSI) ;
            
            lRunMap = max(lRunMap,dXGrammMaskR.*lRIndex);
            
            %hmm
            lNPxLeft    = sum(dXGrammMaskR(:));
            
            if(lNPxLeft == 0)
                bIsDone = true;
                break;
            end
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
            fprintf('SkinPixels left: %9.0f; Run = %5.0f', lNPxLeft,lRIndex);
        end
        
        
        lRIndex = lRIndex + 1;
        %for safety....
        if(lRIndex == size(trajsks,1))
            bIsDone = true;
            break;
        end
        
    end
    %**********************************************************************
    %**********************************************************************
    
    dCurSarSkRR = skinmask;
    lRunMapRR   = skinmask;
    dCurSarSkRR(skinmaskInd) = dCurSarSkR;
    lRunMapRR(skinmaskInd) = lRunMap;
    sarAll = sarNonSkin+dCurSarSkRR.*double(skinmask);
end
%replace NaNvalues in sarmap by Zeros
bNaNVal = isnan(sarAll);
sarAll(bNaNVal) = 0;
%remove the padding which was added at the very beginning:
sarFIN          = removePadding(sarAll,nz);
massAct         = removePadding(massAct,nz);
fprintf('\n\n');
outprop = prop;
outprop.sarNonSkin = sarNonSkin;
outprop.dCurSarSk = dCurSarSk;
end