function OutNdArray = addPadding(inNdArray, lNoOfPaddingVox)
%%**********************************************************************
%
%   adds zero Voxels at the and of a N-D Matrix in each dimension.
%   e.g. a 50x30x20 matrix with lNoOfPaddingVox == 5 becomes a 55x35x25
%   matrix where the last 5 elements in each dimension are 0
%
%
%   INPUT:                                                  [unit]
%   ----------------------------------------------------------------
%   corners:            inNdArray can be an arbitrary N-dim Array
%   lNoOfPaddingVox:    number of zero voxel added in each dim
%  
%
%   OUTPUT:
%   ----------------------------------------------------------------
%   OutNdArray;         OutNdArray size of output matrix is 
%                       [Nx+lNoOfPaddingVox,Ny+lNoOfPaddingVox,Nz+lNoOfPadd
%                       ingVox]
%
%%**********************************************************************

lP          = lNoOfPaddingVox;   %just to make it shorter..
arrsize     = size(inNdArray);

ldim        = max(size(arrsize));

s           ='OutNdArray(';
for lD=1:ldim
   st   = sprintf('1:%d',arrsize(lD));
   if lD==1
        s    = [s,st];
   else
        s    = [s,',',st];
   end
end

s           = [s,') = inNdArray;'];

OutNdArray  = zeros(arrsize+lP);
eval(s);

end 
