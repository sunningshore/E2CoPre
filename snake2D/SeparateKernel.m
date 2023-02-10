function [K1 KN ERR]=SeparateKernel(H)
% This function SEPARATEKERNEL will separate ( do decomposition ) any 
% 2D, 3D or nD kernel into 1D kernels. Ofcourse only a sub-set of Kernels
% are separable such as a Gaussian Kernel, but it will give least-squares
% sollutions for non-separatable kernels.
% 
% Separating a 3D or 4D image filter in to 1D filters will give an large
% speed-up in image filtering with for instance the function imfilter.
%
% [K1 KN ERR]=SeparateKernel(H);
%   
% inputs,
%   H : The 2D, 3D ..., ND kernel
%   
% outputs,
%   K1 : Cell array with the 1D kernels
%   KN : Approximation of the ND input kernel by the 1D kernels
%   ERR : The sum of absolute difference between approximation and input kernel
%
% 
% We first make some structure which contains information about
% the transformation from kernel to 1D kernel array, number of dimensions
% and other stuff
data=InitializeDataStruct(H);

% Make the matrix of c = M * d; 
M=makeMatrix(data);

% Solve c = M * d with least squares
warning('off','MATLAB:rankDeficientMatrix');
par=exp(M\log(abs(data.H(:))));

% Improve the values by solving the remaining difference
KN = Filter1DtoFilterND(par,data);
par2=exp(M\log(abs(KN(:)./data.H(:))));
par=par./par2;

% Change the sign of a 1D filtering value if it decrease the error
par = FilterCorrSign(par,data);

% Split the solution d in separate 1D kernels
K1 = ValueList2Filter1D(par,data);

% Re-add the removed zero rows/planes to the 1D vectors
K1=re_add_zero_rows(data,K1);

% Calculate the approximation of the ND kernel if using the 1D kernels
KN = Filter1DtoFilterND(par,data,K1);

% Calculate the absolute error
ERR =sum(abs(H(:)-KN(:)));

function par = FilterCorrSign(par,data)
Ert=zeros(1,length(par));
ERR=inf; t=0;
par=sign(rand(size(par))-0.5).*par;
while(t<ERR)
    % Calculate the approximation of the ND kernel if using the 1D kernels
    KN = Filter1DtoFilterND(par,data);
    % Calculate the absolute error
    ERR =sum(abs(data.H(:)-KN(:)));
    % Flip the sign of every 1D filter value, and look if the error
    % improves
    for i=1:length(par)
        par2=par; par2(i)=-par2(i);
        KN = Filter1DtoFilterND(par2,data);
        Ert(i) =sum(abs(data.H(:)-KN(:)));
    end
    % Flip the sign of the 1D filter value with the largest improvement
    [t,j]=min(Ert); if(t<ERR), par(j)=-par(j); end
end

function data=InitializeDataStruct(H)
data.sizeHreal=size(H);
data.nreal=ndims(H);
[H,preserve_zeros]=remove_zero_rows(H);
data.H=H;
data.n=ndims(H);
data.preserve_zeros=preserve_zeros;
data.H(H==0)=eps;
data.sizeH=size(data.H);
data.sep_parb=cumsum([1 data.sizeH(1:data.n-1)]);
data.sep_pare=cumsum(data.sizeH);
data.sep_parl=data.sep_pare-data.sep_parb+1;
data.par=(1:numel(H))+1;

function [H,preserve_zeros]=remove_zero_rows(H)
% Remove whole columns/rows/planes with zeros,
% because we know at forehand that they will give a kernel 1D value of 0
% and will otherwise increase the error in the end result.
preserve_zeros=zeros(numel(H),2); pz=0;
sizeH=size(H);
for i=1:ndims(H)
    H2D=reshape(H,size(H,1),[]);
    check_zero=~any(H2D,2);
    if(any(check_zero))
        zero_rows=find(check_zero);
        for j=1:length(zero_rows)
            pz=pz+1;
            preserve_zeros(pz,:)=[i zero_rows(j)];
            sizeH(1)=sizeH(1)-1;
        end
        H2D(check_zero,:)=[];
        H=reshape(H2D,sizeH);
    end
    H=shiftdim(H,1);
    sizeH=circshift(sizeH,[0 -1]);
    H=reshape(H,sizeH);
end
preserve_zeros=preserve_zeros(1:pz,:);

function K1=re_add_zero_rows(data,K1)
% Re-add the 1D kernel values responding to a whole column/row or plane
% of zeros
for i=1:size(data.preserve_zeros,1)
    di=data.preserve_zeros(i,1);
    pos=data.preserve_zeros(i,2);
    if(di>length(K1)), K1{di}=1; end
    val=K1{di};
    val=val(:);
    val=[val(1:pos-1);0;val(pos:end)];
    dim=ones(1,data.nreal); dim(di)=length(val);
    K1{di}=reshape(val,dim);
end

function M=makeMatrix(data)
 M = zeros(numel(data.H),sum(data.sizeH));
 K1 = (1:numel(data.H))';
 for i=1:data.n;
    p=data.par(data.sep_parb(i):data.sep_pare(i)); p=p(:);
    dim=ones(1,data.n); dim(i)=data.sep_parl(i);
    Ki=reshape(p(:),dim);
    dim=data.sizeH; dim(i)=1;
    K2=repmat(Ki,dim)-1;
    M(sub2ind(size(M),K1(:),K2(:)))=1;
 end
 
function Kt = Filter1DtoFilterND(par,data,K1)
if(nargin==2)
 Kt=ones(data.sizeH);
 for i=1:data.n
     p=par(data.sep_parb(i):data.sep_pare(i)); p=p(:);
     dim=ones(1,data.n); dim(i)=data.sep_parl(i);
     Ki=reshape(p(:),dim);
     dim=data.sizeH; dim(i)=1;
     Kt=Kt.*repmat(Ki,dim);
 end
else
  Kt=ones(data.sizeHreal);
  for i=1:data.n
    dim=data.sizeHreal; dim(i)=1;
    Kt=Kt.*repmat(K1{i},dim);
  end
end

function K = ValueList2Filter1D(par,data)
 K=cell(1,data.n);
 for i=1:data.n
     p=par(data.sep_parb(i):data.sep_pare(i)); p=p(:);
     dim=ones(1,data.n); dim(i)=data.sep_parl(i);
     K{i}=reshape(p(:),dim);
 end