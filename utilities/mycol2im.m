function A = mycol2im(X, varargin)
% mycol2im        Rearrange matrix columns into blocks.
%                 This is more flexible variant of Matlab col2im
% The 'inverse' is myim2col.m
%
% A = mycol2im(X, argName, argVal, ...);
%----------------------------------------------------------------------
%  X     a matrix
%  There may be an additional number of input arguments, a struct, a cell 
%        or as pairs: argName, argVal, ...
%        These are mostly as in myim2col. Note imsize
%  'imsize' to give size of buildt image, [M,N]. Note this is not
%        necessarily equal to original size of image. With L = size(X,2) we
%        should have: L = ((M-Mn+Mi)/Mi)*((N-Nn+Ni)/Ni)
%  'o' or 'offset' is ignored (set to be [0,0])
%  't' or 'transform'  as in myim2col
%  's' or 'size' as in myim2col
%  'a' or 'adjust' is ignored here
%  'n' or 'neighborhood' as in myim2col, [Mn,Nn] = size(nei);
%  'i' or 'increment' as in myim2col, [Mi, Ni]
%  'v' or 'verbose' to indicate verboseness
%----------------------------------------------------------------------
% example:
%  A = double(imread('lena.bmp')) - 128;
%  p = struct('t','m79', 's',8, 'i','distinct','v',1);
%  X = myim2col(A, p);
%  Ar = mycol2im(X, p, 'imsize', size(A));
%  disp(num2str(norm(A-Ar)));

%----------------------------------------------------------------------
% Copyright (c) 2009.  Karl Skretting.  All rights reserved.
% University of Stavanger (Stavanger University), Signal Processing Group
% Mail:  karl.skretting@uis.no   Homepage:  http://www.ux.his.no/~karlsk/
% 
% HISTORY:  dd.mm.yyyy
% Ver. 1.0  26.11.2009  Made function
% Ver. 1.1  08.03.2010  KS: made the function 'cleaner'
%----------------------------------------------------------------------

%% check if imwrite should be called 
if (nargin < 1)
    error('myim2col: too few input arguments.');
end

%% default options
M=0; N=0;     % size of reconstructed image
tr = 'none';  % transform
Ms=0; Ns=0;   % size (of transform)
nei = [];     % neighborhood, input or given by tr
Mn=0; Nn=0;   % size of nei (neighborhood)
Mi=0; Ni=0;   % increment
verbose = 0;
lotrho = 0.95;  % these three must be the same in myim2col
eltrho = 0.95;
eltarg2 = 0.7;

%% get the options
nofOptions = nargin-1;
optionNumber = 1;
fieldNumber = 1;
while (optionNumber <= nofOptions)
    if isstruct(varargin{optionNumber})
        sOptions = varargin{optionNumber}; 
        sNames = fieldnames(sOptions);
        opName = sNames{fieldNumber};
        opVal = sOptions.(opName);
        % next option is next field or next (pair of) arguments
        fieldNumber = fieldNumber + 1;  % next field
        if (fieldNumber > numel(sNames)) 
            fieldNumber = 1;
            optionNumber = optionNumber + 1;  % next pair of options
        end
    elseif iscell(varargin{optionNumber})
        sOptions = varargin{optionNumber}; 
        opName = sOptions{fieldNumber};
        opVal = sOptions{fieldNumber+1};
        % next option is next pair in cell or next (pair of) arguments
        fieldNumber = fieldNumber + 2;  % next pair in cell
        if (fieldNumber > numel(sOptions)) 
            fieldNumber = 1;
            optionNumber = optionNumber + 1;  % next pair of options
        end
    else
        opName = varargin{optionNumber};
        opVal = varargin{optionNumber+1};
        optionNumber = optionNumber + 2;  % next pair of options
    end
    % interpret opName and opVal
    if strcmpi(opName,'imsize') 
        if ((numel(opVal) == 1) && isnumeric(opVal))
            M = max(1, floor(opVal));
            % L = ((M-Mn+Mi)/Mi)*((N-Nn+Ni)/Ni)
            N = 0; % to be evaluated later
        elseif ((numel(opVal) == 2) && isnumeric(opVal))
            M = max(1, floor(opVal(1)));
            N = max(1, floor(opVal(2)));
        else
            disp('mycol2im: illegal option imsize. We ignore it.');
        end
    end
    if (strcmpi(opName,'transform') || strcmpi(opName,'t'))
        tr = opVal;
    end
    if (strcmpi(opName,'size') || strcmpi(opName,'s'))
        if ((numel(opVal) == 1) && isnumeric(opVal))
            Ms = max(0, floor(opVal));
            Ns = Ms;
        elseif ((numel(opVal) == 2) && isnumeric(opVal))
            Ms = max(0, floor(opVal(1)));
            Ns = max(0, floor(opVal(2)));
        else
            error('mycol2im: illegal option size. We ignore it.');
        end
    end
    if (strcmpi(opName,'neighborhood') || strcmpi(opName,'n'))
        if ((numel(opVal) <= 2) && isnumeric(opVal))
            nei = ones(opVal);
            [Mn, Nn] = size(nei);
        elseif ((numel(opVal) > 2) && isnumeric(opVal))
            nei = opVal;
            [Mn, Nn] = size(nei);
        else
            error('mycol2im: illegal option neighborhood. We ignore it.');
        end
    end
    if (strcmpi(opName,'increment') || strcmpi(opName,'i'))
        if ((numel(opVal) == 1) && isnumeric(opVal))
            Mi = max(0, floor(opVal));
            Ni = Mi;
        elseif ((numel(opVal) == 2) && isnumeric(opVal))
            Mi = max(0, floor(opVal(1)));
            Ni = max(0, floor(opVal(2)));
        elseif strcmpi(opVal, 'distinct');
            Mi = 0;
            Ni = 0;
        elseif strcmpi(opVal, 'sliding');
            Mi = 1;
            Ni = 1;
        else
            error('mycol2im: illegal option increment. We ignore it.');
        end
    end
    if strcmpi(opName,'verbose') || strcmpi(opName,'v')
        if (islogical(opVal) && opVal); verbose = 1; end;
        if isnumeric(opVal); verbose = opVal(1); end;
    end
end

%% check and display options
if ((Ms == 0) || (Ns == 0))
    if strcmpi(tr,'none')
        Ms = 4; Ns = 4;
    else
        Ms = 8; Ns = 8;
    end
end
if (~iscell(tr) && ((strcmpi(tr,'lot') || strcmpi(tr,'elt'))))
    Ms = Ms - mod(Ms,2);   % Ms and Ns is made even
    Ns = Ns - mod(Ns,2);   
end
if ( iscell(tr) || ... 
        ( ischar(tr) && ... 
          ~strcmpi(tr,'none') && ...
          ~strcmpi(tr,'lot') && ...
          ~strcmpi(tr,'elt') && ...
          ~strcmpi(tr,'dct') )  )
    % a wavelet as in mylwt2
    Ms = max(2,2^floor(log2(Ms)));
    Ns = Ms;
end
if (isnumeric(tr)) % transform is given as a matrix
    Ms = size(tr,2);
    Ns = Ms;
    P = floor(size(tr,1)/Ms);  % overlap factor, should be integer
    if (Ms*P ~= size(tr,1))
        tr = tr(1:(Ms*P));
    end
end

if ((Mn == 0) || (Nn == 0))
    nei = ones(Ms,Ns);
    [Mn, Nn] = size(nei);
end
if ((Mi == 0) || (Ni == 0))
    Mi = Mn; 
    Ni = Nn;
end
% find size of restored image
% according to: L = ((M-Mn+Mi)/Mi)*((N-Nn+Ni)/Ni)
L = size(X,2);
if (M == 0)   
    factors = factor(L);
    M = prod(factors(1:2:end))*Mi + Mn - Mi;
end
if (N == 0)   
    N = (L*Mi/(M-Mn+Mi))*Ni + Nn - Ni;
end
if ~( L == (((M-Mn+Mi)/Mi)*((N-Nn+Ni)/Ni)) ) 
    disp(['mycol2im: Given  image size is ',int2str(M),'x',int2str(N)]);
    disp('Can not make given image size fit data X.');
    factors = factor(L);
    M = prod(factors(1:2:end))*Mi + Mn - Mi;
    N = (L*Mi/(M-Mn+Mi))*Ni + Nn - Ni;
end
            
if verbose
    disp(['mycol2im: Restored  image is ',int2str(M),'x',int2str(N)]);
end
        
%% put the blocks into the image and taking average
A = zeros(M,N);
Ac = zeros(M,N);   % counts
index = find(nei);
I = 1:Mi:(size(A,1)-Mi+1);
J = 1:Ni:(size(A,2)-Ni+1);
if ~(L == numel(I)*numel(J))
    disp('Can still not make image size fit data X.');
end
if verbose
    disp(['neighborhood is ',int2str(numel(index)),' pixels from ',...
        int2str(Mn),'x',int2str(Nn),' block.']);
    disp(['Get ',int2str(L),' columns using ',...
        'increment step ',int2str(Mi),' and ',int2str(Ni),'.']);
end
block = zeros(Mn,Nn);
Lb = size(X,1);
m = Mn-1; 
n = Nn-1;
k = 1;
for j=J
    for i=I
        block(index) = X(:,k);
        A(i:(i+m),j:(j+n)) = A(i:(i+m),j:(j+n)) + block;
        block(index) = ones(Lb,1);
        Ac(i:(i+m),j:(j+n)) = Ac(i:(i+m),j:(j+n)) + block;
        k = k+1;
    end
end
index = find(Ac);
A(index) = A(index)./Ac(index);

%% do the inverse transform
if strcmpi(tr,'none')
    if verbose
        disp('Do no transform.');
    end
elseif strcmpi(tr,'dct')
    if verbose
        disp(['Do IDCT with size ',int2str(Ms),'x',int2str(Ns),'.']);
    end
    if (Ns > 1)
        A = idct( reshape(A',Ns,M*N/Ns) ); % the rows
        A = reshape(A,N,M)';  % back to size MxN
    end
    if (Ms > 1)
        A = idct( reshape(A,Ms,M*N/Ms) );    % the columns (of each block)
        A = reshape(A,M,N);  % back to size MxN
    end
elseif strcmpi(tr,'lot')
    if verbose
        disp(['Do inverse LOT with size ',int2str(Ms),'x',int2str(Ns),'.']);
    end
    if (Ns > 1)
        F = getLOT(Ns, lotrho);
        F1 = F(1:Ns,:);
        F2 = F((Ns+1):(2*Ns),:);
        A = F1 * reshape(A',Ns,M*N/Ns) + ...
            F2 * reshape([A(:,(N-Ns+1):N),A(:,1:(N-Ns))]',Ns,M*N/Ns);
        A = reshape(A,N,M)';
    end
    if (Ms > 1)
        F = getLOT(Ms, lotrho);  % get the synthesis vectors of LOT
        F1 = F(1:Ms,:);
        F2 = F((Ms+1):(2*Ms),:);
        A = F1 * reshape(A,Ms,M*N/Ms) + ...
            F2 * reshape([A((M-Ms+1):M,:);A(1:(M-Ms),:)],Ms,M*N/Ms);
        A = reshape(A,M,N);
    end
elseif strcmpi(tr,'elt')
    if verbose
        disp(['Do inverse ELT with size ',int2str(Ms),'x',int2str(Ns),'.']);
    end
    if (Ns > 1)
        F = reshape(getELT(Ns,eltarg2,eltrho)',Ns,Ns,4); 
        A = F(:,:,1)' * reshape(A',Ns,M*N/Ns) + ...
            F(:,:,2)' * reshape([A(:,(N-Ns+1):N),A(:,1:(N-Ns))]',Ns,M*N/Ns) + ...
            F(:,:,3)' * reshape([A(:,(N-2*Ns+1):N),A(:,1:(N-2*Ns))]',Ns,M*N/Ns) + ...
            F(:,:,4)' * reshape([A(:,(N-3*Ns+1):N),A(:,1:(N-3*Ns))]',Ns,M*N/Ns);
        A = reshape(A,N,M)';
    end
    if (Ms > 1)
        F = reshape(getELT(Ms,eltarg2,eltrho)',Ms,Ms,4); 
        A = F(:,:,1)' * reshape(A,Ms,M*N/Ms) + ...
            F(:,:,2)' * reshape([A((M-Ms+1):M,:);A(1:(M-Ms),:)],Ms,M*N/Ms) + ...
            F(:,:,3)' * reshape([A((M-2*Ms+1):M,:);A(1:(M-2*Ms),:)],Ms,M*N/Ms) + ...
            F(:,:,4)' * reshape([A((M-3*Ms+1):M,:);A(1:(M-3*Ms),:)],Ms,M*N/Ms);
        A = reshape(A,M,N);
    end
elseif isnumeric(tr)   % tr is (Ms*P x Ns)
    if P==1  % no overlap
        F = tr';   % Ms x Ms  (= Ns x Ns), the transposed
        if verbose
            disp(['Do supplied transform with size ',...
                int2str(Ms),'x',int2str(Ms),'.']);
        end
    else   
        F = reshape(tr',Ms,Ms,P);  % Ms x Ms x P  % the transposed
        if verbose
            disp(['Do supplied overlapping transform with size ',...
                int2str(Ms),'x',int2str(Ms),' (P=',int2str(P),').']);
        end
    end
    %
    Y = F(:,:,1)' * reshape(A',Ns,M*N/Ns);
    for p=2:P
        Y = Y + F(:,:,p)' * reshape([A(:,(N-(p-1)*Ns+1):N),A(:,1:(N-(p-1)*Ns))]',Ns,M*N/Ns);
    end
    A = reshape(Y,N,M)';
    % 
    Y = F(:,:,1)' * reshape(A,Ms,M*N/Ms);
    for p=2:P
        Y = Y + F(:,:,p)' * reshape([A((M-(p-1)*Ms+1):M,:);A(1:(M-(p-1)*Ms),:)],Ms,M*N/Ms);
    end
    A = reshape(Y,M,N);
elseif (iscell(tr) || ischar(tr)) % a wavelet as in myilwt2
    if verbose
        if (iscell(tr)); temp = 'by given liftwave'; 
        else temp = tr;
        end
        disp(['Do inverse wavelet ',temp,', size ',int2str(Ms),'x',int2str(Ns),'.']);
    end
    A = myilwt2(A, tr, log2(Ms));
end


return