function res = dictshow(D, varargin)
% dictshow        Show atoms (1D or 2D) for a dictionary 
% 
% General use with an NxK matrix as input returns a handle to the figure
%   h = dictshow(D, argName, argVal, ...);         
% 
% Several special cases may be used as well. These return a handle if the
% atoms are 1D and a matrix (image) if the atoms are 2D, that is represent
% an image patch or a block of coefficients. The first argument should in
% these cases be a struct (containing the dictionary as a field D), or a
% character string containing the name of a mat-file where the dictionary
% is stored. Examples are:
%   Dstruct = load('ex312Jan191930.mat');
%   A = dictshow(Dstruct);     
%   A = dictshow('ex311Feb031626.mat');     
% Special variants where the first argument is 'dct','lot','elt' or 'm79'
% shows the basis images for the corresponding 8x8 coefficient block. Ex:
%   A = dictshow('m79');
%-------------------------------------------------------------------------
% Output argument may be 
%   h           a handle to the figure
%   A           a matrix representing the image of 2D atoms
%   em          an error message (char)
% First input argument may be 
%   D           an NxK matrix representing the dictionary
%   Dstruct     a struct with variables (as in a dictionary mat-file)
%   filename    a string for the name of mat-file 
%   transform   a string, 'dct','lot','elt' or 'm79' (9/7 wavelet)
%               to show the (N1*N2) separable basis images
% There may be an additional number of input arguments, given as a struct,
% a cell or as pairs of the form: argName, argVal, ... 
%   'size'      [N1,N2]: for 2D atoms (N = N1*N2), 
%               N1: for the case where N2==N1 only N1 may be given 
%   'useK'      a list of indexes, plot only the listed atoms
%   'plotsize'  [K1,K2] plot the K, or numel(useK), atoms in a K1xK2 grid
%   'plottype'  0 - make no plot
%               1 - plot as image, the only alternative if N2>1
%               2 - (N2==1), plot atoms (columns of D) as vertical lines
%               3 - (for N2==1), plot atoms as horizontal lines
%               4 - (for N2==1), plot atoms as 'discrete sequence' plot
%-------------------------------------------------------------------------
% Examples: 
%   A = dictshow('dct');
%   A = dictshow('dct','size',[4,6]);  % size of transform [vert, hor]
%   A = dictshow('lot','size',[10,10]); 
%   A = dictshow('ex311Feb031626.mat');  % all atoms
%   K1=12; A = dictshow('dict_sine.mat','plotsize',[K1,ceil(256/K1)]);
%   Dstruct = load('ex312Jan191930.mat');
%   A = dictshow(Dstruct,'useK',[1:20,100:135]);  % selected atoms
%   D = get8x21sine();
%   dictshow(D,'plottype',2);
%   dictshow(D,'plottype',3, 'plotsize',[6,4])
%   dictshow(D,'plottype',4, 'plotsize',[7,3])
%   dictshow(get16x144(), 'plottype',3, 'plotsize',[9,16])
%   dictshow(get16x144(), 'plottype',1, 'size',[4,4], 'plotsize',[9,16]);

%% ---------------------------------------------------------------------
% Copyright (c) 2009.  Karl Skretting.  All rights reserved.
% University of Stavanger.
% Mail:  karl.skretting@uis.no   Homepage:  http://www.ux.uis.no/~karlsk/
% 
% HISTORY:  dd.mm.yyyy
% Ver. 1.0  27.01.2010  KS: function made (as ex324.m)
% Ver. 2.0  31.05.2011  KS: new function made based on ex324.m
% Ver. 3.0  23.09.2011  KS: tried to make the function cleaner without
%                       reducing the flexibility
% ---------------------------------------------------------------------

mfile = 'dictshow';
if (nargin < 1)
    error('dictshow: wrong number of arguments, see help.');
end

%% get the options first
N1 = 8; N2 = 8;  % a N1 x N2 transform, default size
plottype = 1;
titletext = '';

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
    if strcmpi(opName,'size')
        N1 = floor(opVal(1));
        if numel(opVal > 1)
            N2 = floor(opVal(2));
        else
            N2 = N1;
        end
    end
    if strcmpi(opName,'useK')
        useK = floor(opVal);
    end
    if strcmpi(opName,'plotsize')
        K1 = floor(opVal(1));
        if numel(opVal > 1)
            K2 = floor(opVal(2));
        else
            K2 = K1;
        end
    end
    if strcmpi(opName,'plottype')
        plottype = floor(opVal(1));
    end
    if strcmpi(opName,'makeFigure')   % obsolete parameter
        if isnumeric(opVal) 
            plottype = floor(opVal(1)); 
        elseif islogical(opVal) 
            plottype = double(opVal);   % 0 or 1
        else
            disp([mfile,': Option makeFigure should be logical (or numeric).']);
        end
    end
    if strcmpi(opName,'title') 
        if ischar(opVal) 
            titletext = opVal; 
        else
            disp([mfile,': Option title should be a string (char).']);
        end
    end
end

%% check (and process) the first argument D. Take the exceptions first
disp([mfile,': Class of first argument is ',class(D)]);

if (ischar(D) && strcmpi(D,'dct'))
    X = eye(N1*N2);
    A = mycol2im(X, 'imsize',[N1^2, N2^2], 'transform','dct', 'size',[N1,N2]);
    A = (A-min(A(:)))/(max(A(:))-min(A(:)));  % range is 0 to 1
    res = resFromA(A, N1, N2, 1, 0);
    if numel(titletext)==0
        titletext = ['Basis images, or atoms, for ',int2str(N1),'x',int2str(N2),' DCT.'];
    end
elseif (ischar(D) && strcmpi(D,'lot'))
    P = 2;  % overlap factor
    X = setXsparse(N1,N2,P);
    A = mycol2im(X, 'imsize',[N1*(P*N1+P-1), N2*(P*N2+P-1)], ...
        'transform','lot', 'size',[N1,N2]);
    A = (A-min(A(:)))/(max(A(:))-min(A(:)));  % range is 0 to 1
    res = resFromA(A, N1, N2, P, P-1);
    if numel(titletext)==0
        titletext = ['Basis images, or atoms, for ',int2str(N1),'x',int2str(N2),...
            ' LOT, overlap factor P = ',int2str(P),'.'];
    end
elseif (ischar(D) && strcmpi(D,'elt'))
    P = 4;  % overlap factor
    X = setXsparse(N1,N2,P);
    A = mycol2im(X, 'imsize',[N1*(P*N1+P-1), N2*(P*N2+P-1)], ...
        'transform','elt', 'size',[N1,N2]);
    A = (A-min(A(:)))/(max(A(:))-min(A(:)));  % range is 0 to 1
    res = resFromA(A, N1, N2, P, P-1);
    if numel(titletext)==0
        titletext = ['Basis images, or atoms, for ',int2str(N1),'x',int2str(N2),...
            ' ELT, overlap factor P = ',int2str(P),'.'];
    end
elseif (ischar(D) && strcmpi(D,'m79'))
    if (N1 ~= 8) || (N2 ~= 8)
        disp([mfile,': only show 8x8 basis images for 9/7 wavelet.']);
        N1 = 8; N2 = 8;   
    end
    P = 8;  % overlap factor
    offset = 3; 
    X = setXsparse(N1,N2,P);   % 71x71 blocks (each is 8x8)
    A = mycol2im(X, 'imsize',[N1*(P*N1+P-1), N2*(P*N2+P-1)], ...
        'transform','m79', 'size',[N1,N2]);
    A(A==0) = max(A(:));
    A = (A-min(A(:)))/(max(A(:))-min(A(:)));  % range is 0 to 1
    % res = resFromA(A, N1, N2, P, offset); % make all 64, (many duplicates)
    % no duplicates below
    res = ones(N1*P*4,N2*P*4);
    i1 = 1:(P*N1);
    for k1 = [0,4,6,7]
        j1 = (1:(P*N1))+(k1*(P*N1)+offset*N1);
        i2 = 1:(P*N2);
        for k2 = [0,4,6,7]
            if (((k1*k2) == 0) || (k1 == k2))
                j2 = (1:(P*N2))+(k2*(P*N2)+offset*N2);
                res(i1,i2) = A(j1,j2);
            end
            i2 = i2 + (P*N2);
        end
        i1 = i1 + (P*N1);
    end
    if numel(titletext)==0
        titletext = 'Basis images, or atoms, for 3-level 9/7 wavelet.';
    end
%    
% D is filename (char) or a struct
elseif (ischar(D) || isstruct(D))
    % if it is a filename, the file is loaded
    if ischar(D)
        if (exist(D,'file')==2)
            Ds = load(D);
        else
            res = [mfile,': filename ',D,' does not exist.'];
            disp(res);
            return                                % ============ > RETURN
        end
    else
        Ds = D;
    end
    % now Ds should be a struct, check field(s)
    tc = fieldnames(Ds);   % a cell array 
    while (numel(tc)==1) && isstruct( Ds.(tc{1}) )
        Ds = Ds.(tc{1});
        tc = fieldnames(Ds);   % a cell array 
    end
    if ~isfield(Ds,'D')
        res = [mfile,': Can not find a dictionary in given stucture or file.'];
        disp(res);
        return                                    % ============ > RETURN
    end
    %
    if ~(isfield(Ds,'transform')) 
        Ds.transform = 'none';
    end
    N = size(Ds.D,1);
    K = size(Ds.D,2);
    if (N ~= (N1*N2))   % size does not match
        N1 = N;  
        N2 = 1;
    end
    if ~exist('useK','var')
        useK = 1:K;   % all
    end
    if ~exist('K1','var')
        K1 = floor(sqrt(numel(useK)));
        K2 = ceil(numel(useK)/K1);
    end
    % 
    % prepare image A to be shown
    P = 1;
    tt = [int2str(N1),'x',int2str(N2),' blocks in domain ',Ds.transform,...
        ', ',int2str( min(numel(useK),K1*K2) ),' atoms shown.'];
    if (strcmpi(Ds.transform,'lot'));
        P=2;
        tt = [int2str(2*N1),'x',int2str(2*N2),' blocks in domain ',Ds.transform,...
            ', ',int2str(numel(useK)),' atoms shown.'];
    elseif (strcmpi(Ds.transform,'elt'));
        tt = [int2str(4*N1),'x',int2str(4*N2),' blocks in domain ',Ds.transform,...
            ', ',int2str(numel(useK)),' atoms shown.'];
        P=4;
    elseif (strcmpi(Ds.transform,'m79')) && (N1==8) && (N2==8);
        P=8;
        tt = ['51x51 blocks of ',Ds.transform,...
            ' (8x8 overlap), ',int2str(numel(useK)),' atoms shown.'];
    end
    if isfield(Ds,'ResultFile')
        tt = {tt; ['Dictionary file: ',Ds.ResultFile]};
    end
    if isfield(Ds,'resultFile')
        tt = {tt; ['Dictionary file: ',Ds.resultFile]};
    end
    if numel(titletext)==0
        titletext = tt;
    end
    %
    X = setXuseD(Ds.D,useK,K1,K2,P);
    A = mycol2im(X, 'imsize',[N1*(P*K1+P-1), N2*(P*K2+P-1)], ...
        'transform',Ds.transform, 'size',[N1,N2]);
    A(A==0)  = max(A(:));     % why this ???
    A = (A-min(A(:)))/(max(A(:))-min(A(:)));  % range is 0 to 1
    if (P < 8)
        res = ones(K1*P*8+K1+1,K2*8*P+K2+1);
        for k1 = 0:(K1-1)
            i1 = (1:(P*8))+(k1*(P*8+1)+1);
            j1 = (1:(P*8))+(k1*(P*8)+(P-1)*8);
            for k2 = 0:(K2-1)
                i2 = (1:(P*8))+(k2*(P*8+1)+1);
                j2 = (1:(P*8))+(k2*(P*8)+(P-1)*8);
                res(i1,i2) = A(j1,j2);
            end
        end
    else
        res = A;
    end
    D = Ds.D;
%
else
    %  D is numerical, NxK array
    N = size(D,1);
    K = size(D,2);
    if numel(titletext)==0
%         titletext = ['The given matrix D has size ',int2str(N),'x',int2str(K)];
    end
    disp([mfile,': ',titletext]);
end

if (plottype == 1) && ~exist('res','var') && isnumeric(D)    
    % prepare to plot D as image
    if (N ~= (N1*N2))   % size does not match
        N1 = N;  
        N2 = 1;
    end
    if ~exist('useK','var')
        useK = 1:K;   % all
    end
    if ~exist('K1','var')
        K1 = floor(sqrt(numel(useK)));
        K2 = ceil(numel(useK)/K1);
    end
    % 
    % prepare image A to be shown
%     tt = [int2str(N1),'x',int2str(N2),' blocks in ''pixel'' domain ',...
%         ', ',int2str( min(numel(useK),K1*K2) ),' atoms shown.'];
    tt = [int2str(N1),'x',int2str(N2),' block',...
        ', ',int2str( min(numel(useK),K1*K2) ),' atoms'];
    if numel(titletext)==0
        titletext = tt;
    end
    %
    X = setXuseD(D,useK,K1,K2,1);
    A = mycol2im(X, 'imsize',[N1*K1, N2*K2], ...
        'transform','none', 'size',[N1,N2]);
    % A(A==0)  = max(A(:));
    A = (A-min(A(:)))/(max(A(:))-min(A(:)));  % range is 0 to 1
    res = ones(K1*N1+K1+1,K2*N2+K2+1);
    for k1 = 0:(K1-1)
        i1 = (1:N1) + k1*(N1+1) + 1;
        j1 = (1:N1) + k1*N1;
        for k2 = 0:(K2-1)
            i2 = (1:N2) + k2*(N2+1) + 1;
            j2 = (1:N2) + k2*N2;
            res(i1,i2) = A(j1,j2);
        end
    end
end

%% plot, make figure
if (plottype == 1) && exist('res','var') && (numel(res) > 1)
    % A was made as res above, plot it as image
    clf;
    imagesc(res);              
    colormap(gray);
    axis equal
    axis off
    title(titletext);
    if (size(res,1)==256) && (size(res,2)==256)  
        % i.e. (ischar(Dstruct) && strcmpi(Dstruct,'m79'))
        t = {'LLL','LLH','LH','H'};
        ytab = [60, 132, 185, 242];
        for k1 = 1:4
            x = k1*64-30;
            for k2 = 1:4
                y = ytab(k2);
                if ((k1 == 1) || (k2 == 1) || (k1 == k2))
                    len = min(length(t{k1}),length(t{k2}));
                    h = text(x,y,[t{k2}(1:len),' -- ',t{k1}(1:len)]);
                    set(h,'HorizontalAlignment','center');
                end
            end
        end
        h = text(128,220,'Filtering in different levels,');
        set(h,'HorizontalAlignment','center');
        h = text(128,230,'Vertical band -- Horizontal band');
        set(h,'HorizontalAlignment','center');
    end
end

if isnumeric(D)
    N = size(D,1);
    K = size(D,2);
    if ~exist('useK','var')
        useK = 1:K;   % all
    end
    scale = 0.5;
    for k2 = 2:length(useK)
        scale = max(scale, max(D(:,useK(k2-1))-D(:,useK(k2))) );
    end
    scale = floor(20/scale)/20;
end

if (plottype > 1) && isnumeric(D)
    % prepare for the figure
    if nargout > 0
        res = figure(1);
    else
        figure(1);
    end
    clf;
    hold on;
    if ~exist('K1','var')
        K1 = floor(sqrt(numel(useK)));
        K2 = ceil(numel(useK)/K1);
    end
end
if (plottype == 2) && isnumeric(D) % plot D
    xpos=1;
    for k=useK
        left=1;
        right=N;
        while ~D(left,k);
            left=left+1; if left>right; break; end;
        end
        while ~D(right,k);
            right=right-1; if left>right; break; end;
        end
        plot(scale*D(left:right,k)+xpos, left:right, '-b');
        plot(xpos*ones(1,right-left+1), left:right, '-k');
        for i=left:right
            plot([xpos,scale*D(i,k)+xpos],[i,i],'-b');
        end
        xpos = xpos+1;
    end
    % number the vectors
    xpos = 1-0.1;
    ypos = N + N/20;
    if length(useK)<21
        for k=useK
            text(xpos,ypos,int2str(k));
            xpos=xpos+1;
        end
    else
        for k=1:ceil(length(useK)/20):length(useK)
            text(xpos,ypos,int2str(k));
            xpos=xpos+ceil(length(useK)/20);
        end
    end
    ypos=ypos+N/20;
    axis([-1,length(useK)+1,1-N/20,ypos]);
    %
    axis off;
    if numel(titletext)>0
        title(titletext);
    end
    % 
    % t1=['the N-dimension, N=',int2str(N)];
    text(1,ypos,['Vector number   (',mfile,': ',datestr(now),')']);
    set(gca,'ydir','reverse');
    % set(gcf,'PaperType','a4letter');
    % set(gcf,'Position',[300 300 800 600]);
    hold off;
end

if (plottype == 3) && isnumeric(D) % plot D
    dx = N*0.2;
    dy = 1.2;
    for uk = 1:numel(useK)
        k = useK(uk);
        ypos = dy*(K1-ceil((uk-0.5)/K2));
        xpos = (dx+N)*mod(uk-1,K2);
        left=1; right=N;
        plot((left:right)+xpos,scale*D(left:right,k)+ypos,'-b');
        xpos = xpos-dx/2;
        if (numel(useK) < 65)
            H = text(xpos,ypos,int2str(k));
            set(H,'FontSize',8);
        end
    end
    axis off;
    if numel(titletext)>0
        title(titletext);
    end
    % set(gcf,'PaperType','a4letter');
    % set(gcf,'Position',[300 300 800 600]);
    hold off;
end

if (plottype == 4) && isnumeric(D) % plot D
    dx = N*0.2;
    dy = 1.2;
    hold on;
    q1=1;q2=1;       % ypos and xpos
    for uk=1:length(useK)
        k=useK(uk);
        ypos = dy*(K1-q1);
        xpos = (dx+N)*(q2-1);
        left=1;right=N;
        plot((left:right)+xpos,0*(left:right)+ypos,'-b');
        for n=left:right
            plot([n+xpos,n+xpos],[ypos,scale*D(n,k)+ypos],'-b');
            plot(n+xpos,scale*D(n,k)+ypos,'b.');
        end
        xpos=xpos-dx/2;
        H=text(xpos,ypos,int2str(k));
        set(H,'FontSize',8);
        q2=q2+1;
        if q2>K2; q2=1; q1=q1+1; end;
    end
    axis off;
    if numel(titletext)>0
        title(titletext);
    end
    % set(gcf,'PaperType','a4letter');
    % set(gcf,'Position',[100 300 800 300]);
    % set(gcf,'PaperPosition',[0.2500 2.5000 8 3]);
    hold off;
end

return

%% *******************  subfunctions
%
function X = setXuseD(D,useK,K1,K2,P)
i1 = (P*K1+P-1);
i2 = (P*K2+P-1);
X = zeros(size(D,1), i1*i2);
k = 1;
for k1 = 0:(K1-1)
    j1 = (k1+1)*P-1;
    for k2 = 0:(K2-1)
        j2 = (k2+1)*P-1;
        if (k <= numel(useK)) && (useK(k) <= size(D,2)) 
            X( :,j2*i1+j1+1 ) = D(:,useK(k));
        end
        k = k+1;
    end
end
return

function X = setXsparse(K1,K2,P)
% if P=1, this is like X = eye(K1*K2)
i1 = (P*K1+P-1);
i2 = (P*K2+P-1);
X = zeros(K1*K2, i1*i2);
for k1 = 0:(K1-1)
    j1 = (k1+1)*P-1;
    for k2 = 0:(K2-1)
        j2 = (k2+1)*P-1;
        X( k2*K1+k1+1,j2*i1+j1+1 ) = 1;
    end
end
return

function res = resFromA(A, K1, K2, P, offset)
% disp([size(A), K1, K2, P, offset]);
res = ones(K1*P*K1+K1+1,K2*K2*P+K2+1);
for k1 = 0:(K1-1)
    i1 = (1:(P*K1))+(k1*(P*K1+1)+1);
    j1 = (1:(P*K1))+(k1*(P*K1)+offset*K1);
    for k2 = 0:(K2-1)
        i2 = (1:(P*K2))+(k2*(P*K2+1)+1);
        j2 = (1:(P*K2))+(k2*(P*K2)+offset*K2);
        res(i1,i2) = A(j1,j2);
    end
end
return