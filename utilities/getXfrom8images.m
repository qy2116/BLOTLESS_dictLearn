function X = getXfrom8images(varargin)
% getXfrom8images Get (training) vectors generated from 8 (USC) images
%
% X = getXfrom8images(argName, argVal, ...);
%----------------------------------------------------------------------
%        The input arguments (options) may be:
%  'getFixedSet' is 1 (true) or 0 (false), default is 0.
%        If true the random state is set to a fixed value first
%  'transform' or 't' is as in myim2col, default is 'none',
%        Other values are: 'dct', 'lot', 'elt', 'db79', 'm79', ...
%  'myim2colPar' more parameters to use in myim2col, default: 
%        struct('s',[8,8], 'n',[8,8], 'i',[8,8], 'a','none', 'v',0);
%  'noFromEachIm' how many training vectors to get from each image
%        default: 1500
%  'images', a cell of strings with the filenames (images). Default
%        {'elaine.bmp', 'lake.bmp', 'man.bmp', 'crowd.bmp',...
%         'couple.bmp', 'woman1.bmp', 'woman2.bmp', 'baboon.bmp'};
%  'imCat', a cell of strings with the catalog names where the images
%        maybe located. All images should be in the same catalog, and the
%        first catalog where the first image is found is used. Default:
%        {my_matlab_path('USC'), pwd} 
%  'verbose' or 'v' to indicate verboseness, default 0
%----------------------------------------------------------------------
% example:
%  X = getXfrom8images('t','m79');
%  X = getXfrom8images('t','m79', 'getFixedSet',1, 'v',1);

%----------------------------------------------------------------------
% Copyright (c) 2010.  Karl Skretting.  All rights reserved.
% University of Stavanger (Stavanger University), Signal Processing Group
% Mail:  karl.skretting@uis.no   Homepage:  http://www.ux.his.no/~karlsk/
% 
% HISTORY:  dd.mm.yyyy
% Ver. 1.0  19.01.2010  Made function
% Ver. 1.1  11.01.2013  may use my_matlab_path(..) to set imCat
%----------------------------------------------------------------------

%% default options
getFixedSet = 0;
transform = 'none';       
myim2colPar = struct('s',[8,8], 'n',[8,8], 'i',[8,8], 'a','none', 'v',0);
noFromEachIm = 1500;
images = {'elaine.bmp','lake.bmp','man.bmp','crowd.bmp',...
          'couple.bmp','woman1.bmp','woman2.bmp','baboon.bmp'};
% imCat = {my_matlab_path('USC'), pwd};  % {pwd, 'd:/bilder/USC_tiff'};  
imCat = {pwd,'G:/01-Study/08-Temp Papers/L-S method/Code/matlab_code/USCimages_bmp'};
verbose = 0;

%%  get the options
nofOptions = nargin;
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
    if strcmpi(opName,'getFixedSet') 
        if (islogical(opVal) && opVal); getFixedSet = 1; end;
        if isnumeric(opVal); getFixedSet = opVal(1); end;
    end
    if (strcmpi(opName,'transform') || strcmpi(opName,'t'))
        transform = opVal;
    end
    if strcmpi(opName,'myim2colPar') 
        if (iscell(opVal) || isstruct(opVal))
            myim2colPar = opVal;
        else
            error('getXfrom8images: illegal option for myim2colPar, it is ignored.');
        end
    end
    if ( (strcmpi(opName,'noFromEachIm')) && isnumeric(opVal) )
        noFromEachIm = opVal(1);
    end
    if ( (strcmpi(opName,'images')) && iscell(opVal) )
        images = opVal;
    end
    if ( (strcmpi(opName,'images')) && ischar(opVal) )
        images = cell(1,size(opVal,1));
        for i=1:numel(images)
            images{i} = strtrim(opVal(i,:));
        end
    end
    if ( (strcmpi(opName,'imCat')) && iscell(opVal) )
        imCat = opVal;
    end
    if ( (strcmpi(opName,'imCat')) && ischar(opVal) )
        imCat = cell(1,size(opVal,1));
        for i=1:numel(imCat)
            imCat{i} = strtrim(opVal(i,:));
        end
    end
    if strcmpi(opName,'verbose') || strcmpi(opName,'v')
        if (islogical(opVal) && opVal); verbose = 1; end;
        if isnumeric(opVal); verbose = opVal(1); end;
    end
end

%% check where to find the images
for i=1:numel(imCat)
    catalog = imCat{i};
    if exist([catalog,'/',images{1}],'file'); break; end
end
for i=1:numel(images)
    if ~exist([catalog,'/',images{i}],'file'); 
        error(['getXfrom8images: did not find ',catalog,images{i},'.']);
    end
end
% catalog ok    
if verbose
    disp(['getXfrom8images: Use ',int2str(numel(images)),...
        ' images from ',catalog]);
    disp(['Get ',int2str(noFromEachIm),' vectors from each image.']);
end

if (getFixedSet == 1)
    if verbose
        disp('Use a fixed dataset, i.e. set random state to a fixed value.');
    end
    if (exist('RandStream','class') == 8) 
        % make the random generator go to a particular state and back
        defaultStream = RandStream.getDefaultStream;
        savedState = defaultStream.State;
        s = RandStream('mt19937ar', 'Seed', 5489);
        RandStream.setDefaultStream(s);
    else
        rand('twister',5489);
    end
    imOffset = [0,0];
end

%% make the training data
L = numel(images)*noFromEachIm;
X = zeros(prod(myim2colPar.n), L);
if exist('ME','var'); clear('ME'); end;
for i = 1:numel(images)
    imFilename = [catalog,'/',images{i}];
    try
        A = double( imread(imFilename) ) - 128;  % subtract global 'mean'
        if (getFixedSet == 0)  % random offset !
            imOffset = floor(8*rand(1,2));
        end
        Xi = myim2col( A, myim2colPar, 't',transform, 'offset',imOffset );
        I = randperm(size(Xi,2));
        X(:, i:numel(images):L ) = Xi(:, I(1:noFromEachIm) );   % some random vectors
    catch ME
        disp(['getXfrom8images: skip image ',imFilename]);
        disp(['An error occurred: ',ME.message]);
    end
end
if exist('ME','var')
    X = X(:,sum(X)~=0);
end

% reset random state    
if (getFixedSet == 1) 
    if (exist('RandStream','class') == 8) 
        RandStream.setDefaultStream(defaultStream);
        defaultStream.State = savedState;
    else
        rand('twister',sum(100*clock));
    end
end


return