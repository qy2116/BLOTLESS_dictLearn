function [IOut,timecost] = denoiseImage(Image,Param)

[NN1,NN2] = size(Image);
% C = 1.15;
C = 1.15;
bb = 8;
maxNumBlocksToTrainOn = 1000;
sigma = Param.noise;
K = Param.k;


% first, train a dictionary on blocks from the noisy image

if(prod([NN1,NN2]-bb+1)> maxNumBlocksToTrainOn)
    randPermutation =  randperm(prod([NN1,NN2]-bb+1));
    selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);

    blkMatrix = zeros(bb^2,maxNumBlocksToTrainOn);
    for i = 1:maxNumBlocksToTrainOn
        [row,col] = ind2sub(size(Image)-bb+1,selectedBlocks(i));
        currBlock = Image(row:row+bb-1,col:col+bb-1);
        blkMatrix(:,i) = currBlock(:);
    end
else
    blkMatrix = im2col(Image,[bb,bb],'sliding');
end
% load('tdata4.mat');
% blkMatrix = 255*Y_temp(:,1:K);clear Y_temp

param.K = K;
param.I= 1:K;
param.errorGoal = sigma*C;

% make initial dictionary
Pn=ceil(sqrt(K));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end
    DCT(:,k+1)=V/norm(V);
end
DCT=kron(DCT,DCT);

param.initialDictionary = DCT(:,1:param.K );


% param.initialDictionary = randn(bb^2,param.K);
% param.initialDictionary = column_normalize(param.initialDictionary,2);

%reducedc
% blkMatrix = column_normalize(blkMatrix,2);
vecOfMeans = mean(blkMatrix);
blkMatrix(:,2:end) = blkMatrix(:,2:end)-ones(size(blkMatrix(:,2:end),1),1)*vecOfMeans(2:end);
blkMatrix = column_normalize(blkMatrix,2);

param.data = blkMatrix;
% param.Tdata = 5;
param.Edata = 0.01;
param.dictsize = K;
param.iternum = 20;
% param.initdict = DCT(:,1:param.K );
param.memusage = 'high';

time_start = clock;

if strcmp(Param.method, 'KSVD')
    load('KSVD_dic.mat')
    Dictionary = Dksvd;
end

if strcmp(Param.method, 'SimCO')
    load('SimCO_dic.mat')
    Dictionary = D_simco;
end

if strcmp(Param.method, 'PSimCO')
    [Dictionary] = PSimCO(blkMatrix, param);
end

if strcmp(Param.method, 'MOD')
    load('MOD_dic.mat');
    Dictionary = D_mod;
end
if strcmp(Param.method, 'BLOTLESS')
    load('Blotless_dict.mat');
    Dictionary = D_Blot;
end
if strcmp(Param.method, 'None')
    [Dictionary] = param.initialDictionary;
end


time_end = clock;
timecost = etime(time_end,time_start);

% denoise the image using the resulted dictionary
errT = sigma*C;

%blocks = im2col(Image,[NN1,NN2],[bb,bb],'sliding');

blocks = im2col(Image,[bb,bb],'sliding');
idx = 1:size(blocks,2);

% go with jumps of 30000
for jj = 1:30000:size(blocks,2)
  
    jumpSize = min(jj+30000-1,size(blocks,2));
    
    %reduceDC
    vecOfMeans = mean(blocks(:,jj:jumpSize));
    blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
    
    
    Coefs = OMPerr(Dictionary,blocks(:,jj:jumpSize),errT);
   
    %reducedc
    blocks(:,jj:jumpSize)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
    
    %blocks(:,jj:jumpSize)= Dictionary*Coefs ;  
    
end

count = 1;
Weight = zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(Image)-bb+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);        
%     block =reshape(blocks(:,count),[bb,bb]);
    block =reshape(blocks(:,count),[bb,bb]);
    IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;
    Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);
    count = count+1;
end

% IOut = (Image+0.034*sigma*IMout)./(1+0.034*sigma*Weight);
IOut = (Image+0.034*sigma*IMout)./(1+0.034*sigma*Weight);
% IOut = (IMout)./(Weight);


