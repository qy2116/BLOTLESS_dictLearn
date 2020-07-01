% MainTest - denoise an image
% See section V-D2 in our paper.

clear; close all; clc;
rand('state',0);
randn('state',0);

flag_fig = 1; % Draw images or not
Param.k = 128; % number of atoms in the dictionary
Param.noise = 30; % noise level
Param.method = 'BLOTLESS';
ImageName = '6.pgm'; %image source

OriginalImage=im2double(imread(ImageName)); 
OriginalImage = OriginalImage*255;

NoisedImage=OriginalImage+Param.noise*randn(size(OriginalImage));

% Denoise the corrupted image using learned dicitionary from corrupted image 
[DenoisedImage1, timecost1] = denoiseImage(NoisedImage, Param);
if flag_fig
figure;imshow(DenoisedImage1/255,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,92,152])
axis normal;
end

NoisedPSNR = 20*log10(255/sqrt(mean((NoisedImage(:)-OriginalImage(:)).^2)));
DenoisedPSNRblotless = 20*log10(255/sqrt(mean((DenoisedImage1(:)-OriginalImage(:)).^2)));


Param.method = 'SimCO';
[DenoisedImage2, timecost2] = denoiseImage(NoisedImage, Param);
DenoisedPSNRsimco = 20*log10(255/sqrt(mean((DenoisedImage2(:)-OriginalImage(:)).^2)));

if flag_fig
figure;imshow(DenoisedImage2/255,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,92,152])
axis normal;
end

Param.method = 'KSVD';
[DenoisedImage3, timecost3] = denoiseImage(NoisedImage, Param);
DenoisedPSNRksvd = 20*log10(255/sqrt(mean((DenoisedImage3(:)-OriginalImage(:)).^2)));

if flag_fig
figure;imshow(DenoisedImage3/255,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,92,152])
axis normal;
end

Param.method = 'MOD';
[DenoisedImage4, timecost4] = denoiseImage(NoisedImage, Param);
DenoisedPSNRmod = 20*log10(255/sqrt(mean((DenoisedImage4(:)-OriginalImage(:)).^2)));

if flag_fig
figure;imshow(DenoisedImage4/255,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,92,152])
axis normal;
end
