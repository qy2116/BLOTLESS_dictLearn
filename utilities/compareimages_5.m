function compareimages_5(A,ATitle,B,BTitle,C,CTitle,D,DTitle,E,ETitle)
%COMPAREIMAGES   Displays five images side by side with linked axes
%   COMPAREIMAGES(A,B) displays images A and B, where A and B are either
%   grayscale or RGB color images with values in [0,1].  The images are
%   displayed with linked axes for convenient panning and zooming.
%
%   COMPAREIMAGES(A,'A title',B,'B title') specifies titles above the
%   images.
%
%   See also linkaxes.

% Pascal Getreuer 2009

ax(1) = subplot(1,5,1);
hold off

if ndims(A) == 2
    imagesc(A);
    colormap(gray(256));
elseif ndims(A) == 3
    image(min(max(A,0),1));
end

% set(gca,'Units','Normalized','Position',[0,0.1,0.5,0.8]);
axis image
axis off
title(ATitle);
zoom;

ax(2) = subplot(1,5,2);
hold off

if ndims(B) == 2
    imagesc(B);
    colormap(gray(256));
elseif ndims(B) == 3
    image(min(max(B,0),1));
end

% set(gca,'Units','Normalized','Position',[0.2,0.1,0.5,0.8]);
axis image
axis off
title(BTitle);


ax(3) = subplot(1,5,3);
hold off

if ndims(C) == 2
    imagesc(C);
    colormap(gray(256));
elseif ndims(C) == 3
    image(min(max(C,0),1));
end

% set(gca,'Units','Normalized','Position',[0.4,0.1,0.5,0.8]);
axis image
axis off
title(CTitle);

ax(4) = subplot(1,5,4);
hold off

if ndims(D) == 2
    imagesc(D);
    colormap(gray(256));
elseif ndims(D) == 3
    image(min(max(D,0),1));
end

% set(gca,'Units','Normalized','Position',[0.6,0.1,0.5,0.8]);
axis image
axis off
title(DTitle);


ax(5) = subplot(1,5,5);
hold off

if ndims(E) == 2
    imagesc(E);
    colormap(gray(256));
elseif ndims(E) == 3
    image(min(max(E,0),1));
end

% set(gca,'Units','Normalized','Position',[0.8,0.1,0.5,0.8]);
axis image
axis off
title(ETitle);


linkaxes(ax,'xy');
