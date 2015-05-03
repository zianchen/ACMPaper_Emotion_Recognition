function [variety] = example_color_naming(im)
% load the word to color names matrix. The words are a 32x32x32 grid on the sRGB space. 

load('w2c.mat');

% first example
im=double(imread(im));       % load test image

% compute the color name assignment for all pixels in image im:
variety = im2c(im,w2c,-1);                % using im2c(im,w2c,0) is much faster
variety = cell2mat(variety);

%figure(1);
%subplot(1,2,1);imshow(uint8(im));
%subplot(1,2,2);imshow(uint8(out));
%sprintf('out :%f\n',out)

% second example:
%im2=double(imread('opp_color_circle.tif'));     

%out2=im2c(im2,w2c,-1);   

%figure(2);
%subplot(1,2,1);imshow(uint8(im2));
%subplot(1,2,2);imshow(uint8(out2));