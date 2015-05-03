function [RFA] = emphasis(img)
img = imread('C:\Users\chenzian\Documents\MATLAB\Image\happybaby.jpg');
saliency_gbvs = gbvs(img);

% saliency_gbvs.master_map_resized is this saliency map interpolated (bicubic) 
% to the resolution of the original image.

% ?????gbvs???itt???? http://www.vision.caltech.edu/~harel/share/gbvs.php
%?????
%out_itt = ittikochmap(img);
%saliency_gbvs.master_map
%figure;
%subplot(2, 3, 2); imshow(img);
%subplot(2, 3, 4); imshow(saliency_gbvs.master_map);
%subplot(2, 3, 6); imshow(out_itt.master_map);
%show_imgnmap(img,out_gbvs);


%RGB->HSV  ???H??????S?????V?
HSV = rgb2hsv(img);
H = HSV(:, :, 1);
S = HSV(:, :, 2);
V = HSV(:, :, 3);

double_V = double(V);
%?????????
threshold = graythresh(double_V);
Mask = im2bw(double_V,threshold);

%imshow(Mask);
%subplot(2, 1, 2); imshow(BW);
%imwrite(Mask,'/Users/chenzian/863_image/threshold.png');

SM_multiple = saliency_gbvs.master_map_resized.*Mask;
SM_sum = sum(SM_multiple(:));
S_sum = sum(saliency_gbvs.master_map_resized(:));
RFA = SM_sum/S_sum;
end