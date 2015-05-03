function im = resize_image(im, x_mult, y_mult);
%im = resize_image(im, x_mult, y_mult);
%
%  Resizes a 2D image using bilinear interpolation
%
%  INPUTS:
%    x_mult specifies the stretch/shrink in the x direction
%    y_mult specifies the stretch/shrink in the y direction
%
%  This function is useful if you do not have imresize (from 
%  the Matlab Image Processing toolbox).

% Gareth Loy, KTH, 2004

new_width = round(x_mult * size(im,2));
new_height = round(y_mult * size(im,1));
Coords_x = repmat([1:size(im,2)], size(im,1), 1);             % - Construct mx of x pt coordinates
Coords_y = repmat([1:size(im,1)]', 1, size(im,2));            % - Construct mx of y pt coordinates

%new_x_coords = repmat([1: (size(im,2)-1) / (new_width-1) :new_width], new_height, 1) ;
new_x_coords = repmat([1: (size(im,2)-1) / (new_width-1) :size(im,2)], new_height, 1) ;
new_y_coords = repmat([1: (size(im,1)-1) / (new_height-1) :size(im,1)]', 1, new_width) ;

for i = 1:size(im,3)
    im_new(:,:,i) = interp2(Coords_x, Coords_y, im(:,:,i), new_x_coords, new_y_coords); 
end;
im = im_new;
