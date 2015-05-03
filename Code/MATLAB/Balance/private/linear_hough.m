function [Hough_im,max_r,max_ang,r,strength] = linear_hough(im,x,y,ang,weights);
% use linear Hough Transform to vote for symmetry axis candidates
%
% INPUTS:
%   im = the original image
%   x, y  = the coordinates of the sym axis 'particle'
%   ang = the angle of line joining the matching pair, i.e. perpendicular
%         to the sym axis particle's orientation.
%   weights = the symmetry magnitude weights for each sym axis particle.

% determine the maximum radii possible
rad_upper_bound = round(sqrt((size(im,2)/2).^2 + (size(im,1)/2).^2));

% init Hough vote space
Hough_im = zeros(180, 2*rad_upper_bound);

% convert x and y to polar coordinates with origin in centre of image
ang_h = mod(ang,pi);
xx = x - size(im,2)/2;
yy = y - size(im,1)/2;
r = xx .* cos(ang_h) + yy .* sin(ang_h);

% convert angles and radii to points in Hough vote image
ang_hs = floor(180/pi*ang_h)+1;
r_hs = round(r + rad_upper_bound);

% cast votes in Hough vote image
ind_hs = sub2ind(size(Hough_im),ang_hs,r_hs);
%Hough_im = reshape(hist(ind_hs,[1:prod(size(Hough_im))]), size(Hough_im));
Hough_im = reshape(histc_weighted(ind_hs,weights,[1:prod(size(Hough_im))]), size(Hough_im));

th_rad = 5;
[max_r,max_ang,strength,Hough_im] = find_maxes(Hough_im, rad_upper_bound, th_rad);



function [max_r,max_ang,strength,Hough_im] = find_maxes(Hough_im, rad_upper_bound,th_rad)

% blur vote image
g_vec = make_gauss_vec(2*th_rad+1,th_rad/2);
Hough_im = [fliplr(Hough_im(end-th_rad+1:end, :)); Hough_im; fliplr(Hough_im(1:th_rad, :))]; % wrap angular axis
Hough_im = conv2(g_vec,g_vec,Hough_im,'same');
Hough_im = Hough_im(th_rad+1:end-th_rad,:); % unwrap angular axis

% find maxima
some_value = max(Hough_im(:))*0.5;
Hough_im_mod = Hough_im;
i = 1;
rr = 10;
while max(Hough_im_mod(:))>some_value
    strength(i) = max(Hough_im_mod(:));
    [max_ang_hs, max_r_hs] = find(Hough_im_mod == strength(i));
    max_ang(i) = max_ang_hs(1)/180*pi;
    max_r(i) = max_r_hs(1) - rad_upper_bound;
    
    if (max_ang_hs-rr)<1
        Hough_im_mod(end + 1 + max_ang_hs-rr : end, ...
            max(1,max_r_hs-rr) : min(size(Hough_im,2),max_r_hs+rr) ) = 0;
    end
    if (max_ang_hs+rr)>size(Hough_im,1)
        Hough_im_mod(1 : max_ang_hs+rr - size(Hough_im,1), ...
            max(1,max_r_hs-rr) : min(size(Hough_im,2),max_r_hs+rr) ) = 0;
    end 
    
    Hough_im_mod(max(1,max_ang_hs-rr) : min(size(Hough_im,1),max_ang_hs+rr), ...
        max(1,max_r_hs-rr) : min(size(Hough_im,2),max_r_hs+rr) ) = 0;
    i = i+1;
end

