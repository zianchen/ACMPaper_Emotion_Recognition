function max_array = find_maxs(im, rad, thresh)
% Find maximums in image
% FIND_MAXS(IM, RAD, THRESH) returns an array of the x and y locations of
% maximum pts in the image IM.  The maximums are at least RAD a part and
% greater than THRESH in value.
%
% e.g.
%   max_array = find_maxs(im, rad, thresh)

% Gareth Loy, KTH, Stockholm, 2005.

% if there are no maxs above threshold, then exit
if max(im(:))<=thresh
    disp(sprintf('No maxs high enough, threshold = %g.  highest max = %g', ...
        thresh, max(im(:))));
    return
end

% build array of maxs
max_array = [];
max_val = max(im(:));
while max_val>thresh
    
    % find x and y coordinates of maximum pt, and add to array of maxs
    [y, x] = find(im==max_val);
    max_array = [max_array; x(1) y(1)];

    % set value of im to zero within rad of maximum pt
    im(y(1)-rad : y(1)+rad, x(1)-rad : x(1)+rad) = 0;

    % find new maximum value
    max_val = max(im(:));
end
