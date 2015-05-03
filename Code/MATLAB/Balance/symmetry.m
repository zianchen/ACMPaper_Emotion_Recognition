function output = symmetry(imageString, mode, scaling_factor)
%
% The mode parameter controls how the function behaves:
%   mode = 'mirror' ---> Detect bilateral symmetry, output step by step example figures
%   mode = 'mirror results ---> Detect bilateral symmetry, output initial and final figures
%   mode = 'rotational' ---> Detect rotational symmetry, output step by step example figures
%   mode = 'rotational results' ---> Detect rotational symmetry, output initial and final figures
%
% E.g.
% symmetry('test images\m a c s f  butterfly square.jpg','mirror',0.8);
% symmetry('test images\smenzel_cobra.jpg','rotational');

% Gareth Loy, KTH Sweden, 2005.

% ----------------------------------
% INITIALISE PARAMETERS
% ----------------------------------

% Set matlab path
adjust_path;

% Load image. If scaling parameter included, resize image accordingly
if nargin == 2
    scaling_factor = 1;
end
if min(size(imageString))==1
    im = load_image(imageString,scaling_factor);
else
    im = imageString;
end

% angular and radial tolerances for sym particles associated with the same symmetry axis
grouping_thresh = [3/180*pi, 3];

% percentage by which scale can vary within a matched pair
tol = 0.2;


% ----------------------------------
% COMPUTE KEYPOINTS
% ----------------------------------
[im, keyDescriptor, keyVector] = siftV4_mod(im);

% correct keypoint locatoins to be 1 pixel to the right and lower
keyVector(:,[1,2]) = keyVector(:,[1,2]) + 1;

% directly modify SIFT feature vector to represent reflected regions
keyDescriptor_m = mirror_SIFT_descriptor(keyDescriptor);
keyVector_m = keyVector;

% ----------------------------------
% COMPUTE symmetry
% ----------------------------------
switch lower(mode)
    case{'rotational','rotational results'}
        % ROTATIONAL symmetry
        num_matches_per_feature = 4;
        [match_ind,match_ind_m] = match(keyDescriptor, num_matches_per_feature, keyVector);
        [match_ind,match_ind_m] = reject_matches_based_on_scale(...
            tol, keyVector,keyVector_m,match_ind,match_ind_m);rot_centre = find_rotation_centres(keyVector,match_ind,match_ind_m);
        [max_rot_centre, Hough_im,sym_strength] = find_dominant_rot_centres(rot_centre,im);
    %    display_rot_output(im,keyVector,match_ind,match_ind_m,...
    %        rot_centre,Hough_im,max_rot_centre,sym_strength,mode);
        mm = find(sym_strength==max(max(sym_strength)));
        output = [max_rot_centre(mm),max(sym_strength)];
    case{'mirror','mirror results'}
        % REFLECTIVE symmetry
        num_matches_per_feature = 1;
        [match_ind,match_ind_m] = match(keyDescriptor, num_matches_per_feature, ...
            keyVector, keyDescriptor_m, keyVector_m);
        [match_ind,match_ind_m] = reject_matches_based_on_scale(...
            tol, keyVector,keyVector_m,match_ind,match_ind_m);
        [match_ind,match_ind_m,ang,phase_weight] = angular_constraint(keyVector, keyVector_m, ...
            match_ind, match_ind_m);
        % Compute optional distance weighting:
        % dist_weight = distance_weighting(keyVector, keyVector_m, ...
        %    match_ind, match_ind_m);

        % identify sym axes particles
        sym_x = mean([keyVector(match_ind,2), keyVector_m(match_ind_m,2)],2);
        sym_y = mean([keyVector(match_ind,1), keyVector_m(match_ind_m,1)],2);
        [Hough_im, max_r, max_ang, r, sym_strength] = linear_hough(im,sym_x,sym_y,ang,phase_weight);

        for i=1:length(max_r) %min(3,length(max_r))
            ind{i} = assign_to_axis(r,ang,max_r(i),max_ang(i),grouping_thresh);
            [left_ind{i},right_ind{i},pts{i},pts_m{i}] = ...
                group_left_right( keyVector,keyVector_m,match_ind, ...
                match_ind_m,r,ang,max_r,max_ang(i),grouping_thresh,ind{i});
        end

      %  display_output(im,keyVector,ang,sym_x,sym_y,ind,match_ind, match_ind_m,...
      %      Hough_im,keyVector_m,max_ang,max_r,left_ind,right_ind,pts,pts_m,phase_weight,sym_strength,mode);
        output = [length(ind{1}),max(max_r), max(max_ang), max(sym_strength)];
end

%----------------------------------------------------
% CODE SHARED BY BOTH ROTATIONAL AND MIRROR symmetry
%----------------------------------------------------

function adjust_path
% Set Matlab path for access to libraries and image files

% WINDOWS:
addpath(genpath('C:\Users\chenzian\Documents\MATLAB\'))


function im = load_image(imageString,scaling_factor)
% Load image, convert to greyscale and resize according to scaling factor.
%
% INPUTS:
%   imageString    - name of image file
%   scaling_factor - scaling factor, 1 => 100%
%
% OUTPUT:
%   im = image
%
if imageString(end)=='m' % i.e. is pgm image
    im = uint8(loadpgm(imageString));
else
    im = imread(imageString);
end

% convert to grayscale
if size(im,3)>1
    %    im = rgb2gray(im);  % requires image processing toolbox
    im = sum(im,3)/3;
end
im = uint8(resize_image(double(im), scaling_factor, scaling_factor));


function [match_ind,match_ind_m] = match(keyDescriptor, ...
    num_matches_per_feature, keyVector, keyDescriptor_m, keyVector_m)
% Match keypoint pairs.
%
% INPUTS:
%   keyDescriptor           - matrix of descriptor vectors (n x 128)
%   num_matches_per_feature - number of matches per feature
%   keyVector               - matrix of feature vectors (n x 4)
%   keyDescriptor_m         - matrix of mirrored descriptor vectors (n x 128)
%   keyVector_m             - matrix of mirrored feature vectors (n x 4)
%
% OUTPUTS:
%   match_ind, match_ind_m :-
%       indices of matched pairs, i.e. keyDescriptor(match_ind) matches
%       keyDescriptor_m(match_ind_m)
%
% Only need to compute upper triangular portion of distance matrix between
% the point sets since points cannot match with themselves and we don't
% want double matches. This assumes the kps in keyDescriptor and keyDescriptor_m are aligned,
% i.e., keyDescriptor_m(i,:) is the mirrored version of keyDescriptor(i,:)
%
if nargin < 5
    keyDescriptor_m = keyDescriptor;
    keyVector_m = keyVector;
end

% Build distance matrix describing distances between mirrored and
% unmirrored keypoints.
sz = size(keyDescriptor,1);
dist_mx = zeros(sz);
for k1 = 1:size(keyDescriptor,1)

    % compare a single keypoint with all other keypoints
    diff = ones(size(keyDescriptor_m,1)-k1,1) * keyDescriptor(k1,:) - keyDescriptor_m(k1+1:end,:);
    dist_mx(k1,k1+1:end) = sum(diff.^2,2)';
end

% make dist_mx symmetric
dist_mx = dist_mx + dist_mx';

dist_mx(dist_mx==0) = max(dist_mx(:));

match_ind = [];
match_ind_m = [];
for i=1:num_matches_per_feature
    match_ind = [match_ind, 1:size(keyDescriptor,1)];
    [temp, match_ind_m_new] = min(dist_mx,[],2);
    dist_mx(:,match_ind_m_new) = 1e4;
    match_ind_m = [match_ind_m; match_ind_m_new];
end


function [match_ind,match_ind_m] = reject_matches_based_on_scale(...
    tol, keyVector, keyVector_m, match_ind, match_ind_m);
% Reject matches at different scales
%
% INPUTS:
%   tol                     - percentage by which scale can vary within a matched pair
%   keyVector               - matrix of feature vectors (n x 4)
%   keyVector_m             - matrix of mirrored feature vectors (n x 4)
%   match_ind, match_ind_m :-
%       indices of matched pairs, i.e. keyDescriptor(match_ind) matches
%       keyDescriptor_m(match_ind_m)
%
% OUTPUTS:
%   match_ind, match_ind_m :-
%       indices of matched pairs, with matches at across scales removed.
%
scale = keyVector(:,3);
scale_m = keyVector_m(:,3);
scale_sim = (abs(scale(match_ind) - scale_m(match_ind_m))) ./ ...
    max([scale(match_ind), scale_m(match_ind_m)]')' < tol;
match_ind = match_ind(scale_sim);
match_ind_m = match_ind_m(scale_sim);


function dist_weight = distance_weighting(keyVector, keyVector_m, ...
    match_ind, match_ind_m)
% Compute optional distance weighting.  This function is not used at
% present.
%
% INPUTS:
%   keyVector               - matrix of feature vectors (n x 4)
%   keyVector_m             - matrix of mirrored feature vectors (n x 4)
%   match_ind, match_ind_m :-
%       indices of matched pairs, i.e. keyDescriptor(match_ind) matches
%       keyDescriptor_m(match_ind_m)
%
% OUTPUT:
%  dist_weight = distance weighting
%
dist_sq = (keyVector(match_ind,1) - keyVector_m(match_ind_m,1)).^2 + ...
    (keyVector(match_ind,2) - keyVector_m(match_ind_m,2)).^2;
sigma = sqrt(max(dist_sq))/6;
dist_weight = 1/(sigma*sqrt(2*pi)) * exp(-dist_sq/(2*sigma^2));


%----------------------------------
% CODE SPECIFIC TO MIRROR symmetry
%----------------------------------

function keyDescriptor_m = mirror_SIFT_descriptor(keyDescriptor)
% Use lookup table to generate the SIFT descriptor for the mirrored
% version of the image patch.  This will only work for SIFT, and requires
% the precise ordering returned by Lowes SIFT code (version 4).
%
% INPUT:  keyDescriptor    - 128 element SIFT descriptor vector of some
%                            image patch P.
% OUTPUT: keyDescriptor_m  - 128 element SIFT descriptor vector of some
%                            image patch P reflected about the y-axis.

ind(1:8:25) = [97:8:121];
ind(33:8:33+24) = [65:8:65+24];
ind(65:8:65+24) = [33:8:33+24];
ind(97:8:121) = [1:8:25];

ind(2:2+6) = [98+6:-1:98];
ind(10:10+6) = [106+6:-1:106];
ind(18:18+6) = [114+6:-1:114];
ind(26:26+6) = [122+6:-1:122];

ind(34:34+6) = [66+6:-1:66];
ind(42:42+6) = [74+6:-1:74];
ind(50:50+6) = [82+6:-1:82];
ind(58:58+6) = [90+6:-1:90];

ind(66:66+6) = [34+6:-1:34];
ind(74:74+6) = [42+6:-1:42];
ind(82:82+6) = [50+6:-1:50];
ind(90:90+6) = [58+6:-1:58];

ind(98:98+6) = [2+6:-1:2];
ind(106:106+6) = [10+6:-1:10];
ind(114:114+6) = [18+6:-1:18];
ind(122:122+6) = [26+6:-1:26];

keyDescriptor_m = keyDescriptor(:,ind);


function [match_ind,match_ind_m,ang,phase_weight] = angular_constraint(...
    keyVector, keyVector_m, match_ind, match_ind_m);
% Reject matches that are not symmetrically oriented
%
% INPUTS:
%   keyVector               - matrix of feature vectors (n x 4)
%   keyVector_m             - matrix of mirrored feature vectors (n x 4)
%   match_ind, match_ind_m :- indices of matched pairs, i.e.
%                             keyDescriptor(match_ind) matches
%                             keyDescriptor_m(match_ind_m)
%
% OUTPUTS:
%   match_ind, match_ind_m :- indices of matched pairs, after removing
%                             non-symmetric matches
%   ang                     - ang of lines joining feature pairs
%   phase_weight            - phase weight contribution to symmetry
%                             magnitude (n x 1) vector.

% determine orientation of line joining the pair
ang = atan2( (keyVector(match_ind,1) - keyVector_m(match_ind_m,1)), ...
    (keyVector(match_ind,2) - keyVector_m(match_ind_m,2)) );

% flip angle of mirrored points about y-axis.
% necessary because not a right hand coordinate frame.
r = ones(size(keyVector_m(:,4)));
[uxm,uym] = pol2cart(keyVector_m(:,4), r);
uxm = -uxm;
[keyVector(:,4),r] = cart2pol(uxm,uym);
[keyVector_m(:,4),r] = cart2pol(uxm,uym);

% compute angular component of first part of Reisfeld's phase weight function [Reisfeld94]
ang_i = keyVector(match_ind,4);
ang_j = keyVector_m(match_ind_m,4);
phase = ang_i + ang_j - 2*ang;
phase_weight = - cos(phase);
phase_weight(phase_weight < 0) = 0;

% discard pairs that are assymetrically oriented
sym_orient = phase_weight > 0;
match_ind = match_ind(sym_orient);
match_ind_m = match_ind_m(sym_orient);
ang = ang(sym_orient);
phase_weight = phase_weight(sym_orient);


function ind = assign_to_axis(r,ang,max_r,max_ang,grouping_thresh)
% Assign pairs to the dominant axis of symmetry. Tag all pairs that agree
% with the dominant axis, and determine extent of axis in image.
%
% INPUTS:
%   r, ang          - Hough space polar coordinates of axes of symmetry
%   max_r, max_ang  - Hough space polar coordinates of dominant axis of symmetry
%   grouping_thresh - angular and radial tolerances for sym particles
%                     associated with the same symmetry axis.
%
% OUTPUT:
%   ind - indices of pairs assigned to dominant axis of symmetry

% compute angular (d_ang1) radial (d_r1) difference between each pair and
% the dominant axis.
temp = unwrap( [ang, repmat(max_ang,length(ang),1)], [], 2);
d_ang1 = abs(temp(:,1) - temp(:,2));
d_r1 = abs(r - max_r);

% Do the same again, but with the axis oriented 180 degrees in the
% opposite directioncompute angular (d_ang2) radial (d_r2) difference
% between each pair and the dominant axis.
temp = unwrap( [ang+pi, repmat(max_ang,length(ang),1)], [], 2);
d_ang2 = abs(temp(:,1) - temp(:,2));
d_r2 = abs(-r - max_r);

% choose the result that is closest to the dominant axis
[d_ang iind] = min([d_ang1, d_ang2],[],2);
d_r_comb = [d_r1, d_r2];
d_r(iind==1) = d_r1(iind==1);
d_r(iind==2) = d_r2(iind==2);

% Assign pairs that are sufficiently close to the dominant axis to this
% axis.
ind = find((d_ang<grouping_thresh(1)) & (d_r'<grouping_thresh(2)));


function display_output(im,keyVector,ang,sym_x,sym_y,... % first 5 parameters
    ind,match_ind,match_ind_m,Hough_im,keyVector_m,...
    max_ang,max_r,left_ind,right_ind,pts,...
    pts_m,phase_weight,sym_strength,mode);
% Disply output for bilateral symmetry detection

color = hsv(256);
% color_cycle_rate = 50;

% DISPLAY OUTPUT
im_dark = uint8(round(double(im)*0.55));

switch lower(mode)
    case{'mirror'}
        for i = 1:length(ind)
            col = [1 1 1];  % col = color(floor((i-1)/length(ind)*size(color,1))+1, :);
            for j = 1:length(ind{i})
            end
        end
        disp(sprintf('%g reflective matches',length(sym_x)));
        disp(sprintf('%g symmetric features',length(ind{1})));
 
    case{'mirror results'}
        clf;
        display_output_6(im_dark,max_r,max_ang,ind,sym_x,sym_y,left_ind,right_ind,pts,pts_m,sym_strength);
end
%set(gcf,'Color',[1 1 1])
%format_figure(im);

% % draw scale of features
% sc = 2;
% for i = 1:length(ind)
%     h = plot(keyVector(match_ind(ind(i)),2), keyVector(match_ind(ind(i)),1),...
%         'ro','MarkerSize',sc*keyVector(match_ind(ind(i)),3));
%     h = plot(keyVector_m(match_ind_m(ind(i)),2), keyVector_m(match_ind_m(ind(i)),1),...
%         'ro','MarkerSize',sc*keyVector_m(match_ind_m(ind(i)),3));
% end;


function display_output_6(im_dark,max_r,max_ang,ind,sym_x,sym_y,left_ind,right_ind,pts,pts_m,sym_strength);
% display matched pairs of mirrored keypoints associated with dominant axis

%figure(6);
imshow(im_dark);
hold on;
color = colormap(hsv); colormap(gray);
centre = size(im_dark)/2;
tol = 5; %10;
if mean(size(im_dark)) < 400
    linewidth = mean(size(im_dark))*0.005; %*0.01; %20;%10;
    markersize = mean(size(im_dark))*0.035; %*0.075; %20;%10;
else
    linewidth = mean(size(im_dark))*0.0075; %*0.015; %20;%10;
    markersize = mean(size(im_dark))*0.05; %*0.1; %20;%10;
end
for i = 1:length(max_r)

    if isempty(ind{i}), continue, end

    %col = color(floor((i-1)/length(ind)*size(color,1))+1, :);
    shade = sym_strength(i)/max(sym_strength); %max([sym_strength(i)/max(sym_strength), min([sym_strength(i)/0.3, 1])]);
    col = [1 1 1] * shade;
    % c = mod(color_cycle_rate*i,size(color,1))+1;

    % construct set of points describing current symmetry axis and draw
    % axis.
    if abs(diff(unwrap([max_ang(i),pi]))) < (pi/4)
        Y = [min(sym_y(ind{i})) - tol : ...
            (max(sym_y(ind{i})) - min(sym_y(ind{i})))/100 : ...
            max(sym_y(ind{i})) + tol];
        X = ((-Y + centre(1))*sin(max_ang(i)) + max_r(i)) / cos(max_ang(i)) + centre(2);
    else
        X = [min(sym_x(ind{i})) - tol : ...
            (max(sym_x(ind{i})) - min(sym_x(ind{i})))/100 : ...
            max(sym_x(ind{i})) + tol];
        Y = ((-X + centre(2))*cos(max_ang(i)) + max_r(i)) / sin(max_ang(i)) + centre(1);
    end

    p_ind = find((X>=(min(sym_x(ind{i}))-tol)) & (X<=(max(sym_x(ind{i}))+tol)) & ...
        (Y>=(min(sym_y(ind{i}))-tol)) & (Y<=(max(sym_y(ind{i}))+tol)));
    plot(X(p_ind),Y(p_ind),'y-','LineWidth', linewidth,'Color',col);

    % draw left and right point sets (currently uses the same colour for left and right sets)
    plot(pts{i}(right_ind{i},1), pts{i}(right_ind{i},2), 'r.','Color',col,'MarkerSize',markersize);
    plot(pts_m{i}(left_ind{i},1), pts_m{i}(left_ind{i},2), 'r.','Color',col,'MarkerSize',markersize);
    plot(pts{i}(left_ind{i},1), pts{i}(left_ind{i},2), 'g.','Color', col,'MarkerSize',markersize);
    plot(pts_m{i}(right_ind{i},1), pts_m{i}(right_ind{i},2), 'g.','Color', col,'MarkerSize',markersize);
end
%disp([gcf, sym_strength]');
%disp(sprintf('Figure %g: ', gcf));
disp(sprintf('Symmetry Magnitude = %g  ', sym_strength));
hold off
drawnow


function [left_ind,right_ind,pts,pts_m] = group_left_right(...
    keyVector,keyVector_m,match_ind,match_ind_m,r,ang,max_r,max_ang,...
    grouping_thresh,ind,i)
% Group left and right points of symmetric feature pairs into two distinct
% sets.
%
% INPUTS:
%   keyVector,keyVector_m,match_ind,match_ind_m,r,ang,max_r,max_ang,
%   grouping_thresh
%
% OUTPUTS:
%   ind,left_ind,right_ind,pts,pts_m
%

pts = [keyVector(match_ind(ind),2), keyVector(match_ind(ind),1)];
pts_m = [keyVector_m(match_ind_m(ind),2), keyVector_m(match_ind_m(ind),1)];

% divide matched points into two sets about symmetry axis
norm_vec = [cos(max_ang) sin(max_ang)];
left_ind = find(((pts - pts_m) * norm_vec')>0);
right_ind = find(((pts - pts_m) * norm_vec')<0);


%--------------------------------------
% CODE SPECIFIC TO ROTATIONAL symmetry
%--------------------------------------

function [max_rot_centre, Hough_im, sym_strength] = ...
    find_dominant_rot_centres(rot_centre,im)
% Find dominant centres of rotational symetry.
%
% INPUTS:
%   rot_centre - array of all rotation centres
%   im         - input image (only required to give image dimensions)
%
% OUTPUTS:
%   max_rot_centre - dominant rotation centre
%   Hough_im       - Hough image, result of Hough transform in x,y space
%   sym_strength   - symmetry magnitude of dominant rotation centre.

max_rot_centre = [];

% discard centres of rotation that lie outside the image
r1 = find(rot_centre(:,1)>=1 & rot_centre(:,1)<size(im,2) & ...
    rot_centre(:,2)>=1 & rot_centre(:,2)<size(im,1));
rot_centre = floor(rot_centre(r1,:))+1;

% cast votes in vote image
ind_hs = sub2ind(size(im),rot_centre(:,2),rot_centre(:,1));
Hough_im = reshape(hist(ind_hs,(1:numel(im))), size(im));

% blur vote image
th_rad = 5;
g_vec = make_gauss_vec(2*th_rad+1,th_rad/2);
Hough_im = conv2(g_vec,g_vec,Hough_im,'same');

% find maxima
some_value = max([4*max(g_vec(:))^2, max(Hough_im(:))*0.5]);
Hough_im_mod = Hough_im;
i = 1;
rr = 10;
if max(Hough_im_mod(:))<=some_value
    disp('No centres of rotation found');
    max_rot_centre = [0,0];
    sym_strength = 0;
    return;
end
while max(Hough_im_mod(:))>some_value
    sym_strength(i) = max(Hough_im_mod(:));
    [i2, i1] = find(Hough_im_mod == sym_strength(i));
    max_rot_centre(i,2) = i2(1);
    max_rot_centre(i,1) = i1(1);
    Hough_im_mod( max(1,max_rot_centre(i,2)-rr) : ...
        min(size(Hough_im,1),max_rot_centre(i,2)+rr ), ...
        max(1,max_rot_centre(i,1)-rr) : ...
        min(size(Hough_im,2),max_rot_centre(i,1)+rr ) ) = 0;
    i = i+1;
end


function rot_centre = find_rotation_centres(keyVector,match_ind,match_ind_m)
% Find centres of rotation.
%
% INPUTS:
%   keyVector      - matrix of feature vectors (n x 4)
%   keyVector_m    - matrix of mirrored feature vectors (n x 4)
%   match_ind and  - indices of matched pairs i.e. keyDescriptor(match_ind)
%     match_ind_m    matches keyDescriptor_m(match_ind_m)
%
% OUTPUT:
%   rot_centre  - matrix of x,y coordinates of centres of rotation (n x 2)
%
% Refer to ECCV paper for diagram and description of formulas

% beta is the fixed angle each pt vector makes with the radial
% line joining it to the centre of rotation
beta = (mod(keyVector(match_ind_m,4) - pi - keyVector(match_ind,4) + pi, 2*pi)  - pi)/2;

% gamma is the angle of the line joining the two pts
gamma = atan2( (keyVector(match_ind_m,1) - keyVector(match_ind,1)), ...
    (keyVector(match_ind_m,2) - keyVector(match_ind,2)) );

% d is the distance between the two pts
d = sqrt( (keyVector(match_ind,2) - keyVector(match_ind_m,2)).^2 + ...
    (keyVector(match_ind,1) - keyVector(match_ind_m,1)).^2 );

% r is the length of the radial line
r = 0.5*d .* sqrt(1 + (tan(beta)).^2);

% finally determine the centre of rotation:
rot_centre = [keyVector(match_ind,2) + (r .* cos(beta+gamma)), ...
    keyVector(match_ind,1) + (r .* sin(beta+gamma))];


function display_rot_output(im,keyVector,match_ind,match_ind_m, ...
    rot_centre, Hough_im, max_rot_centre, sym_strength, mode)
% Display output for rotational symmetry

im_dark = uint8(round(double(im)*0.55));

switch lower(mode)
    case{'rotational'}
        figure(2);
        imshow(im_dark);
        hold on;
        %plot(keyVector(:,2), keyVector(:,1), 'w.');
        draw_keypoints_fast(im,keyVector);
        hold off;
        format_figure(im);
        print('-djpeg','results/rot_eg2');
        disp(sprintf('%g keypoints detected',size(keyVector,1)));

        % draw centres of rotation
        figure(3);
        imshow(im_dark);
        axis image;
        hold on,
        plot(rot_centre(:,1), rot_centre(:,2), 'wx', 'MarkerSize', 10, 'LineWidth', 2);
        hold off;
        format_figure(im);
        print('-djpeg','results/rot_eg3');

        figure(4);
end

imshow(im_dark);
hold on;
%plot(keyVector(:,2), keyVector(:,1), 'r.');
grouping_sq_dist = 25;
pointsize = 1;%5;
centresize = 7; %20
if mean(size(im_dark,1)) == 234
    linewidth = mean(size(im_dark))*0.005; %*0.01; %20;%10;
    markersize = mean(size(im_dark))*0.035; %*0.075; %20;%10;
else
    linewidth = mean(size(im_dark))*0.0075; %*0.015; %20;%10;
    markersize = mean(size(im_dark))*0.05; %*0.1; %20;%10;
%     linewidth = mean(size(im_dark))*0.0075; %20;%10;
%     markersize = mean(size(im_dark))*0.075; %20;%10;
% else
%     linewidth = mean(size(im_dark))*0.010; %20;%10;
%     markersize = mean(size(im_dark))*0.1; %20;%10;
end

for i=1:min(3,size(max_rot_centre,1))

    % tag all sym-axis-particles that agree with the current centre of
    % rotation

    dist_sq = (max_rot_centre(i,1) - rot_centre(:,1)).^2 + ...
        (max_rot_centre(i,2) - rot_centre(:,2)).^2;
    ind{i} = find(dist_sq<grouping_sq_dist);
    pts{i} = [keyVector(match_ind(ind{i}),2), keyVector(match_ind(ind{i}),1)];
    pts_m{i} = [keyVector(match_ind_m(ind{i}),2), keyVector(match_ind_m(ind{i}),1)];

    % plot(keyVector(match_ind(ind{i}),2), keyVector(match_ind(ind{i}),1), 'g*');
    % plot(keyVector(match_ind_m(ind{i}),2), keyVector(match_ind_m(ind{i}),1), 'g*');

    %     plot([keyVector(match_ind(ind{i}),2), keyVector(match_ind_m(ind{i}),2)]', ...
    %         [keyVector(match_ind(ind{i}),1), keyVector(match_ind_m(ind{i}),1)]', 'g-');

    c_map = colormap('hsv');
    for j = 1:size(ind{i})
        colour = c_map(mod(7*j,256),:);
        colour = [1 1 1];
        k = match_ind(ind{i}(j));
        k_m = match_ind_m(ind{i}(j));
        pt1 = [keyVector(k,2), keyVector(k,1)];
        pt2 = [keyVector(k_m,2), keyVector(k_m,1)];
        d_theta = pi/200;
        [h,th] = draw_arc(d_theta,rot_centre(ind{i}(j),:),pt1,pt2);
        set(h,'Color',colour,'LineWidth',linewidth);
        ang_of_rotation(j) = abs(th(2)-th(1));
        % plot centres of rotation of individual pairs
        %plot(rot_centre(ind{i}(j),1), rot_centre(ind{i}(j),2), 'gx','Color',colour);

        plot(keyVector(k,2), keyVector(k,1), 'g.','Color',colour,'LineWidth',pointsize,'MarkerSize',markersize);
        plot(keyVector(k_m,2), keyVector(k_m,1), 'g.','Color',colour,'LineWidth',pointsize,'MarkerSize',markersize);
        %         plot(keyVector(k,2), keyVector(k,1), 'go','Color',colour,'LineWidth',pointsize);
        %         plot(keyVector(k_m,2), keyVector(k_m,1), 'go','Color',colour,'LineWidth',pointsize);
    end
    h = plot(max_rot_centre(i,1), max_rot_centre(i,2), 'wx');
    %set(h,'LineWidth',centresize);
    set(h,'MarkerSize',10*centresize,'LineWidth',2*linewidth);
    disp(sprintf('%g rotationally symmetric matches, symmetry strength = %g',length(ind{i}), sym_strength(i)));

    fig = gcf;
    h = figure(10+i);
    set(h,'Color',[1 1 1]);
    histogram = hist(ang_of_rotation*180/pi,[1:180]); axis([0 180 0 10]);
    order = quantify_rot_order(histogram,mode);
    switch lower(mode)
        case{'rotational'}
            print('-dpng','results/rot_eg6');
    end
    figure(fig);
end
colormap(gray);
hold off;
set(gcf,'Color',[1 1 1])
drawnow;
format_figure(im);
switch lower(mode)
    case{'rotational'}
        print('-djpeg','results/rot_eg4');
end


function format_figure(im);
% Format single image figure for printing for inclusion in paper
%
axis image;  % makes spacing
axis([0 size(im,2) 0 size(im,1)]);
set(gcf,'PaperSize',[size(im,2), size(im,1)]);
set(gcf,'Position',[0,0,size(im,2),size(im,1)]);
set(gcf,'PaperPosition',[1,1,size(im,2)/50,size(im,1)/50]);
set(gca,'Position',[0,0,1,1]);
box off; axis off;
set(gcf,'InvertHardcopy','off');


function order = quantify_rot_order(histogram,mode);
% Quantify the order of rotational symmetry
%
for order = 2:12
    g = make_gauss_vec(40,10);
    g = g/max(g);
    g = ones(1,20);
    rot = 360/order;
    pts = [rot:rot:180];
    x = zeros(1,180);
    x(round(pts)) = 1;
    x = conv2(x,g,'same');
    offset_pts = [rot/2:rot:180];
    x_offset = zeros(1,180);
    x_offset(round(offset_pts)) = 1;
    x_offset = conv2(x_offset,g,'same');

    switch lower(mode)
        case{'rotational'}
            if order == 10
                %plot(x);
                h = bar(x);
                set(h,'FaceColor',0.8*[1 1 1],'EdgeColor',0.8*[1 1 1]);
                hold on;
                h = bar(-x_offset);
                set(h,'FaceColor',0.8*[1 1 1],'EdgeColor',0.8*[1 1 1]);
                h = bar(histogram);
                set(h,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
                set(gca,'FontSize',16);
                set(gcf,'PaperPosition',[1,1,400/50,200/50]);
                axis([0 180 -1.5 2.5]);
                print('-dpng',sprintf('results/rot_eg5'));
                hold off;
            end
    end
    response(order) = (sum(x.*histogram) - sum(x_offset.*histogram))/order;
end
response(response<0) = 0;
response = response/max(response(:));
%subplot(6,2,1);
plot([2:length(response)], response(2:end),'k-','Linewidth',3);
set(gca,'FontSize',16);
set(gcf,'PaperPosition',[1,1,400/50,200/50]);
axis tight
[temp,order] = max(response);


function [h,th] = draw_arc(d_theta,rot_centre,pt1,pt2);
% Draws and arc from pt1 to pt1 with centre at rot_centre
% Returns:
%  h  = handle to the plotted arc
%  th = the angles determined specifying the pts positions on the arc
%       the angle of rotation = abs(th(2)-th(1))

% find radius
r = sqrt(sum((rot_centre - pt1).^2));

% find theta at pt1 and pt2
th(1) = atan2(pt1(2) - rot_centre(2), pt1(1) - rot_centre(1));
th(2) = atan2(pt2(2) - rot_centre(2), pt2(1) - rot_centre(1));
th = unwrap(th);

% make list of thetas
th_list = [min(th):d_theta:max(th), max(th)];

% compute list of x,y points
x_list = rot_centre(1) + r * cos(th_list);
y_list = rot_centre(2) + r * sin(th_list);

% draw lines between points
h =line(x_list,y_list);


