function draw_keypoints_fast(im,keyVector,colour);
% Draw image with keypoints.
%
% From code by 
% D. Alvaro & J.J. Guerrero, Universidad de Zaragoza, Nov 2003.  
% Vectorised by Gareth Loy, KTH, 2005.
%
% e.g.
%   draw_keypoints(im,keyVector);

if nargin==2
    colour = [1 1 1];
end
    
[r, c] = size(im); 
% Draw a unitary horizontal arrow, transformed according to
% keypoint parameters.
TransformLine([r c], keyVector, 0.0, 0.0, 1.0, 0.0, colour);
TransformLine([r c], keyVector, 0.85, 0.1, 1.0, 0.0, colour);
TransformLine([r c], keyVector, 0.85, -0.1, 1.0, 0.0, colour);



% Draw the given line in the im, but first translate, rotate, and
% scale according to the keypoint parameters.
% Function adapted from 'lowe_keys\source\main.c'
%
% Parameters:
%   Arrays:
%    im = [rows columns] of im
%    keypoint = [subpixel_row subpixel_column scale orientation]
%
%   Scalars:
%    x1, y1; begining of vector
%    x2, y2; ending of vector
function TransformLine(im, keypoint, x1, y1, x2, y2, colour)

% The scaling of the unit length arrow is set to half the width
% of the region used to compute the keypoint descriptor.
len = 6 * keypoint(:,3);

% Rotate the keypoints by 'ori' = keypoint(4)
s = sin(keypoint(:,4));
c = cos(keypoint(:,4));

% Apply transform
r1 = keypoint(:,1) - len .* ( c*y1 + s*x1);
c1 = keypoint(:,2) + len .* (-s*y1 + c*x1);
r2 = keypoint(:,1) - len .* ( c*y2 + s*x2);
c2 = keypoint(:,2) + len .* (-s*y2 + c*x2);

% Draw line
line([c1'; c2'], [r1'; r2'], 'Color', colour,'LineWidth',1);
