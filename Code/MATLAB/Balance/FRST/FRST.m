function max_grad = FRST(im,radii,thresh,alpha,orient,A);
% FRST   compute Fast Radial Symmetry Transform (FRST).
%     Computes the Fast Radial Symmetry Transform (FRST) as described in 
%       
%       Loy & Zelinsky (2003), 
%       Fast Radial Symmetry for Detecting Points of Interest,
%       IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%       August 2003.
%
%     This code returns the combined absolute magnitudes of the dark and 
%     light symmetry components which is well suited for generic interest 
%     point detection (this varies slightly from the combined output presented
%     in the papers where dark symmetry is shown as negative and light as 
%     positive values).
%
%     To restrict this code to only detect dark or light comment out
%     the 3 lines suffixed with "% DISABLE FOR DARK/LIGHT".  With these
%     lines disabled dark symmetry is computed for positive radii
%     parameter values and light symmetry for negative radii values.  (The
%     resulting dark and light symmetries can be combined to generate the
%     combined output presented in the papers.)
%
%     Parameters:
%       im     = input image (gray scale)
%       radii  = vector of radii at which to compute transform
%       thresh = optional gradient threshold (beta in ECCV and PAMI papers), 
%                default 0.05.
%       A      = optional set of 1D separable component of 2D blurring kernel, default
%                A{i} = gaussian size 2*radii with sd radii/2.
%       alpha  = optional radial strictness parameter, default 2.
%       orient = optional flag (1 or 0) for computing 'orientation symmetry', i.e. only
%                using gradient orientations and not magnitude, see PAMI paper,
%                default 0.
%
%     Examples:
%       S = FRST(im, [2 4 7 10]);
%       S = FRST(double(imread('g.jpg')), 5);
%
%
% release 1.0.
% Gareth Loy, KTH, Stockholm, 2005.


if size(im,3)==3
   im = rgb2gray(im);
end
% set default values for undefined parameters
if nargin<=5,
    for i = 1:length(radii)
        A{i} = make_gauss_vec(round(2*radii(i)), 0.5*radii(i));
    end
end
if nargin<=4, orient = 0; end
if nargin<=3, alpha = 2; end
if nargin<=2, thresh = 0.05; end

% determien gradient magnitude and unit gradient
grad = gradient(im);
Mag = sqrt(grad.x.^2 + grad.y.^2);
unit_Grad.x = grad.x./(Mag + 1e-5);
unit_Grad.y = grad.y./(Mag + 1e-5);

% set edge mag to zero to avoid edge effects
Mag([1:2, end-2:end], :) = 0;
Mag(:, [1:2, end-2:end]) = 0;

% pad with zeros to accommodate voting outside the image
frame_r = max(abs(radii));
unit_Grad.x = pad(unit_Grad.x, frame_r, frame_r);
unit_Grad.y = pad(unit_Grad.y, frame_r, frame_r);
Mag = pad(Mag, frame_r, frame_r);

[Coords_x, Coords_y] = meshgrid(1:size(unit_Grad.x,2), 1:size(unit_Grad.x,1));

% only consider significant gradient elements
max_grad = max(Mag(:))
ind = find(Mag > thresh*255);
num = length(ind);

affected_pix_x = zeros(length(radii)*num,1);
affected_pix_y = zeros(length(radii)*num,1);
affected_pix_xn = zeros(length(radii)*num,1);
affected_pix_yn = zeros(length(radii)*num,1);
affecting_ind = repmat(ind,length(radii),1);

% if there are significant gradient elements present
if ~isempty(ind)
    S = 0;
    for i = 1:length(radii)

        % compute affected pixels and indices of affecting gradient
        % elements                
        for r = 1:length(radii)           
            affected_pix_x((r-1)*num+1:r*num) = ...
                round(Coords_x(ind) + unit_Grad.x(ind)*radii(r));
            affected_pix_y((r-1)*num+1:r*num) = ...
                round(Coords_y(ind) + unit_Grad.y(ind)*radii(r));
        end

        % TO ONLY COMPUTE DARK OR LIGHT SYMMETRY, COMMENT OUT FROM HERE
        % ------------------------
        for r = 1:length(radii)                             % DISABLE FOR DARK/LIGHT
            affected_pix_xn((r-1)*num+1:r*num) = ...
                round(Coords_x(ind) - unit_Grad.x(ind)*radii(r));% DISABLE FOR DARK/LIGHT
            affected_pix_yn((r-1)*num+1:r*num) = ...
                round(Coords_y(ind) - unit_Grad.y(ind)*radii(r));% DISABLE FOR DARK/LIGHT
        end                                                 % DISABLE FOR DARK/LIGHT
        affecting_ind = [affecting_ind; affecting_ind];     % DISABLE FOR DARK/LIGHT
        affected_pix_x = [affected_pix_x; affected_pix_xn]; % DISABLE FOR DARK/LIGHT 
        affected_pix_y = [affected_pix_y; affected_pix_yn]; % DISABLE FOR DARK/LIGHT
        % ----------------------------
        % TO ONLY COMPUTE DARK OR LIGHT SYMMETRY, COMMENT OUT ABOVE
        
        affected_ind = (affected_pix_x-1)*size(unit_Grad.x,1) + affected_pix_y;

        % compute orientation image,
        hist_1D = histc(affected_ind, [1:size(unit_Grad.x,1)*size(unit_Grad.x,2)]);
        Oa_1D = abs(hist_1D).^(alpha);
        Oa = reshape(Oa_1D, size(unit_Grad.x,1), size(unit_Grad.x,2));

        if orient~=1
            % compute magnitude image, and combine with orientation image
            hist_1D_M = histc_weighted(affected_ind, Mag(affecting_ind(:)), ...
                [1:size(unit_Grad.x,1)*size(unit_Grad.x,2)]);
            Ma = reshape(hist_1D_M, size(unit_Grad.x,1), size(unit_Grad.x,2));
            Oa = Oa.*Ma;
        end

        % remove frame of zeros
        Oa = unpad(Oa, frame_r, frame_r);

        % Use a separable kernal to blur result
        S = S + conv2(A{i}, A{i}', Oa, 'same');
    end
else
    % return an image of zeros
    warning('IND is empty, returning an image of zeros');
    S = zeros(size(unit_Grad.x)-2*frame_r);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE FUNCTIONS BELOW ARE GENERIC LIBRARY FUNCTIONS. THEY ARE INCLUDED IN 
% THIS FILE HERE AS THEY ARE REQUIRED BY THE FRST FUNCTION, HOWEVER, THEY 
% ARE USEFUL ON THEIR OWN AND AS SUCH CAN BE MOVED TO INDIVIDUAL FILES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function grad = gradient(I)
% grad = GRADIENT(I)
%
% determines the gradient of I returns the result as a vector field.
% 
% INPUT:  I = grayscale image
%
% OUTPUT: grad.x = x gradient image
%         grad.y = y gradient image
% 
% e.g.    grad = gradient(I);
% 

% Gareth Loy, KTH, Stockholm, 2003-2005.

if size(I,3) > 1
    error('gradient.m: unable to compute gradient of colour images.');
end

% g = makeGaussVec(5,1); 
% h = [-1; -2; 0; 2; 1];
% inv_step_edge_response = 0.4466; 

g = [1; 1; 1];
h = [-1; 0; 1];
inv_step_edge_response = 0.1768;

Ix = conv2(I,h','same');            
Iy = conv2(I,h,'same');            
Ix = conv2(Ix,g,'same');            
Iy = conv2(Iy,g','same');            
% Ix = conv2(g,g,Ix,'same');            
% Iy = conv2(g,g,Iy,'same');            
grad.x = Ix*inv_step_edge_response;
grad.y = Iy*inv_step_edge_response;


function mx_new = pad(mx,x,y,val)
%
% Pads a mx with a frame x wide and y high consisting of the value val.
% Default value for 'val' is 0.
%
% E.g. to pad a matrix mx with a border x_rad and y_rad of zeros:
%  pad(mx,x_rad,y_rad,0);
%  pad(mx,x_rad,y_rad);
%

% Gareth Loy, KTH, Stockholm, 2003-2005.

if nargin < 4
    val = 0;
end
mx_new = repmat(val,[2*[y,x]+[size(mx,1),size(mx,2)], size(mx,3)]);
mx_new(y+1:end-y,x+1:end-x,:) = mx;


function mx = unpad(mx,x,y)
%
% Unpads matrix mx by removing the outer x-wide and y-high frame, i.e.
% "cropping" the matrix.
%
% E.g. to pad a matrix mx with a border x_rad and y_rad of zeros:
%  unpad(max,x_rad,y_rad);
%

% Gareth Loy, KTH, Stockholm, 2003-2005.

mx = mx(y+1:end-y, x+1:end-x, :);


function [gauss_vec,d_gauss_vec] = make_gauss_vec(n,sd)
% gauss_vec = make_gauss_vec(n,sd)
%
% returns a normalised 1D Gaussian mask with
%  - n elements, and 
%  - sd standard deviation.
% 
% [gauss_vec,d_gauss_vec] = make_gauss_vec(n,sd)
%  also returns the derivative of the gaussian 
%

% Gareth Loy, KTH, Stockholm, 2003-2005.

for i=1:n
    dist = i-(n+1)/2;
    gauss_vec(i) = exp(-0.5*(dist/sd)^2)  / (sqrt(2*pi)*sd);
end; %for
gauss_vec = gauss_vec/(sum(sum(gauss_vec)));

if nargout==2
    for i=1:n
        dist = i-(n+1)/2;
        d_gauss_vec(i) = -2*dist/(2*sd)/(sqrt(2*pi)*sd) * exp(-0.5*(dist/sd)^2);
    end
    d_gauss_vec = d_gauss_vec/(sum(sum(gauss_vec)));
end
