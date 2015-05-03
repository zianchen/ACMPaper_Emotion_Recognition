function [gauss_vec,d_gauss_vec] = make_gauss_vec(n,sd)
% gauss_vec = make_gauss_vec(n,sd)
%
% returns a normalised 1D Gaussian mask with
%  - n elements, and 
%  - sd standard deviation.
%
% 
% [gauss_vec,d_gauss_vec] = make_gauss_vec(n,sd)
% Also returns the derivative of the gaussian 
%
i=[1:n];
dist = i-(n+1)/2;
gauss_vec(i) = exp(-0.5*(dist/sd).^2)  / (sqrt(2*pi).*sd);

gauss_vec = gauss_vec/(sum(sum(gauss_vec)));

if nargout==2
    i=[1:n]
    dist = i-(n+1)/2;
    d_gauss_vec(i) = -2*dist/(2*sd)/(sqrt(2*pi).*sd) .* exp(-0.5*(dist/sd).^2);
    d_gauss_vec = d_gauss_vec/(sum(sum(gauss_vec)));
end