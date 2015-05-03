function [RG,AGTx,AGTy,AGIx,AGIy] = get_Gradation(img)

I = imread(img);
HSV = rgb2hsv(I);
V = HSV(:, :, 3);
[m,n] = size(V);
[Gx, Gy] = imgradientxy(V);
[Dx,Dy,Lx,Ly] = getMatrix(V,Gx,Gy);
RG=0;
AGTx = 0;
AGTy = 0;
AGIx = 0;
AGIy = 0;
for i=1:m
    for j=1:n
        RG = RG + Dx(i,j)/(Lx(i,j)+0.1) + Dy(i,j)/(Ly(i,j)+0.1);
        AGTx = AGTx + Dx(i,j);
        AGTy = AGTy + Dy(i,j);
        AGIx = AGIx + Lx(i,j);
        AGIy = AGIy + Ly(i,j);
    end
end

