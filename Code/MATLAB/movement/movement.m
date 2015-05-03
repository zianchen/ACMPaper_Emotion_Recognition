%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Matlab implementation of the following paper:                          %
%  "What Are We Looking For: Towards Statistical Modeling of Saccadic Eye %
%  Movements and Visual Saliency", pp. 1552-1559, CVPR(2012).             %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gaze_vector] = movement(img)
%% Load image and parameters
OriImg = im2double(img);
%OriImg = im2double(imread('C:\Users\chenzian\Documents\MATLAB\Image\happybaby.jpg'));
% '27.jpg'-- short rage attention shifts
% '28.jpg'-- long rage attention shifts
[Hei Wid Dim] = size(OriImg);
Bsize = 5;        % size of the patch
RStep = 1;        % sampling step length, 1 means dense sampling.
numSaccades = 10; % number of eye-movements
Rsize = floor(Bsize/2);
FilterDim = Bsize*Bsize*Dim;
bins = 255;       % bin number for non-parametric probability estimation.
DHei = 60;        % resized height
DWid = 80;        % resized width

Image = imresize(OriImg, [DHei DWid]);
% for grey image
if Dim == 1 
    Img(:,:,1) = OriImg;
    Img(:,:,2) = OriImg;
    Img(:,:,3) = OriImg;
    OriImg = Img;
end
SaliencyMap = zeros(DHei,DWid);
%% Patch based representation
PatchNum = 0;
for m= Rsize+1:RStep: DHei-Rsize
    for n= Rsize+1:RStep: DWid-Rsize
         PatchNum = PatchNum+1;
         patch=Image(m-Rsize: m+Rsize, n-Rsize: n+Rsize,:);  
         mixedsig(:,PatchNum)=reshape(patch,FilterDim,1);   
    end
end
%% Super Gaussian Component Analysis
[ppsig, A, W] = SGCA(mixedsig, numSaccades);

%% get distribution of gaze vector -Anson
W = sum(W,2)
gaze_vector = W'
[FeaDim,SampleNum] = size(ppsig);
PatchNum = 0;
%% Obtain all Response Maps
RM = zeros(DHei,DWid,FeaDim);
for k= Rsize+1:RStep: DHei-Rsize
    for l= Rsize+1:RStep: DWid-Rsize
         PatchNum = PatchNum+1;
         patch=Image(k-Rsize: k+Rsize, l-Rsize: l+Rsize,:);  
         RM(k,l,:) = abs(ppsig(:,PatchNum));  
    end
end
%% Display eye-movements and saliency maps
ScanImg = OriImg;
for m = 1:FeaDim
    %% Gaze Selection
    DisImg = OriImg; 
    RMC = abs(RM(:,:,m));
    RMC(1:Bsize,:) = 0;
    RMC(:,1:Bsize) = 0;
    RMC(DHei-Bsize:DHei,:) = 0;
    RMC(:, DWid-Bsize:DWid) = 0;
    RMC = imresize(RMC,[Hei Wid],'box');
    RMC = RMC/max(RMC(:));
    [Tx Ty] = find(RMC == max(RMC(:))); 
    if (m>1)
        ScanImg=bitmapplot([Px,mean(Tx)], [Py,mean(Ty)],ScanImg,struct('LineWidth',3,'Color',[1 0 0 1]));
    end
    Px =  floor(mean(Tx));
    Py =  floor(mean(Ty));
    Scale = floor((Hei/DHei)*Bsize);
    
%    DisImg(Px-Scale:Px+Scale,Py-Scale:Py+Scale,1)=1;
    %% Saliency Map
    ppsig(m,:)=ppsig(m,:)-min(min(ppsig(m,:)));
    ppsig(m,:)=ppsig(m,:)/(max(max(ppsig(m,:))));
    histo = hist(ppsig(m,:),bins);   
    ppsig(m,:)= -log(histo(round(ppsig(m,:).*(bins-1)+1))./sum(histo)+0.000001);
    ppsig((ppsig<0))=0;
    PatchNum = 0;
    for k= Rsize+1:RStep: DHei-Rsize
        for l= Rsize+1:RStep: DWid-Rsize
            PatchNum = PatchNum+1;
            SaliencyMap(k,l) = SaliencyMap(k,l)+ppsig(m,PatchNum);
        end
    end   
    SMap = filter2(fspecial('gaussian',8,3),SaliencyMap);
    SMap = SMap/max(SMap(:));
    SMap = imresize(SMap.^4,[Hei Wid],'bilinear');   
    
    subplot(2,2,1);
    imshow(DisImg);
    title('Input Image');
    subplot(2,2,2);
    imshow(RMC);
    title('Current Response Map');
    subplot(2,2,3);
    imshow(ScanImg);
    title('Scan Path');
    subplot(2,2,4);
    imshow(SMap);
    title('Current Saliency Map');
    
    pause(0.1);
end  
