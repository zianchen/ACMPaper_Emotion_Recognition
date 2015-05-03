addpath('C:\Users\chenzian\Documents\MATLAB\Balance\siftDemoV4')
%traversal all images
path = 'C:\Users\chenzian\Documents\MATLAB\Image\emotion_test\sadness\';  
fileExt = '*.jpg';  
files = dir(fullfile(path,fileExt));  
len = size(files,1);  
fid = fopen('C:\Users\chenzian\Documents\MATLAB\result.txt','w');
for i=1:len  
fileName = strcat(path,files(i,1).name), 
img = imread(fileName);

%define the label
fprintf(fid,'7 ');
%generate balance vector
out1 = symmetry(img,'mirror',0.8);
fprintf(fid,'1:%f 2:%f 3:%f 4:%f ',out1);
out2 = symmetry(img,'rotational');
fprintf(fid,'5:%f 6:%f ',out2);
S = FRST(double(img), 5);
fprintf(fid,'7:%f ',S);
%generate movement vector
[gaze_vector] = movement(img);
fprintf(fid,'8:%f 9:%f 10:%f 11:%f 12:%f 13:%f 14:%f 15:%f 16:%f 17:%f ',gaze_vector);
%generate emphasis vector
[RFA] = emphasis(img);
fprintf(fid,'18:%f ',RFA);
%generate gradation vector
[RG,AGTx,AGTy,AGIx,AGIy] = get_Gradation(fileName);
fprintf(fid,'19:%f 20:%f 21:%f 22:%f 23:%f ',RG,AGTx,AGTy,AGIx,AGIy);
%generate variety vector
[variety] = get_Variety(fileName);
fprintf(fid,'24:%f 25:%f 26:%f 27:%f 28:%f 29:%f 30:%f 31:%f 32:%f 33:%f 34:%f\n',variety);
end;  
fclose(fid);













