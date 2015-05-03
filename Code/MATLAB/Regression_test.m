addpath('C:\Users\sliu5\Desktop\first\Balance\siftDemoV4')
[t4,t5,t6,valence,arousal,len] = Getdata();
fid = fopen('C:\Users\sliu5\Desktop\first\valenceT.txt','w');
fid2 = fopen('C:\Users\sliu5\Desktop\first\arousalT.txt','w');
path = 'C:\Users\sliu5\Desktop\first\re\retest\';
for i=1:len
%modify file name    
fileExt = strcat(num2str(i),'.jpg');   
fileName = strcat(path,fileExt); 
sprintf('file:%c',fileName);
img = imread(fileName);

% write valence value
fprintf(fid,'%f ',valence(i));
fprintf(fid2,'%f ',arousal(i));
%generate balance vector
out1 = symmetry(img,'mirror',0.8);
fprintf(fid,'1:%f 2:%f 3:%f 4:%f',out1);
fprintf(fid2,'1:%f 2:%f 3:%f 4:%f',out1);
out2 = symmetry(img,'rotational');
fprintf(fid,'5:%f 6:%f',out2);
fprintf(fid2,'5:%f 6:%f',out2);
S = FRST(double(img), 5);
fprintf(fid,'7:%f ',S);
fprintf(fid2,'7:%f ',S);
%generate movement vector
[gaze_vector] = movement(img);
fprintf(fid,'8:%f 9:%f 10:%f 11:%f 12:%f 13:%f 14:%f 15:%f 16:%f 17:%f ',gaze_vector);
fprintf(fid2,'8:%f 9:%f 10:%f 11:%f 12:%f 13:%f 14:%f 15:%f 16:%f 17:%f ',gaze_vector);
%generate emphasis vector
[RFA] = emphasis(img);
fprintf(fid,'18:%f ',RFA);
fprintf(fid2,'18:%f ',RFA);
%generate gradation vector
[RG,AGTx,AGTy,AGIx,AGIy] = get_Gradation(fileName);
fprintf(fid,'19:%f 20:%f 21:%f 22:%f 23:%f ',RG,AGTx,AGTy,AGIx,AGIy);
fprintf(fid2,'19:%f 20:%f 21:%f 22:%f 23:%f ',RG,AGTx,AGTy,AGIx,AGIy);
%generate variety vector
[variety] = get_Variety(fileName);
fprintf(fid,'24:%f 25:%f 26:%f 27:%f 28:%f 29:%f 30:%f 31:%f 32:%f 33:%f 34:%f\n',variety);
fprintf(fid2,'24:%f 25:%f 26:%f 27:%f 28:%f 29:%f 30:%f 31:%f 32:%f 33:%f 34:%f\n',variety);
sprintf('-------------%d-----is----done!!!',i)
end;  
fclose(fid);
fclose(fid2);
