function [a,pixmax] = loadpgm(filename)

fid = fopen(filename,'r');
magic = fgetl(fid);
imsize = fscanf(fid,'%d',2);
pixmax = fscanf(fid,'%d',1);

a=zeros(imsize(2),imsize(1));

status=fseek(fid,1,'cof');

if magic=='P5'
  for(i=1:imsize(2))
    a(i,:) = transpose(fread(fid,imsize(1)));
  end
else
  a = transpose(fscanf(fid,'%d',[imsize(1) imsize(2)]));
end

fclose(fid);

