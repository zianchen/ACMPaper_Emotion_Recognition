function get_symmetry(img,file)

out1 = symmetry(img,'mirror',0.8);
fprintf(file,'%f ',out1);
out2 = symmetry(img,'rotational');
fprintf(file,'%f ',out2);
out3 = FRST(double(img), 5);
fprintf(file,'%f ',out3);
end