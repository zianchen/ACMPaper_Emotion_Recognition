figure(1);
symmetry('test images\m a c s f  butterfly square.jpg','mirror',0.8);
print('-djpeg','results/mirror_result1');

% make results for rotational symmetry
figure(1);
symmetry('test images\smenzel_jensen.jpg','rotational results');
print('-djpeg','results/rot_result1');

figure(2);
symmetry('test images\cards.jpg','rotational results');
print('-djpeg','results/rot_result2');

figure(3);
symmetry('test images\olivander.jpg','rotational results',0.5);
print('-djpeg','results/rot_result3');

figure(4);
symmetry('test images\creativity_1.jpg','rotational results',0.5);
print('-djpeg','results/rot_result4');

figure(5);
symmetry('test images\gregw.jpg','rotational results',0.7);
print('-djpeg','results/rot_result5');
 
symmetry('test images\smenzel_cobra.jpg','rotational');


% make results for mirror symmetry
figure(8);
symmetry('test images\ze1 two butterflys.jpg','mirror results');
print('-djpeg','results/mirror_result8');

figure(5);
symmetry('test images\Leo_L30_lamp.jpg','mirror results');
print('-djpeg','results/mirror_result5');

figure(4);
symmetry('test images\Stuart_Maxwell_cheetah.jpg','mirror results',0.75);
print('-djpeg','results/mirror_result4');

figure(2);
symmetry('test images\smenzel.jpg','mirror results');
print('-djpeg','results/mirror_result2');

figure(3);
symmetry('F:\sym images non-commercial\BioID face database\images\BioID_1132.pgm','mirror results');
print('-djpeg','results/mirror_resultBioID');

figure(6);
symmetry('test images\david_martin_helicopter.jpg','mirror results',0.66);
print('-djpeg','results/mirror_result6');

figure(7);
symmetry('test images\elfintech.jpg','mirror results',0.95);
print('-djpeg','results/mirror_result7');







