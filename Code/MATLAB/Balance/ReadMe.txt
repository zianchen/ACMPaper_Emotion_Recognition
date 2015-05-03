This folder contains Matlab code for running symmetry detection algorithm presented in

Loy and Eklundh, Detecting Symmetry and Symmetric Constellations of Features, ECCV, May, 2006.

Type 'help symmetry' at the Matlab command prompt to get help on this function.  Type 'make_results' at the Matlab command prompt to generate results from the images in the 'test images' directory.  Results are stored in the local 'results' folder.


Some additional points:

The current set-up runs under Windows XP.  To run on other operating systems it will be necessary to modify the call to the binary file that generates the SIFT features, see SIFTDemoV4/siftV4_mod.m lines 46-49.

The code expects the local directory structure provided. 

The code requires the ímshow.m function from the Matlab Image Processing Toolbox. 

David Lowe's siftDemoV4 has been included in this bundle, see licence and terms of use in that folder.

The test images provided were obtained from Flickr and are licenced under Creative Commons, see above ECCV publication for appropriate attributions if you wish to use these images in a publication.


