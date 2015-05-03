function [ppsig, A, W] = SGCA(mixedsig, numSaccades)
if length (size (mixedsig)) > 2,
  error ('Input data can not have more than two dimensions.');
end
if any (any (isnan (mixedsig))),
  error ('Input data contains NaN''s.');
end
if ~isa (mixedsig, 'double')
  fprintf ('Warning: converting input data into regular (double) precision.\n');
  mixedsig = double (mixedsig);
end

% Remove the mean and check the data
[mixedsig, mixedmean] = remmean(mixedsig);

[Dim, NumOfSample] = size(mixedsig); 
% Calculating PCA
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';
verbose           = 'off';
[E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);
% Whitening
[whitesig, whiteningMatrix, dewhiteningMatrix] = whitenv ...
						     (mixedsig, E, D, verbose);
% Projection Persuit
[A, W] = PJPersuit(whitesig,  whiteningMatrix, dewhiteningMatrix, numSaccades); 
% Check for valid return
ppsig = W*mixedsig;
