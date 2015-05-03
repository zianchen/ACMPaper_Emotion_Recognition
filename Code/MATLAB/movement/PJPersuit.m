function [A, W] = PJPersuit(X, whiteningMatrix, dewhiteningMatrix, numSaccades)
[vectorSize, numSamples] = size(X);
%% initialize the main parameters
maxNumIterations = 500;  %Iteration Limit
epsilon = 0.0001;        %Parameter for testing iteration stop.
if numSaccades>vectorSize
    numOfPP =vectorSize; %Number of Projection Persuit.
else
    numOfPP = numSaccades; 
end
b_verbose = 0;           %Open Process indicator
failureLimit = 5;        %How many times do we try for convergence until we give up.

%% Check the data
if ~isreal(X)
  error('Input has an imaginary part.');
end
if b_verbose, fprintf('Starting Projection Pursuit...\n'); end
B = zeros(vectorSize);
round = 1; 
numFailures = 0;
while round<=numOfPP
    % Show the progress...
    if b_verbose, fprintf('Projection Vector %d ', round); end
    % Take a random initial vector of lenght 1 and orthogonalize it with respect to the other vectors.
    w = randn (vectorSize, 1);    
    w = w - B * B' * w;
    w = w / norm(w); 
    wOld = zeros(size(w));  
    % This is the fixed-point iteration loop.
    i = 1;
    while i <= maxNumIterations  
      % Project the vector into the space orthogonal to the space
      % spanned by the earlier found basis vectors. 
      w = w - B * B' * w;
      % Normalize the new w.
      w = w / norm(w);
      % If the iteration is not converged in N times, stop the process.
      if i == maxNumIterations
          round = round - 1;
          numFailures = numFailures + 1;
          %Failure Test
          if numFailures > failureLimit     
              if round == 0
                  A=[];
                  W=[];
              end
              return;
          end
          break;
      end
      % Show the progress...
      if b_verbose, fprintf('.'); end;   
      % Test for termination condition. The algorithm is converged if the direction of w and wOld is the same.
      if norm(w - wOld) < epsilon
          numFailures = 0;
          % Save the projection vector
          B(:, round) = w;
          % Calculate the de-whitened projection vector so it can be used for the orginal patch data.
          A(:,round) = dewhiteningMatrix * w;
          % Calculate corresponding filter.
          W(round,:) = w' * whiteningMatrix;
	      break;
      end      
      wOld = w;     
      w = (X * ((X' * w) .^ 3)) / numSamples;
      i = i + 1;
    end
    round = round + 1;
end 
if b_verbose, fprintf('Done.\n'); end 

%% Again, security test...
if ~isreal(A)
  A = real(A);
  W = real(W);
end

