%   Optimal Fast Numerical Quantization-based Sampling
% INPUTS
%   uOpt(t): optimal input to be quantized
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation
%      (initial guess)
%   N: number of sampling instants in [0, t_f]
%   t_f: final instant
%   numSteps: number of integration steps (should be much larger than N)
%   allT: should be linspace(0,t_f,numSteps+1);
% INPUTS
%     A B n m Phi Gamma x_0 t_f Q R S Kinf L N tauK
% OUTPUTS
%   K0: cost matric so that the cost is J(N) = x_0'*K0*x_0
%   Uvec: Uvec(:,k) is the optimal k-th input
%   minCost: cost with this sampling method
%   tK: tK(i) is the i-th sampling instant (vector of size N+1)
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation


% the optimal input WARNING: this is true only when S=Kinf
for i=1:numSteps+1,
  allF(i) = uOpt(allT(i));
endfor

tKfromTauK;

iter=0;
while (iter<=10)     % STOP CONDITION (1): too many iterations

  %disp("------------------------------------------------------");
  %disp(iter);
  tauK = diff(tK);

  % computing Uvec as the average of the opt continuous time
  computeAvgs;
  for k=1:N-1,
    % index in allT such that allT(curK)=tK(k+1)
    curK = floor(tK(k+1)/t_f*numSteps+1);
    if norm(Uvec(k)-allF(curK)) < norm(Uvec(k+1)-allF(curK))
      % should increase curK
      while norm(Uvec(k)-allF(curK)) < norm(Uvec(k+1)-allF(curK))
	curK = curK+1;
	if curK>numSteps+1
	  disp("there are problems with the algorithm");
	endif
      endwhile
    else
      % should decrease curK
      while norm(Uvec(k)-allF(curK)) >= norm(Uvec(k+1)-allF(curK))
	curK = curK-1;
	if curK <= 0
	  disp("there are problems with the algorithm");
	endif
      endwhile
      curK = curK+1;
    endif
    % computing the intercept between the curve of the optimal u and 
    % plane of equidistance between Uvec(k) and Uvec(k+1)

    %scalarMid = .5*(dot(Uvec(k),Uvec(k))+dot(Uvec(k+1),Uvec(k+1)));
    %scalar0 = dot(Uvec(k+1)-Uvec(k),allF(curK-1));
    %scalar1 = dot(Uvec(k+1)-Uvec(k),allF(curK));
    %lambda = (scalarMid-scalar0)/(scalar1-scalar0);
    %nextTk = allT(curK-1)+lambda*(allT(curK)-allT(curK-1));
    %tK(k+1) = nextTk;

    tK(k+1) = allT(curK);
  endfor
  iter = iter+1;
endwhile

% the quantized input
%[Uvec] = computeAvgs(A,B,L,x_0,tauK);

computeDiscDyn;                 % the discrete dynamics matrices
computeDiscWeights;             % the discrete-time weights matrices
computeDiscRiccati;             % the solution of the discrete Riccati Eqs
computeOptUk;                   % the system dynamics and optimal input
minCost = x_0'*K0*x_0;

clear allF