%   Optimal Periodic Sampling
% INPUTS
%   n: dimension of the state
%   m: dimension of the input
%   t_f: final instant
%   N: number of sampling instants in [0, t_f]
%   Phi(t), Gamma(t): discretization over interval of lenght t
%   Q, R: weight to the state and to the input
%   S: weight to the final state x(t_f)
%   x_0: initial state
% OUTPUTS
%   K0: cost matric so that the cost is J(N) = x_0'*K0*x_0
%   Uvec: Uvec(:,k) is the optimal k-th input
%   minCost: cost with this sampling method
%   tK: tK(i) is the i-th sampling instant (vector of size N+1)
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation

% periodic samples
tK = linspace(0,t_f,N+1);
tauK = diff(tK);

% compute the matrixes of the dynamics and store them as columns
computeDiscDyn;

% compute the cost matrices for discrete-time
computeDiscWeights;

% computes optimal feedback
computeDiscRiccati;

 % compute the optimal input sequence and the corresponding state
computeOptUk;
  
minCost = x_0'*K0*x_0;

clear Xvec Avec Bvec QDvec RDvec PDvec Kvec
