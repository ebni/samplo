%   Computes the discrete dynamics from the continuous one
% INPUTS
%   n: dimension of the state
%   m: dimension of the input
%   N: number of sampling instants in [0, t_f]
%   Phi(t), Gamma(t): discretization over interval of lenght t
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation
% OUTPUTS
%   Avec, Bvec: discrete-time dynamics (see the code below)

Avec = zeros(n,n,N);
Bvec = zeros(n,m,N);
for k=1:N
  Avec(:,:,k) = Phi(tauK(k));
  Bvec(:,:,k) = Gamma(tauK(k));
end
clear k
