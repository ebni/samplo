%   Computes the sampling instants from the intersample separations
% INPUTS
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation
%      (initial guess)
%   N: number of sampling instants in [0, t_f]
% OUTPUTS
%   tK: tK(i) is the i-th sampling instant (vector of size N+1)


tK = zeros(1,N+1);       % init
for j=1:N
  tK(j+1) = tK(j)+tauK(j);
end
