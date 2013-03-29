%   Compute the average of the optima input over the sampling intervals
% INPUTS
%   uOpt(t): optimal input
%   tK: tK(i) is the i-th sampling instant (vector of size N+1)
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation
%   N: number of sampling instants in [0, t_f]
% OUTPUTS
%   Uvec: Uvec(:,k) is the average of uOpt over [tK(k), tK(k+1)]

Uvec = zeros(1,N);
for k=1:N
  if (tauK(k) > 0)
    Uvec(k) = quad(uOpt, tK(k), tK(k+1))/tauK(k);
  endif
  if (tauK(k) == 0)
    Uvec(k) = uOpt(tK(k));
  endif
  if (tauK(k) < 0)
    disp("computeAvgs: WARNING: integrating on a negative interval")
    Uvec(k) = uOpt(tK(k));
  endif
end
clear k