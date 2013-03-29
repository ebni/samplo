%   Optimal Sampling based on density
% INPUTS
%   N: number of sampling instants in [0, t_f]
%   t_f: final instant
%   numSteps: number of integration steps (should be much larger than N)
%   allT: should be linspace(0,t_f,numSteps+1);
%   uDot: should be diff(uOpt(allT))./diff(allT)
%   mA: sampling according to |uDot|^mA
% OUTPUTS
%   K0: cost matric so that the cost is J(N) = x_0'*K0*x_0
%   Uvec: Uvec(:,k) is the optimal k-th input
%   minCost: cost with this sampling method
%   tK: tK(i) is the i-th sampling instant (vector of size N+1)
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation

uDotM = abs(uDot).^mA;

% Integration of uDotM by hand
intDens = zeros(1,numSteps+1);   %init, length is numSteps+1
for i=2:numSteps+1
  % rectangle rule
  intDens(i) = intDens(i-1)+uDotM(i-1)*(allT(i)-allT(i-1));
end
intDens = intDens./intDens(numSteps);     % normalize

densToSamp = linspace(0,1,N+1);
% find the sampling instants by computing the inverse if intDens by hand
ind = 1;
tK = zeros(1,N+1);    % init
for i=2:N,
  while (intDens(ind+1) < densToSamp(i))
    ind = ind+1;
  end
  lll = (densToSamp(i)-intDens(ind))/(intDens(ind+1)-intDens(ind));
  tK(i) = allT(ind)+lll*(allT(ind+1)-allT(ind));
end
tK(N+1) = t_f;
tauK = diff(tK);

computeDiscDyn;                 % the discrete dynamics matrices
computeDiscWeights;             % the discrete-time weights matrices
computeDiscRiccati;             % the solution of the discrete Riccati Eqs
computeOptUk;                   % the system dynamics and optimal input
minCost = x_0'*K0*x_0;

clear uDotM intDens densToSamp i ind lll