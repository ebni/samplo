%   Compute the optimal discrete-time Riccati matrixes
% INPUTS
%   n: dimension of the state
%   m: dimension of the input
%   N: number of sampling instants in [0, t_f]
%   Avec, Bvec: discrete-time dynamics (see the code below)
%   QDvec, RDvec, PDvec: distrete-time weighting matrices
%      (see computeDiscWeights.m for more details)
%   S: weight to the final state x(t_f)
% OUTPUTS
%   Kvec: sequence of discrete-time Riccati matrices
%   K0: Kvec(:,:,1) so that the cost is J(N) = x_0'*K0*x_0

Kvec=zeros(n,n,N+1);    % init
Kvec(:,:,N+1) = S;
for k=(N+1-(1:N))
  QD=QDvec(:,:,k);
  RD=RDvec(:,:,k);
  PD=PDvec(:,:,k);
  AD=Avec(:,:,k);
  BD=Bvec(:,:,k);
  KD=Kvec(:,:,k+1);
  curK = QD+AD'*KD*AD;
  curH = PD+BD'*KD*AD;
  curK = curK-curH'*inv(RD+BD'*KD*BD)*curH;
  Kvec(:,:,k) = curK;
  %sort(eig(curK))
end
K0 = Kvec(:,:,1);
clear k QD RD PD AD BD KD curK curH