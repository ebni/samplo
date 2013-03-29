%   Compute the optimal input sequence and the corresponding state
% INPUTS
%   n: dimension of the state
%   m: dimension of the input
%   N: number of sampling instants in [0, t_f]
%   x_0: initial state
%   Avec, Bvec: discrete-time dynamics (see the code below)
%   QDvec, RDvec, PDvec: distrete-time weighting matrices
%      (see computeDiscWeights.m for more details)
%   Kvec: sequence of discrete-time Riccati matrices
% OUTPUTS
%   Uvec: Uvec(:,k) is the optimal k-th input
%   Xvec: Xvec(:,k) is the k-th state is optimal input is applied

% Cvec = zeros(1,N+1);
% Cvec(1) = x_0'*Kvec(:,:,1)*x_0;
Uvec = zeros(m,N);
Xvec = zeros(n,N+1);
Xvec(:,1) = x_0;
for k=1:N
  RD=RDvec(:,:,k);
  PD=PDvec(:,:,k);
  AD=Avec(:,:,k);
  BD=Bvec(:,:,k);
  KD=Kvec(:,:,k+1);
  
  Uvec(:,k) = -inv(RD+BD'*KD*BD)*(PD+BD'*KD*AD)*Xvec(:,k);
  Xvec(:,k+1) = AD*Xvec(:,k)+BD*Uvec(:,k);

%  Cvec(k+1) = Xvec(:,k+1)'*KD*Xvec(:,k+1);
end
% minCost = Cvec(1);
clear RD PD AD BD KD k