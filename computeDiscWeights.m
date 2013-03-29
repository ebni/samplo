%   Compute the weights of discrete-time
% INPUTS
%   n: dimension of the state
%   m: dimension of the input
%   N: number of sampling instants in [0, t_f]
%   Phi(t), Gamma(t): discretization over interval of lenght t
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation
%   Q, R: weight to the state and to the input
% OUTPUTS
%   QDvec, RDvec, PDvec: distrete-time weighting matrices so that
%       the cost of the discrete-time action is
%       J(N) = sum_k x(k)'*QDvec(:,:,k)*x(k)
%            + 2*x(k)'*PDvec(:,:,k)*u(k)+u(k)'*RDvec(:,:,k)*u(k)
%            + x(N)'*S*x(N)

QDvec = zeros(n,n,N);
RDvec = zeros(m,m,N);
PDvec = zeros(m,n,N);
for k=1:N,
  QDvec(:,:,k) = dQ(tauK(k));
  RDvec(:,:,k) = dR(tauK(k));
  PDvec(:,:,k) = dP(tauK(k));
end
