% INPUTS
%   A: system dynamics is \dot x = A x + u
%   Q: weight to the state (the weight to the input is assumed R=1)

% OUTPUTS
%   A, B: system dynamics is \dot x = Ax+Bu
%   n: dimension of the state
%   m: dimension of the input
%   Phi(t), Gamma(t): discretization over interval of lenght t
%   x_0: initial state
%   t_f: final instant
%   Q, R: weight to the state and to the input
%   dQ(t), dR(t), dQ(t): discretization of the state, input, and cross-term
%       weights over an interval of length t
%   Kinf: solution of ARE
%   S: weight to the final state x(t_f), assumed S=Kinf
%   Acl: closed loop matrix, assuming optimal state feedback
%   expAcl(t): analytical expression of expm(Acl*t)
%   uOpt(t): optimal input, assuming x_0 = [1;0]

% parameters of the problem
% A = 1;
% Q = 8;

% system dynamics
A = 1;
B = 1;
[n,m] = size(B);

if (A==0)
  Phi = inline("1","x");
  Gamma = inline("x*B","x");
else
  Phi = inline("exp(A*x)","x");
  Gamma = inline("(exp(A*x)-1)*B/A","x");
endif

% boundary conditions
x_0 = 1;
t_f =1;

% init weight matrices
Q = 8;
R = 1;
if (A==0)
  dQ = @(x) Q*x;
  dR = @(x) R*x+x*x*x*Q*B*B/3;
  dP = @(x) .5*x*x*Q*B;
else
  dQ = @(x) .5/A*(exp(2*A*x)-1)*Q;
  dR = @(x) R*x+.5*B*B*Q/(A*A*A)*(exp(2*A*x)-4*exp(A*x)+2*A*x+3);
  dP = @(x) .5*Q*B/(A*A)*(exp(2*A*x)-2*exp(A*x)+1);
endif

% solution of ARE
Kinf = A+sqrt(A*A+Q);

% the weight to the final state is set equal to ARE
S=Kinf;  % this way the optimal feedback is constant

% closed loop matrix
Acl = -sqrt(A*A+Q);
expAcl = @(t) exp(-Acl.*t);

% setting the optimal input
uOpt = @(t) -(A+sqrt(A*A+Q)).*exp(-sqrt(A*A+Q)*t);
