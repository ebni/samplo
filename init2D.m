%   Initialize the plant
% INPUTS
%   w: pulsation of the system
%   q: weight to the state (the weight to the input is assumed R=1)

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
% w = 25;     % frequency of resonance  w\neq 0
% q = 100;     % state weight

% system dynamics
A = [0 -w; w 0];
B = [1; 0];
[n,m] = size(B);

Phi = @(x) [cos(w*x) -sin(w*x); sin(w*x) cos(w*x)];
Gamma = @(x) 1/w*[sin(w*x); (1-cos(w*x))];

% boundary conditions
x_0 = [1; 0];    % WARNING: if changed then uOpt must be changed
t_f =1;

% init weight matrices
Q = q*eye(2);
R = 1;
dQ = @(x) q*x*eye(n);
dR = @(x) -2*q*sin(w*x)/w^3+2*q*x/w^2+x;
dP = @(x) [q/w^2-q*cos(w*x)/w^2; q*sin(w*x)/w^2-q*x/w];

% auxiliary variables
v = sqrt(w*w+q);
a = sqrt((v+3*w)*(v-w));
b = sqrt((v-3*w)*(v+w));
aa = 8*w^2/q;

% solution of ARE
Kinf = [a v-w; v-w a*v/w];

% the weight to the final state is set equal to ARE
S=Kinf;  % this way the optimal feedback is constant

% closed loop matrix
Acl = [-a -v; w 0];
expAcl = @(t) exp(-a.*t./2).*[cosh(b.*t./2)-a/b.*sinh(b.*t./2) -2*v/b*sinh(b.*t/2); 2*w/b*sinh(b.*t/2) cosh(b.*t/2)+a/b*sinh(b.*t/2)];    % OK

% setting the optimal input
if (q > 8*w*w)
  % if b is real <==> q > 8*w^2
  d = log((1+sqrt(1-aa))/sqrt(aa));
  uOpt = @(t) 2*sqrt(2)*exp(-a.*t./2).*sqrt((v-w)/(v-3*w))*w.*sinh(b.*t./2-d); % OK
elseif (q < 8*w*w)
  % if b is imaginary <==> q < 8*w^2
  b = abs(b);
  d = atan(sqrt(aa-1));
  uOpt = @(t) 2*sqrt(2)*exp(-a.*t./2).*sqrt((v-w)/(3*w-v))*w.*sin(b.*t./2-d);  % OK
else
  % if b == 0 <==> q = 8*w^2
  uOpt = @(t) (4*w*w.*t-2*w*sqrt(3)).*exp(-a.*t./2);  % OK
end

clear a aa b d v