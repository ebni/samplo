% general settings
N = 500;                            % number of samples
numSteps = N*1000;                  % steps of the function
vecAlpha = linspace(0,1,100);

% cases to be tested
w = 5;
q = 1;
init2D;
allT = linspace(0,t_f,numSteps+1);
% costCont = x_0'*Kinf*x_0;
% uDot = diff(uOpt(allT))./diff(allT);    % needed in optDensM
% for iter=1:length(vecAlpha)
%   mA = vecAlpha(iter);
%   optDensM;
%   cMa = N*N*(minCost-costCont)/(costCont*t_f*t_f);
%   vec1(iter) = cMa;
% end
figure(1);
plot(allT,uOpt(allT),'r');
hold on
figure(2);
semilogy(vecAlpha,vec1,'r');
hold on

q = 10;
init2D;
% costCont = x_0'*Kinf*x_0;
% uDot = diff(uOpt(allT))./diff(allT);    % needed in optDensM
% for iter=1:length(vecAlpha)
%   mA = vecAlpha(iter);
%   optDensM;
%   cMa = N*N*(minCost-costCont)/(costCont*t_f*t_f);
%   vec2(iter) = cMa;
% end
figure(1);
plot(allT,uOpt(allT),'k');
figure(2);
semilogy(vecAlpha,vec2,'k');

q = 100;
init2D;
% costCont = x_0'*Kinf*x_0;
% uDot = diff(uOpt(allT))./diff(allT);    % needed in optDensM
% for iter=1:length(vecAlpha)
%   mA = vecAlpha(iter);
%   optDensM;
%   cMa = N*N*(minCost-costCont)/(costCont*t_f*t_f);
%   vec3(iter) = cMa;
% end
figure(1);
plot(allT,uOpt(allT),'b');
figure(2);
semilogy(vecAlpha,vec3,'b');

w = 25;
q = 1;
init2D;
% costCont = x_0'*Kinf*x_0;
% uDot = diff(uOpt(allT))./diff(allT);    % needed in optDensM
% for iter=1:length(vecAlpha)
%   mA = vecAlpha(iter);
%   optDensM;
%   cMa = N*N*(minCost-costCont)/(costCont*t_f*t_f);
%   vec4(iter) = cMa;
% end
figure(1);
plot(allT,uOpt(allT),'r--');
hold on
figure(2);
semilogy(vecAlpha,vec4,'r--');
hold on

q = 10;
init2D;
% costCont = x_0'*Kinf*x_0;
% uDot = diff(uOpt(allT))./diff(allT);    % needed in optDensM
% for iter=1:length(vecAlpha)
%   mA = vecAlpha(iter);
%   optDensM;
%   cMa = N*N*(minCost-costCont)/(costCont*t_f*t_f);
%   vec5(iter) = cMa;
% end
figure(1);
plot(allT,uOpt(allT),'k--');
figure(2);
semilogy(vecAlpha,vec5,'k--');

q = 100;
init2D;
% costCont = x_0'*Kinf*x_0;
% uDot = diff(uOpt(allT))./diff(allT);    % needed in optDensM
% for iter=1:length(vecAlpha)
%   mA = vecAlpha(iter);
%   optDensM;
%   cMa = N*N*(minCost-costCont)/(costCont*t_f*t_f);
%   vec6(iter) = cMa;
% end
figure(1);
plot(allT,uOpt(allT),'b--');
figure(2);
semilogy(vecAlpha,vec6,'b--');

clear AD Avec BD Bvec Cvec K0 KD Kvec PD PDvec QD QDvec RD RDvec Uvec Xvec cMa costCont curH curK densToSamp i ind intDens iter k mA minCost lll tk tauK uDot uDotM uOpt
save('CvsAlphaN500.mat');
