% general settings
w = 25;
q = 100;
init2D;
costCont = x_0'*Kinf*x_0;
numSteps = 1000000;                  % steps of the function
allT = linspace(0,t_f,numSteps+1);
uDot = diff(uOpt(allT))./diff(allT);    % needed in sampleUmA

vecAlpha = linspace(0,1,100);

N = 512;                            % number of samples
for iter=1:length(vecAlpha)
 mA = vecAlpha(iter);
 sampleUmA;
 vec1(iter) = cMa;
end
figure(1);
semilogy(vecAlpha,vec1,'b');
hold on

N = 256;                            % number of samples
for iter=1:length(vecAlpha)
  mA = vecAlpha(iter);
  sampleUmA;
  vec2(iter) = cMa;
end
figure(1);
semilogy(vecAlpha,vec2,'r');

N = 128;                            % number of samples
for iter=1:length(vecAlpha)
  mA = vecAlpha(iter);
  sampleUmA;
  vec3(iter) = cMa;
end
figure(1);
semilogy(vecAlpha,vec3,'k');

N = 64;                            % number of samples
for iter=1:length(vecAlpha)
  mA = vecAlpha(iter);
  sampleUmA;
  vec4(iter) = cMa;
end
figure(1);
semilogy(vecAlpha,vec4,'b--');

N = 32;                            % number of samples
for iter=1:length(vecAlpha)
  mA = vecAlpha(iter);
  sampleUmA;
  vec5(iter) = cMa;
end
figure(1);
semilogy(vecAlpha,vec5,'r--');

N = 16;                            % number of samples
for iter=1:length(vecAlpha)
  mA = vecAlpha(iter);
  sampleUmA;
  vec6(iter) = cMa;
end
figure(1);
semilogy(vecAlpha,vec6,'k--');

N = 8;                            % number of samples
for iter=1:length(vecAlpha)
  mA = vecAlpha(iter);
  sampleUmA;
  vec7(iter) = cMa;
end
figure(1);
semilogy(vecAlpha,vec7,'g');

N = 4;                            % number of samples
for iter=1:length(vecAlpha)
  mA = vecAlpha(iter);
  sampleUmA;
  vec8(iter) = cMa;
end
figure(1);
semilogy(vecAlpha,vec8,'m');

clear AD Avec BD Bvec Cvec K0 KD Kvec PD PDvec QD QDvec RD RDvec Uvec Xvec cMa costCont curH curK densToSamp i ind intDens iter k mA minCost lll tk tauK uDot uDotM uOpt
save('CvsAlphaW25q100.mat');
