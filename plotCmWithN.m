% general settings
Nvec = [8 16 32 64 128 256 512];

% first-order case to be tested
A = 1;
Q = 8;
init1D

% second-order case to be tested
%w = 5;
%q = 100;
%init2D;

%numSteps = 1000000;
numSteps = 10000;
allT = linspace(0,t_f,numSteps+1);
uDot = diff(uOpt(allT))./diff(allT);

% initialize the data struct to store data
discSolVec = [];
myDiscSol = struct('N', 0, 'tauK', [], 'Uvec', [], 'cost', 0, 'K0', zeros(n,n), 'method', '  ', 'coef', 0);
mySolution = myDiscSol;

mySolution.N = Inf;
mySolution.K0 = Kinf;
mySolution.cost = x_0'*Kinf*x_0;
mySolution.method = 'inf';
costInf = mySolution.cost;

discSolVec = [discSolVec mySolution];
cPerVec = [];
cDLsVec = []
cQ23Vec = [];
cQntVec = [];
cNumVec = [];
for N=Nvec,
  % provide initial guess for tauK
  tK = linspace(0,t_f,N+1);
  tauK = diff(tK);

  % Periodic sampling
  optPeriodic;
  mySolution.N = N;
  mySolution.tauK = tauK;
  mySolution.Uvec = Uvec;
  mySolution.K0 = K0;
  mySolution.cost = minCost;
  mySolution.method = 'per';
  % cPer = N*N*(minCost-costCont)/costCont/t_f/t_f;
  mySolution.coef = (minCost-costInf)/costInf*N*N/t_f/t_f;
  discSolVec = [discSolVec mySolution];
  cPerVec = [cPerVec mySolution.coef];

  % Deterministic Lebesgue sampling
  mA = 1;
  optDensM;
  mySolution.N = N;
  mySolution.tauK = tauK;
  mySolution.Uvec = Uvec;
  mySolution.K0 = Kvec(:,:,1);
  mySolution.cost = minCost;
  mySolution.method = 'dLs';
  mySolution.coef = (minCost-costInf)/costInf*N*N/t_f/t_f;
  discSolVec = [discSolVec mySolution];
  cDLsVec = [cDLsVec mySolution.coef];
  
  % Quantization-based sampling (based on the asymptotic density)
  mA = 2/3;
  optDensM;
  mySolution.N = N;
  mySolution.tauK = tauK;
  mySolution.Uvec = Uvec;
  mySolution.K0 = Kvec(:,:,1);
  mySolution.cost = minCost;
  mySolution.method = 'q23';
  mySolution.coef = (minCost-costInf)/costInf*N*N/t_f/t_f;
  discSolVec = [discSolVec mySolution];
  cQ23Vec = [cQ23Vec mySolution.coef];
  
  % Quantization-based sampling (based on Newton's method)
  optQuanta;
%  exitMsg
  mySolution.N = N;
  mySolution.tauK = tauK;
  mySolution.Uvec = Uvec;
  mySolution.K0 = Kvec(:,:,1);
  mySolution.cost = minCost;
  mySolution.method = 'qnt';
  mySolution.coef = (minCost-costInf)/costInf*N*N/t_f/t_f;
  discSolVec = [discSolVec mySolution];
  cQntVec = [cQntVec mySolution.coef];
  
  % Choosing the best available initial guess for the numerical algorithm
  bestCost = +inf;
  bestInd = 0;
  for i=1:length(discSolVec),
    curSol = discSolVec(i);
    if (curSol.N<=N && curSol.cost<bestCost)
      bestInd = i;
    end
  end
  if (bestInd == 0)
    tK = linspace(0,t_f,N+1);
    tauK = diff(tK);   % just in case it can find anything better
  else
    curSol = discSolVec(bestInd);
    if (curSol.N = N)
      tauK = curSol.tauK;
    else  % curSol.N =< N
      tauK = curSol.tauK;
      tKfromTauK;
      c = ceil(N/curSol.N);
      curTk = zeros(1,c*curSol.N+1);
      for i=1:curSol.N
	curTk((1+(i-1)*c):(i*c+1)) = linspace(tK(i),tK(i+1),c+1);
      end
      % truncating curTk, it can be done better
      curTk = [curTk(1:N) t_f];
      tauK = diff(curTk);
    end
  end
  optNumeric;
%  exitMsg
  mySolution.N = N;
  mySolution.tauK = tauK;
  mySolution.Uvec = Uvec;
  mySolution.K0 = Kvec(:,:,1);
  mySolution.cost = minCost;
  mySolution.method = 'num';
  mySolution.coef = (minCost-costInf)/costInf*N*N/t_f/t_f;
  discSolVec = [discSolVec mySolution];
  cNumVec = [cNumVec mySolution.coef];
end

semilogy(cPerVec,'o-');
hold on
semilogy(cDLsVec,'+-');
semilogy(cQ23Vec,'x-');
semilogy(cQntVec,'kx-');
semilogy(cNumVec,'ko-');
set(gca,"xlim",[0.5 7.5]);
%set(gca,"ylim",[0.5 8]);
print("CvsN.fig","-dfig");
hold off
