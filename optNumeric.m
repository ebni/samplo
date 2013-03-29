%   Optimal Sampling (with gradient descent)
% INPUTS
%   n: dimension of the state
%   m: dimension of the input
%   t_f: final instant
%   N: number of sampling instants in [0, t_f]
%   Phi(t), Gamma(t): discretization over interval of lenght t
%   Q, R: weight to the state and to the input
%   S: weight to the final state x(t_f)
%   x_0: initial state
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation
%      (initial guess)
% OUTPUTS
%   K0: cost matric so that the cost is J(N) = x_0'*K0*x_0
%   Uvec: Uvec(:,k) is the optimal k-th input
%   minCost: cost with this sampling method
%   tK: tK(i) is the i-th sampling instant (vector of size N+1)
%   tauK: tauK(i), with i=1,...,N, is the i-th intersample separation

iter=0;
minCost = +inf;

lScale = 1;
factorScale = 2;

% MAGIC NUMBER here
minTau = t_f/(N*1000);    % cannot smaller than this

normGvec = [];
distVec = [];
while (1)
  %disp("------------------------------------------------------");
  %disp(iter);
  computeDiscDyn;
  computeDiscWeights;
  computeDiscRiccati;

  %disp("minCost");
  %disp(minCost);
  
  curCost = x_0'*K0*x_0;
  %disp("curCost");
  %disp(curCost);
  if (curCost - minCost < 0)  % take a step along the -gradient
    %disp("good step");
    % update the minimum
    minSol = tauK;
    % MAGIC NUMBER here
%    if ((minCost-curCost)/minCost < 1e-10)    % however if the improvement is small we can stop
 %     exitMsg = "lowErr";
  %    break
   % endif
    minCost = curCost;

    computeDiscRiccatiDiff;

    % compute the gradient at the current solution
    gradCost = zeros(1,N);
    for k=1:N
      gradCost(k) = x_0'*DiffKvec(:,:,1,k)*x_0;
    endfor

    % STOP CONDITION (1): gradient orthogonal to constraint
    cosAngle = sum(gradCost)/(sqrt(N)*norm(gradCost));
    if (1-cosAngle < 1e-12)
      exitMsg = "gradOrtho";
      break;
    endif

    %disp("tauK");
    %disp(tauK);
    %disp("gradCost");
    %disp(gradCost);

    % the projection matrix over sum(tauK) = constant
    %H=null(conProj);
    proj = eye(N)-1/N*ones(N,N);

    % the projected gradient
    gradCostPr = gradCost*proj;
    normGrad = norm(gradCostPr);   % and its norm
    %disp("normGrad");
    %disp(normGrad);
    normGvec = [normGvec normGrad];

    % STOP CONDITION (2): projected gradient close to zero
    if (normGrad < 1e-11)
      exitMsg = "smallGrad";
      break
    endif

    %disp("gradCostPr");
    %disp(gradCostPr);

    % compute the longest step along -gradCostPr
    mayHit = setdiff((gradCostPr >0).*(1:N),[0]);

    %disp("mayHit");
    %disp(mayHit);
    if (!isempty(mayHit))
      allScale = (tauK(mayHit)-ones(1,length(mayHit))*minTau)./(gradCostPr(mayHit));

      %disp("allScale");
      %disp(allScale);

      maxScale = min(allScale);
    else
      maxScale = inf;
      disp("WARNING: this should never be possible because feasible \
	  region is bounded")
    endif

    % compute the step length
    if (maxScale<lScale/normGrad)
      %disp("too long");
      stepL = maxScale;
      lScale = lScale/factorScale;  % next time take shorter step
    else
      %disp("short enough");
      stepL = lScale/normGrad;
      lScale = lScale*factorScale;  % next time take longer step
    endif
    %disp("stepL");
    %disp(stepL);

    % finally, take a step
    newTauK = tauK-gradCostPr*stepL;
    newTauK(N) = t_f-sum(newTauK(1:N-1));  % in case the sum is no longer t_f
    %disp("dist");
    %disp(dist);

    % STOP CONDITION (3): new solution is very similar to old one
    dist = norm(newTauK-tauK);
    distVec = [distVec dist];
    if (dist < 1e-11)
      exitMsg = "smallStep";
      break
    endif
    tauK = newTauK;
  else  % go back on your steps
    %disp("bad step");
    % step in between is a good idea
    % MAGIC NUMBER here
    % tauK = .5*(tauK+minSol);
    tauK = 0.2*tauK+0.8*minSol;
    
    lScale = lScale/factorScale*0.8;  % next time take shorter step
  endif
  
  % STOP CONDITION (4): too many iterations
  if (iter > 10)
    exitMsg = "maxIter";
    break
  endif
  iter = iter+1;

  %disp("leqIsEq");
  %disp(leqIsEq);
endwhile

tauK = minSol;
computeDiscDyn;                 % the discrete dynamics matrices
computeDiscWeights;             % the discrete-time weights matrices
computeDiscRiccati;             % the solution of the discrete Riccati Eqs
computeOptUk;                   % the system dynamics and optimal input
minCost = x_0'*K0*x_0;
