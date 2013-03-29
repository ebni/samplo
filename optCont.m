% Compute the continuous dynamics of optimal input
% INPUTS
%     A B n m x_0 t_f Q R S Kinf
% OUTPUTS
%     uFun xFun K0 costCont allT

numSteps = floor(max(abs(eig(A)))*100);   % related to the fastest eig
allT = linspace(0,t_f,numSteps+1);   % instants of continuous-time simulation


if (S == Kinf)
  % if S=Kinf, the optimal continuous-time feedback is constant
  L = -inv(R)*B'*Kinf;

  % the optimal state, input
  for i=1:length(allT)
    xFun(:,i) = expm((A+B*L)*allT(i))*x_0;
    uFun(:,i) = L*xFun(:,i);
  endfor

  % the cost matrix at 0
  K0 = Kinf;
else
  Tric=t_f-linspace(0,t_f,100);  % time goes backward
  [KcontVec, ISTATE, MSG] = lsode(@(x,t) riccatiDifEq(x, A, B, Q, R), 
				  mat2col(S), 
				  Tric);
  if (ISTATE != 2)
    MSG
  endif
    % compute the system dynamics
  [xFun, istate, msg] = lsode (@(x,t) (A-B*inv(R)*B'*col2mat(interp1(Tric,KcontVec,t)',n))*x,
			    x_0, allT);
  xFun = xFun';
  uFun = zeros(m,length(allT));
  for i=1:length(allT),
    uFun(:,i) = -inv(R)*B'*col2mat(interp1(Tric,KcontVec,allT(i))',n)*xFun(:,i);
  endfor
  K0 = col2mat(interp1(Tric,KcontVec,0)',n);
endif

costCont = x_0'*K0*x_0;
