%   Derivatives of Riccati matrices over the intersample separation
% INPUTS
%   Avec, Bvec: discrete-time dynamics (see the code below)
%   QDvec, RDvec, PDvec: distrete-time weighting matrices
%      (see computeDiscWeights.m for more details)
%   Kvec: sequence of discrete-time Riccati matrices
% OUTPUTS
%   DiffKvec: Diff(Kvec(:,:,k,h) is the derivative
%      of Kvec(:,:,k) w.r.t tau_h

DiffKvec=zeros(n,n,N+1,N);    % init
for k=(N+1-(1:N))
  QD=QDvec(:,:,k);
  RD=RDvec(:,:,k);
  PD=PDvec(:,:,k);
  AD=Avec(:,:,k);
  BD=Bvec(:,:,k);
  KD=Kvec(:,:,k+1);
  for h=1:N
    if (h<k)
      continue
    endif
    WD=QD+AD'*KD*AD;
    TD=RD+BD'*KD*BD;
    HD=PD+BD'*KD*AD;
    if (h==k)
      DiffWD = AD'*Q*AD+A'*(WD-QD)+(WD-QD)*A;
      DiffTD = R+BD'*Q*BD+B'*(HD-PD)'+(HD-PD)*B;
      DiffHD = BD'*Q*AD+B'*AD'*KD*AD+BD'*KD*AD*A;
    else % h>k
      DiffWD = AD'*DiffKvec(:,:,k+1,h)*AD;
      DiffTD = BD'*DiffKvec(:,:,k+1,h)*BD;
      DiffHD = BD'*DiffKvec(:,:,k+1,h)*AD;
    endif
    DiffKvec(:,:,k,h) = DiffWD-DiffHD'*inv(TD)*HD+HD'*inv(TD)*DiffTD*inv(TD)*HD-HD'*inv(TD)*DiffHD;
    %sort(eig(curK))
  endfor
endfor

clear k QD RD PD AD BD KD h WD TD HD DiffWD DiffTD DiffHD
