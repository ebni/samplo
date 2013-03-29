% Plotting sampled-time control input
% INPUTS
%     N tauK Uvec
% OUTPUTS
%     

tK = zeros(1,N+1);
for j=1:N
  tK(j+1) = sum(tauK(1:j));
endfor
hold on
for k=1:N,
  plot([tK(k) tK(k+1)],[Uvec(:,k) Uvec(:,k)],'k');
endfor

%  print(strcat("plotOpt",num2str(N),"TF2.fig"),"-dfig");
