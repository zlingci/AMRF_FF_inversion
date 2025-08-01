function  [AMRF_est,MRFd_est] = RT_matrix_1D(model,F,L,t_num,ts_num,G)
 
F_num = length(F);
L_num = length(L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating AMRF_est
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = G*model;
AMRF_est = zeros(F_num,t_num);
for i=1:1:F_num
    AMRF_est(i,1:t_num) = data((i-1)*t_num+1:i*t_num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating MRFd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 MRFd_est = zeros(L_num,ts_num);
for i=1:1:L_num
    MRFd_est(i,1:ts_num) = model((i-1)*ts_num+1:i*ts_num);
end

end
