function  [AMRF_est,MRFd_est] = RT_matrix_2D(model,F1,F2,L,t_num,ts_num,G,tri_num)
 
F_num = length(F1(:,1));
L_num = length(L(:,1));
FL = zeros(F_num,L_num);
for i=1:1:F_num
    for j=1:1:L_num
        FL(i,j) = F1(i,j)*L(j,1)+ F2(i,j)*L(j,2);
    end
end
point_num = 2*round((ts_num-1)/(tri_num + 1)) + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating MRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = G*model;
AMRF_est = zeros(F_num,t_num);
for i=1:1:F_num
    AMRF_est(i,1:t_num) = data((i-1)*t_num+1:i*t_num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate MRFd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gnum_col = tri_num;
gnum_row = ts_num;
G_m = zeros(gnum_row,gnum_col);
for i=1:1:gnum_row
    for j=1:1:gnum_col
        if  i <=  j*(point_num-1)/2 + 1 && i >= (j-1)*(point_num-1)/2 + 1
            G_m(i,j) = (i-(j-1)*(point_num-1)/2-1)/((point_num-1)/2);
        elseif  i >  j*(point_num-1)/2 + 1 && i <= (j+1)*(point_num-1)/2 + 1     
            G_m(i,j) = 1-(i-j*(point_num-1)/2-1)/((point_num-1)/2);
        else
            G_m(i,j) = 0;
        end
    end
end
MRFd_est = zeros(L_num,ts_num);
for i=1:1:L_num
    MRFd_est(i,:) = G_m*model((i-1)*tri_num+1:i*tri_num);
end

