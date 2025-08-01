function  [AMRF,G,Gw,data] = LSQ_G_1D(STFs,F,L,dL,t_num,ts_num,dt,dts,Vr,weight,damp)

F_num = length(F);
L_num = length(L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate G matrix and data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% G1 and data1 (G_m, G_1) RT %%%%%%%%%%%%%%%%%%%%%%%%% 
gnum_col = ts_num*L_num;
gnum_row = t_num*F_num;
G_m = zeros(gnum_row,gnum_col);
i = 0;
for F_index=1:1:F_num
    for t_index1=1:1:t_num
        i = i + 1;
        j = 0;
        for L_index=1:1:L_num
            for t_index2=1:1:ts_num
                j = j + 1;
                if t_index1 == round(((t_index2-1)*dts + abs(L(L_index))/Vr+F(F_index)*L(L_index))/dt) + 1 
                    G_m(i,j) = dL;
                else
                    G_m(i,j) = 0;
                end
            end
        end
    end
end
data_m = zeros(gnum_row,1);
for i=1:1:F_num
    for j=1:1:t_num
        data_m((i-1)*t_num+j) = STFs(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gnum_row1 = 0;
for i=1:1:F_num
    if round((max(L(:)/Vr+F(i)*L(:)) + (ts_num-1)*dts)/dt) + 1 >= 1 + t_num
        gnum_row1 = gnum_row1 + round((max(L(:)/Vr+F(i)*L(:)) + (ts_num-1)*dts)/dt) + 1 - t_num;
    elseif round(min(L(:)/Vr+F(i)*L(:))/dt) + 1 <= 1
        gnum_row1 = gnum_row1 + 2 - (round(min(L(:)/Vr+F(i)*L(:))/dt) + 1);
    end
end
G_1 = zeros(gnum_row1,gnum_col);
i = 0;
for F_index=1:1:F_num
    if  round((max(L(:)/Vr+F(F_index)*L(:)) + (ts_num-1)*dts)/dt) + 1  >= 1 + t_num
        for t_index1=t_num+1:1:round((max(L(:)/Vr+F(F_index)*L(:)) + (ts_num-1)*dts)/dt) + 1
            i = i + 1;
            j = 0;
            for L_index=1:1:L_num
                for t_index2= 1:1:ts_num
                    j = j + 1;
                    if t_index1 == round(((t_index2-1)*dts + (abs(L(L_index))/Vr+F(F_index)*L(L_index)))/dt) + 1 
                        G_1(i,j) = dL;
                    else
                        G_1(i,j) = 0;
                    end
                end
            end
        end
    elseif round(min(L(:)/Vr+F(F_index)*L(:))/dt) + 1 <= 1
        for t_index1=round(min(L(:)/Vr+F(F_index)*L(:))/dt)+1:1:1
            i = i + 1;
            j = 0;
            for L_index=1:1:L_num
                for t_index2= 1:1:ts_num
                    j = j + 1;
                    if t_index1 == round(((t_index2-1)*dts + (abs(L(L_index))/Vr+F(F_index)*L(L_index)))/dt) + 1 
                        G_1(i,j) = dL;
                    else
                        G_1(i,j) = 0;
                    end
                end
            end
        end
    end
end
data_1 = zeros(gnum_row1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% G2 and data2 (moment) %%%%%%%%%%%%%%%%%%%%%%%%
G2 = zeros(F_num,gnum_col);
for i=1:1:F_num
    for t_index1=1:1:t_num
        j = t_index1 + (i-1)*t_num;
        G2(i,:) = G2(i,:) + G_m(j,:)*dt;
    end
    j = 0; 
    for F_index=1:1:F_num
        if round((max(L(:)/Vr+F(F_index)*L(:)) + (ts_num-1)*dts)/dt) + 1 - t_num < 1
            continue;
        end       
        for t_index1=t_num+1:1:round((max(L(:)/Vr+F(F_index)*L(:)) + (ts_num-1)*dts)/dt) + 1
            j = j + 1;
            if F_index == i
                G2(i,:) = G2(i,:) + G_1(j,:)*dt;
            end
        end
    end
end
data2 = zeros(F_num,1);
for i=1:1:F_num
    data2(i) = sum(STFs(i,:))*dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G3 and data3 (tempotal smooth) %%%%%%%%%%%% 
G3 = zeros(gnum_col,gnum_col);
for j=1:1:gnum_col
    L_index = floor(j/ts_num) + 1;
    t_index2 = j - (L_index-1)*ts_num;
    if t_index2 == 0
        t_index2 = ts_num;
    end
    if t_index2 > 1 && t_index2 < ts_num
        G3(j,j) = -2;
        G3(j,j+1) = 1;
        G3(j,j-1) = 1;
    elseif t_index2 == 1
        G3(j,j) = 1e3;
    else
        G3(j,j) = 1;
        G3(j,j-1) = -1;
    end
end
data3 = zeros(gnum_col,1);
%%%%%%%%%%%%%%%%%% G4 and data4(spatial smooth) %%%%%%%%%%%%%%%%%%%%%%%%%
G4 = zeros(L_num,gnum_col);
for i=1:1:L_num
    if i==1 
        G4(i,(i - 1)*ts_num+1:i*ts_num) = -1;
        G4(i,(i + 1 -1)*ts_num+1:(i + 1)*ts_num) = 1;
    elseif i==L_num
        G4(i,(i - 1)*ts_num+1:i*ts_num) = -1;
        G4(i,(i - 1 -1)*ts_num+1:(i - 1)*ts_num) = 1;
    else
        G4(i,(i - 1 -1)*ts_num+1:(i - 1)*ts_num) = 1;
        G4(i,(i + 1 -1)*ts_num+1:(i + 1)*ts_num) = 1;
        G4(i,(i - 1)*ts_num+1:i*ts_num) = -2;
    end
end
data4 = zeros(L_num,1);
% conbine G and data
G = zeros(gnum_row+gnum_row1+F_num+gnum_col+L_num,gnum_col);
G(1:gnum_row,:) = G_m;
G(gnum_row+1:gnum_row+gnum_row1,:) = G_1;
G(gnum_row+gnum_row1+1:gnum_row+gnum_row1+F_num,:) = G2;
G(gnum_row+gnum_row1+F_num+1:gnum_row+gnum_row1+F_num+gnum_col,:) = G3;
G(gnum_row+gnum_row1+F_num+gnum_col+1:gnum_row+gnum_row1+F_num+gnum_col+L_num,:) = G4;
data = zeros(gnum_row+gnum_row1+F_num+gnum_col+L_num,1);
data(1:gnum_row) = data_m;
data(gnum_row+1:gnum_row+gnum_row1) = data_1;
data(gnum_row+gnum_row1+1:gnum_row+gnum_row1+F_num) = data2;
data(gnum_row+gnum_row1+F_num+1:gnum_row+gnum_row1+F_num+gnum_col) = data3;
data(gnum_row+gnum_row1+F_num+gnum_col+1:gnum_row+gnum_row1+F_num+gnum_col+L_num) = data4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate weight matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gw = zeros(gnum_row+gnum_row1+F_num+gnum_col+L_num,gnum_row+gnum_row1+F_num+gnum_col+L_num);
NORM = norm(G_m)^2;
NORM1 = norm(G_1)^2;
NORM2 = norm(G2)^2;
NORM3 = norm(G3)^2;
NORM4 = norm(G4)^2;
for i=1:1:gnum_row+gnum_row1+F_num+gnum_col+L_num
    if i <= gnum_row
        F_index = floor(i/t_num) + 1;
        t_index1 = i - (F_index-1)*t_num;
        if t_index1 == 0
            F_index = F_index - 1;
        end
        Gw(i,i) = weight(F_index);
    elseif i <= gnum_row+gnum_row1 && i > gnum_row
        Gw(i,i) = 1e3*NORM/NORM1;
    elseif i <= gnum_row+gnum_row1+F_num && i > gnum_row+gnum_row1
        Gw(i,i) = damp(1)*NORM/NORM2;
    elseif i <= gnum_row+gnum_row1+F_num+gnum_col && i > gnum_row+gnum_row1+F_num
        Gw(i,i) = damp(2)*(dt/dts)^2*NORM/NORM3;
    elseif  i <= gnum_row+gnum_row1+F_num+gnum_col+L_num && i > gnum_row+gnum_row1+F_num+gnum_col
        Gw(i,i) = damp(3)*NORM/NORM4; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate AMRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_num1 = fix(0.15*t_num);
t_num1 = 0;
t_num2 = 0;
% for i=1:1:F_num
%     if fix(max(L(:)/Vr+F(i)*L(:))/dt) + ts_num -t_num > t_num2
%         t_num2 = fix(max(L(:)/Vr+F(i)*L(:))/dt)+ts_num-t_num;
%     end
% end
AMRF = zeros(F_num,t_num+t_num1+t_num2);
% k = 0;
for i=1:1:F_num
    AMRF(i,t_num1+1:t_num1+t_num) = data_m((i-1)*t_num+1:i*t_num);
%     if fix(max(L(:)/Vr+F(i)*L(:))/dt) + ts_num -t_num >= 1
%         MRF(i,t_num1+t_num+1:t_num1+t_num+fix(max(L(:)/Vr+F(i)*L(:))/dt)+ts_num-t_num) = data_1(k+1:k+fix(max(L(:)/Vr+F(i)*L(:))/dt)+ts_num-t_num);  
%         k = k + fix(max(L(:)/Vr+F(i)*L(:))/dt) + ts_num -t_num;
%     end
end
