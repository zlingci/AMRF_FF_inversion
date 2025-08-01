function  MRF = RT_stack_normal_2D(MRFd,L,dx,dy,dt,dts,Vr)

L_num = length(L(:,1));
t_num = round(((length(MRFd(1,:))-1)*dts + max(abs(L(:,3)))/Vr)/dt) + 1;
ts_num = length(MRFd(1,:));

MRF = zeros(t_num,1);
%for j=1:1:length(MRFd(1,:))
%    for k=1:1:L_num
%        MRF_index = round(((j-1)*dts + abs(L(k,3))/Vr)/dt) + 1;
%        MRF(MRF_index) = MRF(MRF_index) + MRFd(k,j)*dx*dy;
%    end
%end
for t_index1=1:1:t_num
    for k=1:1:L_num
        for j=1:1:ts_num
            if round(((t_index1-1)*dt - abs(L(k,3)/Vr))/dts) + 1 == j
                MRF(t_index1) = MRF(t_index1) + MRFd(k,j)*dx*dy;
            end
        end
    end
end
