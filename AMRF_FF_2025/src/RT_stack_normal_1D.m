function  MRF = RT_stack_normal_1D(MRFd,L,dL,dt,dts,Vr)

L_num = length(L);

t_num = round(((length(MRFd(1,:))-1)*dts + max(abs(L))/Vr)/dt) + 1;
MRF = zeros(1,t_num);  
for j=1:1:length(MRFd(1,:))
    for k=1:1:L_num
        MRF_index = round(((j-1)*dts + abs(L(k))/Vr)/dt) + 1;
        MRF(MRF_index) = MRF(MRF_index) + MRFd(k,j)*dL;
    end
end
 