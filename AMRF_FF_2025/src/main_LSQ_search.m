clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input parameter %%%%%%%%%%%%%%%%%%%%%%%5%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputf = strcat('LSQ_search.txt');
strminputf = fopen(inputf,'r');
inputfiles = textscan(strminputf,'%s %*[^\n]',1);
inputfiles = char(string(inputfiles));
outputfiles = textscan(strminputf,'%s %*[^\n]',1);
outputfiles = char(string(outputfiles));
weightf = textscan(strminputf,'%s %*[^\n]',1);
weightf = char(string(weightf));
state_data = textscan(strminputf,'%f %*[^\n]',1);
state_data = cell2mat(state_data);
DT = textscan(strminputf,'%f %f %*[^\n]',1);
DT = cell2mat(DT);
dt = DT(1);
dts = DT(2);
duration = textscan(strminputf,'%f %*[^\n]',1);
duration = cell2mat(duration);
V = textscan(strminputf,'%f %f %f %*[^\n]',1);
V = cell2mat(V);
Vp = V(1);
Vs = V(2);
density = V(3);
PLANE = textscan(strminputf,'%f %f %*[^\n]',1);
PLANE = cell2mat(PLANE);
strike = PLANE(1);
dip = PLANE(2);
VR = textscan(strminputf,'%f %f %f %*[^\n]',1);
VR = cell2mat(VR);
RANK = textscan(strminputf,'%f %f %f %f %*[^\n]',1);
RANK = cell2mat(RANK);
Moment = textscan(strminputf,'%f %f %*[^\n]',1);
Moment = cell2mat(Moment);
ts_max = textscan(strminputf,'%f %*[^\n]',1);
ts_max = cell2mat(ts_max);
state_model = textscan(strminputf,'%f %*[^\n]',1);
state_model = cell2mat(state_model);
damp = textscan(strminputf,'%f %f %f %*[^\n]',1);
damp = cell2mat(damp);
abs_dir = textscan(strminputf,'%s %*[^\n]',1);
abs_dir = char(string(abs_dir));
fclose(strminputf);
inputfiles = [abs_dir '/' inputfiles];
outputfiles = [abs_dir '/' outputfiles];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get STFs, dirctivity parameter (F), model setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting Moment %%%%%%%%%%%%%%%%%%%%%%%%%%%
if Moment(2) ~= 0.0
    M01 = Moment(1);
    M02 = Moment(2);
else
    M02 = fix((Moment(1)+10.73)*3/2) - 7;
    M01 = 10^((Moment(1)+10.73)*3/2 - fix((Moment(1)+10.73)*3/2));
end
M0 = M01*10^M02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting normalized STFs %%%%%%%%%%%%%%%%%%%%
% weightf = strcat(inputfiles,'/',weightf);
[egfname,sname,cha,dist,az,takeoff,weight] = textread(weightf,'%s %s %s %f %f %f %f','headerlines',1);
F_num = length(az);
if state_data == 1
    for i=1:1:F_num
        sacf = strcat(inputfiles,'/',char(egfname(i)),'.',char(sname(i)),'.',char(cha(i)),'.out');
        [hd,data] = irdsac(sacf);
        t_start = round((hd.a - hd.b)/hd.delta) + 1;
        if duration==0
            t_end = round((hd.e - hd.b)/hd.delta) + 1;
        else
            t_end = round((hd.a + duration - hd.b)/hd.delta) + 1;
        end
        if i==1
            if dt == 0
                dt = hd.delta;
            end
            if dt ~= hd.delta  
                t_num = round((t_end - t_start)*hd.delta/dt) + 1;
                t = 0:dt:(t_num-1)*dt;
            else
                t_num = t_end - t_start + 1;
                t = 0:hd.delta:(t_num-1)*hd.delta;
            end
            STFs = zeros(length(az),t_num);
        end
        if t_end > length(data)
            data_temp = zeros(t_end-t_start+1,1);
            data_temp(1:length(data)-t_start+1) = data(t_start:end);
        else
            data_temp = data(t_start:t_end);
        end
        t_t1 = round((hd.t1 - hd.b)/hd.delta) + 1;
        if t_t1 <= 0
            t_t1 = t_end;
        end
        for k=1:1:t_end-t_start+1
            if data_temp(k)*10^6 <= 0 || k + t_start - 1 >= t_t1 || k == 1
                data_temp(k) = 0;
            end
        end
        if abs(dt - hd.delta) > 10^-6
            STFs(i,:) = spline(0:hd.delta:(t_end - t_start)*hd.delta,data_temp,t);   
            for k=1:1:length(STFs(i,:))
                if STFs(i,k) < 0
                     STFs(i,k) = 0;
                end
            end
        else
            STFs(i,:) = data_temp;
        end
    end
    for i=1:1:F_num
        STFs(i,:) = STFs(i,:)*M01/(sum(STFs(i,:))*dt);
    end
else
    MRFf = strcat(inputfiles,'/MRF.txt');
    data = load(MRFf);
    STFs = zeros(F_num,length(data(:,1)));
    for i=1:1:F_num
        STFs(i,:) = data(:,i+1);
    end
    for i=1:1:F_num
        STFs(i,:) = STFs(i,:)*M01/(sum(STFs(i,:))*dt);
    end
    t = data(:,1);
    t_num = length(t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Get dirctivity parameter (F) %%%%%%%%%%%%%%%%%%
rank = RANK(1):RANK(3):RANK(2);
Vr = VR(1):VR(3):VR(2);
rank_num = length(rank);
Vr_num = length(Vr);
ph = zeros(F_num,1);
pv = zeros(F_num,1);
for i=1:1:F_num
    if strcmp(char(cha(i)),'zp') || strcmp(char(cha(i)),'rp')
        ph(i) = sin(takeoff(i)*pi/180)/Vp;
        pv(i) = cos(takeoff(i)*pi/180)/Vp;
    else
        ph(i) = sin(takeoff(i)*pi/180)/Vs;
        pv(i) = cos(takeoff(i)*pi/180)/Vs;
    end
end
if exist(outputfiles,'dir') == 0 
    mkdir(outputfiles)
end
VRf = strcat(outputfiles,'/VR.txt');
strmVRf = fopen(VRf,'w');
for m=1:1:rank_num
    F = zeros(F_num,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%% Get dirctivity parameter (F) %%%%%%%%%%%%%%%%%% 
    for i=1:1:F_num
        F(i) = -ph(i)*cos((strike-az(i))*pi/180)*cos(rank(m)*pi/180) - ph(i)*sin((strike-az(i))*pi/180)*sin(rank(m)*pi/180)*sin((90-dip)*pi/180) + pv(i)*sin(rank(m)*pi/180)*cos((90-dip)*pi/180);
    end
    for n=1:1:Vr_num
        %%%%%%%%%%%%%%%%%%%%%%%% Constructing model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L_e = t(end)/(1.2/Vr(n) + max(abs(F)));
        if state_model == 1 
            L_num = 31;
            dL = L_e/(L_num - 1);
            L = zeros(L_num,1);
            for i=1:1:L_num
                L(i) = (i-1)*dL;
            end
        else
            L_s = -1*L_e;
            L_num = 61;
            dL = L_e*2/(L_num - 1);
            L = zeros(L_num,1);
            for i=1:1:L_num
                L(i) = (i - 1 - (L_num-1)/2)*dL;
            end 
        end
        ts_num = round(0.3*L_e/Vr(n)/dts) + 1;
        ts = 0:dts:(ts_num-1)*dts;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%% 1D Inversion on the fault plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Constructing data, and G %%%%%%%%%%%%%%%%%%%%%%% 
        [AMRF,G,Gw,data] = LSQ_G_1D(STFs,F,L,dL,t_num,ts_num,dt,dts,Vr(n),weight,damp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NNLSQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        G1 = sqrt(Gw)*G;
        data1 = sqrt(Gw)*data;
        [model,resnorm,residual,exitflag,output,lambda] = lsqnonneg(G1,data1,optimset('Display','notify'));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Synthetic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [AMRF_est,MRFd_est] = RT_matrix_1D(model,F,L,t_num,ts_num,G);
        AMRF_est_normal = RT_stack_normal_1D(MRFd_est,L,dL,dt,dts,Vr(n));
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output Variance reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        M0_est = sum(AMRF_est_normal)*dt*10^M02;
        M0_res = (1-abs(M0_est - M0)/M0)*100;
        VRall = 0;
        powerall = 0;
        for i=1:1:F_num
            powerall = powerall + sum(AMRF(i,:).*AMRF(i,:))*weight(i);
            VRall = VRall + sum((AMRF_est(i,:)-AMRF(i,:)).*(AMRF_est(i,:)-AMRF(i,:)))*weight(i);
        end
        VRall = 100*(1 - VRall/powerall);
        fprintf(strmVRf,'%-.3f %-.3f %.3f %-.6f\n',rank(m),Vr(n),VRall,M0_res);
    end
end
fclose(strmVRf);





