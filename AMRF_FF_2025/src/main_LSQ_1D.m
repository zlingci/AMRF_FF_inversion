clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputf = strcat('LSQ_1D.txt');
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
Rupture = textscan(strminputf,'%f %f %*[^\n]',1);
Rupture = cell2mat(Rupture);
Vr = Rupture(1);
rank = Rupture(2);
Moment = textscan(strminputf,'%f %f %*[^\n]',1);
Moment = cell2mat(Moment);
R = textscan(strminputf,'%f %f %f %f %*[^\n]',1);
R = cell2mat(R);
dL = R(4);
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
F = zeros(F_num,1);
for i=1:1:F_num
    F(i) = -ph(i)*cos((strike-az(i))*pi/180)*cos(rank*pi/180) - ph(i)*sin((strike-az(i))*pi/180)*sin(rank*pi/180)*sin((90-dip)*pi/180) + pv(i)*sin(rank*pi/180)*cos((90-dip)*pi/180);
end
%%%%%%%%%%%%%%%%%%%%%%%% Constructing model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dts == 0
    dts = dt;
end
if R(1) == 0 && R(2) == 0
    L_e = t(end)/(1.2/Vr + max(abs(F)));
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
    ts_num = round(0.3*L_e/Vr/dts) + 1;
    ts = 0:dts:(ts_num-1)*dts;
else
    L_0 = R(3);
    L_s = R(1) + L_0;
    L_e = R(2) + L_0;
    L_num = ceil((L_e - L_0)/dL) + ceil((L_0 - L_s)/dL) + 1;
    L = zeros(L_num,1);
    for i=1:1:L_num
        L(i) = (i - 1 - ceil((L_0 - L_s)/dL))*dL;
    end
    if ts_max == 0
        ts_max = 0.3*t(end)/(1.2+max(F)*Vr);
    end
    ts_num = ceil(ts_max/dts) + 1;
    ts = 0:dts:dts*(ts_num - 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D Inversion on the fault plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Constructing data, and G %%%%%%%%%%%%%%%%%%%%%%% 
[AMRF,G,Gw,data] = LSQ_G_1D(STFs,F,L,dL,t_num,ts_num,dt,dts,Vr,weight,damp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NNLSQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
G1 = sqrt(Gw)*G;
data1 = sqrt(Gw)*data;
[model,resnorm,residual,exitflag,output,lambda] = lsqnonneg(G1,data1,optimset('Display','notify'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Synthetic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AMRF_est,MRFd_est] = RT_matrix_1D(model,F,L,t_num,ts_num,G);
AMRF_est_normal = RT_stack_normal_1D(MRFd_est,L,dL,dt,dts,Vr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% slip displacement of AMRF %%%%%%%%%%%%%%%%%%%%%
D0_est = zeros(L_num,2);
for i=1:1:L_num
    D0_est(i,1) = L(i);
    D0_est(i,2) = sum(MRFd_est(i,:))*dts/(density*Vs^2);
end
M0_est = sum(AMRF_est_normal)*dt*10^M02;
%%%%%%%%%%%%%%%%%%%%%%% Moment time histories of AMRF %%%%%%%%%%%%%%%%%%%%%
moment_AMRF = zeros(F_num,length(AMRF(1,:)));
for i=1:1:F_num
    for j=1:1:length(AMRF(1,:))
        if j==1
            moment_AMRF(i,j) =  AMRF(i,j)*dt;
        else
            moment_AMRF(i,j) =  moment_AMRF(i,j-1) + AMRF(i,j)*dt; 
        end
    end
end
%%%%%%%%%%%%%%%%%%%% Moment time histories of AMRF_est %%%%%%%%%%%%%%%%%%%%
moment_AMRF_est = zeros(F_num,length(AMRF_est(1,:)));
for i=1:1:F_num
    for j=1:1:length(AMRF_est(1,:))
        if j==1
            moment_AMRF_est(i,j) =  AMRF_est(i,j)*dt;
        else
            moment_AMRF_est(i,j) =  moment_AMRF_est(i,j-1) + AMRF_est(i,j)*dt;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% Fitting variance reduction %%%%%%%%%%%%%%%%%%%%%%%
VR = zeros(F_num,1);
powerall = 0;
VRall = 0;
for i=1:1:F_num
    power = sum(AMRF(i,:).*AMRF(i,:));
    VR(i) = 100*(1-sum((AMRF_est(i,:)-AMRF(i,:)).*(AMRF_est(i,:)-AMRF(i,:)))/power);
    powerall = powerall + power*weight(i);
    VRall =  VRall + sum((AMRF_est(i,:)-AMRF(i,:)).*(AMRF_est(i,:)-AMRF(i,:)))*weight(i);
end
VRall = 100*(1 - VRall/powerall);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRFd_est, AMRF_est_normal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AMRF_est, moment_MRF_est %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AMRF, moment_AMRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fitting_variance reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(outputfiles,'dir') == 0 
    mkdir(outputfiles)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D0_est %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D0_estf = strcat(outputfiles,'/D0_est.txt');
strmD0_estf = fopen(D0_estf,'w'); 
for i=1:1:L_num
    fprintf(strmD0_estf,'%-.6f %-.6f\n',D0_est(i,1),D0_est(i,2));
end
fclose(strmD0_estf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRFd_est %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MRFd_estf = strcat(outputfiles,'/MRFd_est.txt');
strmMRFd_estf = fopen(MRFd_estf,'w'); 
for i=1:1:length(MRFd_est(1,:))+2
    fprintf(strmMRFd_estf,'%-.3f ',dts*(i-2));
    for j=1:1:L_num
        if i==1 || i == length(MRFd_est(1,:))+2
            fprintf(strmMRFd_estf,'%-.6f ',0);
        else
            fprintf(strmMRFd_estf,'%-.6f ',MRFd_est(j,i-1));
        end
    end
    fprintf(strmMRFd_estf,'\n');
end
fclose(strmMRFd_estf);
%%%%%%%%%%%%%%%%%%%% AMRF_est_normal, AMRF, and AMRF_est %%%%%%%%%%%%%%%%%
AMRF_est_normalf = strcat(outputfiles,'/AMRF_est_normal.txt');
strmAMRF_est_normalf = fopen(AMRF_est_normalf,'w'); 
for j=1:1:length(AMRF_est_normal)
    fprintf(strmAMRF_est_normalf,'%-.3f %-.6f\n',dt*(j-1),AMRF_est_normal(j));
end
fclose(strmAMRF_est_normalf);
AMRFf = strcat(outputfiles,'/AMRF.txt');
strmAMRFf = fopen(AMRFf,'w'); 
for i=1:1:length(AMRF(1,:))
    fprintf(strmAMRFf,'%-.3f ',dt*(i-1));
    for j=1:1:F_num
        fprintf(strmAMRFf,'%-.6f ',AMRF(j,i));
    end
    fprintf(strmAMRFf,'\n');
end
fclose(strmAMRFf);
AMRF_estf = strcat(outputfiles,'/AMRF_est.txt');
strmAMRF_estf = fopen(AMRF_estf,'w'); 
for i=1:1:length(AMRF_est(1,:))
    fprintf(strmAMRF_estf,'%-.3f ',dt*(i-1));
    for j=1:1:F_num
        fprintf(strmAMRF_estf,'%-.6f ',AMRF_est(j,i));
    end
    fprintf(strmAMRF_estf,'\n');
end
fclose(strmAMRF_estf);

%%%%%%%%%%%%%% moment_AMRF and moment_AMRF_est %%%%%%%%%%%%%%%%%%%%%%%%%%%%
moment_AMRFf = strcat(outputfiles,'/moment_AMRF.txt');
strmmoment_AMRFf = fopen(moment_AMRFf,'w'); 
for i=1:1:length(moment_AMRF(1,:))
    fprintf(strmmoment_AMRFf,'%-.3f ',dt*(i-1));
    for j=1:1:F_num
        fprintf(strmmoment_AMRFf,'%-.6f ',moment_AMRF(j,i));
    end
    fprintf(strmmoment_AMRFf,'\n');
end
fclose(strmmoment_AMRFf);
moment_AMRF_estf = strcat(outputfiles,'/moment_AMRF_est.txt');
strmmoment_AMRF_estf = fopen(moment_AMRF_estf,'w'); 
for i=1:1:length(moment_AMRF_est(1,:))
    fprintf(strmmoment_AMRF_estf,'%-.3f ',dt*(i-1));
    for j=1:1:F_num
        fprintf(strmmoment_AMRF_estf,'%-.6f ',moment_AMRF_est(j,i));
    end
    fprintf(strmmoment_AMRF_estf,'\n');
end
fclose(strmmoment_AMRF_estf);
%%%%%%%%%%%%%%%%%%%%%%%% Fitting_variance reduction %%%%%%%%%%%%%%%%%%%%%%%
VRf = strcat(outputfiles,'/VR.txt');
strmVRf = fopen(VRf,'w'); 
for i=1:1:length(az)
    fprintf(strmVRf,'%s %s %s %-.3f %-.3f %-.3f %-.3f\n',char(egfname(i)),char(sname(i)),char(cha(i)),az(i),takeoff(i),weight(i),VR(i));
end
fclose(strmVRf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% other parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameterf = strcat(outputfiles,'/parameter.txt');
strmparameterf = fopen(parameterf,'w');
fprintf(strmparameterf,'%s %-.3f\n','Variance_reduction',VRall);
fprintf(strmparameterf,'%s %-.3e\n','Moment',M0);
fprintf(strmparameterf,'%s %-.3e\n','Moment_est',M0_est);
fprintf(strmparameterf,'%s %-.3e\n','Moment_scale',10^M02);
fclose(strmparameterf);


