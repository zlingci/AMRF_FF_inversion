clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% input parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputf = 'LSQ_2D.txt';
strminputf = fopen(inputf,'r');
inputfiles = textscan(strminputf,'%s %*[^\n]',1);
inputfiles = char(string(inputfiles));
outputfiles = textscan(strminputf,'%s %*[^\n]',1);
outputfiles = char(string(outputfiles));
if exist(outputfiles,'dir') == 0 
    mkdir(outputfiles)
end
weightf = textscan(strminputf,'%s %*[^\n]',1);
weightf = char(string(weightf));
modelf = textscan(strminputf,'%s %*[^\n]',1);
modelf = char(string(modelf));
state_data = textscan(strminputf,'%f %*[^\n]',1);
state_data = cell2mat(state_data);
DT = textscan(strminputf,'%f %f %*[^\n]',1);
DT = cell2mat(DT);
dt = DT(1);
dts = DT(2);
duration = textscan(strminputf,'%f %*[^\n]',1);
duration = cell2mat(duration);
PLANE = textscan(strminputf,'%f %f %f %*[^\n]',1);
PLANE = cell2mat(PLANE);
strike = PLANE(1);
dip = PLANE(2);
rake = PLANE(3);
Vr = textscan(strminputf,'%f %*[^\n]',1);
Vr = cell2mat(Vr);
Moment = textscan(strminputf,'%f %f %*[^\n]',1);
Moment = cell2mat(Moment);
Ly_0 = textscan(strminputf,'%f %*[^\n]',1);
Ly_0 = cell2mat(Ly_0);
GRID_X = textscan(strminputf,'%f %f %f %*[^\n]',1);
GRID_X = cell2mat(GRID_X);
Lx_s = GRID_X(1);
Lx_e = GRID_X(2);
dx = GRID_X(3);
GRID_Y = textscan(strminputf,'%f %f %f %*[^\n]',1);
GRID_Y = cell2mat(GRID_Y);
Ly_s = GRID_Y(1);
Ly_e = GRID_Y(2);
dy = GRID_Y(3);
TRIANGLE = textscan(strminputf,'%f %f %*[^\n]',1);
TRIANGLE = cell2mat(TRIANGLE);
tri_num = TRIANGLE(1);
point_num = TRIANGLE(2);
damp = textscan(strminputf,'%f %f %f %*[^\n]',1);
damp = cell2mat(damp);
freq_cut = textscan(strminputf,'%f %f %*[^\n]',1);
freq_cut = cell2mat(freq_cut);
abs_dir = textscan(strminputf,'%s %*[^\n]',1);
abs_dir = char(string(abs_dir));
fclose(strminputf);
Duration_level = 0.99;
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
Mw = 2/3*log10(M01*10^(M02+7)) - 10.73;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting normalized STFs %%%%%%%%%%%%%%%%%%%%
% weightf = strcat(inputfiles,'/',weightf);
[egfname,sname,cha,dist,az,weight] = textread(weightf,'%s %s %s %f %f %*f %f','headerlines',1);
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
%%%%%%%%%%%%%%%%%%%%%%%% Constructing model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx_num = round((Lx_e - Lx_s)/dx) + 1;
Ly_num = round((Ly_e - Ly_s)/dy) + 1;
L_num = Lx_num*Ly_num;
Lx = zeros(Lx_num,1);
for i=1:1:Lx_num
    Lx(i) = Lx_s + (i-1)*dx;
end
Ly = zeros(Ly_num,1);
for i=1:1:Ly_num
    Ly(i) = Ly_s + (i-1)*dy;
end
L = zeros(L_num,3);
Dep = zeros(L_num,1);
for i=1:1:Lx_num
    for j=1:1:Ly_num
        L(j+(i-1)*Ly_num,1) = Lx(i);
        L(j+(i-1)*Ly_num,2) = Ly(j);
        Dep(j+(i-1)*Ly_num,1) = Ly(j)*cos((90-dip)*pi/180) + Ly_0;
        L(j+(i-1)*Ly_num,3) = sqrt(Lx(i)^2 + Ly(j)^2);
    end
end
if point_num == 0
     ts_num = fix(0.3*t(end)/dts);
     point_num = round((ts_num-1)/(tri_num + 1))*2 + 1;
     ts_num = (point_num-1)/2*(tri_num + 1) + 1;
else
    Ts_max = (point_num - 1)*dts;
    ts_max =  Ts_max/2*(tri_num + 1);
    ts_num = round(ts_max/dts) + 1;
end
%%%%%%%%%%%%%%%%%%%%%%% Getting dirctivity parameter %%%%%%%%%%%%%%%%%%%%%
modelf1 = strcat(abs_dir,'/inp/',modelf,'.txt');
[depth_m,Vp_m,Vs_m,density_m] = textread(modelf1,'%f %f %f %f %*f %*f');
Vp = zeros(L_num,1);
Vs = zeros(L_num,1);
density = zeros(L_num,1);
for i=1:1:L_num
    for j=1:1:length(depth_m)
        if Dep(i) <= depth_m(j)
            Vp(i) = Vp_m(j);
            Vs(i) = Vs_m(j);
            density(i) = density_m(j);
            break;
        end
    end
end
tempf1 = strcat(outputfiles,'/dep.txt');
tempf2 = strcat(outputfiles,'/dist.txt');
tempf3 = strcat(outputfiles,'/takeoff_p.txt');
tempf4 = strcat(outputfiles,'/takeoff_s.txt');
strmtempf1 = fopen(tempf1,'w');
fprintf(strmtempf1,'%s\n',strcat(inputfiles,'/',modelf,'.nd'));
% for i=1:1:Ly_num
%     fprintf(strmtempf1,'%-.3f\n',Dep(i));
% end
fprintf(strmtempf1,'%-.3f\n',Ly_0);
fclose(strmtempf1);
strmtempf2 = fopen(tempf2,'w');
for i=1:1:F_num
    fprintf(strmtempf2,'%-.3f\n',dist(i));
end
system('../bin/gettakeoff.sh');
takeoff_p = load(tempf3);
takeoff_s = load(tempf4);
ph = zeros(F_num,L_num);
pv = zeros(F_num,L_num);
for i=1:1:F_num
    for j=1:1:L_num
        if strcmp(char(cha(i)),'zp') || strcmp(char(cha(i)),'rp')
            ph(i,j) = sin(takeoff_p(i,1)*pi/180)/Vp(j);
            pv(i,j) = cos(takeoff_p(i,1)*pi/180)/Vp(j);
        else
            ph(i,j) = sin(takeoff_s(i,1)*pi/180)/Vs(j);
            pv(i,j) = cos(takeoff_s(i,1)*pi/180)/Vs(j);
        end
    end
end
% for i=1:1:F_num
%     for j=1:1:L_num
%         index_temp = mod(j,Ly_num);
%         if index_temp == 0
%             index_temp = Ly_num;
%         end
%         if strcmp(char(cha(i)),'zp') || strcmp(char(cha(i)),'rp')
%             ph(i,j) = sin(takeoff_p(i,index_temp)*pi/180)/Vp(j);
%             pv(i,j) = cos(takeoff_p(i,index_temp)*pi/180)/Vp(j);
%         else
%             ph(i,j) = sin(takeoff_s(i,index_temp)*pi/180)/Vs(j);
%             pv(i,j) = cos(takeoff_s(i,index_temp)*pi/180)/Vs(j);
%         end
%     end
% end
F1 = zeros(F_num,L_num);
F2 = zeros(F_num,L_num);
for i=1:1:F_num
    for j=1:1:L_num
        F1(i,j) = -ph(i,j)*cos((strike-az(i))*pi/180);
        F2(i,j) = ph(i,j)*sin((strike-az(i))*pi/180)*sin((90-dip)*pi/180) - pv(i,j)*cos((90-dip)*pi/180);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D Inversion on the fault plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AMRF,G,Gw,data] = LSQ_G_2D(STFs,F1,F2,L,dx,dy,Lx_num,Ly_num,t_num,ts_num,dt,dts,Vr,tri_num,weight,damp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NNLSQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
G1 = sqrt(Gw)*G;
data1 = sqrt(Gw)*data;
[model,resnorm,residual,exitflag,output,lambda] = lsqnonneg(G1,data1,optimset('Display','notify'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Synthetic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AMRF_est,MRFd_est]= RT_matrix_2D(model,F1,F2,L,t_num,ts_num,G,tri_num);
AMRF_est_normal = RT_stack_normal_2D(MRFd_est,L,dx,dy,dt,dts,Vr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% slip displacement ands monent of AMRF %%%%%%%%%%%%%%%
mu = density.*(Vs.^2)*10^9;
lamba = density.*(Vp.^2)*10^9 - 2*mu;
D0_est = zeros(L_num,5);
for i=1:1:L_num
    D0_est(i,1) = L(i,1);
    D0_est(i,2) = L(i,2);
    D0_est(i,3) = sum(MRFd_est(i,:))*dts*10^(M02-6)/mu(i);
    D0_est(i,4) = sum(MRFd_est(i,:))*dts*10^M02*dx*dy;
    D0_est(i,5) = rake;
end
M0_est = sum(AMRF_est_normal)*dt*10^M02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ave mu and lamba, density, Vs, Vp,%%%%%%%%%%%%%%%%%%%%%%%%
mu0 = sum(D0_est(:,4).*mu(:))/sum(D0_est(:,4));
lamba0 = sum(D0_est(:,4).*lamba(:))/sum(D0_est(:,4));
density0 = sum(D0_est(:,4).*density(:))/sum(D0_est(:,4));
Vs0 = sum(D0_est(:,4).*Vs(:))/sum(D0_est(:,4));
Vp0 = sum(D0_est(:,4).*Vp(:))/sum(D0_est(:,4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Centroid slip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dep_centroid = sum(D0_est(:,3).*Dep(:))/sum(D0_est(:,3));
Lx_centroid = sum(D0_est(:,3).*D0_est(:,1))/sum(D0_est(:,3));
Ly_centroid = sum(D0_est(:,3).*D0_est(:,2))/sum(D0_est(:,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Centroid time,Duration,Skewness,Kurtosis %%%%%%%%%%%%%%%%%%%%%
time = 0:dt:dt*(length(AMRF_est_normal)-1);
time_centroid = sum(AMRF_est_normal(:).*time(:))/sum(AMRF_est_normal);
for i=1:1:length(AMRF_est_normal)
    if sum(AMRF_est_normal(1:i))/sum(AMRF_est_normal) > Duration_level
        Duration = time(i);
        break;
    end
end
Cmoment_2 = sum(((time(:) - time_centroid).^2).*AMRF_est_normal(:))/sum(AMRF_est_normal); 
Cmoment_3 = sum(((time(:) - time_centroid).^3).*AMRF_est_normal(:))/sum(AMRF_est_normal);
Cmoment_4 = sum(((time(:) - time_centroid).^4).*AMRF_est_normal(:))/sum(AMRF_est_normal);
Skew = Cmoment_3/(Cmoment_2)^(3/2);
Kurtosis = Cmoment_4/(Cmoment_2)^2 - 3;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRFd_est, AMRF_est_normal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AMRF_est, moment_MRF_est %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AMRF, moment_AMRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fitting_variance reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D0_est %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D0_estf = strcat(outputfiles,'/D0_est.txt');
strmD0_estf = fopen(D0_estf,'w'); 
fprintf(strmD0_estf,'%s %s\n','hypocenter','grid_number');
fprintf(strmD0_estf,'%-.3f %d\n',Ly_0,L_num);
fprintf(strmD0_estf,'%s %s %s %s %s\n','X','Y','M0','Rake','slip_m');
for i=1:1:L_num
    fprintf(strmD0_estf,'%-.3f %-.3f %-.3e %-.2f %-.3e\n',D0_est(i,1),D0_est(i,2),D0_est(i,4),D0_est(i,5),D0_est(i,3));
end
fclose(strmD0_estf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRFd_est %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MRFd_estf = strcat(outputfiles,'/time_slip.txt');
strmMRFd_estf = fopen(MRFd_estf,'w'); 
fprintf(strmMRFd_estf,'%s %s %s %s %s %s %s %s %s\n','X','Y','Depth','rake','slip','mo','starttime','length','dt');
fprintf(strmMRFd_estf,'%s %s %s\n','time','slip','mo');
fprintf(strmMRFd_estf,'%d\n',L_num);
for i=1:1:L_num
    fprintf(strmMRFd_estf,'%s %-.3f %-.3f %-.3f %-.2f %-.3e %-.3e %-.3f %d %-.3f\n','#',D0_est(i,1),D0_est(i,2),Dep(i),D0_est(i,5),D0_est(i,3),D0_est(i,4),L(i,3)/Vr,length(MRFd_est(1,:)),dts);
    for j=1:1:length(MRFd_est(1,:)) 
        fprintf(strmMRFd_estf,'%-.3f %-.3e %-.3e\n',dts*(j-1), MRFd_est(i,j)*10^(M02-6)/mu(i),MRFd_est(i,j)*10^M02*dx*dy);
    end
end
fclose(strmMRFd_estf);
%%%%%%%%%%%%%%%%%%%% AMRF_est_normal, AMRF, and AMRF_est %%%%%%%%%%%%%%%%%
AMRF_est_normalf = strcat(outputfiles,'/AMRF_est_normal.txt');
strmAMRF_est_normalf = fopen(AMRF_est_normalf,'w'); 
for i=1:1:length(AMRF_est_normal)
    fprintf(strmAMRF_est_normalf,'%-.3f %-.3e\n',(i-1)*dt,AMRF_est_normal(i)*10^M02);
end
fclose(strmAMRF_est_normalf);
AMRFf = strcat(outputfiles,'/AMRF.txt');
strmAMRFf = fopen(AMRFf,'w'); 
fprintf(strmAMRFf,'%s %s %s %s %s\n','dist','azimuth','mo','length','dt');
fprintf(strmAMRFf,'%s %s\n','time','mo');
fprintf(strmAMRFf,'%d\n',F_num);
for i=1:1:F_num
    fprintf(strmAMRFf,'%s %-.3f %-.3f %-.3e %d %-.3f\n','#',dist(i),az(i),M0,length(AMRF(1,:)),dt);
    for j=1:1:length(AMRF(i,:))
        fprintf(strmAMRFf,'%-.3f %-.3e\n',(j-1)*dt,AMRF(i,j)*10^M02);
    end
end
fclose(strmAMRFf);
AMRF_estf = strcat(outputfiles,'/AMRF_est.txt');
strmAMRF_estf = fopen(AMRF_estf,'w'); 
fprintf(strmAMRF_estf,'%s %s %s %s %s\n','dist','azimuth','mo','length','dt');
fprintf(strmAMRF_estf,'%s %s\n','time','mo');
fprintf(strmAMRF_estf,'%d\n',F_num);
for i=1:1:F_num
    fprintf(strmAMRF_estf,'%s %-.3f %-.3f %-.3e %d %-.3f\n','#',dist(i),az(i),M0_est,length(AMRF_est(1,:)),dt);
    for j=1:1:length(AMRF_est(i,:))
        fprintf(strmAMRF_estf,'%-.3f %-.3e\n',(j-1)*dt,AMRF_est(i,j)*10^M02);
    end
end
fclose(strmAMRF_estf);
%%%%%%%%%%%%%% moment_AMRF and moment_AMRF_est %%%%%%%%%%%%%%%%%%%%%%%%%%%%
moment_AMRFf = strcat(outputfiles,'/moment_AMRF.txt');
strmmoment_AMRFf = fopen(moment_AMRFf,'w'); 
fprintf(strmmoment_AMRFf,'%s %s %s %s %s\n','dist','azimuth','mo','length','dt');
fprintf(strmmoment_AMRFf,'%s %s\n','time','mo');
fprintf(strmmoment_AMRFf,'%d\n',F_num);
for i=1:1:F_num
    fprintf(strmmoment_AMRFf,'%s %-.3f %-.3f %-.3e %d %-.3f\n','#',dist(i),az(i),M0,length(moment_AMRF(1,:)),dt);
    for j=1:1:length(moment_AMRF(i,:))
        fprintf(strmmoment_AMRFf,'%-.3f %-.3e\n',(j-1)*dt,moment_AMRF(i,j)*10^M02);
    end
end
fclose(strmmoment_AMRFf);
moment_AMRF_estf = strcat(outputfiles,'/moment_AMRF_est.txt');
strmmoment_AMRF_estf = fopen(moment_AMRF_estf,'w'); 
fprintf(strmmoment_AMRF_estf,'%s %s %s %s %s\n','dist','azimuth','mo','length','dt');
fprintf(strmmoment_AMRF_estf,'%s %s\n','time','mo');
fprintf(strmmoment_AMRF_estf,'%d\n',F_num);
for i=1:1:F_num
    fprintf(strmmoment_AMRF_estf,'%s %-.3f %-.3f %-.3e %d %-.3f\n','#',dist(i),az(i),M0_est,length(moment_AMRF_est(1,:)),dt);
    for j=1:1:length(moment_AMRF_est(i,:))
        fprintf(strmmoment_AMRF_estf,'%-.3f %-.3e\n',(j-1)*dt,moment_AMRF_est(i,j)*10^M02);
    end
end
fclose(strmmoment_AMRF_estf);
%%%%%%%%%%%%%%%%%%%%%%%% Fitting_variance reduction %%%%%%%%%%%%%%%%%%%%%%%
VRf = strcat(outputfiles,'/VR.txt');
strmVRf = fopen(VRf,'w'); 
for i=1:1:length(az)
    fprintf(strmVRf,'%s %s %s %-.3f %-.3f %-.3f %-.3f\n',char(egfname(i)),char(sname(i)),char(cha(i)),az(i),weight(i),VR(i),F1(i,1));
end
fclose(strmVRf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% other parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameterf = strcat(outputfiles,'/parameter.txt');
strmparameterf = fopen(parameterf,'w');
fprintf(strmparameterf,'%s %-.3f\n','Variance_reduction',VRall);
fprintf(strmparameterf,'%s %-.3e\n','Moment',M0);
fprintf(strmparameterf,'%s %-.3e\n','Moment',M0_est);
fprintf(strmparameterf,'%s %-.3e\n','Moment_scale',10^M02);
fprintf(strmparameterf,'%s %-.3f\n','Centroid_depth',Dep_centroid);
fprintf(strmparameterf,'%s %-.3f\n','Centroid_Lx',Lx_centroid);
fprintf(strmparameterf,'%s %-.3f\n','Centroid_Ly',Ly_centroid);
fprintf(strmparameterf,'%s %-.3f\n','Centroid_time',time_centroid);
fprintf(strmparameterf,'%s %-.3e\n','ave_mu',mu0);
fprintf(strmparameterf,'%s %-.3e\n','ave_lamba',lamba0);
fprintf(strmparameterf,'%s %-.4f\n','ave_Vp',Vp0);
fprintf(strmparameterf,'%s %-.4f\n','ave_Vs',Vs0);
fprintf(strmparameterf,'%s %-.4f\n','ave_density',density0);
fprintf(strmparameterf,'%s %-.3f\n','Duration',Duration);
fclose(strmparameterf);

