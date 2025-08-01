clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameter
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
inputfile = strcat(outputfiles,'/VR.txt');
oufputfile1 = strcat(outputfiles,'/ellipse.txt');
oufputfile2 = strcat(outputfiles,'/mean_var.txt');
oufputfile3 = strcat(outputfiles,'/count_drection.txt');
oufputfile4 = strcat(outputfiles,'/count_speed.txt');
oufputfile5 = strcat(outputfiles,'/VR_use.txt');
oufputfile6 = strcat(outputfiles,'/normal_speed.txt');
oufputfile7 = strcat(outputfiles,'/normal_direction.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shift = RANK(4);
data = load(inputfile);
X = data(:,1);
Y = data(:,2);
for i=1:1:length(X)
    if X(i) < shift 
        X(i) = X(i) + 360;
    end
end
VAR = data(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 0;
for i=1:1:length(X)
    if VAR(i)>= max(VAR) - 0.1*max(VAR)
        k = k + 1;
    end
end
X1 = zeros(k,1);        
Y1 = zeros(k,1);
x_index = shift:RANK(3):360+shift;
y_index = VR(1):VR(3):VR(2);
count_x = zeros(length(x_index),1);
count_y = zeros(length(y_index),1);
k = 0;
for i=1:1:length(X)
    if VAR(i) >= max(VAR) - 0.1*max(VAR)
        k = k + 1;
        X1(k) = X(i);
        Y1(k) = Y(i);
        for j=1:1:length(x_index)
            if abs(x_index(j) - X(i)) < 10^-3 
                count_x(j) = count_x(j) + 1;
            end
        end
        for j=1:1:length(y_index)
            if abs(y_index(j) - Y(i)) < 10^-3 
                count_y(j) = count_y(j) + 1;
            end
        end
    end
end
mea = [mean(X1);mean(Y1)];
covar = cov(X1,Y1);        
[vector,value] = eig(covar);
c = 2.447747; % 95% 
a = sqrt(value(1,1)*c);
b = sqrt(value(2,2)*c);
N = 100;
dx = 2*a/(N-1);
xa = zeros(2*N,1);
xb = zeros(2*N,1); 
for i=1:1:N
    xa(i) = -a + dx*(i-1);
    xb(i) = b/a*sqrt(a^2-xa(i)^2); 
end
for i=N+1:1:2*N
    xa(i) = a - dx*(i-N-1);
    xb(i) = -b/a*sqrt(a^2-xa(i)^2); 
end
Xa = zeros(2*N,1);
Xb = zeros(2*N,1); 
for i=1:1:2*N
    temp = vector*[xa(i);xb(i)];
    Xa(i) = temp(1)+mea(1);
    Xb(i) = temp(2)+mea(2);
end
strmoufputfile1 = fopen(oufputfile1,'w');
for i=1:1:2*N
    fprintf(strmoufputfile1,'%-.6f %-.6f\n',Xa(i),Xb(i));
end
fclose(strmoufputfile1);
strmoufputfile2 = fopen(oufputfile2,'w');
fprintf(strmoufputfile2,'%-.2f %-.2f %-.3f %-.3f\n',mea(1),mea(2),sqrt(value(2,2)),sqrt(value(1,1)));
fclose(strmoufputfile2);
strmoufputfile3 = fopen(oufputfile3,'w');
for i=1:1:length(x_index)
    fprintf(strmoufputfile3,'%-.3f %d\n',x_index(i),count_x(i));
end
fclose(strmoufputfile3);
strmoufputfile4 = fopen(oufputfile4,'w');
for i=1:1:length(y_index)
    fprintf(strmoufputfile4,'%-.3f %d\n',y_index(i),count_y(i));
end
fclose(strmoufputfile4);
strmoufputfile5 = fopen(oufputfile5,'w');
for i=1:1:length(X1)
    fprintf(strmoufputfile5,'%.1f %.1f\n',X1(i),Y1(i));
end
fclose(strmoufputfile5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = VR(1):0.01:VR(2);
b = shift:0.5:360+shift;
a1 = zeros(length(a),1);
b1 = zeros(length(b),1);
for i=1:1:length(a)
    a1(i) = 1/(sqrt(2*pi)*sqrt(value(1,1)))*exp(-(a(i)-mea(2))^2/2/value(1,1));
end
strmoufputfile6 = fopen(oufputfile6,'w');
for i=1:1:length(a)
    fprintf(strmoufputfile6,'%f %f\n',a(i),a1(i));
end
fclose(strmoufputfile6);
for i=1:1:length(b)
    b1(i) = 1/(sqrt(2*pi)*sqrt(value(2,2)))*exp(-(b(i)-mea(1))^2/2/value(2,2));
end
strmoufputfile7 = fopen(oufputfile7,'w');
for i=1:1:length(b)
    fprintf(strmoufputfile7,'%f %f\n',b(i),b1(i)*100);
end
fclose(strmoufputfile7);
% figure(1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = min(X)-5/2:5:max(X)+5/2;
% y = min(Y)-0.1/2:0.1:max(Y)+0.1/2;
% [X2,Y2] = meshgrid(x,y);
% Z2 = griddata(X-5/2,Y-0.1/2,res,X2,Y2);
% pcolor(X2,Y2,Z2);
% axis([min(X)-5/2 max(X)+5/2 min(Y)-0.1/2 max(Y)+0.1/2]); 
% colormap('jet');
% c = colorbar('FontSize',14);
% c.Label.String = 'Fitting Vriance Reduction';
% caxis([max(max(res))-5,max(max(res))]);
% hold on;
% plot(Xa,Xb,'r--','LineWidth',3);
% hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gca,'Fontname', 'Times New Roman','fontsize',16,'LineWidth',2);
% xlabel('Rupture direction','Fontname','Times New Roman','fontsize',20);
% ylabel('Rupture speed (km/s)','Fontname','Times New Roman','fontsize',20);
% set(gca,'Position',[0.10 0.13 0.78 0.80]);
% outfiguref = strcat('Direction_Speed_Variance_reduction');
% set(gcf,'unit','centimeters','position',[10 5 30 15]);
% saveas(gcf,outfiguref, 'jpg');
