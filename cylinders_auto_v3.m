clear all; close all; clc

number = 4;

%number = 1 -> Model: CC. Fit: 1.
%number = 2 -> Model: EC. Fit: 1.
%number = 3 -> Model: CC. Fit: 2.
%number = 4 -> Model: EC. Fit: 2.

phi_vec_sta =   [10.79 46.54 157.83 141.86]-90;
theta_vec_sta = [-40.80 -50.64 -54.13 -47.72];
psi_vec_sta = [90 104.18 90 90];
radio_vec_sta = [0.4234 0.7753 0.1368 0.1638];
delta_sta = [1 0.6874 1 0.6937];
y0_vec_sta = [0.3987 0.6506 -0.0931 -0.0908];

phi_vec_ms =    [93.37 77.67 156.03 130.36]-90;
theta_vec_ms =  [-42.19 -39.62 9.24 23.73];
psi_vec_ms = [90 90 90 53.0812];
radio_vec_ms = [0.0851 0.0830 0.0436 0.0610];
delta_ms = [1 1.3661 1 0.7499];
y0_vec_ms = [0.0454 0.059 -0.0403 -0.0403];

n1 = 400;
h1 = 0.2;
n2 = 200;
h2 = 0.1;

phi_cc_sta = phi_vec_sta(number); %STA
theta_cc_sta = theta_vec_sta(number); %STA
psi_cc_sta = psi_vec_sta(number);

phi_cc_ms = phi_vec_ms(number); %MESSENGER
theta_cc_ms = theta_vec_ms(number); %MESSENGER
psi_cc_ms = psi_vec_ms(number);


%Reading files
C0 = csvread('275_C0.csv');
C1 = csvread('275_C1.csv');
C2 = csvread('275_C2.csv');

%Ejes globales
figure('position',[10 0 800 1000])
%{
factor = 0.5;
grosor = 1.5;
plot3([-1*factor 0],[0 0],[0 0],'k','LineWidth',grosor); hold on
xGSE_text = text(-1*factor,0,0,'\bf x^{GSE}','HorizontalAlignment','left','FontSize',15);
set(xGSE_text,'Color','k')
plot3([0 0],[-1*factor 0],[0 0],'k','LineWidth',grosor); hold on
yGSE_text = text(0,-1*factor,0,'\bf y^{GSE}','HorizontalAlignment','left','FontSize',15);
set(yGSE_text,'Color','k')
plot3([0 0],[0 0],[0 1*factor],'k','LineWidth',grosor); hold on
zGSE_text = text(0,0,1*factor,'\bf z^{GSE}','HorizontalAlignment','left','FontSize',15);
set(zGSE_text,'Color','k')
%}
%Rotation Angles
ang_c = pi/180;
Long = (156.759-180)*ang_c; %-8.841 -(180-180.27)
Lati = -2*ang_c; 
Tilt = 36.335*ang_c; %59.814

%Rotation Matrix
Long_Mat = @(x)([1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)]);
Lati_Mat = @(x)([cos(x) 0 sin(x); 0 1 0; -sin(x) 0 cos(x)]);
Tilt_Mat = @(x)([cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1]);
R_model = @(phi,theta,psi) ([cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) -cos(theta)*sin(phi) cos(phi)*sin(psi)+cos(psi)*sin(theta)*sin(phi);
                     cos(psi)*sin(phi)+cos(phi)*sin(theta)*sin(psi) cos(phi)*cos(theta) sin(phi)*sin(psi)-cos(psi)*cos(phi)*sin(theta);
                    -cos(theta)*sin(psi) sin(theta) cos(theta)*cos(psi)]);


Long_Rot = Long_Mat(Long);
Lati_Rot = Lati_Mat(Lati);
Tilt_Rot = Tilt_Mat(Tilt);

%Cloud points array
conv = 695700/1.496e8;
Nube = Tilt_Rot*Lati_Rot*Long_Rot*[C0,C1,C2]'*conv;
Nube = Nube';

%CME direction
dir_vec = [-1.4 0 0]';
Long_Dir = Tilt_Mat(Long);
Lati_Dir = Lati_Mat(Lati);
Tilt_Dir = Long_Mat(-Tilt);

dir_vec = Tilt_Dir*Lati_Dir*Long_Dir*dir_vec; %
plot3([0 dir_vec(1)],[0 dir_vec(2)],[0 dir_vec(3)],'k','LineWidth',1.5); hold on
%plot3(dir_vec(1),dir_vec(2),dir_vec(3),'Marker','<','MarkerEdgeColor','k','MarkerFaceColor','k'); hold on

R = [0.02 0];
n = 20;
cyl_color = 'k';
closed = 1;
lines = 0;
fact0 = 1.025;

dir_cme1 = [dir_vec(1) dir_vec(2) dir_vec(3)];
dir_cme2 = [dir_vec(1) dir_vec(2) dir_vec(3)]*fact0;
[Cone0,EndPlate01,EndPlate02] = Cone(dir_cme1,dir_cme2,R,n,cyl_color,closed,lines);

%Spacecraft positions
sta_x = -0.789586;
sta_y = 0.558722;
sta_z = -0.00207215;

stb_x = -0.761043;
stb_y = -0.686539;
stb_z = -0.00299437;

ms_x = -0.3149;
ms_y = 0.0983;
ms_z = 0.0398;

%3D scatter plot
p_size = 100;

nube_x = -Nube(:,3);
nube_y = Nube(:,2);
nube_z = Nube(:,1);

plot3(nube_x,nube_y,nube_z,'-','Color',[120,120,120]/255,'LineWidth',0.25); hold on; grid on
plot3(nube_x,nube_y,nube_z,'.k','LineWidth',0.25); hold on; grid on
scatter3(0,0,0,400,'y','o','MarkerFaceColor','y'); hold on; grid on
scatter3(sta_x,sta_y,sta_z,p_size,'r','o','MarkerFaceColor','r'); hold on
scatter3(stb_x,stb_y,stb_z,p_size,'b','o','MarkerFaceColor','b'); hold on
scatter3(ms_x,ms_y,ms_z,p_size,'m','o','MarkerFaceColor','m'); hold on
scatter3(1,0,0,p_size,'g','o','MarkerFaceColor','g'); hold on

%Text 
%{
sta_text = text(sta_x-0.1,sta_y+0.1,sta_z+0.1,'\bf STA','HorizontalAlignment','left','FontSize',15);
set(sta_text,'Color','k')
stb_text = text(stb_x-0.1,stb_y+0.1,stb_z+0.1,'\bf STB','HorizontalAlignment','left','FontSize',15);
set(stb_text,'Color','k')
ms_text = text(ms_x-0.1,ms_y+0.1,ms_z+0.1,'\bf MS','HorizontalAlignment','left','FontSize',15);
set(ms_text,'Color','k')
%}
sun_text = text(0.1,0.1,0.1,'\bf Sun','HorizontalAlignment','left','FontSize',15);
set(sun_text,'Color','k')
earth_text = text(1.1,0.1,0.1,'\bf Earth','HorizontalAlignment','left','FontSize',15);
set(earth_text,'Color','k')


mult_sta = 1.45;
mult_ms = 4.3;

plot3([0 sta_x*mult_sta],[0 sta_y*mult_sta],[0 sta_z*mult_sta],'r','LineWidth',1.5); hold on
plot3([0 stb_x],[0 stb_y],[0 stb_z],'b','LineWidth',1.5); hold on
plot3([0 mult_ms*ms_x],[0 mult_ms*ms_y],[0 mult_ms*ms_z],'m','LineWidth',1.5); hold on
plot3([0 1],[0 0],[0 0],'g','LineWidth',1.5); hold on
%{
R = [0.02 0];
n = 20;
cyl_color = 'r';
closed = 1;
lines = 0;
fact0 = 1.025;

dir_cme1 = [dir_vec(1) dir_vec(2) dir_vec(3)];
dir_cme2 = [dir_vec(1) dir_vec(2) dir_vec(3)]*fact0;
[Cone0,EndPlate01,EndPlate02] = Cone(dir_cme1,dir_cme2,R,n,cyl_color,closed,lines);
%}

xlabel('\fontsize{14} x axis (AU)')
ylabel('\fontsize{14} y axis (AU)')
zlabel('\fontsize{14} z axis (AU)')
%{
xlim([-1.5 0.5])
ylim([-1 1])
zlim([-1 1])
%}
xlim([-1.5 1])
ylim([-1 1.5])
zlim([-1 1.5])

% Tal y como esta ahora:
% 1 ROTA THETA. 2 ROTA PHI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%-------- CC Fit 1 STEREO A --------%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1 = radio_vec_sta(number); %0.1368;

a_sta=delta_sta(number)*radio_vec_sta(number); % horizontal radius
b_sta=radio_vec_sta(number); % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
t=linspace(-pi,pi,n1);
x_sta=x0+a_sta*cos(t);
y_sta = ones(1,length(x_sta))*h1/2;
y2_sta = h1/2;
z_sta=y0+b_sta*sin(t);
xcc_sta = [x_sta,x_sta];
ycc_sta = zeros(1,n1*2);
for m = 1:n1*2
    if mod(m,2) == 1
        ycc_sta(m) = y2_sta;
    else
        ycc_sta(m) = -y2_sta;
    end
end
zcc_sta = [z_sta,z_sta];
cc_sta1 = [xcc_sta;ycc_sta;zcc_sta];

fact_sta = 0.45;
cc_sta_axis = [0 1 0;0 -1 0]'*fact_sta;

%Apply Rotations:
local_axis = [0 0 0; 1 0 0; 0 1 0; 0 0 1]';
local_axis = R_model((phi_cc_sta-(90-atan(-sta_x/sta_y)/ang_c))*ang_c,theta_cc_sta*ang_c,psi_cc_sta*ang_c)*local_axis;
cc_sta2 = R_model((phi_cc_sta-(90-atan(-sta_x/sta_y)/ang_c))*ang_c,theta_cc_sta*ang_c,psi_cc_sta*ang_c)*cc_sta1; 
cc_sta_axis = R_model((phi_cc_sta-(90-atan(-sta_x/sta_y)/ang_c))*ang_c,theta_cc_sta*ang_c,psi_cc_sta*ang_c)*cc_sta_axis;%

%view([0,0,90])

plot3(cc_sta2(1,:)+sta_x,cc_sta2(2,:)+sta_y,cc_sta2(3,:)+sta_z+y0_vec_sta(number),'-r'); hold on
plot3(cc_sta2(1,:)+sta_x,cc_sta2(2,:)+sta_y,cc_sta2(3,:)+sta_z+y0_vec_sta(number),'.k','LineWidth',0.5); hold on
plot3([cc_sta_axis(1,1) cc_sta_axis(1,2)]+sta_x,[cc_sta_axis(2,1) cc_sta_axis(2,2)]+sta_y,[cc_sta_axis(3,1) cc_sta_axis(3,2)]+sta_z+y0_vec_sta(number),'--k','LineWidth',2.5); hold on

%Axis
factor = 0.3;
grosor = 1.5;
aux_axis_sta = [-1 0 0; 0 -1 0; 0 0 1];
axis_sta = Tilt_Mat((270-atan(-sta_x/sta_y)/ang_c)*ang_c)*aux_axis_sta;%

plot3([axis_sta(1,1)*factor 0]+sta_x,[axis_sta(1,2)*factor 0]+sta_y,[axis_sta(1,3)*factor 0]+sta_z+0.001,'k','LineWidth',grosor); hold on
plot3([axis_sta(2,1)*factor 0]+sta_x,[axis_sta(2,2)*factor 0]+sta_y,[axis_sta(2,3)*factor 0]+sta_z,'k','LineWidth',grosor); hold on
plot3([axis_sta(3,1)*factor 0]+sta_x,[axis_sta(3,2)*factor 0]+sta_y,[axis_sta(3,3)*factor 0]+sta_z,'k','LineWidth',grosor); hold on
%{
plot3([local_axis(1,1) local_axis(1,2)]+sta_x,[local_axis(2,1) local_axis(2,2)]+sta_y,[local_axis(3,1) local_axis(3,2)]+sta_z,'-c','LineWidth',2); hold on
plot3([local_axis(1,1) local_axis(1,3)]+sta_x,[local_axis(2,1) local_axis(2,3)]+sta_y,[local_axis(3,1) local_axis(3,3)]+sta_z,'-m','LineWidth',2); hold on
plot3([local_axis(1,1) local_axis(1,4)]+sta_x,[local_axis(2,1) local_axis(2,4)]+sta_y,[local_axis(3,1) local_axis(3,4)]+sta_z,'-k','LineWidth',2); hold on
%}
R = [0.02 0];
n = 20;
cyl_color = 'k';
closed = 1;
lines = 0;
fact1 = 0.3;
fact2 = 0.34;

x_sta1 = [axis_sta(1,1)*fact1+sta_x axis_sta(1,2)*fact1+sta_y axis_sta(1,3)*fact1+sta_z];
x_sta2 = [axis_sta(1,1)*fact2+sta_x axis_sta(1,2)*fact2+sta_y axis_sta(1,3)*fact2+sta_z];
[Cone4,EndPlate7,EndPlate8] = Cone(x_sta1,x_sta2,R,n,cyl_color,closed,lines);
y_sta1 = [axis_sta(2,1)*fact1+sta_x axis_sta(2,2)*fact1+sta_y axis_sta(2,3)*fact1+sta_z];
y_sta2 = [axis_sta(2,1)*fact2+sta_x axis_sta(2,2)*fact2+sta_y axis_sta(2,3)*fact2+sta_z];
[Cone5,EndPlate9,EndPlate10] = Cone(y_sta1,y_sta2,R,n,cyl_color,closed,lines);
z_sta1 = [axis_sta(3,1)*fact1+sta_x axis_sta(3,2)*fact1+sta_y axis_sta(3,3)*fact1+sta_z];
z_sta2 = [axis_sta(3,1)*fact2+sta_x axis_sta(3,2)*fact2+sta_y axis_sta(3,3)*fact2+sta_z];
[Cone6,EndPlate11,EndPlate12] = Cone(z_sta1,z_sta2,R,n,cyl_color,closed,lines);

%plot3(axis_sta(1,1)*factor+sta_x,axis_sta(1,2)*factor+sta_y,axis_sta(1,3)*factor+sta_z,'Marker','>','MarkerEdgeColor','k','MarkerFaceColor','k'); hold on
%plot3(axis_sta(2,1)*factor+sta_x,axis_sta(2,2)*factor+sta_y,axis_sta(2,3)*factor+sta_z,'Marker','<','MarkerEdgeColor','k','MarkerFaceColor','k'); hold on
%plot3(axis_sta(3,1)*factor+sta_x,axis_sta(3,2)*factor+sta_y,axis_sta(3,3)*factor+sta_z,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','k'); hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%-------- CC Fit 1 MESSENGER --------%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r2 = radio_vec_ms(number); %0.1368;

a_ms=delta_ms(number)*radio_vec_ms(number); % horizontal radius
b_ms=radio_vec_ms(number); % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
t=linspace(-pi,pi,n2);
x_ms=x0+a_ms*cos(t);
y_ms = ones(1,length(x_ms))*h2/2;
y2_ms = h2/2;
z_ms=y0+b_ms*sin(t);
xcc_ms = [x_ms,x_ms];
ycc_ms = zeros(1,n2*2);
for m = 1:n2*2
    if mod(m,2) == 1
        ycc_ms(m) = y2_ms;
    else
        ycc_ms(m) = -y2_ms;
    end
end
zcc_ms = [z_ms,z_ms];
cc_ms = [xcc_ms;ycc_ms;zcc_ms];

fact_ms = 0.25;
cc_ms_axis = [0 1 0;0 -1 0]'*fact_ms;

%Apply Rotations:

local_axis = [0 0 0; 1 0 0; 0 1 0; 0 0 1]';
local_axis = R_model((phi_cc_ms-(90-atan(-ms_x/ms_y)/ang_c))*ang_c,theta_cc_ms*ang_c,psi_cc_ms*ang_c)*local_axis;
cc_ms = R_model((phi_cc_ms-(90-atan(-ms_x/ms_y)/ang_c))*ang_c,theta_cc_ms*ang_c,psi_cc_ms*ang_c)*cc_ms; 
cc_ms_axis = R_model((phi_cc_ms-(90-atan(-ms_x/ms_y)/ang_c))*ang_c,theta_cc_ms*ang_c,psi_cc_ms*ang_c)*cc_ms_axis;%

plot3(cc_ms(1,:)+ms_x,cc_ms(2,:)+ms_y,cc_ms(3,:)+ms_z+y0_vec_ms(number),'-m'); hold on
plot3(cc_ms(1,:)+ms_x,cc_ms(2,:)+ms_y,cc_ms(3,:)+ms_z+y0_vec_ms(number),'.k'); hold on
plot3([cc_ms_axis(1,1) cc_ms_axis(1,2)]+ms_x,[cc_ms_axis(2,1) cc_ms_axis(2,2)]+ms_y,[cc_ms_axis(3,1) cc_ms_axis(3,2)]+ms_z+y0_vec_ms(number),'--k','LineWidth',2.5); hold on

%Axis
factor = 0.2;
grosor = 1.5;
aux_axis_ms = [-1 0 0; 0 -1 0; 0 0 1];
axis_ms = Tilt_Mat((270-atan(-ms_x/ms_y)/ang_c)*ang_c)*aux_axis_ms;

plot3([axis_ms(1,1)*factor 0]+ms_x,[axis_ms(1,2)*factor 0]+ms_y,[axis_ms(1,3)*factor 0]+ms_z+0.001,'k','LineWidth',grosor); hold on
plot3([axis_ms(2,1)*factor 0]+ms_x,[axis_ms(2,2)*factor 0]+ms_y,[axis_ms(2,3)*factor 0]+ms_z,'k','LineWidth',grosor); hold on
plot3([axis_ms(3,1)*factor 0]+ms_x,[axis_ms(3,2)*factor 0]+ms_y,[axis_ms(3,3)*factor 0]+ms_z,'k','LineWidth',grosor); hold on
%{
plot3([local_axis(1,1) local_axis(1,2)]+ms_x,[local_axis(2,1) local_axis(2,2)]+ms_y,[local_axis(3,1) local_axis(3,2)]+ms_z,'-c','LineWidth',2); hold on
plot3([local_axis(1,1) local_axis(1,3)]+ms_x,[local_axis(2,1) local_axis(2,3)]+ms_y,[local_axis(3,1) local_axis(3,3)]+ms_z,'-m','LineWidth',2); hold on
plot3([local_axis(1,1) local_axis(1,4)]+ms_x,[local_axis(2,1) local_axis(2,4)]+ms_y,[local_axis(3,1) local_axis(3,4)]+ms_z,'-k','LineWidth',2); hold on
%}
R = [0.02 0];
n = 20;
cyl_color = 'k';
closed = 1;
lines = 0;
fact1 = 0.2;
fact2 = 0.24;

x_ms1 = [axis_ms(1,1)*fact1+ms_x axis_ms(1,2)*fact1+ms_y axis_ms(1,3)*fact1+ms_z];
x_ms2 = [axis_ms(1,1)*fact2+ms_x axis_ms(1,2)*fact2+ms_y axis_ms(1,3)*fact2+ms_z];
[Cone1,EndPlate1,EndPlate2] = Cone(x_ms1,x_ms2,R,n,cyl_color,closed,lines);
y_ms1 = [axis_ms(2,1)*fact1+ms_x axis_ms(2,2)*fact1+ms_y axis_ms(2,3)*fact1+ms_z];
y_ms2 = [axis_ms(2,1)*fact2+ms_x axis_ms(2,2)*fact2+ms_y axis_ms(2,3)*fact2+ms_z];
[Cone2,EndPlate3,EndPlate4] = Cone(y_ms1,y_ms2,R,n,cyl_color,closed,lines);
z_ms1 = [axis_ms(3,1)*fact1+ms_x axis_ms(3,2)*fact1+ms_y axis_ms(3,3)*fact1+ms_z];
z_ms2 = [axis_ms(3,1)*fact2+ms_x axis_ms(3,2)*fact2+ms_y axis_ms(3,3)*fact2+ms_z];
[Cone3,EndPlate5,EndPlate6] = Cone(z_ms1,z_ms2,R,n,cyl_color,closed,lines);

if number == 1
    model = 'CC';
    fit = '1';
elseif number == 2
    model = 'EC';
    fit = '1';
elseif number == 3
    model = 'CC';
    fit = '2';
else
    model = 'EC';
    fit = '2';
end

title({['\fontsize{16} Model: ',model,'. Fit: ',fit,'.'],
       ['CME parameters -> Long=156.759 Stony, Lat=-2, Tilt=36.335'],
       ['STA Model -> \phi=',num2str(phi_cc_sta),', \theta=',num2str(theta_cc_sta),', \psi=',num2str(psi_cc_sta),', R=',num2str(radio_vec_sta(number)),' AU.'],
       ['MES Model -> \phi=',num2str(phi_cc_ms),', \theta=',num2str(theta_cc_ms),', \psi=',num2str(psi_cc_ms),', R=',num2str(radio_vec_ms(number)),' AU.']},'fontsize',14)

x = reshape(nube_x,[20 47]);
y = reshape(nube_y,[20 47]);
z = reshape(nube_z,[20 47]);
x = [x;x(1,:)];
y = [y;y(1,:)];
z = [z;z(1,:)];
C = sqrt(x.^2+y.^2+z.^2);
h = surf(x,y,z,C);
colormap('jet')
%shading interp
set(h,'edgecolor','none','facealpha',0.4)




