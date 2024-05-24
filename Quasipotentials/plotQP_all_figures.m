% plot QP, vector field components and norm of vector field
clear all; clc; clf;

u=load("Qpot.txt");
vfx=load("VFx.txt");
vfy=load("VFy.txt");
xx=load("Xaxis.txt");
yy=load("Yaxis.txt");
path = load('Instanton.txt');

nvf=vfy.^2+vfx.^2;

ind = find(u > 1e5);
u(ind) = NaN;


%%
%PLOT QP with basin boundary
figure(1);
contourf(xx,yy,u,[0 0.005 0.01  0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07]);
c = colorbar;
c.Ticks = [0 0.005 0.01  0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07];

xlabel('\boldmath $S_{\rm{Nor.}}$','Interpreter','latex','FontSize',12);
ylabel('\boldmath $S_{\rm{Trop.}}$','Interpreter','latex','FontSize',12);
title('FAMOUS B 2xCO2 from ON (N = 1024, K = 22)')
hold on;
%plot point of stable on state
SNi = 35.3722;
STi = 36.3670;
plot(SNi,STi,'w*','LineWidth',1.2);
text(35.5,36.4,'On','Color','white');
%plot point of stable off state
SNend = 33.0155;
STend = 36.4995;
grey = [0.8 0.8 0.8];
plot(SNend, STend,'w*','LineWidth',1.2);
text(32.9,36.9,'Off','Color','white');
%plot position of saddle
SNsaddle = 33.7996;
STsaddle = 43.4262;
plot(SNsaddle, STsaddle,'w*','LineWidth',1.2);
text(33.9,43.5,'Saddle','Color','white');


%Calculate basin boundary at H=0 on QP
t = 5000:-1:0;
H0 = hosingfunc0(t);
Dt = -1;
N = length(t);
L = [0,0;0,0]; %no noise
p = FamousB2xCO2Parameters;
Sinput4 = [0.0338, 0.0434262];%just above saddle point
SxEM4 = ThreeboxEMvectors(L,N,Dt,t,Sinput4,p,H0);
SxEMpsu4 = SxEM4*1000;
Sinput5 = [0.0338,0.0433]; %just below saddle point
SxEM5 = ThreeboxEMvectors(L,N,Dt,t,Sinput5,p,H0);
SxEMpsu5 = SxEM5*1000;

%plot basin boundary
plot(SxEMpsu4(1,:),SxEMpsu4(2,:),'w','LineWidth',1.2);
plot(SxEMpsu5(1,:),SxEMpsu5(2,:),'w','LineWidth',1.2);
xlim([32 39]);
ylim([35 45]);

%adding instanton to plot
x_path = path(:,1);
y_path = path(:,2);
plot(x_path,y_path,'r','LineWidth',1.2);

%print('AMOC-QP_N_1024_K_22_FAMOUSB2xco2_on_SADDLE_BOUNDARY.png','-dpng');

%%
figure(2);
subplot(2,2,1);
contourf(xx,yy,u);
colorbar;
xlabel("S_N");
ylabel("S_T");
title('Quasipotential');
hold on;
%adding instanton to plot
x_path = path(:,1);
y_path = path(:,2);
plot(x_path,y_path,'r');


subplot(2,2,2);
contourf(xx,yy,log(nvf));
colorbar;
xlabel("S_N");
ylabel("S_T");
title('Log Norm of VF');

subplot(2,2,3);
contourf(xx,yy,vfx);
colorbar;
xlabel("S_N");
ylabel("S_T");
title('VF-x');


subplot(2,2,4);
contourf(xx,yy,vfy);
colorbar;
xlabel("S_N");
ylabel("S_T");
title('VF-y');
%print('AMOC-QP-portrait4.png','-dpng');

%%
figure(3);
s = surf(xx,yy,u,'FaceAlpha',0.7);
%colorbar;
s.EdgeColor = 'none';
hold on;
xlabel('$\mathbf{S_{Nor.}}$','Interpreter','latex','FontSize',12);
ylabel('$\mathbf{S_{Trop.}}$','Interpreter','latex','FontSize',12);
zlabel("Quasipotential value")

%print('AMOC_QP_N_1024_K_22_FAMOUSB2xco2_on_SURFACE.png','-dpng');

%%
function [h] = hosingfunc0(time)
%input t is time steps being used
%output h is in m^3 s^-1

len = length(time);
h = zeros(1,len);

end
%%
%%
function SxEM = ThreeboxEMvectors(L,N,Dt,t,Sinput,p,H);
%Function takes inputs and uses threeboxfunction.m to find the
%differentials then uses Euler-Maruyama method to solve those differentials
%at each time step

%Inputs- L is matrix of noise parameters, N is size of time vector 
%,Dt is time step of method, t is time series being used, 
%Sinput is starting values in raw salinity, p is parameter file in SI units, 
%H is hosing in m^3 s^-1
%Output- SxEM is salinity in SI units 

SxEM = zeros(3,N); %preallocate array
% randn('state',100);
len = length(t);
rng('default');
r1 = normrnd(0,1,[1,len]);%first noise vector
r2 = normrnd(0,1,[1,len]);%second noise vector
SN = Sinput(1);%raw salinities
ST = Sinput(2);
SS = p.Ss;        
SIP = p.Sip;
SB = p.Sb;
S0 = p.S0;

for j = 1:N
    Winc1 = sqrt(Dt)*r1(j);
    Winc2 = sqrt(Dt)*r2(j);
    S = [SN,ST]; %initial values which get updated
    dSdt = threeboxfunction(S,p,H(j));
    
    SN = SN + dSdt(1,:)*Dt + L(1,1)*Winc1 + L(1,2)*Winc2;
    SxEM(1,j) = SN; %read out result to array
    ST = ST + dSdt(2,:)*Dt + L(2,1)*Winc1 + L(2,2)*Winc2;
    SxEM(2,j) = ST;
    SIP = (p.C-p.Vn.*SN - p.Vt.*ST - p.Vs.*p.Ss - p.Vb.*p.Sb)/p.Vip;
    SxEM(3,j) = SIP;
end
end
%%
%% Wood three box model
%
function [dydt,q] = threeboxfunction(Sinput,p,H)
% input is raw salinities, as is output
% dydt(1) is S_N and dydt(2) is S_T
dydt = [0;0];

SN=Sinput(1);%inputs
ST=Sinput(2);
SB=p.Sb;
SS=p.Ss;
S0=p.S0;
% include p.dV change in Vn/Vt

% seconds per year
Y = 3.15e7;

SIP = (p.C-(p.Vn+p.dV).*SN - (p.Vt-p.dV).*ST - p.Vs.*SS - p.Vb.*SB)/p.Vip;

% switch depending on sign of q:
q = (p.lambda.*(p.alpha.*(p.Ts-p.T0) +p.beta.*(SN-SS)))...
    ./(1+ p.lambda.*p.alpha.*p.mu);
aq =abs(q);

%H is hosing given as input
FN = p.Fn + p.An*H;
FT = p.Ft + p.At*H;

if q>=0
    dydt(1) = (aq*(ST-SN)+p.Kn*(ST-SN)-FN*S0).*(Y/(p.Vn+p.dV));
    dydt(2) = (aq*(p.gamma*SS+(1-p.gamma).*SIP-ST) ...
        +p.Ks*(SS-ST)+p.Kn*(SN-ST)-FT*S0).*(Y/(p.Vt-p.dV));
else
    dydt(1) = (aq*(SB-SN)+p.Kn*(ST-SN)-FN*S0).*(Y/(p.Vn+p.dV));
    dydt(2) = (aq*(SN-ST)+p.Ks*(SS-ST)+p.Kn*(SN-ST)-FT*S0).*(Y/(p.Vt-p.dV));
end

end
%%
function p = FamousB2xCO2Parameters()
%FAMOUS B 2xCO2 Parameters from (Wood et. al 209)

p.Vn = 3.683e16; p.Vt = 5.418e16; p.Vs = 6.097e16; p.Vip = 14.86e16; p.Vb = 99.25e16;
p.Fn = 0.486e6; p.Ft = -0.997e6; p.Fs = 1.265e6; p.Fip = -0.754e6; p.Fb= 0;
p.Sn = 0.034912; p.St = 0.035435; p.Ss = 0.034427; p.Sip = 0.034668; p.Sb = 0.034538; p.S0 = 0.035;
%Salinities from FAMOUS B (Alkhayuon et.al.2019)
p.Ts = 7.919; p.T0 = 3.87;
p.Kn = 1.762e6; p.Ks = 1.872e6; p.Kip = 99.977e6;
p.alpha = 0.12; p.beta = 790; p.eta = 33.264e6;
%alpha and beta from text in section 2.2
p.lambda = 1.62e7; p.gamma = 0.36; p.mu = 22.0e-8;
p.C = p.Vn*p.Sn + p.Vt*p.St + p.Vs*p.Ss + p.Vip*p.Sip + p.Vb*p.Sb;
p.An = 0.131; p.At = 0.696; p.As = -0.263; p.Aip = -0.564;

p.dV = 0;
end
