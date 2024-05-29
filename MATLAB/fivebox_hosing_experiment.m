%Three box stochastic + hosing
%RRC
clear all; clc;

%define problem parameters
p = FamousAParameters();
Dt = 1; %delta t
T = 25000;
Ti = 0;
t = [Ti:Dt:T-1];
SB = p.Sb*1000; SS = p.Ss*1000; SIP = p.Sip*1000; %convert constants to psu
%N, T, S, IP, B
Sinput = [0.0353278, 0.0357409, p.Ss, p.Sip, p.Sb];
L= [0,0;0,0]; %no noise

N = floor((T-Ti)./Dt);

%TIPS
%time,Twait,Tpert,Trise,Tfall,height
H2 = hosingfunc3(t);
SxEM2 = FiveboxEMvectors(L,N,Dt,t,Sinput,p,H2);
SxEMpsu2 = SxEM2*1000;

%% FIGURES
grayColor = [.4 .4 .4];
printfigs = false;
q2 = p.lambda*(p.alpha*(p.Ts - p.T0) + p.beta*(SxEMpsu2(1,:) - p.Ss*1e3))./(1 + p.lambda*p.alpha*p.mu);

% PANEL A
subplot(4,4,[1,2,3,4])
hold on;
plot(t(:),H2(:)*1e-6,'b','LineWidth',1.2);
xlabel('Time (years)');
ylabel({'Freshwater';'Forcing (Sv)'});
text(0.005, 0.9,'(a)','Units','normalized');
set(gca,'LineWidth',1.2,'FontSize',12);
pbaspect([15 2 1])
box off;
hold off;

%PANEL C
subplot(4,4,[9,16])
%plot q on bif diagram
hold on;
plot(H2(1,:)*1e-6,q2*1e-9,'b','LineWidth',1.2);
xlabel('Freshwater Forcing (Sv)');
ylabel('AMOC strength (Sv)');
%overlay bifurcation diagram
SNH = readmatrix("SN_H_RRC_AMOC_5box_FAMOUSA_1xCO2.dat"); %column 1 H, column 2 S_N
STH = readmatrix("ST_H_RRC_AMOC_5box_FAMOUSA_1xCO2.dat"); %column 1 H, column 2 S_T
q_FAMOUSA = (p.lambda*(p.alpha*(p.Ts - p.T0) + p.beta*((SNH(:,2)*1e3) - p.Ss*1e3)))./(1 + p.lambda*p.alpha*p.mu);
H_FAMOUSA = SNH(:,1);
%upper branch
plot(H_FAMOUSA(1:373),q_FAMOUSA(1:373)*1e-9,'k','LineWidth',1.2);
%hopf branch
plot(H_FAMOUSA(374:382),q_FAMOUSA(374:382)*1e-9,'--k','LineWidth',1.2);
%saddle
plot(H_FAMOUSA(383:664),q_FAMOUSA(383:664)*1e-9,'--k','LineWidth',1.2);
%lower branch
plot(H_FAMOUSA(665:1000),q_FAMOUSA(665:1000)*1e-9,'k','LineWidth',1.2);
plot(H_FAMOUSA(1002:1622),q_FAMOUSA(1002:1622)*1e-9,'k','LineWidth',1.2);
%upper branch pt 2
plot(H_FAMOUSA(1624:2622),q_FAMOUSA(1624:2622)*1e-9,'k','LineWidth',1.2);
xlim([-0.5 1]);
plot(H_FAMOUSA(373),q_FAMOUSA(373)*1e-9,'k*','LineWidth',1.2); %position of hopf
plot(H_FAMOUSA(382),q_FAMOUSA(382)*1e-9,'k*','LineWidth',1.2); %position of upper fold
plot(H_FAMOUSA(664),q_FAMOUSA(664)*1e-9,'k*','LineWidth',1.2); %position of upper fold
set(gca,'LineWidth',1.2,'FontSize',12);
box off;
text(0.005, 0.95,'(c)','Units','normalized');
hold off;

%PANEL B
%plot SN against time
subplot(4,4,[5,6,7,8])
hold on;
plot(t(:),SxEMpsu2(1,:),'LineWidth',1.2); %SN
plot(t(:),SxEMpsu2(2,:),'LineWidth',1.2); %ST
plot(t(:),SxEMpsu2(3,:),'LineWidth',1.2); %S
plot(t(:),SxEMpsu2(4,:),'LineWidth',1.2); %S
plot(t(:),SxEMpsu2(5,:),'LineWidth',1.2); %S
legend('SN','ST','SS','SIP','SB');
xlabel('Time (years)');
ylabel({'Box';'salinities (psu)'});
set(gca,'LineWidth',1.2,'FontSize',12);
text(0.005, 0.9,'(b)','Units','normalized');
pbaspect([15 2 1])

hold off;
%%
function SxEM = FiveboxEMvectors(L,N,Dt,t,Sinput,p,H);
%Function takes inputs and uses threeboxfunction.m to find the
%differentials then uses Euler-Maruyama method to solve those differentials
%at each time step

%Inputs- L is matrix of noise parameters, N is size of time vector 
%,Dt is time step of method, t is time series being used, 
%Sinput is starting values in raw salinity, p is parameter file in SI units, 
%H is hosing in m^3 s^-1
%Output- SxEM is salinity in SI units 

SxEM = zeros(5,N); %preallocate array
len = length(t);
rng('default');
r1 = normrnd(0,1,[1,len]);%first noise vector
r2 = normrnd(0,1,[1,len]);%second noise vector
SN = Sinput(1);%raw salinities
ST = Sinput(2);
SS = Sinput(3);
SIP = Sinput(4);
SB = Sinput(5);
S0 = p.S0;

for j = 1:N
    Winc1 = sqrt(Dt)*r1(j);
    Winc2 = sqrt(Dt)*r2(j);
    S = [SN, ST, SS, SIP, SB]; %initial values which get updated
    dSdt = fiveboxfunction(S,p,H(j)); %!!!%
    
    SN = SN + dSdt(1,:)*Dt + L(1,1)*Winc1 + L(1,2)*Winc2;
    SxEM(1,j) = SN; %read out result to array
    ST = ST + dSdt(2,:)*Dt + L(2,1)*Winc1 + L(2,2)*Winc2;
    SxEM(2,j) = ST;
    SS = SS + dSdt(3,:)*Dt;
    SxEM(3,j) = SS;
    SIP = SIP + dSdt(4,:)*Dt;
    SxEM(4,j) = SIP;
    SB = SB + dSdt(5,:)*Dt;
    SxEM(5,j) = SB;
end
end
%%
function p = FamousAParameters()
%FAMOUS A parameters from Wood et. al. 2019 1xCO2
%initial salinities from FAMOUS B (not given in paper),
%but shouldnt affect results
%all parameters in SI units

p.Vn = 3.683e16; p.Vt = 5.151e16; p.Vs = 10.28e16; p.Vip = 21.29e16; p.Vb = 88.12e16;
p.Fn = 0.375e6; p.Ft = -0.723e6; p.Fs = 1.014e6; p.Fip = -0.666e6; p.Fb = 0;
p.Sn = 0.034912; p.St = 0.035435; p.Ss = 0.034427; p.Sip = 0.034668; p.Sb = 0.034538; p.S0 = 0.0350;
%Salinities from FAMOUS B (Alkhayuon et.al.2019)
p.Ts = 5.571; p.T0 = 3.26;
p.Kn = 5.439e6; p.Ks = 1.880e6; p.Kip = 89.778e6;
p.alpha = 0.12; p.beta = 790; p.eta = 66.061e6;
%alpha and beta from text in section 2.2
p.lambda = 2.66e7; p.gamma = 0.58; p.mu = 7e-8;
p.C = p.Vn*p.Sn + p.Vt*p.St + p.Vs*p.Ss + p.Vip*p.Sip + p.Vb*p.Sb;
p.An = 0.194; p.At = 0.597; p.As = -0.226; p.Aip = -0.565;

% used for modifications
p.offset=0;
p.dV = 0;
p.qscale =0;
end

%%
function [h] = hosingfunc0(time)
%input t is time steps being used
%output h is in m^3 s^-1

len = length(time);
h = zeros(1,len);

end

%%
function [h] = hosingfunc(time,Twait,Tpert,Trise,Tfall,height)
%input t is time steps being used
%output h is in m^3 s^-1

raterise = height./Trise;
ratefall = -(height-0.1)./Tfall;
tst0 = time>0;
tst1 = time<Twait;
tst2 = time<Twait+Trise;
tst3 = time<Twait+Trise+Tpert;
tst4 = time<Twait+Trise+Tpert+Tfall;
tst5 = time<Twait+Trise+Tpert+Tfall+5000;
h = 0.*tst0.*tst1.*time + (raterise.*(time-Twait)).*tst0.*(1 - tst1).*tst2 + height.*tst0.*(1 - tst2).*tst3 + ((ratefall.*(time-(Twait+Trise+Tpert+Tfall))+0.1).*tst0.*(1 - tst3).*tst4) + 0.1.*tst0.*(1-tst4).*tst5;% + 0.2*tst0.*tst4.*(1-tst5);
h = 1e6*h;
end
%%
function [h] = hosingfunc3(time)
%input t is time steps being used
%output h is in m^3 s^-1
%'Slow' Hosing version of H11, increases to 1Sv at 10000 years, decreases
%at same rate to 25000 years.

tst0 = time>0;
tst = time<10000;
    h = (1e-4.*time).*tst.*tst0 + (1-1e-4 .*(time-10000)).*(1-tst);
    h = h*1e6;
end

%% Wood five box model
%
function [dydt,q] = fiveboxfunction(Sinput,p,H)
% input is raw salinities, as is output
% dydt(1) is S_N and dydt(2) is S_T
dydt = zeros(5,1);

SN=Sinput(1);%inputs
ST=Sinput(2);
SS = Sinput(3);
SIP = Sinput(4);
SB = Sinput(5);
S0=p.S0;

% seconds per year
Y = 3.15e7;

% switch depending on sign of q:
q = (p.lambda.*(p.alpha.*(p.Ts-p.T0) +p.beta.*(SN-SS)))...
    ./(1+ p.lambda.*p.alpha.*p.mu);
aq =abs(q);

%H is hosing given as input
FN = p.Fn + p.An*H;
FT = p.Ft + p.At*H;
FS = p.Fs + p.As*H;
FIP = -FN - FT - FS; %for conservation 

if q>=0 
    dydt(1) = (aq*(ST-SN)+p.Kn*(ST-SN)-FN*S0).*(Y/p.Vn);
    dydt(2) = (aq*(p.gamma*SS+(1-p.gamma).*SIP-ST) ...
        +p.Ks*(SS-ST)+p.Kn*(SN-ST)-FT*S0).*(Y/p.Vt);
    dydt(3) = ((p.gamma*Y./p.Vs)*aq*(SB - SS) + (p.Kip*Y./p.Vs)*(SIP - SS) + (p.Ks*Y./p.Vs)*(ST - SS)...
        +(p.eta*Y./p.Vs)*(SB - SS) -(FS*Y./p.Vs)*S0);
    dydt(4) = ((1-p.gamma)*aq*(SB - SIP) +p.Kip*(SS - SIP)...
        - FIP*S0)*(Y./p.Vip);
    dydt(5) = (aq*(SN - SB) + p.eta*(SS - SB))*(Y./p.Vb);
else
    dydt(1) = (aq*(SB-SN)+p.Kn*(ST-SN)-FN*S0).*(Y/p.Vn);
    dydt(2) = (aq*(SN-ST)+p.Ks*(SS-ST)+p.Kn*(SN-ST)-FT*S0).*(Y/p.Vt);
    dydt(3) = ((p.gamma*Y./p.Vs)*aq*(ST - SS)+(p.Kip*Y./p.Vs)*(SIP - SS) +(p.Ks*Y./p.Vs)*(ST - SS)...
        +(p.eta*Y./p.Vs)*(SB - SS) - (FS*Y./p.Vs)*S0);
    dydt(4) = ((1-p.gamma)*aq*(ST - SIP) +p.Kip*(SS - SIP)...
        -FIP*S0).*(Y./p.Vip);
    dydt(5) = (p.gamma*aq*SS + (1-p.gamma)*aq*SIP - aq*SB...
        +p.eta*(SS - SB)).*(Y./p.Vb);
end

end

