%plot bifurcation diagram three box model HadGEM3MM calibration (optimised
%F values)
%RRC Dec 2023
clc; clear all; clf;

p = HadGEM3MMParameters();
Dt = 1; %delta t
T = 3500;
Ti = 0;
t = [Ti:Dt:T-1];
Sinput = [0.0353278, 0.0357409];
L= [0,0;0,0]; %no noise
SB = p.Sb*1000; SS = p.Ss*1000;%convert constants to psu
N = floor((T-Ti)./Dt);

SNH = readmatrix("SN_H_PR_RRC_AMOC_3box_HadGEM3MM.dat");
STH = readmatrix("ST_H_PR_RRC_AMOC_3box_HadGEM3MM.dat");

%%

%upper branch
SNH_upper = SNH(1:246,:);
STH_upper = STH(1:246,:);
q_upper = p.lambda*(p.alpha*(p.Ts - p.T0) + p.beta*(((SNH_upper(:,2)/100 + p.S0) - p.Ss)*1e3))./(1 + p.lambda*p.alpha*p.mu);

%unstable branch
SNH_saddle = SNH(247:512,:);
STH_saddle = STH(247:512,:);
q_saddle = p.lambda*(p.alpha*(p.Ts - p.T0) + p.beta*(((SNH_saddle(:,2)/100 + p.S0) - p.Ss)*1e3))./(1 + p.lambda*p.alpha*p.mu);

%lower branch
SNH_lower = SNH(513:1550,:);
STH_lower = STH(513:1550,:);
q_lower = p.lambda*(p.alpha*(p.Ts - p.T0) + p.beta*(((SNH_lower(:,2)/100 + p.S0) - p.Ss)*1e3))./(1 + p.lambda*p.alpha*p.mu);

printfigs=false;

figure(1);

% Basin Boundary at H = 0.1
H6 = hosingfunc(t,0, 5000,1,1,0.1);
t6 = 5000:-1:0;
Dt6 = -1;
Sinput6 = [(-0.0608/100)+p.S0, (0.15/100)+p.S0];%just above saddle point
SxEM6 = ThreeboxEMvectors(L,N,Dt6,t6,Sinput6,p,H6);
SxEMpsu6 = SxEM6*1000;
Sinput7 = [(-0.0609/100)+p.S0,(0.15/100)+p.S0]; %just below saddle point
SxEM7 = ThreeboxEMvectors(L,N,Dt6,t6,Sinput7,p,H6);
SxEMpsu7 = SxEM7*1000;
hold on;

grayColor = [.7 .7 .7];
X_basin = [SxEMpsu6(1,:) SxEMpsu7(1,:)];
Y_basin = [SxEMpsu6(2,:) SxEMpsu7(2,:)];
a = area(X_basin,Y_basin,'FaceColor',grayColor);
a.FaceAlpha = 0.5;
a.LineWidth = 1.3;
%ON
plot(((0.0193071/100)+p.S0)*1000, ((0.0816398/100)+p.S0)*1000,'k*','LineWidth',1.2)
text(35.2,35.9,'On');
%OFF
plot(((-0.115287/100)+p.S0)*1000,((0.0429525/100)+p.S0)*1000,'k*','LineWidth',1.2);
text(33.75,35.55,'Off');
%SADDLE
plot(((-0.0618019/100)+p.S0)*1000,((0.150674/100)+p.S0)*1000,'k*','LineWidth',1.2)
text(34.45,36.47,'Saddle');
ht = text(34.5,35.5,{'Basin boundary'});
set(ht,'Rotation',-60);
xlabel('Northern box Salinity (psu)');
ylabel('Tropical box Salinity (psu)');
set(gca,'LineWidth',1.2,'FontSize',12);

%plot example trajectories
Dt = 1;
t0 = 0:1:10000;
%time,Twait,Tpert,Trise,Tfall,height
H0 = hosingfunc(t0,0, 10000,1,1,0.1);
%TO ON
Sinput8 = [0.0345, 0.037];
SxEM8 = ThreeboxEMvectors(L,N,Dt,t0,Sinput8,p,H0);
SxEMpsu8 = SxEM8*1000;
plot(SxEMpsu8(1,:),SxEMpsu8(2,:),'k');
Sinput9 = [0.035, 0.037];
SxEM9 = ThreeboxEMvectors(L,N,Dt,t0,Sinput9,p,H0);
SxEMpsu9 = SxEM9*1000;
plot(SxEMpsu9(1,:),SxEMpsu9(2,:),'k');
Sinput10 = [0.0349, 0.035];
SxEM10 = ThreeboxEMvectors(L,N,Dt,t0,Sinput10,p,H0);
SxEMpsu10 = SxEM10*1000;
plot(SxEMpsu10(1,:),SxEMpsu10(2,:),'k');
Sinput11 = [0.0354, 0.035];
SxEM11 = ThreeboxEMvectors(L,N,Dt,t0,Sinput11,p,H0);
SxEMpsu11 = SxEM11*1000;
plot(SxEMpsu11(1,:),SxEMpsu11(2,:),'k');
%TO OFF
Sinput12 = [0.0336, 0.0371];
SxEM12 = ThreeboxEMvectors(L,N,Dt,t0,Sinput12,p,H0);
SxEMpsu12 = SxEM12*1000;
plot(SxEMpsu12(1,:),SxEMpsu12(2,:),'k');
Sinput13 = [0.0342, 0.037];
SxEM13 = ThreeboxEMvectors(L,N,Dt,t0,Sinput13,p,H0);
SxEMpsu13 = SxEM13*1000;
plot(SxEMpsu13(1,:),SxEMpsu13(2,:),'k');
Sinput14 = [0.0345, 0.035];
SxEM14 = ThreeboxEMvectors(L,N,Dt,t0,Sinput14,p,H0);
SxEMpsu14 = SxEM14*1000;
plot(SxEMpsu14(1,:),SxEMpsu14(2,:),'k');
Sinput15 = [0.0336, 0.035];
SxEM15 = ThreeboxEMvectors(L,N,Dt,t0,Sinput15,p,H0);
SxEMpsu15 = SxEM15*1000;
plot(SxEMpsu15(1,:),SxEMpsu15(2,:),'k');

xlim([33.5 35.6]);
ylim([35 37]);
%%
%write out txt file of SN,ST,q for bif diagram in order
SN = [SNH_upper(:,2);fliplr(SNH_saddle(:,2));SNH_lower(:,2)];
ST = [STH_upper(:,2);fliplr(STH_saddle(:,2));STH_lower(:,2)];
q = [q_upper*1e-9;fliplr(q_saddle*1e-9);q_lower*1e-9];
H = [SNH_upper(:,1);fliplr(SNH_saddle(:,1));SNH_lower(:,1)];
%writematrix([SN,ST,q,H],'bif_HadGEM3MM_data.txt');

%%
function p = HadGEM3MMParameters()
%HadGEM3-GC3.1-MM parameters from Jon Baker 1xCO2
%all parameters in SI units
%F values optimised

p.Vn = 4.192e16; p.Vt = 4.192e16; p.Vs = 13.26e16; p.Vip = 16.95e16; p.Vb = 96.76e16;
p.Fn = 0.2799e6; p.Ft = -0.7920e6; p.Fs = 1.485e6; p.Fip = -0.9729e6; p.Fb = 0;
p.Sn = 0.034912; p.St = 0.035435; p.Ss = 0.034427; p.Sip = 0.0345; p.Sb = 0.034538; p.S0 = 0.0350;

p.Ts = 5.349; p.T0 = 4.514;
p.Kn = 4.73e6; p.Ks = 7.68e6; p.Kip = 123.49e6;
%these still from FAMOUS A (no update)
p.alpha = 0.12; p.beta = 790; p.eta = 66.061e6;
%alpha and beta from text in section 2.2
%gamma still from FAMOUS A, no update
p.lambda = 2.328e7; p.gamma = 0.58; p.mu = 0;
p.C = 4.6994e16;
p.An = 0.9841; p.At = -0.1853; p.As = -0.1742; p.Aip = -0.6245; 
% used for modifications
p.offset=0;
p.dV = 0;
p.qscale =0;
end

%%
function savefigure_RRC(name,number,printfigs)
% Save figure(number) as f-name-number.pdf/.png if printfigs==true
%
ffname = sprintf('Figures//fig-%s-%i', name, number);
h=figure(number);
set(h,'Units','centimeters');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
if(printfigs==true)
    %savefig(sprintf('%s.fig',ffname));
    print(sprintf('%s.pdf',ffname),'-dpdf','-bestfit');
    %print(sprintf('%s.png',ffname),'-dpng');
end
end
%%
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

