%% Tipping Probabilities for HadGEM3MM
%RRC Oct 2023
%02/10/23 adjust to extract time of tipping to then later calculate escape
%times

%Estimate the probability of the three box model tipping with a prescribed
%noise using MC simulations
%Noise level and hosing function can be manually changed

clear all; clc; close all;
%% define problem parameters
p = HadGEM3MMParameters();
Dt = 1; 

%HadGEM3MM
%updated for HadGEM3MM decadal estimates of noise
sigma = [-27.1641 , -26.9686 , -0.7315 , -0.3197, -0.0393]; %HadGEM3MM results
Q = covmat(sigma); %is positive definite
noise_matrix = chol(Q,'lower'); %HadGEM3MM estimated nosie
noise_matrix = 5*noise_matrix;

%optimise freshwater flux parameters 
p.Fn = p.Fno + sigma(4)*1e6; 
p.Ft = p.Fto + sigma(5)*1e6; %try also perfecting Ft
p.Fs = p.Fso - (sigma(4)+sigma(5))*1e6;

%% Reconstruct the model with noise estimate
t = [-200:1:1000]; %allow for 200 years spin up of box model
n = size(t,2);
Sinput = [p.Sn p.St];
H = u03_r70(t);
I = eye(2);

%realisation of noise seed number, rel 
%total number of realisations, rels
rels = 1000;
%for loop to go trhough different seed numbers for realisations
%initialise array Sx
Sx = zeros(4,n,rels);
%initialise array of tipping (1 is tipped, 0 is not tipped)
tip = zeros(1,rels);
%initialise array of tipping times
tip_time = zeros(1,rels);

for rel = 0:rels

%salinities *1000 to convert from raw salinities to psu
    Sx(:,:,rel+1) = 1000*ThreeboxEMvectors(noise_matrix,I,n,Dt,t,Sinput,p,H,rel);%psu
    q = Sx(4,200:end,rel+1)/1e9; %vector of q in Sv

      B = any(q < 5); %Does it tip at all? 1 if yes, 0 if no.
    tip(1,rel+1) = B;
    
%find the time of the tipping
k = find(q<5,1); %find index of first time the condition is filled
if k>0
    tip_time(:,rel+1) = k;
else
    tip_time(:,rel+1) = 0;
end
    
end
    % plot(t,Sx(1,:,6),'b-.');%Northern
% hold on;
    % plot(t,Sx(2,:,6),'r-.');%Tropical
    % plot(t,Sx(3,:,6))%IP
figure;
%plot of q
%plot shaded area where model is being forced
xlim = [0 1000];
ylim = [2 18];
Y = [0 20 20 0];
X = [0 0 70 70];
patch(X,Y,'k','FaceAlpha',.2);
hold on;
    plot(t(201:500),Sx(4,201:500,6)/1e9,'LineWidth',1.2)%q /1000 to get back to SI units
    plot(t(201:500),Sx(4,201:500,25)/1e9,'LineWidth',1.2) %tips
    plot(t(201:500),Sx(4,201:500,54)/1e9,'LineWidth',1.2)
    plot(t(201:500),Sx(4,201:500,77)/1e9,'LineWidth',1.2)
    plot(t(201:500),Sx(4,201:500,90)/1e9,'LineWidth',1.2)
    plot(t(201:500),Sx(4,201:500,668)/1e9,'LineWidth',1.2) %tips
    plot(t(201:500),Sx(4,201:500,201)/1e9,'LineWidth',1.2)
% legend('$S_{N}\ HadGEM3-GC3.1-LL$','$S_{T}\ HadGEM3-GC3.1-LL$','$S_{N}\ \rm{box\ model\ with\ estimated\ noise}$','$S_{T}\ \rm{box\ model\ with\ estimated\ noise}$','interpreter','latex','FontSize',12);
set(gca,'LineWidth',1.2,'FontSize',12);
xlabel('\boldmath $\rm{Time}\ (years)$','Interpreter','latex','FontSize',12);
ylabel('\boldmath $q (\rm{Sv})$','Interpreter','latex','fontsize',12);
yline(5);
%text(0.02, 0.95,'(a)','Units','normalized','FontSize',16);
%title('u03 70yrs');
printfigs = false;
savefigure_RRC('Tipping_trajectories_5_HadGEM3MM_noise_u03_r70_300_years',1,printfigs)


%% Tipping Probability
%number of non-zero elements in 'tip' (number of realisations that tipped)
tipped = nnz(tip);
prob = tipped/(rels+1)

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


%% Euler Maruyama integration with given noise
%
function SxEM = ThreeboxEMvectors(Linput,Q,N,Dt,t,Sinput,p,H,rel);

SxEM = zeros(4,N); %preallocate array
% randn('state',100);
len = length(t);
rng(rel);
r1 = normrnd(0,1,[1,len]);%first noise vector
r2 = normrnd(0,1,[1,len]);
SN = Sinput(1);
ST = Sinput(2);
SS = p.Ss;
SIP = p.Sip;
SB = p.Sb;
S0 = p.S0;

% Cholesky factorization such that SQ*SQ'=Q
SQ = chol(Q, 'lower');
% following also seems to work..
% SQ=sqrtm(Q);

% SQ = Q;

for j = 1:N
    %L=L*q
    Winc1 = sqrt(Dt)*r1(j);
    Winc2 = sqrt(Dt)*r2(j);
    S = [SN,ST]; %updated each iteration
    [dSdt,q] = threeboxfunction(S,p,H(j));
%     qSv = q./ p.qscale;
    
    L = (1+p.qscale*(q<0))*Linput;
%     L = Linput;
    ML = L*SQ;
    % the output of threeboxfunction is not in psu
    SN = SN + dSdt(1,:)*Dt + ML(1,1)*Winc1 + ML(1,2)*Winc2;
    SxEM(1,j) = SN; % read out result to array
    ST = ST + dSdt(2,:)*Dt + ML(2,1)*Winc1 + ML(2,2)*Winc2;
    SxEM(2,j) = ST; % this is now iterating S :)
    SIP = (p.C-(p.Vn+p.dV).*SN - (p.Vt-p.dV).*ST - p.Vs.*p.Ss - p.Vb.*p.Sb)/p.Vip;
    SxEM(3,j) = SIP;
    SxEM(4,j) = q;
    % update q here and do L*q at beginning of loop
end

end

%% HadGEM3MMparameter set
%
function p = HadGEM3MMParameters()

%
p.Vn = 4.192e16; p.Vt = 4.192e16; p.Vs = 13.26e16; p.Vip = 16.95e16; p.Vb = 96.76e16;
p.Fn = 0.5996e6; p.Ft = -0.7527e6; p.Fs = 1.126e6; p.Fip = -0.6728e6; p.Fb = 0;
p.Fno = p.Fn; p.Fso = p.Fs; p.Fto = p.Ft;
p.Sn = 0.03488; p.St = 0.03563; p.Ss = 0.03444; p.Sip = 0.0345; p.Sb = 0.03475; p.S0 = 0.0350;

%
p.Ts = 5.349; p.T0 = 4.514;
p.Kn = 4.73e6; p.Ks = 7.68e6; p.Kip = 123.49e6;
%these still from FAMOUS A (no update)
p.alpha = 0.12; p.beta = 790; p.eta = 66.061e6;
%gamma still from FAMOUS A, no update
p.lambda = 2.328e7; p.gamma = 0.58; p.mu = 0;
p.C = p.Vn*p.Sn + p.Vt*p.St + p.Vs*p.Ss + p.Vip*p.Sip + p.Vb*p.Sb;
p.An = 0.9841; p.At = -0.1853; p.As = -0.1742; p.Aip = -0.6245;

% used for modifications
p.offset=0;
p.dV = 0;
p.qscale =0;


end
%% parameterize covariance matrix
%
function QM=covmat(sigma)

% always produces pos definite matrix
var1=exp(sigma(1));
var2=exp(sigma(2));
cov=exp((sigma(1)+sigma(2))/2)*tanh(sigma(3));
QM = [var1,cov;cov,var2];%
end

%%
function [h] = hosingfunc0(time)
%input t is time steps being used
%output h is in m^3 s^-1

len = length(time);
h = zeros(1,len);

end

%%
 function [h] = u03_hos(time)
%input t is time steps being used
%output h is in Sv

Tpert = 1000; %length of perturbation
Tpert_start = 0; %start perturbation after spin up
tst0 = time>0;
tst1 = time<Tpert_start;
tst2 = time<Tpert_start + Tpert;
step = 0.3;%Sv
h = 0.*tst0.*tst1.*time + step.*(1-tst1).*tst2 + 0.*tst0.*(1-tst2).*time;
h = 1e6*h;
end
%%
 function [h] = u03_r100(time)
%input t is time steps being used
%output h is in m^3 s^-1

Tpert = 100; %length of perturbation
Tpert_start = 0; %start perturbation after spin up
tst0 = time>0;
tst1 = time<Tpert_start;
tst2 = time<Tpert_start + Tpert;
step = 0.3;%Sv
h = 0.*tst0.*tst1.*time + step.*(1-tst1).*tst2 + 0.*tst0.*(1-tst2).*time;
h = 1e6*h;
end
%%
function [h] = u03_r70(time)
%input t is time steps being used
%output h is in m^3 s^-1

Tpert = 70; %length of perturbation
Tpert_start = 0; %start perturbation after spin up
tst0 = time>0;
tst1 = time<Tpert_start;
tst2 = time<Tpert_start + Tpert;
step = 0.3;%Sv
h = 0.*tst0.*tst1.*time + step.*(1-tst1).*tst2 + 0.*tst0.*(1-tst2).*time;
h = 1e6*h;
end
%%
 function [h] = u03_r50(time)
%input t is time steps being used
%output h is in m^3 s^-1

Tpert = 50; %length of perturbation
Tpert_start = 0; %start perturbation after spin up
tst0 = time>0;
tst1 = time<Tpert_start;
tst2 = time<Tpert_start + Tpert;
step = 0.3;%Sv
h = 0.*tst0.*tst1.*time + step.*(1-tst1).*tst2 + 0.*tst0.*(1-tst2).*time;
h = 1e6*h;
 end
%%
 function [h] = u03_r20(time)
%input t is time steps being used
%output h is in m^3 s^-1

Tpert = 20; %length of perturbation
Tpert_start = 0; %start perturbation after spin up
tst0 = time>0;
tst1 = time<Tpert_start;
tst2 = time<Tpert_start + Tpert;
step = 0.3;%Sv
h = 0.*tst0.*tst1.*time + step.*(1-tst1).*tst2 + 0.*tst0.*(1-tst2).*time;
h = 1e6*h;
 end
%%
 function [h] = u01_hos(time)
%input t is time steps being used
%output h is in Sv

Tpert = 1000; %length of perturbation
Tpert_start = 0; %start perturbation after spin up
tst0 = time>0;
tst1 = time<Tpert_start;
tst2 = time<Tpert_start + Tpert;
step = 0.1;%Sv
h = 0.*tst0.*tst1.*time + step.*(1-tst1).*tst2 + 0.*tst0.*(1-tst2).*time;
h = 1e6*h;
 end
%%
 function [h] = u01_r100(time)
%input t is time steps being used
%output h is in m^3 s^-1

Tpert = 100; %length of perturbation
Tpert_start = 0; %start perturbation after spin up
tst0 = time>0;
tst1 = time<Tpert_start;
tst2 = time<Tpert_start + Tpert;
step = 0.1;%Sv
h = 0.*tst0.*tst1.*time + step.*(1-tst1).*tst2 + 0.*tst0.*(1-tst2).*time;
h = 1e6*h;
end
%%
 function [h] = u02_hos(time)
%input t is time steps being used
%output h is in Sv

Tpert = 1000; %length of perturbation
Tpert_start = 0; %start perturbation after spin up
tst0 = time>0;
tst1 = time<Tpert_start;
tst2 = time<Tpert_start + Tpert;
step = 0.2;%Sv
h = 0.*tst0.*tst1.*time + step.*(1-tst1).*tst2 + 0.*tst0.*(1-tst2).*time;
h = 1e6*h;
 end
