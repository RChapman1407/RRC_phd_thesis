%% doPlotData - plot and analyse VAR models for salinities.
%% DECADAL AVERAGES!
%% including standard error calculations and Q analysis
%% modified 29th Feb 2024 PA and RC
%
close all
clear all

%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%%
printfigs=false;

%% Load timeseries
TSS(1).odata=load('HadGEM3LL.txt');
TSS(1).name='HadGEM3LL';
TSS(2).odata=load('HadGEM3MM.txt');
TSS(2).name='HadGEM3MM';
TSS(3).odata=load('CanESM5.txt');
TSS(3).name='CanESM5';
TSS(4).odata=load('MPI-ESM1-2-LR.txt');
TSS(4).name='MPI-ESM1-2-LR';

%% make decadal averages

for i=1:4
    tav=10;
    temp=movmean(TSS(i).odata,tav,1);
    sampletimes=1:tav:size(TSS(i).odata,1);
    temp2=temp(sampletimes,:);
    TSS(i).data=temp2;
end

%% Plot S_N and S_T timeseries

figure(1);
f1=gcf();
f1.Position(3:4)=[400 600];

for i=1:4
    subplot(4,2,2*i-1)
    plotSNST(TSS(i).data,TSS(i).name);
    subplot(4,2,2*i)
    plotTNTT(TSS(i).data,TSS(i).name);
end
hold off;

savefigure('compare-salins-dec-v2',1,printfigs);

keyboard;

figure(2);
f1=gcf();
f1.Position(3:4)=[550 400];

subplot(2,2,1)
[EM(1),EMSE(1),EMX(1),EMXSE(1)]=plotTSS(TSS(1).data,TSS(1).name);
legend('$S_T$','$S_N$','$T_T$','$T_N$',...
    'Location','east','Interpreter','latex');
subplot(2,2,2)
[EM(2),EMSE(2),EMX(2),EMXSE(2)]=plotTSS(TSS(2).data,TSS(2).name);
subplot(2,2,3)
[EM(3),EMSE(3),EMX(3),EMXSE(3)]=plotTSS(TSS(3).data,TSS(3).name);
subplot(2,2,4)
[EM(4),EMSE(4),EMX(4),EMXSE(4)]=plotTSS(TSS(4).data,TSS(4).name);
hold off;

savefigure('compare-salins-dec',1,printfigs);

for i=1:4
    %% Estimate equilibrium values from VAR fit
    EstMdl=EM(i);
    EstMdlSE=EMSE(i);

    Phi=cell2mat(EstMdl.AR);
    PhiSE=cell2mat(EstMdlSE.AR);
    AR11(i)=Phi(1,1);
    AR12(i)=Phi(1,2);
    AR21(i)=Phi(2,1);
    AR22(i)=Phi(2,2);
    AR11SE(i)=PhiSE(1,1);
    AR12SE(i)=PhiSE(1,2);
    AR21SE(i)=PhiSE(2,1);
    AR22SE(i)=PhiSE(2,2);
    cc=EstMdl.Constant;
    ccse=EstMdlSE.Constant;
    C1(i)=cc(1);
    C2(i)=cc(2);
    C1SE(i)=ccse(1);
    C2SE(i)=ccse(2);
    ES=inv(eye(2)-Phi)*cc;
    ESN(1,i)=ES(1);
    ESN(2,i)=ES(2);
    CV=EstMdl.Covariance;
    CV11(i)=CV(1,1);
    CV22(i)=CV(2,2);
    CV12(i)=CV(1,2);
    CV21(i)=CV(2,1);

    %% Estimate equilibrium values from VARX fit
    EstMdlX=EMX(i);
    EstMdlXSE=EMXSE(i);

    PhiX=cell2mat(EstMdlX.AR);
    PhiXSE=cell2mat(EstMdlXSE.AR);
    ARX11(i)=PhiX(1,1);
    ARX12(i)=PhiX(1,2);
    ARX21(i)=PhiX(2,1);
    ARX22(i)=PhiX(2,2);
    ARX11SE(i)=PhiXSE(1,1);
    ARX12SE(i)=PhiXSE(1,2);
    ARX21SE(i)=PhiXSE(2,1);
    ARX22SE(i)=PhiXSE(2,2);
    ccX=EstMdlX.Constant;
    ccXse=EstMdlXSE.Constant;
    CX1(i)=ccX(1);
    CX2(i)=ccX(2);
    CX1SE(i)=ccXse(1);
    CX2SE(i)=ccXse(2);
    ESX=inv(eye(2)-PhiX)*ccX;
    ESNX(1,i)=ESX(1);
    ESNX(2,i)=ESX(2);
    CVX=EstMdlX.Covariance;
    CVX11(i)=CVX(1,1);
    CVX22(i)=CVX(2,2);
    CVX12(i)=CVX(1,2);
    CVX21(i)=CVX(2,1);
    beta=EstMdlX.Beta;
    betaSE=EstMdlXSE.Beta;
    BX11(i)=beta(1,1);
    BX22(i)=beta(2,2);
    BX11SE(i)=betaSE(1,1);
    BX22SE(i)=betaSE(2,2);


    % output to command window
    TSS(i).name
    % VAR
    summarize(EstMdl);
    % VARX
    summarize(EstMdlX);
end

%% Plot bar charts showing coefficients for var fit for salinities
figure(3);
f2=gcf();
f2.InnerPosition(3:4)=[750 500];
set(f2,'Units','Inches');
pos = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

subplot(2,5,1);
bar(AR11);
hold on
errorbar(AR11,AR11SE,'.');
hold off
title('$A_{11}$','Interpreter','latex');
subplot(2,5,2);
bar(AR22);
hold on
errorbar(AR22,AR22SE,'.');
hold off
title('$A_{22}$','Interpreter','latex');
subplot(2,5,3);
bar(AR12);
hold on
errorbar(AR12,AR12SE,'.');
hold off
title('$A_{12}$','Interpreter','latex');
subplot(2,5,4);
bar(AR21);
hold on
errorbar(AR21,AR21SE,'.');
hold off
title('$A_{21}$','Interpreter','latex');
subplot(2,5,5);
bar(C1);
hold on
errorbar(C1,C1SE,'.');
hold off
title('$c_{1}$','Interpreter','latex');

subplot(2,5,6);
bar(CV11);
title('$CV_{11}$','Interpreter','latex');
subplot(2,5,7);
bar(CV22);
title('$CV_{22}$','Interpreter','latex');
subplot(2,5,8);
bar(CV12);
title('$CV_{12}$','Interpreter','latex');
subplot(2,5,9);
bar(CV21);
title('$CV_{21}$','Interpreter','latex');
subplot(2,5,10);
bar(C2);
hold on
errorbar(C2,C2SE,'.');
hold off
title('$c_{2}$','Interpreter','latex');

savefigure('compare-var-fits-dec',2,printfigs);

%% Plot bar charts showing coefficients for var X fit for salinities
figure(4);
f3=gcf();
f3.InnerPosition(3:4)=[900 500];
set(f3,'Units','Inches');
pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


subplot(2,6,1);
bar(ARX11);
hold on
errorbar(ARX11,ARX11SE,'.');
hold off
title("$A'_{11}$",'Interpreter','latex');
subplot(2,6,2);
bar(ARX22);
hold on
errorbar(ARX22,ARX22SE,'.');
hold off
title("$A'_{22}$",'Interpreter','latex');
subplot(2,6,3);
bar(ARX12);
hold on
errorbar(ARX12,ARX12SE,'.');
hold off
title("$A'_{12}$",'Interpreter','latex');
subplot(2,6,4);
bar(ARX21);
hold on
errorbar(ARX21,ARX21SE,'.');
hold off
title("$A'_{21}$",'Interpreter','latex');
subplot(2,6,5);
bar(CX1);
hold on
errorbar(CX1,CX1SE,'.');
hold off
title("$c'_{1}$",'Interpreter','latex');
subplot(2,6,6);
bar(BX11);
hold on
errorbar(BX11,BX11SE,'.');
hold off
title("$\beta_{11}$",'Interpreter','latex');

subplot(2,6,7);
bar(CVX11);
title("$CV'_{11}$",'Interpreter','latex');
subplot(2,6,8);
bar(CVX22);
title("$CV'_{22}$",'Interpreter','latex');
subplot(2,6,9);
bar(CVX12);
title("$CV'_{12}$",'Interpreter','latex');
subplot(2,6,10);
bar(CVX21);
title("$CV'_{21}$",'Interpreter','latex');
subplot(2,6,11);
bar(CX2);
hold on
errorbar(CX2,CX2SE,'.');
hold off
title("$c'_{2}$",'Interpreter','latex');
subplot(2,6,12);
bar(BX22);
hold on
errorbar(BX22,BX22SE,'.');
hold off
title("$\beta_{22}$",'Interpreter','latex');

savefigure('compare-varx-fits-dec',3,printfigs);

%% plot and fit Q timeseries

figure(5);
f1=gcf();
f1.Position(3:4)=[550 400];

subplot(2,1,1)
[EMQ(1),EMQSE(1),logL(1),EMQS(1),EMQSSE(1),logLS(1),EMQST(1),EMQSTSE(1),logLST(1)]=plotQ(TSS(1).data,TSS(1).name);
subplot(2,1,2)
[EMQ(2),EMQSE(2),logL(2),EMQS(2),EMQSSE(2),logLS(2),EMQST(2),EMQSTSE(2),logLST(2)]=plotQ(TSS(2).data,TSS(2).name);


savefigure('compare-Q-dec',4,printfigs);

%%
TSS(1).name
% model with q only
summarize(EMQ(1))
% model with q and S
summarize(EMQS(1))
% model with q, S and T
summarize(EMQST(1))


%%
TSS(2).name
% model with q only
summarize(EMQ(2))
% model with q and S
summarize(EMQS(2))
% model with q, S and T
summarize(EMQST(2))

%% Plot bar charts showing coefficients for Q fits
for j=1:2
    figure(5+j);
    f5=gcf();
    f5.InnerPosition(3:4)=[600 400];
    set(f5,'Units','Inches');
    pos = get(f5,'Position');
    set(f5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    AR=cell2mat(EMQ(j).AR);
    ARSE=cell2mat(EMQSE(j).AR);
    C=EMQ(j).Constant;
    CSE=EMQSE(j).Constant;
    CV=EMQ(j).Covariance;
    YY1= [AR C CV 0 0 0 0 ];
    YYSE1 = [ARSE CSE 0 0 0 0 0];

    AR=cell2mat(EMQS(j).AR);
    ARSE=cell2mat(EMQSSE(j).AR);
    C=EMQS(j).Constant;
    CSE=EMQSSE(j).Constant;
    CV=EMQS(j).Covariance;
    betas=EMQS(j).Beta;
    betasSE=EMQSSE(j).Beta;
    YY2= [AR C CV betas(1) betas(2) 0 0 ];
    YYSE2= [ARSE CSE 0 betasSE(1) betasSE(2) 0 0 ];

    AR=cell2mat(EMQST(j).AR);
    ARSE=cell2mat(EMQSTSE(j).AR);
    C=EMQST(j).Constant;
    CSE=EMQSTSE(j).Constant;
    betas=EMQST(j).Beta;
    betasSE=EMQSTSE(j).Beta;
    CV=EMQST(j).Covariance;
    YY3= [AR C CV betas(1) betas(2) betas(3) betas(4) ];
    YYSE3= [ARSE CSE 0 betasSE(1) betasSE(2) betasSE(3) betasSE(4) ];

    XX=categorical({'-','S','ST'});
    subplot(2,4,1)
    bar(XX,[YY1(1) YY2(1) YY3(1)]);
    hold on
    errorbar([YY1(1) YY2(1) YY3(1)],[YYSE1(1) YYSE2(1) YYSE3(1)],'.');
    hold off
    title("$A$",'Interpreter','latex');

    subplot(2,4,2)
    bar(XX,[YY1(2) YY2(2) YY3(2)]);
    hold on
    errorbar([YY1(2) YY2(2) YY3(2)],[YYSE1(2) YYSE2(2) YYSE3(2)],'.');
    hold off
    title("$c$",'Interpreter','latex');

    subplot(2,4,3)
    bar([YY1(3) YY2(3) YY3(3)]);
    title("$CV$",'Interpreter','latex');
    
    subplot(2,4,4)
    AIC1=aicbic(logL(j),2);
    AIC2=aicbic(logLS(j),4);
    AIC3=aicbic(logLST(j),6);
    
    bar([AIC1 AIC2 AIC3]);
    title("$AIC$",'Interpreter','latex');

    subplot(2,4,5)
    bar(XX,[YY1(4) YY2(4) YY3(4)]);
    hold on
    errorbar([YY1(4) YY2(4) YY3(4)],[YYSE1(4) YYSE2(4) YYSE3(4)],'.');
    hold off
    
    title("$\beta_1$",'Interpreter','latex');

    subplot(2,4,6)
    bar(XX,[YY1(5) YY2(5) YY3(5)]);
    hold on
    errorbar([YY1(5) YY2(5) YY3(5)],[YYSE1(5) YYSE2(5) YYSE3(5)],'.');
    hold off
    
    title("$\beta_2$",'Interpreter','latex');

    subplot(2,4,7)
    bar(XX,[YY1(6) YY2(6) YY3(6)]);
    hold on
    errorbar([YY1(6) YY2(6) YY3(6)],[YYSE1(6) YYSE2(6) YYSE3(6)],'.');
    hold off
    title("$\beta_3$",'Interpreter','latex');

    subplot(2,4,8)
    bar(XX,[YY1(7) YY2(7) YY3(7)]);
    hold on
    errorbar([YY1(7) YY2(7) YY3(7)],[YYSE1(7) YYSE2(7) YYSE3(7)],'.');
    hold off
    title("$\beta_4$",'Interpreter','latex');

    savefigure('compare-q-fits-dec',4+j,printfigs);
end


keyboard

%%

%%

%% plot S time series only
%
function plotSNST(TS,Name)
T=TS(:,1);
SN=TS(:,2);
ST=TS(:,3);
TN=TS(:,4);
TT=TS(:,5);
mSN=mean(SN);
mST=mean(ST);
%cSNT=cov(SN,ST);

hold on;
plot(SN-mSN,ST-mST,'r');
xlim([-0.02 0.02]);
ylim([-0.02 0.02]);
title(Name,'Interpreter','latex');
xlabel('$S_N$ anomaly (psu)','Interpreter','latex');
ylabel('$S_T$ anomaly (psu)','Interpreter','latex');
end

%% plot T time series only
%
function plotTNTT(TS,Name)
T=TS(:,1);
SN=TS(:,2);
ST=TS(:,3);
TN=TS(:,4);
TT=TS(:,5);
mTN=mean(TN);
mTT=mean(TT);
%cSNT=cov(SN,ST);

hold on;
plot(TN-mTN,TT-mTT,'r');
xlim([-0.1 0.1]);
ylim([-0.1 0.1]);
title(Name,'Interpreter','latex');
xlabel('$T_N$ anomaly (K)','Interpreter','latex');
ylabel('$T_T$ anomaly (K)','Interpreter','latex');
end

%% plot T and S time series, estimate models
%
function [EstMdl,EstSE,EstMdlX,EstXSE]=plotTSS(TS,Name)

T=TS(:,1);
SN=TS(:,2);
ST=TS(:,3);
TN=TS(:,4);
TT=TS(:,5);
mSN=mean(SN);
mST=mean(ST);
cSNT=cov(SN,ST);

yyaxis left
hold on;
plot(T,ST,'r');
plot(T,SN,'b');
ylim([34.75 36]);
title(Name,'Interpreter','latex');
xlabel('t (yr)','Interpreter','latex');
ylabel('$S_X$ (psu)','Interpreter','latex');
yyaxis right
plot(T,TT,'k');
plot(T,TN,'g');
ylim([0.0 15.5]);
ylabel('$T_X$ (C)','Interpreter','latex');

% VAR model of salinities
Mdl=varm(2,1);
[EstMdl,EstSE,logL,E]=estimate(Mdl,[SN ST]);

% VAR model of salinities forced by temps
MdlX=varm(2,1);
[EstMdlX,EstXSE,logLX,EX]=estimate(MdlX,[SN ST],'X',[TN,TT]);

end

%% plot Q time series, estimate models
%
function [EstMdl,EstSE,logL,EstMdlS,EstSSE,logLS,EstMdlST,EstSTSE,logLST]=plotQ(TS,Name)

T=TS(:,1);
SN=TS(:,2);
ST=TS(:,3);
TN=TS(:,4);
TT=TS(:,5);
Q=TS(:,6);
mSN=mean(SN);
mST=mean(ST);
cSNT=cov(SN,ST);

hold on;
plot(T,Q,'r');
title(Name,'Interpreter','latex');
xlabel('t (yr)','Interpreter','latex');
ylabel('$q$ (Sv)','Interpreter','latex');
hold off;

% AR model of Q
Mdl=varm(1,1);
[EstMdl,EstSE,logL,E]=estimate(Mdl,[Q]);

% AR model of Q forced by salinities only
MdlS=varm(1,1);
[EstMdlS,EstSSE,logLS,ES]=estimate(MdlS,[Q],'X',[SN ST]);

% AR model of Q forced by salinities and temps
MdlST=varm(1,1);
[EstMdlST,EstSTSE,logLST,EST]=estimate(MdlST,[Q],'X',[SN ST TN TT]);

end
