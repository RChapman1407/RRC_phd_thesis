%% doPlotData - plot and analyse VAR models for salinities.
%% DECADAL AVERAGES!
%% including standard error calculations and Q analysis
%% modified 29th Feb 2024
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
f1.Position(3:4)=[900 400];

for i=1:4
    subplot(2,4,i)
    %plot salinties in top row
    plotSNST(TSS(i).data,TSS(i).name);
    subplot(2,4,i+4)
    %plot temperatures
    plotTNTT(TSS(i).data,TSS(i).name);
end
hold off;

savefigure('compare_salins_temp_RRC',1,printfigs);

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
plot(SN-mSN,ST-mST,'LineWidth',1.2);
xlim([-0.02 0.02]);
ylim([-0.02 0.02]);
title(Name,'Interpreter','latex');
xlabel('$S_{\rm{Nor.}}$ anomaly (psu)','Interpreter','latex');
ylabel('$S_{\rm{Trop.}}$ anomaly (psu)','Interpreter','latex');
set(gca,'FontSize',11);
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
plot(TN-mTN,TT-mTT,'LineWidth',1.2);
xlim([-0.1 0.1]);
ylim([-0.1 0.1]);
% title(Name,'Interpreter','latex');
xlabel('$T_{\rm{Nor.}}$ anomaly (K)','Interpreter','latex');
ylabel('$T_{\rm{Trop.}}$ anomaly (K)','Interpreter','latex');
set(gca,'FontSize',11);
end

%%
function savefigure(name,number,printfigs)
% Save figure(number) as f-name-number.pdf/.png if printfigs==true
%
ffname = sprintf('..//Figures//fig-%s-%i', name, number);
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