% HadGEM3MM_run_parameters.m
% script to define run parameters

%updated for HadGEM3MM decadal estimates of noise
par.sigma = [-27.1641 , -26.9686 , -0.7315 , -0.3197, -0.0393]; %HadGEM3MM results
Q = covmat(par.sigma); %is positive definite
par.noise_matrix = chol(Q,'lower'); %HadGEM3MM estimated nosie
% par.noise_matrix = zeros(2); 

%realisation of noise seed number, rel 
%total number of realisations, rels
par.rels = 1000;
par.H = @hosingfunc0;
par.name="default";


%% default parameter set
pars(1)=par;
pars(1).noise_matrix = zeros(2);
pars(1).name="hosingfunc0 noise zero";
%%
pars(2)=par;
pars(2).noise_matrix = [0.3369e-5, 0 ; -0.0551e-5, 0.4006e-5];
pars(2).name="hosingfunc0 noise HadGEM3LL";
%%
pars(3)=par;
%default noise
pars(3).name="hosingfunc0 noise HadGEM3MM";
%%
pars(4)=par;
pars(4).noise_matrix = [0.386e-5, 0 ; -0.1579e-5, 0.2792e-5];
pars(4).name="hosingfunc0 noise CanESM5";
%%
pars(5)=par;
pars(5).noise_matrix = [0.8686e-5, 0 ; -0.0443e-5, 0.3543e-5];
pars(5).name="hosingfunc0 noise MPI";
%% 
pars(6)=par;
pars(6).noise_matrix = 5*par.noise_matrix;
pars(6).name="hosingfunc0 noise 5";
%%
pars(7)=par;
pars(7).noise_matrix = 10*par.noise_matrix;
pars(7).name="hosingfunc0 noise 10";
%%
pars(8) = par;
pars(8).noise_matrix = 15*par.noise_matrix;
pars(8).name = "hosingfunc0 noise 15"
%%
pars(9) = par;
pars(9).noise_matrix = 20*par.noise_matrix;
pars(9).name = "hosingfunc0 noise 20"
%%
pars(10) = par;
pars(10).noise_matrix = zeros(2);
pars(10).H = @u01_r100;
pars(10).name = "u01_r100 noise zero"
%%
pars(11) = par;
pars(11).noise_matrix = [0.3369e-5, 0 ; -0.0551e-5, 0.4006e-5];
pars(11).H = @u01_r100;
pars(11).name = "u01_r100 noise HadGEM3LL"
%%
pars(12) = par;
%default noise
pars(12).H = @u01_r100;
pars(12).name = "u01_r100 noise HadGEM3MM"
%%
pars(13) = par;
pars(13).noise_matrix = [0.386e-5, 0 ; -0.1579e-5, 0.2792e-5];
pars(13).H = @u01_r100;
pars(13).name = "u01_r100 noise CanESM5"
%%
pars(14) = par;
pars(14).noise_matrix = [0.8686e-5, 0 ; -0.0443e-5, 0.3543e-5];
pars(14).H = @u01_r100;
pars(14).name = "u01_r100 noise MPI"
%%
pars(15) = par;
pars(15).noise_matrix = 5*par.noise_matrix;
pars(15).H = @u01_r100;
pars(15).name = "u01_r100 noise 5"
%%
pars(16) = par;
pars(16).noise_matrix = 10*par.noise_matrix;
pars(16).H = @u01_r100;
pars(16).name = "u01_r100 noise 10"
%%
pars(17) = par;
pars(17).noise_matrix = 15*par.noise_matrix;
pars(17).H = @u01_r100;
pars(17).name = "u01_r100 noise 15"
%%
pars(18) = par;
pars(18).noise_matrix = 20*par.noise_matrix;
pars(18).H = @u01_r100;
pars(18).name = "u01_r100 noise 20"
%%
pars(19) = par;
pars(19).noise_matrix = zeros(2);
pars(19).H = @u01_hos;
pars(19).name = "u01_hos noise zero"
%%
pars(20) = par;
pars(20).noise_matrix = [0.3369e-5, 0 ; -0.0551e-5, 0.4006e-5];
pars(20).H = @u01_hos;
pars(20).name = "u01_hos noise HadGEM3LL"
%%
pars(21) = par;
%default noise
pars(21).H = @u01_hos;
pars(21).name = "u01_hos noise HadGEM3MM"
%%
pars(22) = par;
pars(22).noise_matrix = [0.386e-5, 0 ; -0.1579e-5, 0.2792e-5];
pars(22).H = @u01_hos;
pars(22).name = "u01_hos noise CanESM5"
%%
pars(23) = par;
pars(23).noise_matrix = [0.8686e-5, 0 ; -0.0443e-5, 0.3543e-5];
pars(23).H = @u01_hos;
pars(23).name = "u01_hos noise MPI"
%%
pars(24) = par;
pars(24).noise_matrix = 5*par.noise_matrix;
pars(24).H = @u01_hos;
pars(24).name = "u01_hos noise 5"
%%
pars(25) = par;
pars(25).noise_matrix = 10*par.noise_matrix;
pars(25).H = @u01_hos;
pars(25).name = "u01_hos noise 10"
%%
pars(26) = par;
pars(26).noise_matrix = 15*par.noise_matrix;
pars(26).H = @u01_hos;
pars(26).name = "u01_hos noise 15"
%%
pars(27) = par;
pars(27).noise_matrix = 20*par.noise_matrix;
pars(27).H = @u01_hos;
pars(27).name = "u01_hos noise 20"
%%
pars(28) = par;
pars(28).noise_matrix = zeros(2);
pars(28).H = @u03_r20;
pars(28).name = "u03_r20 noise zero"
%%
pars(29) = par;
pars(29).noise_matrix = [0.3369e-5, 0 ; -0.551e-5, 0.4006e-5];
pars(29).H = @u03_r20;
pars(29).name = "u03_r20 noise HadGEM3LL"
%%
pars(30) = par;
% default noise
pars(30).H = @u03_r20;
pars(30).name = "u03_r20 noise HadGEM3MM"
%%
pars(31) = par;
pars(31).noise_matrix = [0.386e-5, 0 ; -0.1579e-5, 0.2792e-5];
pars(31).H = @u03_r20;
pars(31).name = "u03_r20 noise CanESM5"
%%
pars(32) = par;
pars(32).noise_matrix = [0.8686e-5, 0 ; -0.0443e-5, 0.3543e-5];
pars(32).H = @u03_r20;
pars(32).name = "u03_r20 noise MPI"
%%
pars(33) = par;
pars(33).noise_matrix = 5*par.noise_matrix;
pars(33).H = @u03_r20;
pars(33).name = "u03_r20 noise 5"
%%
pars(34) = par;
pars(34).noise_matrix = 10*par.noise_matrix;
pars(34).H = @u03_r20;
pars(34).name = "u03_r20 noise 10"
%%
pars(35) = par;
pars(35).noise_matrix = 15*par.noise_matrix;
pars(35).H = @u03_r20;
pars(35).name = "u03_r20 noise 15"
%%
pars(36) = par;
pars(36).noise_matrix = 20*par.noise_matrix;
pars(36).H = @u03_r20;
pars(36).name = "u03_r20 noise 20"
%%
pars(37) = par;
pars(37).noise_matrix = zeros(2);
pars(37).H = @u03_r50;
pars(37).name = "u03_r50 noise zero"
%%
pars(38) = par;
pars(38).noise_matrix = [0.3369e-5, 0 ; -0.551e-5, 0.4006e-5];
pars(38).H = @u03_r50;
pars(38).name = "u03_r50 noise HadGEM3LL"
%%
pars(39) = par;
% default nosie
pars(39).H = @u03_r50;
pars(39).name = "u03_r50 noise HadGEM3MM"
%%
pars(40) = par;
pars(40).noise_matrix = [0.386e-5, 0 ; -0.1579e-5, 0.2792e-5];
pars(40).H = @u03_r50;
pars(40).name = "u03_r50 noise CanESM5"
%%
pars(41) = par;
pars(41).noise_matrix = [0.8686e-5, 0 ; -0.0443e-5, 0.3543e-5];
pars(41).H = @u03_r50;
pars(41).name = "u03_r50 noise MPI"
%%
pars(42) = par;
pars(42).noise_matrix = 5*par.noise_matrix;
pars(42).H = @u03_r50;
pars(42).name = "u03_r50 noise 5"
%%
pars(43) = par;
pars(43).noise_matrix = 10*par.noise_matrix;
pars(43).H = @u03_r50;
pars(43).name = "u03_r50 noise 10"
%%
pars(44) = par;
pars(44).noise_matrix = 15*par.noise_matrix;
pars(44).H = @u03_r50;
pars(44).name = "u03_r50 noise 15"
%%
pars(45) = par;
pars(45).noise_matrix = 20*par.noise_matrix;
pars(45).H = @u03_r50;
pars(45).name = "u03_r50 noise 20"
%%
pars(46) = par;
pars(46).noise_matrix = zeros(2);
pars(46).H = @u03_r70;
pars(46).name = "u03_r70 noise zero"
%%
pars(47) = par;
pars(47).noise_matrix = [0.3369e-5, 0 ; -0.551e-5, 0.4006e-5];
pars(47).H = @u03_r70;
pars(47).name = "u03_r70 noise HadGEM3LL"
%%
pars(48) = par;
% default noise
pars(48).H = @u03_r70;
pars(48).name = "u03_r70 noise HadGEM3MM"
%%
pars(49) = par;
pars(49).noise_matrix = [0.386e-5, 0 ; -0.1579e-5, 0.2792e-5];
pars(49).H = @u03_r70;
pars(49).name = "u03_r70 noise CanESM5"
%%
pars(50) = par;
pars(50).noise_matrix = [0.8686e-5, 0 ; -0.0443e-5, 0.3543e-5];
pars(50).H = @u03_r70;
pars(50).name = "u03_r70 noise MPI"
%%
pars(51) = par;
pars(51).noise_matrix = 5*par.noise_matrix;
pars(51).H = @u03_r70;
pars(51).name = "u03_r70 noise 5"
%%
pars(52) = par;
pars(52).noise_matrix = 10*par.noise_matrix;
pars(52).H = @u03_r70;
pars(52).name = "u03_r70 noise 10"
%%
pars(53) = par;
pars(53).noise_matrix = 15*par.noise_matrix;
pars(53).H = @u03_r70;
pars(53).name = "u03_r70 noise 15"
%%
pars(54) = par;
pars(54).noise_matrix = 20*par.noise_matrix;
pars(54).H = @u03_r70;
pars(54).name = "u03_r70 noise 20"
%%
pars(55) = par;
pars(55).noise_matrix = zeros(2);
pars(55).H = @u03_r100;
pars(55).name = "u03_r100 noise zero"
%%
pars(56) = par;
pars(56).noise_matrix = [0.3369e-5, 0 ; -0.551e-5, 0.4006e-5];
pars(56).H = @u03_r100;
pars(56).name = "u03_r100 noise HadGEM3LL"
%%
pars(57) = par;
%default noise
pars(57).H = @u03_r100;
pars(57).name = "u03_r100 noise HadGEM3MM"
%%
pars(58) = par;
pars(58).noise_matrix = [0.386e-5, 0 ; -0.1579e-5, 0.2792e-5];
pars(58).H = @u03_r100;
pars(58).name = "u03_r100 noise CanESM5"
%%
pars(59) = par;
pars(59).noise_matrix = [0.8686e-5, 0 ; -0.0443e-5, 0.3543e-5];
pars(59).H = @u03_r100;
pars(59).name = "u03_r100 noise MPI"
%%
pars(60) = par;
pars(60).noise_matrix = 5*par.noise_matrix;
pars(60).H = @u03_r100;
pars(60).name = "u03_r100 noise 5"
%%
pars(61) = par;
pars(61).noise_matrix = 10*par.noise_matrix;
pars(61).H = @u03_r100;
pars(61).name = "u03_r100 noise 10"
%%
pars(62) = par;
pars(62).noise_matrix = 15*par.noise_matrix;
pars(62).H = @u03_r100;
pars(62).name = "u03_r100 noise 15"
%%
pars(63) = par;
pars(63).noise_matrix = 20*par.noise_matrix;
pars(63).H = @u03_r100;
pars(63).name = "u03_r100 noise 20"
%%
pars(64) = par;
pars(64).noise_matrix = zeros(2);
pars(64).H = @u03_hos;
pars(64).name = "u03_hos noise zero"
%%
pars(65) = par;
pars(65).noise_matrix = [0.3369e-5, 0 ; -0.551e-5, 0.4006e-5];
pars(65).H = @u03_hos;
pars(65).name = "u03_hos noise HadGEM3LL"
%%
pars(66) = par;
%default noise
pars(66).H = @u03_hos;
pars(66).name = "u03_hos noise HadGEM3MM"
%%
pars(67) = par;
pars(67).noise_matrix = [0.386e-5, 0 ; -0.1579e-5, 0.2792e-5];
pars(67).H = @u03_hos;
pars(67).name = "u03_hos noise CanESM5"
%%
pars(68) = par;
pars(68).noise_matrix = [0.8686e-5, 0 ; -0.0443e-5, 0.3543e-5];
pars(68).H = @u03_hos;
pars(68).name = "u03_hos noise MPI"
%%
pars(69) = par;
pars(69).noise_matrix = 5*par.noise_matrix;
pars(69).H = @u03_hos;
pars(69).name = "u03_hos noise 5"
%%
pars(70) = par;
pars(70).noise_matrix = 10*par.noise_matrix;
pars(70).H = @u03_hos;
pars(70).name = "u03_hos noise 10"
%%
pars(71) = par;
pars(71).noise_matrix = 15*par.noise_matrix;
pars(71).H = @u03_hos;
pars(71).name = "u03_hos noise 15"
%%
pars(72) = par;
pars(72).noise_matrix = 20*par.noise_matrix;
pars(72).H = @u03_hos;
pars(72).name = "u03_hos noise 20"

%% number of parameter sets to consider
nruns=size(pars,2);


% pars(2).H=@u03_r100;

v
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
%% parameterize covariance matrix
%
function QM=covmat(sigma)

% always produces pos definite matrix
var1=exp(sigma(1));
var2=exp(sigma(2));
cov=exp((sigma(1)+sigma(2))/2)*tanh(sigma(3));
QM = [var1,cov;cov,var2];%
end