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
end