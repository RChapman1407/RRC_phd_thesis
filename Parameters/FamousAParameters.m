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
end