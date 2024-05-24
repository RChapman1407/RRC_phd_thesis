function p = FamousBParameters()
%FAMOUS B parameters 1xCO2

p.Vn = 3.261e16; p.Vt = 7.777e16; p.Vs = 8.897e16; p.Vip = 2.2020e17; p.Vb = 8.6490e17;
p.Fn =  3.84e5; p.Ft = -7.23e5; p.Fs = 1.078e6; p.Fip = -7.38e5; p.Fb = 0;
p.Sn = 0.034912; p.St = 0.035435; p.Ss = 0.034427; p.Sip = 0.034668; p.Sb = 0.034538; p.S0 = 0.035;
p.Ts = 4.773; p.T0 = 2.65; p.Tn = 6.679;
p.Kn = 5.456e6; p.Ks = 5.447e6; p.Kip = 96.817e6;
p.alpha = 0.12; p.beta = 790; p.eta = 74.492e6;
p.lambda = 2.79e7; p.gamma = 0.39; p.mu = 5.5e-8;
p.C = p.Vn*p.Sn + p.Vt*p.St + p.Vs*p.Ss + p.Vip*p.Sip + p.Vb*p.Sb;
%Hosing weightings for five box model
p.An = 0.07; p.At = 0.752; p.As = -0.257; p.Aip = -0.565;
end