function p = HadGEM3MMParameters()
%HadGEM3-GC3.1-MM parameters from Jon Baker 1xCO2
%all parameters in SI units

p.Vn = 4.192e16; p.Vt = 4.192e16; p.Vs = 13.26e16; p.Vip = 16.95e16; p.Vb = 96.76e16;
p.Fn = 0.5996e6; p.Ft = -0.7527e6; p.Fs = 1.126e6; p.Fip = -0.6728e6; p.Fb = 0;
p.Sn = 0.03488; p.St = 0.03563; p.Ss = 0.03444; p.Sip = 0.0345; p.Sb = 0.03475; p.S0 = 0.0350;
p.Ts = 5.349; p.T0 = 4.514;
p.Kn = 4.73e6; p.Ks = 7.68e6; p.Kip = 123.49e6;
%these still from FAMOUS A (no update)
p.alpha = 0.12; p.beta = 790; p.eta = 66.061e6;
%alpha and beta from text in section 2.2
%gamma still from FAMOUS A, no update
p.lambda = 2.328e7; p.gamma = 0.58; p.mu = 0;
p.C = p.Vn*p.Sn + p.Vt*p.St + p.Vs*p.Ss + p.Vip*p.Sip + p.Vb*p.Sb;
p.An = 0.9802; p.At = -0.1752; p.As = -0.2045; p.Aip = -0.6005;
end