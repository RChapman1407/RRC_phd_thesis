function p = HadGEM3LLParameters()
%HadGEM3LLL parameters from Jon Baker (Met Office), Including Med. sea
%all parameters in SI units

p.Vn = 4.633e16; p.Vt = 4.358e16; p.Vs = 11.13e16; p.Vip = 18.46e16; p.Vb = 94.98e16;
p.Fn = 0.5296e6; p.Ft = -0.8022e6; p.Fs = 0.9088e6; p.Fip = -0.8053e6; p.Fb = 0;
p.Sn = 0.03492; p.St = 0.03581; p.Ss = 0.03432; p.Sip = 0.03438; p.Sb = 0.03476; p.S0 = 0.0350;
%Salinities from FAMOUS B (Alkhayuon et.al.2019)
p.Ts = 5.831; p.T0 = 5.33;
p.Kn = 6.58e6; p.Ks = 2.68e6; p.Kip = 106.09e6;
%these still from FAMOUS A (no update)
p.alpha = 0.12; p.beta = 790; p.eta = 66.061e6;
%alpha and beta from text in section 2.2
%gamma still from FAMOUS A, no update
p.lambda = 1.696e7; p.gamma = 0.58; p.mu = 0;
p.C = p.Vn*p.Sn + p.Vt*p.St + p.Vs*p.Ss + p.Vip*p.Sip + p.Vb*p.Sb;

p.An = 0.9802; p.At = -0.1752; p.As = -0.2045; p.Aip = -0.6005;
end