function [h] = hosingfunc(time,Trise,Tpert,Tfall,H0,Hpert,H1)
%input t is time steps being used
%Trise,Tpert,Tfall,H0,Hpert and H1 as described in Section 3.2, Figure 3.3
%and Equation 3.2.2 of Thesis
%output h is in Sverdrups

raterise = (Hpert - H0)./Trise;
ratefall = -(Hpert - H1)./Tfall;

tst0 = time<0;
tst1 = time>0;
tst2 = time<Trise;
tst3 = time<Trise+Tpert;
tst4 = time<Trise+Tpert+Tfall;

h = H0.*tst0.*time + (raterise.*(time)).*tst1.*tst2 + Hpert.*tst1.*(1-tst2).*tst3 + (ratefall.*(time-Trise-Tpert) + Hpert).*tst1.*(1-tst3).*tst4 + H1.*tst1.*(1-tst4);

h = 1e6*h; %in Sv
end