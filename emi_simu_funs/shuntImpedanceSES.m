function out = shuntImpedanceSES(f,Rout,RhoCoat,ThickCoat,RelPermCoat,L)

e0 = 8.8542e-12;

w = 2*pi*f;

sigma = 1/RhoCoat;

e = e0*RelPermCoat;

out = (1/(2*pi*L))*(1/(sigma + 1i*w*e))*log((Rout + ThickCoat)/Rout);

end
