function out = remove_complex_poles(pol)
%%%%%%%%%%%%%%%%%%%%%%%%%%% REMOVING COMPLEX POLES
% Check https://doi.org/10.1007/s00202-019-00807-8
DELTA=.01;
cp_poles=imag(pol)>0;
cp_poles_conj=imag(pol)<0;
pol(cp_poles)=-(sqrt(real(pol(cp_poles)).^2+imag(pol(cp_poles)).^2)+DELTA); %Forcing poles to be real
pol(cp_poles_conj)=-(sqrt(real(pol(cp_poles_conj)).^2+imag(pol(cp_poles_conj)).^2)-DELTA);

out=pol;
end