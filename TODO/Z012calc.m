function  Z012 = Z012calc(Z)

A=exp(1i*2*pi/3);
F=[1 1 1; 1 A^2 A; 1 A A^2];
mZ012=(F\Z)*F;
% Z = F*mZ012*inv(F) para converter de seq -> fase
format short
Z012=diag(mZ012);
re=real(Z012);
im=imag(Z012);
re(re<eps) = 0;
im(im<eps) = 0;
Z012=complex(re,im);


end