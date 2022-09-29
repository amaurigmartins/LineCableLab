function Y012 = Y012calc(Y)
A=exp(1i*2*pi/3);
F=[1 1 1; 1 A^2 A; 1 A A^2];
mY012=(F\Y)*F;
format shortEng
Y012=diag(mY012);
re=real(Y012);
im=imag(Y012);
re(re<eps) = 0;
im(im<eps) = 0;
Y012=complex(re,im);

end