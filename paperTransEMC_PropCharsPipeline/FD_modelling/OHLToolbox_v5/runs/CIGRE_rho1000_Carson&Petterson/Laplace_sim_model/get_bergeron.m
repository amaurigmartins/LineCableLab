format long;
f=750000;
%z=[Z(19,1:3); Z(19,4:6); Z(19,7:9)];
%y=[Y(19,1:3); Y(19,4:6); Y(19,7:9)];
%T=[Ti(19,1:3); Ti(19,4:6); Ti(19,7:9)];
%T=[0.59434182-1i*0.00311810 0.70710678 -0.40948887-1i*0.00280747; 0.54171728+1i*0.00620774 0 0.81523949-1i*0.00255889; 0.59434182-1i*0.00311810 -0.70710678 -0.40948887-1i*0.00280747];

z=[259.592675942351e-003 + 7.36470122572616e+000i	250.002675127155e-003 + 1.53551794054434e+000i	234.771882555142e-003 + 951.080569368971e-003i;
250.002675127155e-003 + 1.53551794054434e+000i	259.592675942351e-003 + 7.36470122572616e+000i	250.002675127155e-003 + 1.53551794054434e+000i;
234.771882555142e-003 + 951.080569368971e-003i	250.002675127155e-003 + 1.53551794054434e+000i	259.592675942351e-003 + 7.36470122572616e+000i;];

y=[0 + 36.9613797805203e-006i	0 - 5.44952137027947e-006i	0 - 2.12693755445726e-006i;
0 - 5.44952137027947e-006i	0 + 37.6424533755553e-006i	0 - 5.44952137027947e-006i;
0 - 2.12693755445726e-006i	0 - 5.44952137027947e-006i	0 + 36.9613797805203e-006i;];

[T,D]=eig(y*z);

%Zmode=inv(transpose(inv(T)))*z*T;

%Zmode=transpose(T)*z*T;
Zmode=transpose(real(T))*z*real(T);
%Ymode=inv(T)*y*transpose(inv(T));
Ymode=inv(real(T))*y*transpose(inv(real(T)));

R=[real(Zmode(1,1)) real(Zmode(2,2)) real(Zmode(3,3))]*1000

%Zs=[sqrt(abs(imag(Zmode(1,1)))/abs(imag(Ymode(1,1)))) sqrt(abs(imag(Zmode(2,2)))/abs(imag(Ymode(2,2)))) sqrt(abs(imag(Zmode(3,3)))/abs(imag(Ymode(3,3))))]
%Zs=[sqrt(Zmode(1,1)/Ymode(1,1)) sqrt(Zmode(2,2)/Ymode(2,2)) sqrt(Zmode(3,3)/Ymode(3,3))]
%Zs=sqrt(Zmode/Ymode)
Zs=sqrt(imag(Zmode)/imag(Ymode)) % lossless

%velocity=[(2*pi*f)/(sqrt(abs(imag(Zmode(1,1)))*abs(imag(Ymode(1,1))))) (2*pi*f)/(sqrt(abs(imag(Zmode(2,2)))*abs(imag(Ymode(2,2)))))]
velocity=[(2*pi*f)/(sqrt(imag(Zmode(1,1))*imag(Ymode(1,1)))) (2*pi*f)/(sqrt(imag(Zmode(2,2))*imag(Ymode(2,2)))) (2*pi*f)/(sqrt(imag(Zmode(3,3))*imag(Ymode(3,3))))]/1000 % lossless

T