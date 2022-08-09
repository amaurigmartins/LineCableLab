% zyph_param.m
function [Z,Y]=zyph_param(Ztmp,Ytmp,ord)

Z_real=zeros(ord,ord);
Z_imag=zeros(ord,ord);
Y_real=zeros(ord,ord);
Y_imag=zeros(ord,ord);
%%%%%%%%%%% Z Part%%%%%%%%%%%%%%
% Real Part
for m=1:1:ord
    for n=1:1:ord
        Z_real(m,n)=Ztmp(m,2*n-1);
    end
end

% Imaginary Part
for m=1:1:ord
    for n=1:1:ord
        Z_imag(m,n)=Ztmp(m,2*n);
    end
end

Z=Z_real+1i*Z_imag;
%Z=Z/1000; % Ώστε να προκύψει σε Ω/m - Μόνο όταν μελετάμε Εναέριες Γραμμές Μεταφοράς!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Y Part%%%%%%%%%%%%%%
% Real Part
for m=1:1:ord
    for n=1:1:ord
        Y_real(m,n)=Ytmp(m,2*n-1);
    end
end

% Imaginary Part
for m=1:1:ord
    for n=1:1:ord
        Y_imag(m,n)=Ytmp(m,2*n);
    end
end

Y=Y_real+1i*Y_imag;
%Y=Y/1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
