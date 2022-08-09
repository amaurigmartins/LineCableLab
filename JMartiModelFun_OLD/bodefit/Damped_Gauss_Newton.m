%================================================
% DAMPED GAUSS NEWTON METHOD
% Authors: Eduardo Salvador Bañuelos Cabral
% José Alberto Gutiérrez Robles
% Bjørn Gustavsen
%===============================================
% Inputs
% --- Np, number of poles
% --- Ns, Number of samples
% --- Fs, function to be fitted
% --- s, complex frequency (rad/s)
% --- X, initial poles, zeros and constant term
% --- ite,iterations
% Outputs
% --- Ks, constant term
% --- Rs, residues
% --- Ps, poles
% --- Ff, objective function (deviation)

function [Ks,Rs,Ps,Ff] = Damped_Gauss_Newton(Np,Ns,Fs,s,X,ite)
a = ones(Ns,1); % Vector column
J1 = zeros(Ns,Np); % Matrix size (J1)
J2 = zeros(Ns,Np); % Matrix size (J2)
Jn = zeros(2*Ns,2*Np+1); % Matrix size (Jn)
At = zeros(2*Np+1,2*Np+1); % Matrix size (At)
en = zeros(2*Ns,1); % Vector size (en)
epn = zeros(2*Ns,1); % Vector size (epn)
Ff = zeros(ite,1); % Vector size (Ff)
Euclidian = zeros(1,1+2*Np); % Vector size
for ki = 1:ite % Damped Gauss Newton methodology
    R = X(1:Np); % Residues
    P = X(Np+1:2*Np); % Poles
    K = X(2*Np+1); % Constant term
    Fa = zeros(1,Ns); % Set the approximation to the function
    for k = 1:length(P)
        Fa = Fa + (R(k)./(s.' - P(k))).';
    end
    Fa = Fa + K; % Construct the approximation to the function
    error = (Fa.' - Fs.'); % Deviation
    for n = 1:1:Np % Loop to construct the Jacobian
        J1(1:Ns,n) = 1./(s-P(n));
        J2(1:Ns,n) = R(n)./((s-P(n)).^2);
    end
    J = [J1 J2 a]; % Jacobian
    [Xmax Ymax] = size(J); % Matrix size (J)
    Jr = real(J); % Real part of vector J
    Ji = imag(J); % Imaginary part of vector J
    er = real(error); % Real part of vector error
    ei = imag(error); % Imaginary part of vector error
    km = 1; % Counter
    for k = 2:2:2*Xmax % Interleaved
        Jn(k-1,:) = Jr(km,:);
        Jn(k,:) = Ji(km,:);
        en(k-1,1) = er(km);
        en(k,1) = ei(km);
        km = km+1;
    end
    F = ((norm(en,2)^2)); % Objective function
    Ff(ki,1) = F; % Storage the objective function must tend to zero
    [Q,R] = qr(Jn); % Matrix Q and R of Jn
    Jn = R; % It updates matrix Jn
    Gf = (Jn.'*Q.')*en; % Gradient (QR decomposition)
    Hess = Jn.'*Jn; % Hessian (QR decomposition)
    for col = 1:Ymax % Euclidian norm
        Euclidian(col) = norm(Hess(:,col),2);
        At(:,col) = Hess(:,col)./Euclidian(col);
    end
    h = (At)\-Gf; % Solution for the system (Ax = b)
    h = h./Euclidian.'; % Real solution
    stop = 1; % Variable to stop the next loop
    al = 1; % Variable to weigh the approximation
    
    while (stop == 1)
        Xp = X + al*h; % New coefficients (without updating x)
        Rp = Xp(1:Np); % Residues
        Pp = Xp(Np+1:2*Np); % Poles
        Kp = Xp(2*Np+1); % Constant term
        Fap = zeros(1,Ns); % Set the new approximation of the function
        for k = 1:length(P)
            Fap = Fap + (Rp(k)./(s.' - Pp(k))).';
        end
        Fap = Fap + Kp; %New approximation of the function
        ep = (Fap.' - Fs.'); % Deviation
        epr = real(ep); % Real part of vector ep
        epi = imag(ep); % Imaginary part of vector ep
        km = 1; % Counter
        for k = 2:2:2*Xmax % Interleaved
            epn(k-1,1) = epr(km);
            epn(k,1) = epi(km);
            km = km+1;
        end
        Fp = ((norm(epn,2)^2)); % Objective function
        % Updating process
        Erel = 1e-4*al*h.'*Gf; % Relative error to stop the while process
        if (Fp < F + Erel)
            X = X + al*h; % Final approximation
            stop = 0; % Go out the while
        else
            al = al*0.9; % Weigh to update X
            if (al < 1e-30)
                stop = 0;
            else
                stop = 1;
            end
        end
    end
end
Rs = X(1:Np); % Residues
Ps = X(Np+1:2*Np); % Poles
Ks = X(2*Np+1); % Constant term