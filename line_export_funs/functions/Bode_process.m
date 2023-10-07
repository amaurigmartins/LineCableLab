%================================================
% BODE PROCESS PROGRAM
% Authors: Eduardo Salvador Bañuelos Cabral
% José Alberto Gutiérrez Robles
% Bjørn Gustavsen
%================================================
% Inputs
% --- Fs, function to be fitted
% --- f, frequency (Hz)
% --- Ns, number of samples
% --- tol, tolerance
% Outputs
% --- P, Poles
% --- Z, Zeros
% --- k, Constant term
function [P,Z,k] = Bode_process(Fs,f,Ns,tol)
P = []; Z = []; % Initialize P and Z like empty matrices
Fsdb = 20*log10(abs(Fs)); % Function in decibels
k0db = Fsdb(1); % k0 in decibels
k0 = 10.^(k0db./20); % k0 in magnitude
Fsdb1 = ones(1,Ns).*k0db; % Constant term
% Plot of the function in decibels
% figure(99)
% semilogx(f,Fsdb,'k',f,Fsdb1,'r--','LineWidth',2)
% ylabel('Decibels'), xlabel('Frequency [Hz]'), hold on
% legend('Data','Bode')
c = 1; h = 1; r = 0; % Counters
Fsdb2 = zeros(1,Ns); % Initialize Fsdb2
while (r < Ns) % Bode process algorithm
    Fsfitdb = Fsdb1 + Fsdb2;
%     semilogx(f,Fsfitdb,'r--','LineWidth',2)
%     pause(0.5)
    
    for r = 1:1:Ns
        error = abs((Fsdb(r))-(Fsfitdb(r))); % Deviation
        if error >= (tol) % Tolerance
            % New zero
            if (Fsfitdb(r)< Fsdb(r)); Z(h)=f(r); h=h+1; break; end
            % New pole
            if (Fsfitdb(r)> Fsdb(r)); P(c)=f(r); c=c+1; break; end
        end
    end

    
    % Contribution of each zero in decibels
    Arg1 = 0;
    for kp = 1:length(Z)
        Arg1 = Arg1 + 20.*log10(abs(1+1i*f/Z(kp)));
    end
    % Contribution of each pole in decibels
    Arg2 = 0;
    for kp = 1:length(P)
        Arg2 = Arg2 - 20.*log10(abs(1+1i*f/P(kp)));
    end
    Fsdb2 = Arg1 + Arg2;
end

% Poles and zeros
P = -P*2*pi;
Z = -Z*2*pi;
% Construct the constant term
Num=1;for k=1:length(P);Num=Num*P(k);end
Den=1;for k=1:length(Z);Den=Den*Z(k);end
k = abs(Num)*k0/abs(Den);
% Plots
% figure(2)
% semilogx(f,Fsdb,'k',f,Fsfitdb,'--r','LineWidth',2);
% xlabel('Frequency [Hz]'); ylabel('Decibels')
% legend('Data','Bode')