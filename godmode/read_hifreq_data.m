clear all
close all

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [9, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data_above = readtable("HI_Aboveground_MetalGPR.csv", opts);
data_under = readtable("HI_Underground_MetalGPR.csv", opts);

%% Convert to output type
data_above = table2array(data_above);
data_under = table2array(data_under);

%% Clear temporary variables
clear opts

%% Extract distances and voltages
f=[1 50 500 1000 10000 100000 500000 1000000 10000000 100000000];
start_rows=find(data_above(:,4)==1);
end_rows=find(data_above(:,4)==2000);

d=data_above(start_rows(1):end_rows(1),4);

nfreq=length(f);

for k=1:nfreq
    v_real_above=data_above(start_rows(k):end_rows(k),12);
    v_imag_above=data_above(start_rows(k):end_rows(k),13);
    v_complx_above=complex(v_real_above,v_imag_above);

    v_real_under=data_under(start_rows(k):end_rows(k),12);
    v_imag_under=data_under(start_rows(k):end_rows(k),13);
    v_complx_under=complex(v_real_under,v_imag_under);

    figure;hold all
    plot(d(1:end-1),abs(v_complx_above(1:end-1)),'-');
    plot(d(1:end-1),abs(v_complx_under(1:end-1)),':');
    legend('Aboveground','Underground')
    xlabel('Distance along pipeline [m]')
    ylabel('Induced voltage [V]')
    title(sprintf('Frequency = %d Hz',f(k)))
    grid on
    box on
%     pause() []

end





for k=1:nfreq

    v_real_above=data_above(end_rows(k)-1,12);
    v_imag_above=data_above(end_rows(k)-1,13);
    v_complx_above_atend(k)=complex(v_real_above,v_imag_above);


    v_real_under=data_under(end_rows(k)-1,12);
    v_imag_under=data_under(end_rows(k)-1,13);
    v_complx_under_atend(k)=complex(v_real_under,v_imag_under);

end

figure;
semilogx(f,abs(v_complx_above_atend(1:end)),'-');hold all
semilogx(f,abs(v_complx_under_atend(1:end)),':');
legend('Aboveground','Underground')
xlabel('Frequency [Hz]')
ylabel('Induced voltage [V]')
title('Response at pipeline end')
grid on
box on