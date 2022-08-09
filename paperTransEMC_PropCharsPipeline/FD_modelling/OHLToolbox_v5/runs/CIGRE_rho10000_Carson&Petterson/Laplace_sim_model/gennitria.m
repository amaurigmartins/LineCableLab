t_sim=3.27675e-3; % Simulation time (sec)
amp=1;

e=20000000; % Sampling rate e=10000000;
data_t_sim=0:1/e:t_sim; % Vector of time

% alpha=14600; % 0.1/50 μs
% beta=86933087.6;

% alpha=288000; % 1.2/5 μs
% beta=1241666.667;

% alpha=14600; % 1.2/50 μs
% beta=2466666.667;

% alpha=3500; % 1.2/200 μs
% beta=2625000;

% alpha=291.9991617; % 250/2500 μs
% beta=16406.80018;

% d_emp_fun=amp*((exp(-alpha*data_t_sim))-(exp(-beta*data_t_sim)));
% d_emp_fun=transpose(d_emp_fun);

vo=[zeros(1,32768),amp*ones(1,32768)]; % Η μοναδιαία "παύλα"
vo=transpose(vo);

% f=1000; % Frequency
% vo=amp*(sin(2*pi*f*data_t_sim)); % Το συνημίτονο
% vo=transpose(vo);

plot(data_t_sim,vo);

csvwrite('Step.csv',vo);