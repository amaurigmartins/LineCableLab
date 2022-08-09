function [g,Z,Y,Ti,freq,num]=pul_input_fun(file)

load(file)
g=g_dis;
Ti=Ti_dis;
Z=Z_dis;
Y=Y_dis;
freq=f;
num=length(freq);