close all
clear all
clc
fload='wise_cg.mat';
load(fload)

line_length = 100;

Z=eval(Zmod_src);
Y=eval(Ymod_src);

f=f';

fname=['E:\Users\Amauri\Documents\EMTP\EMTP_test\' fload];

save(fname,'f','Z','Y','line_length')