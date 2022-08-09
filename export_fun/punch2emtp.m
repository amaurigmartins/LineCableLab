% close all
% clear all
% clc
% fload='wise_cg.mat';
% load(fload)
% % 
% line_length = 100;

Z=eval(Zmod_src);
Y=eval(Ymod_src);

if size(f,1) ~= 1
f=f';
end

fname=fullfile(currPath, ['emtp_export__' jobid]);

save(fname,'f','Z','Y','line_length')