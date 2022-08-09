clear
close all

load('Mode1phase.mat')

[out] = punchJMartiCard(NP, polesZc, resZc, polesA1, resA1, Ti);   
        
% This is for debug purpose:
% fid = fopen('tmp_lcc.atp','wt');
% fprintf(fid, out);
% fclose(fid);        
        
        
        
        
        
        
