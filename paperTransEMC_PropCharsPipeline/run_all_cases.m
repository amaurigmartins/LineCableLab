close all
clear
clc

i=1;
allcases(i).rho=100;
allcases(i).eps=15;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=0;

i=i+1;
allcases(i).rho=1000;
allcases(i).eps=5;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=0;

i=i+1;
allcases(i).rho=5000;
allcases(i).eps=3;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=0;

i=i+1;
allcases(i).rho=104.8;
allcases(i).eps=15;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=1;

i=i+1;
allcases(i).rho=1097;
allcases(i).eps=5;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=1;

i=i+1;
allcases(i).rho=5790.9;
allcases(i).eps=3;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=1;

i=i+1;
allcases(i).rho=100;
allcases(i).eps=15;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=9;

i=i+1;
allcases(i).rho=1000;
allcases(i).eps=5;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=9;

i=i+1;
allcases(i).rho=5000;
allcases(i).eps=3;
allcases(i).Z='Ztot_Carson';
allcases(i).Y='Ytot_Imag';
allcases(i).FD=9;

i=i+1;
allcases(i).rho=100;
allcases(i).eps=15;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=0;

i=i+1;
allcases(i).rho=1000;
allcases(i).eps=5;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=0;

i=i+1;
allcases(i).rho=5000;
allcases(i).eps=3;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=0;

i=i+1;
allcases(i).rho=104.8;
allcases(i).eps=15;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=1;

i=i+1;
allcases(i).rho=1097;
allcases(i).eps=5;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=1;

i=i+1;
allcases(i).rho=5790.9;
allcases(i).eps=3;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=1;

i=i+1;
allcases(i).rho=100;
allcases(i).eps=15;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=9;

i=i+1;
allcases(i).rho=1000;
allcases(i).eps=5;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=9;

i=i+1;
allcases(i).rho=5000;
allcases(i).eps=3;
allcases(i).Z='Ztot_Wise';
allcases(i).Y='Ytot_Wise';
allcases(i).FD=9;

numcases=i;

for i=1:numcases
    soil_rho=allcases(i).rho;%resistivity of earth
    soil_eps=allcases(i).eps; %this is not used if FD_flag != 0
    Zmod_src=allcases(i).Z;
    Ymod_src=allcases(i).Y;
    FD_flag=allcases(i).FD;
    
    clc
    fprintf('Now running case %d...\n', i);
    OHTL_run_v6
end