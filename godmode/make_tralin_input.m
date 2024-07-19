function [] = make_tralin_input(foldername, jobid, f, Geom, soilFD, k)


% Base values in TRALIN
rho0=1.7241e-8;
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
ins0=1.0e-14;


tr_id=sprintf('%s_f%1.0f',jobid,f);
fname=['tr_' tr_id '.f05'];
out{1,1} =     sprintf('TRALIN');
out{end+1,1} = sprintf('TEXT,MODULE,%s','LineCableLab run');
out{end+1,1} = sprintf('OPTIONS');
out{end+1,1} = sprintf('UNITS,METRIC');
out{end+1,1} = sprintf('RUN-IDENTIFICATION,%s',tr_id);
out{end+1,1} = sprintf('SEQUENCE,ON');
out{end+1,1} = sprintf('MULTILAYER,ON');
out{end+1,1} = sprintf('CONDUCTANCE,ON');
out{end+1,1} = sprintf('!KEEP_CIRCUIT_MODE');
out{end+1,1} = sprintf('PARAMETERS');
out{end+1,1} = sprintf('BASE-VALUES');
out{end+1,1} = sprintf('ACCURACY,1e-7');
out{end+1,1} = sprintf('BESSEL');
out{end+1,1} = sprintf('TERMS,300');
out{end+1,1} = sprintf('FREQUENCY,%6.6f',f);
out{end+1,1} = sprintf('INTEGRATION,AUTO-ADJUST,9');
out{end+1,1} = sprintf('STEP,1e-6');
out{end+1,1} = sprintf('UPPER-LIMIT,5.');
out{end+1,1} = sprintf('SERIES-TERMS,300');
out{end+1,1} = sprintf('SOIL-TYPE');
if isfield(soilFD,'layer')
    % will handle this when I handle this
else
    out{end+1,1} = sprintf('UNIFORM,%6.6f,%6.6f,%6.6f',1/soilFD.sigma_g_total(k),soilFD.m_g/mu0,soilFD.erg_total(k));
end
out{end+1,1} = sprintf('SYSTEM');
for i=1:size(Geom,1)
    out{end+1,1} = sprintf('GROUP,PH-%d,%6.6f,%6.6f',Geom(i,1),Geom(i,2),Geom(i,3));
    if isnan(Geom(i,8)); ext_rad=Geom(i,5); else; ext_rad=Geom(i,8); end
    if isnan(Geom(i,10)); epsr_coat=1; else; epsr_coat=Geom(i,10); end
    out{end+1,1} = sprintf('CABLE,CA-%d,%6.6f',Geom(i,1),ext_rad);
    out{end+1,1} = sprintf('CORE,CO-%d,%d,%6.6f,%6.6f,%6.6f,%6.6f,%6.6f,%6.6f',Geom(i,1), ...
        Geom(i,1),Geom(i,5),Geom(i,4),Geom(i,6)/rho0,Geom(i,7),0,epsr_coat);
end
out{end+1,1} = sprintf('ENDPROGRAM');
outstr=strjoin(out, '\n') ;
fid = fopen(fullfile(foldername,fname),'wt');
fprintf(fid, outstr);
fclose(fid);

end