%% Basic flags

% dbstop if error
SESPATH='C:\Program Files (x86)\SES Software\16.0';

RunSolver=true;
PostProcess=true;

%% This is where the fun begins

foldername=fullfile(currPath,'tralin_files');

if RunSolver
    if isfolder(foldername)
        rmdir(foldername, 's');
    end
    mkdir(foldername)
    for k=1:freq_siz
        make_tralin_input(foldername, jobid, f(k), Geom, soilFD, k);
        fname=sprintf('tr_%s_f%d.f05',jobid,k);
        cmd=sprintf('"%s\\SesBat.exe" "%s" /Run /NonBlocking /Shutdown', SESPATH, fullfile(foldername,fname));
        if k==1;pause(5);else; pause(3);end
        system(cmd);
    end
    display(sprintf('All cases done!')) %#ok<*DSPS>
end

if PostProcess
    for k=1:freq_siz
        display(sprintf('Now post-processing frequency sample %d...',k)) %#ok<*DSPS>
        fname=sprintf('tr_%s_f%d.F09',jobid,k);
        [Z_tr(:,:,k), Y_tr(:,:,k),~] = parse_tralin_file(fullfile(foldername,fname), Nph);
    end
    display(sprintf('All files done!')) %#ok<*DSPS>
    
    o=findElementByField(allZ_pul, 'VarName', 'Ztot_tr');
    allZ_pul(o).VarName='Ztot_tr';
    allZ_pul(o).Label='TRALIN';
    allZ_pul(o).Values=Z_tr;
    
    o=findElementByField(allY_pul, 'VarName', 'Ytot_tr');
    allY_pul(o).VarName='Ytot_tr';
    allY_pul(o).Label='TRALIN';
    allY_pul(o).Values=Y_tr;
    
    clear Z_tr Y_tr outputs_list fname foldername run_list k
end




