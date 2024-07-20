%% Basic flags

% dbstop if error
SESPATH='C:\Program Files (x86)\SES Software\16.0';

RUNME=true;

%% This is where the fun begins

if RUNME
    foldername=fullfile(currPath,'tralin_files');

    if isfolder(foldername)
        rmdir(foldername, 's');
    end
    mkdir(foldername)

    for k=1:freq_siz
        make_tralin_input(foldername, jobid, f(k), Geom, soilFD, k);
    end
    run_list=dir(fullfile(foldername,'*.f05'));
    for k=1:numel(run_list)
        fname=run_list(k).name;
        cmd=sprintf('"%s\\SesBat.exe" "%s" /Run /NonBlocking /Shutdown', SESPATH, fullfile(foldername,fname));
        if k==1;pause(5);else; pause(3);end
        system(cmd);
    end
    display(sprintf('All cases done!',k)) %#ok<DSPSY>

    outputs_list=dir(fullfile(foldername,'*.F09'));
    for k=1:numel(outputs_list)
        display(sprintf('Now post-processing frequency sample %d...',k)) %#ok<DSPSY>
        fname=outputs_list(k).name;
        [Z_tr(:,:,k), Y_tr(:,:,k),~] = parse_tralin_file(fullfile(foldername,fname), Nph);
    end

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
