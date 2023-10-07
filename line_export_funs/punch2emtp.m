function [] = punch2emtp(f,Z,Y,line_length,jobid,currPath)

if size(f,1) ~= 1
    f=f';
end

if length(f)<40
    warning('LineCable_Data MAT-file not created. You need at least 40 log-spaced frequency samples to prevent this becoming a wordly-recognized failure. Increase your sample count or maybe reconsider your life choices.');
    return;
end

fname=fullfile(currPath, ['emtp_export__' jobid]);
save(fname,'f','Z','Y','line_length')

end