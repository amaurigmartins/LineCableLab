Z=eval([jobid '_data.Z']);
Y=eval([jobid '_data.Y']);

if size(f,1) ~= 1
f=f';
end

fname=fullfile(currPath, ['emtp_export__' jobid]);

save(fname,'f','Z','Y','line_length')