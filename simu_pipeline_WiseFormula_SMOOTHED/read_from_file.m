function [out] = read_from_file(fname,what)
% Zch_mod, Ych_mod, g_dis, a_dis, vel_dis
load([fname '.mat']);

out = eval([fname '_data.' what]);

end

