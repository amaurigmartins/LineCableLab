function [freqs] = implogspace(dec_start,dec_end, P)
%IMPLOGSPACE Improved logspace

if nargin < 3
    P=40;
end

% Initialize an empty array for frequencies
freqs = [];

% Loop through each decade and generate the points
for dec = dec_start:dec_end-1
    new_freqs = logspace(dec, dec+1, P+1);
    freqs = [freqs, new_freqs(1:end-1)]; % Exclude last point to avoid overlap
end

% Add the final point to complete the sequence
freqs = [freqs, 10^dec_end];

end

