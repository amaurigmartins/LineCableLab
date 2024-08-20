function hex = rgb2hex(rgb, opts)
%

%   Copyright 2023 The MathWorks, Inc.

arguments
    rgb {mustBeA(rgb,["double","single","uint8","uint16"]), mustBeNonnegative}
    opts.Shorthand (1,1) logical
end

rgbType = class(rgb);
sz = size(rgb);

% Validate RGB shape and reshape MxNx3 to a vector of RGB triplets.
if sz(end) ~= 3 && ~isequal(sz,[0 0])
    error(message('MATLAB:graphics:validatecolor:InvalidRGBMatrixShape'))
end

% Change to a vector of rgb triplets.
rgb = reshape(rgb,[],3);

if ~isempty(rgb) && ismember(rgbType, {'double' 'single'}) && (max(rgb,[],"all") > 1)
    error(message('MATLAB:graphics:validatecolor:OutOfRange'));
end

% Scale values to be on the [0,255] range.
if ismember(rgbType, {'double', 'single'})
    rgb = round(rgb * 255);
elseif isequal(rgbType, 'uint16')
    scale = 255 / 65535;
    rgb = rgb * scale;
end

% Determine whether to use shorthand 3-digit representation or standard 6-digit.
% Convert to hexadecimal representation.
if isfield(opts, 'Shorthand') && opts.Shorthand
    % Shorthand represents multiples of '#11', which is 17 in decimal.
    % Divide each rgb value by 17 and round to the closes multiple of 17.
    % dec2hex to convert each scaled rgb value to its corresponding 1-digit
    % hex value. reshape each triple of hex values to a 3-digit hex value.
    hex(:,2:4) = reshape(dec2hex(round(rgb/17)')', 3, [])';
else
    % dec2hex to convert each rgb value to its corresponding 2-digit hex
    % value. reshape each triplet of hex values to a 6-digit hex value.
    hex(:,2:7) = reshape(dec2hex(rgb',2)', 6, [])';
end
hex(:,1) = '#';
if ~isempty(hex)
    hex = string(cellstr(hex));
else
    hex = strings(0); % Create an empty string to prevent the empty char from becoming "".
end
if numel(sz) > 2
    % For NDx3, change back to an ND matrix of hex codes.
    hex = reshape(hex, sz(1:end-1));
end
end