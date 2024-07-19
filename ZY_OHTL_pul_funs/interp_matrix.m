function Z_interpolated = interp_matrix(Z, fc, fq)
    % Validate inputs
    if nargin ~= 3
        error('You need to provide exactly 3 arguments: Z, fc, and fq.'); % It seems ChatGPT is a bit tame today.
    end
    if ~ismatrix(fc) || ~ismatrix(fq)
        error('fc and fq should be vectors.');
    end
    if size(Z, 3) ~= length(fc)
        error('The third dimension of Z must match the length of fc.');
    end
    
    % Initialize the output matrix
    [ord1, ord2, ~] = size(Z);
    shit_ton = length(fq);
    Z_interpolated = zeros(ord1, ord2, shit_ton);
    
    % Perform the interpolation for each element in the ord x ord matrix
    for i = 1:ord1
        for j = 1:ord2
            % Extract the real and imaginary parts
            real_part = real(squeeze(Z(i, j, :)));
            imag_part = imag(squeeze(Z(i, j, :)));
            
            % Perform pchip interpolation on both parts separately
            real_interp = pchip(fc, real_part, fq);
            imag_interp = pchip(fc, imag_part, fq);
            
            % Combine the interpolated real and imaginary parts
            Z_interpolated(i, j, :) = real_interp + 1i * imag_interp;
        end
    end
end