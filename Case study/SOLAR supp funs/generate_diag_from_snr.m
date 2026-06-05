function d = generate_diag_from_snr(target_SNR, n, p, sigma2, q, d_min) % vary between 0.01, 0.05, 0.2

% ------------------------------------------------------------
% Generate singular values d1 >= ... >= dq such that:
%
%   SNR_X = (sum d_k^2) / (n p sigma2) = target_SNR
%
% Inputs:
%   target_SNR : desired matrix SNR
%   n, p       : dimensions
%   sigma2     : noise variance
%   q          : rank
%   d_min      : desired smallest singular value (lower bound)
%
% Output:
%   d          : q x 1 vector of diagonal entries
% ------------------------------------------------------------

% Target signal energy
E_target = target_SNR * n * p * sigma2;

% Template decreasing profile (linear)
a = (q:-1:1)';      

% Normalize template to match energy
c = sqrt(E_target / sum(a.^2));
d = c * a;

% Enforce minimum singular value if needed
if d(end) < d_min
    
    % Increase energy so smallest singular value equals d_min
    scale = d_min / d(end);
    d = d * scale;
    
    warning('Energy increased to enforce minimum singular value.');
end

end