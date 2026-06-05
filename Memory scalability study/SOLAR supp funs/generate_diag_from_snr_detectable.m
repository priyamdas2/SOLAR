function [d_true, info] = generate_diag_from_snr_detectable(target_SNR_X, n, p, sigma2, q_true, alpha_min, alpha_max)
% Generate singular values d_true with:
%  (i) detectable "shape" relative to noise bulk edge
% (ii) exact Frobenius SNR_X matching target_SNR_X

sigma = sqrt(sigma2);

% Noise bulk edge (rough scale for top noise singular value)
bulk_edge = sigma * (sqrt(n) + sqrt(p));

% --- Step A: detectable shape (unscaled)
if nargin < 6 || isempty(alpha_min), alpha_min = 1.15; end
if nargin < 7 || isempty(alpha_max), alpha_max = alpha_min + 0.35; end

alpha = linspace(alpha_max, alpha_min, q_true)';   % decreasing, last > 1
s = bulk_edge * alpha;                              % shape singular values

% --- Step B: scale to hit target SNR_X exactly
S2 = sum(s.^2);
c  = sqrt( target_SNR_X * (n*p*sigma2) / S2 );

d_true = c * s;

% info
info.bulk_edge = bulk_edge;
info.alpha     = alpha;
info.scale_c   = c;
info.SNR_check = sum(d_true.^2) / (n*p*sigma2);
info.d_min     = min(d_true);
info.d_max     = max(d_true);
end