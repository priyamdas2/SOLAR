%% ============================================================
%  Generate ultra-high-dimensional simulated X and Y only
%
%  Saves:
%      Data/X_p_<p>.mat   single precision X
%      Data/Y_p_<p>.csv   response vector Y
%
%  IMPORTANT:
%  This code generates X blockwise to avoid forming full
%  X_signal or full noise matrices in memory.
% ============================================================

clear;
clc;
rng(123);

fprintf('MATLAB mem at start:\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Settings
% ============================================================

n = 1000;       % change manually: 500, 1000, 2000
p = 10^6;       % change manually if needed

q_true = 3;

target_SNR_X = 0.05;
target_SNR_w = 0.2;

sigma2_true = 0.05^2;

block_size = 50000;   % safe for 32 GB; reduce to 25000 if needed

%% ============================================================
%  Output folder
% ============================================================

data_dir = 'Data';

if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

%% ============================================================
%  Generate latent factors and response
% ============================================================

fprintf('\nGenerating latent factors...\n');

[d_true, info] = generate_diag_from_snr_detectable( ...
    target_SNR_X, n, p, sigma2_true, q_true, 1.15, 1.5);

D_true = diag(d_true);

[H_true, ~] = qr(randn(n, q_true), 0);
[V_true, ~] = qr(randn(p, q_true), 0);

beta_true = 5 * randn(q_true, 1);

tau2_true = (norm(beta_true, 2)^2) / (n * target_SNR_w);

w_signal = H_true * beta_true;
Y = w_signal + sqrt(tau2_true) * randn(n, 1);

fprintf('Generated Y.\n');

%% ============================================================
%  Generate X blockwise in single precision
% ============================================================

fprintf('\nAllocating X as single precision: n = %d, p = %d\n', n, p);

X = zeros(n, p, 'single');

fprintf('Generating X blockwise...\n');

gen_tic = tic;

for j0 = 1:block_size:p
    
    j1 = min(p, j0 + block_size - 1);
    jj = j0:j1;
    
    signal_block = single(H_true * D_true * V_true(jj, :)');
    noise_block  = sqrt(single(sigma2_true)) * randn(n, numel(jj), 'single');
    
    X(:, jj) = signal_block + noise_block;
    
    clear signal_block noise_block
    
    fprintf('  Generated columns %d to %d of %d\n', j0, j1, p);
end

gen_time_sec = toc(gen_tic);

fprintf('Finished generating X in %.2f seconds.\n', gen_time_sec);

fprintf('\nClass of X: %s\n', class(X));
fprintf('Approximate X storage in memory: %.2f GB\n', numel(X) * 4 / 1024^3);

%% ============================================================
%  Quick diagnostics without storing extra large matrices
% ============================================================

designed_SNR_X = sum(d_true.^2) / (n * p * sigma2_true);
designed_SNR_w = norm(w_signal, 2)^2 / (n * tau2_true);

fprintf('\nTarget SNR_X = %.4f | Designed SNR_X = %.4f\n', ...
    target_SNR_X, designed_SNR_X);

fprintf('Target SNR_w = %.4f | Designed SNR_w = %.4f\n', ...
    target_SNR_w, designed_SNR_w);

%% ============================================================
%  Save X and Y only
% ============================================================

fprintf('\nSaving X and Y only...\n');

x_file = fullfile(data_dir, sprintf('X_n_%d_p_%d.mat', n, p));
y_file = fullfile(data_dir, sprintf('Y_n_%d_p_%d.csv', n, p));

save(x_file, 'X', '-v7.3');

Y_tbl = table(Y, 'VariableNames', {'Y'});
writetable(Y_tbl, y_file);

fprintf('Saved X: %s\n', x_file);
fprintf('Saved Y: %s\n', y_file);

fprintf('\nMATLAB mem after saving:\n');
try
    disp(memory);
catch
end

fprintf('\nDone.\n');

%% ============================================================
%  Local function
% ============================================================

function [d_true, info] = generate_diag_from_snr_detectable( ...
    target_SNR_X, n, p, sigma2, q_true, alpha_min, alpha_max)

sigma = sqrt(sigma2);

bulk_edge = sigma * (sqrt(n) + sqrt(p));

if nargin < 6 || isempty(alpha_min)
    alpha_min = 1.15;
end

if nargin < 7 || isempty(alpha_max)
    alpha_max = alpha_min + 0.35;
end

alpha = linspace(alpha_max, alpha_min, q_true)';

s = bulk_edge * alpha;

S2 = sum(s.^2);

c = sqrt(target_SNR_X * (n * p * sigma2) / S2);

d_true = c * s;

info.bulk_edge = bulk_edge;
info.alpha = alpha;
info.scale_c = c;
info.SNR_check = sum(d_true.^2) / (n * p * sigma2);
info.d_min = min(d_true);
info.d_max = max(d_true);

end