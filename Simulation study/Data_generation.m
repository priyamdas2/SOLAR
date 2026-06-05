clear all;
rng(123);

fprintf('MATLAB mem (Windows only):\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Settings
% ============================================================
% (n,p) scenarios:  (100, 1000),  (100, 100000), (100, 1000000)
%                   (1000, 100), (1000, 10000), (1000, 1000000)

n = 1000;        % 100, 1000
p = 100000;       % 1000, 10000, 100000, 1000000

q_true = 3;        % 3, 5

target_SNR_X = 0.05;   % 0.01, 0.05, 0.2
target_SNR_w = 0.2;    % 0.05, 0.2, 0.5

sigma2_true = 0.05^2;

% Generate singular values with:
% (i) detectable shape relative to the noise bulk edge
% (ii) exact Frobenius SNR_X matching target_SNR_X
[d_true, info] = generate_diag_from_snr_detectable( ...
    target_SNR_X, n, p, sigma2_true, q_true, 1.15, 1.5);

D_true = diag(d_true);

% Supervised regression coefficients on the true latent factors
beta_true = 5 * randn(q_true, 1);

% Choose tau2_true so that the supervised signal-to-noise ratio matches target_SNR_w
tau2_true = (norm(beta_true, 2)^2) / (n * target_SNR_w);

N_rep = 10;   % <<<<< NUMBER OF REPLICATIONS

%% ============================================================
%  Output folder
% ============================================================

data_dir = 'Data';
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

%% ============================================================
%  Storage for generated data summaries
% ============================================================

snrX_store = zeros(N_rep, 1);
snrw_store = zeros(N_rep, 1);
xnorm_store = zeros(N_rep, 1);
wnorm_store = zeros(N_rep, 1);

%% ============================================================
%  Replication loop
% ============================================================

for rep = 1:N_rep

    fprintf('\n=== Replication %d / %d ===\n', rep, N_rep);

    %% ============================================================
    %  Generate data
    % ============================================================

    % Generate true left and right orthonormal latent factor matrices
    [H_true, ~] = qr(randn(n, q_true), 0);
    [V_true, ~] = qr(randn(p, q_true), 0);

    % Generate low-rank signal matrix
    X_signal = H_true * D_true * V_true';

    % Generate observed high-dimensional predictor matrix
    X = X_signal + sqrt(sigma2_true) * randn(n, p);

    Xnorm2 = norm(X, 'fro')^2;
    Xnorm  = sqrt(Xnorm2);

    % Generate supervised response from the same latent factors
    w_signal = H_true * beta_true;
    w = w_signal + sqrt(tau2_true) * randn(n, 1);

    %% ============================================================
    %  Store data-generation diagnostics
    % ============================================================

    snrX_store(rep) = norm(X_signal, 'fro')^2 / (n * p * sigma2_true);
    snrw_store(rep) = norm(w_signal, 2)^2 / (n * tau2_true);

    xnorm_store(rep) = Xnorm;
    wnorm_store(rep) = norm(w, 2);

    fprintf('Target SNR_X = %.4f | Designed SNR_X = %.4f\n', ...
        target_SNR_X, snrX_store(rep));

    fprintf('Target SNR_w = %.4f | Designed SNR_w = %.4f\n', ...
        target_SNR_w, snrw_store(rep));

        %% ============================================================
    %  Save generated data
    % ============================================================

    simdata.X = X;
    simdata.w = w;

    simdata.H_true = H_true;
    simdata.V_true = V_true;
    simdata.D_true = D_true;
    simdata.beta_true = beta_true;

    simdata.sigma2_true = sigma2_true;
    simdata.tau2_true = tau2_true;

    snrX_str = strrep(sprintf('%.3f', target_SNR_X), '.', 'p');
    snrw_str = strrep(sprintf('%.3f', target_SNR_w), '.', 'p');

    fname = sprintf( ...
        ['Simulated_data_snrX_%s_snrw_%s_qTrue_%d_' ...
         'n_%d_p_%d_rep_%02d.mat'], ...
        snrX_str, snrw_str, q_true, n, p, rep);

    save(fullfile(data_dir, fname), 'simdata', '-v7.3');

end

%% ============================================================
%  SUMMARY OF GENERATED DATA
% ============================================================

fprintf('\n========================================================\n');
fprintf('Generated data summary over %d runs\n', N_rep);
fprintf('========================================================\n');

fprintf('SNR_X: mean=%.4f (sd=%.4f)\n', ...
    mean(snrX_store), std(snrX_store));

fprintf('SNR_w: mean=%.4f (sd=%.4f)\n', ...
    mean(snrw_store), std(snrw_store));

fprintf('||X||_F: mean=%.4f (sd=%.4f)\n', ...
    mean(xnorm_store), std(xnorm_store));

fprintf('||w||_2: mean=%.4f (sd=%.4f)\n', ...
    mean(wnorm_store), std(wnorm_store));

fprintf('========================================================\n');

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