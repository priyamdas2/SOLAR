%% ============================================================
%  SOLAR scalability run with in-memory data generation
%
%  Setting:
%      n = 1000
%      p = 10^7
%
%  Workflow:
%      1. Generate X and Y in memory
%      2. Clear all generation-only objects except X and Y
%      3. Start SOLAR timing/memory tracking
%      4. Save computation summary only
%
%  Saves:
%      Output/Comp_summary_n_1000_p_10000000.csv
%
%  Notes:
%  - X is stored as single precision.
%  - No full n x p double copy is intentionally created.
%  - MATLAB memory tracking may return NaN on non-Windows systems.
% ============================================================

clearvars;
clc;
rng(123);

fprintf('\n============================================================\n');
fprintf('SOLAR scalability run: generate data in memory, then fit SOLAR\n');
fprintf('============================================================\n');

%% ============================================================
%  Settings
% ============================================================

n = 1000;
p = 10^7;

q_true = 3;

target_SNR_X = 0.05;
target_SNR_w = 0.2;

sigma2_true = 0.05^2;

block_size = 50000;

out_dir  = 'Output';
supp_dir = 'SOLAR supp funs';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if ~exist(supp_dir, 'dir')
    error('SOLAR support-function folder not found: %s', supp_dir);
end

addpath(supp_dir);

fprintf('\nProblem size:\n');
fprintf('n = %d\n', n);
fprintf('p = %d\n', p);
fprintf('q_true = %d\n', q_true);
fprintf('Target SNR_X = %.4f\n', target_SNR_X);
fprintf('Target SNR_w = %.4f\n', target_SNR_w);
fprintf('Block size = %d columns\n', block_size);

fprintf('\nMachine information:\n');
try
    m0 = memory;
    fprintf('Physical RAM: %.2f GB\n', ...
        m0.PhysicalMemory.Total / 1024^3);
    fprintf('Available RAM at start: %.2f GB\n', ...
        m0.PhysicalMemory.Available / 1024^3);
    fprintf('MATLAB memory used at start: %.2f GB\n', ...
        m0.MemUsedMATLAB / 1024^3);
catch
    fprintf('memory command not available on this system.\n');
end

fprintf('MATLAB max computational threads: %d\n', maxNumCompThreads);

%% ============================================================
%  Generate latent factors and response
% ============================================================

fprintf('\n============================================================\n');
fprintf('Generating simulated data in memory\n');
fprintf('============================================================\n');

gen_tic = tic;

[d_true, info] = generate_diag_from_snr_detectable( ...
    target_SNR_X, n, p, sigma2_true, q_true, 1.15, 1.5);

D_gen = single(diag(d_true));

% Left latent factors: small, safe to QR in double
[H_tmp, ~] = qr(randn(n, q_true), 0);
H_gen = single(H_tmp);
clear H_tmp

% Right latent factors: memory-conscious orthonormalization
fprintf('\nGenerating V factor matrix as single precision...\n');

V_raw = randn(p, q_true, 'single');
Gv = double(V_raw' * V_raw);
Rv = chol(Gv);

V_gen = V_raw / single(Rv);

clear V_raw Gv Rv

beta_true = single(5 * randn(q_true, 1));

tau2_true = double(norm(beta_true, 2)^2) / (n * target_SNR_w);

w_signal = double(H_gen * beta_true);
Y = w_signal + sqrt(tau2_true) * randn(n, 1);

w_raw = Y;

fprintf('Generated Y.\n');

%% ============================================================
%  Generate X blockwise in single precision
% ============================================================

fprintf('\nAllocating X as single precision: n = %d, p = %d\n', n, p);

X = zeros(n, p, 'single');

fprintf('Approximate X storage in memory: %.2f GB\n', ...
    numel(X) * 4 / 1024^3);

fprintf('\nGenerating X blockwise...\n');

for j0 = 1:block_size:p
    
    j1 = min(p, j0 + block_size - 1);
    jj = j0:j1;
    
    signal_block = H_gen * D_gen * V_gen(jj, :)';
    noise_block  = sqrt(single(sigma2_true)) * ...
        randn(n, numel(jj), 'single');
    
    X(:, jj) = signal_block + noise_block;
    
    clear signal_block noise_block
    
    fprintf('  Generated columns %d to %d of %d\n', j0, j1, p);
end

generation_time_sec = toc(gen_tic);

fprintf('\nFinished generating X and Y in %.2f seconds.\n', generation_time_sec);
fprintf('Class of X: %s\n', class(X));

designed_SNR_X = sum(d_true.^2) / (n * p * sigma2_true);
designed_SNR_w = norm(w_signal, 2)^2 / (n * tau2_true);

fprintf('\nTarget SNR_X = %.4f | Designed SNR_X = %.4f\n', ...
    target_SNR_X, designed_SNR_X);

fprintf('Target SNR_w = %.4f | Designed SNR_w = %.4f\n', ...
    target_SNR_w, designed_SNR_w);

%% ============================================================
%  Clear generation-only objects before SOLAR benchmark
% ============================================================

fprintf('\nClearing generation-only objects before SOLAR fitting...\n');

clear H_gen V_gen D_gen beta_true d_true info w_signal
clear designed_SNR_X designed_SNR_w

fprintf('Only X, Y/w_raw, and scalar settings are retained for SOLAR.\n');

fprintf('\nMemory after data generation and cleanup:\n');
try
    m1 = memory;
    fprintf('MATLAB memory used: %.2f GB\n', ...
        m1.MemUsedMATLAB / 1024^3);
    fprintf('Available RAM: %.2f GB\n', ...
        m1.PhysicalMemory.Available / 1024^3);
catch
    fprintf('memory command not available on this system.\n');
end

%% ============================================================
%  Start SOLAR timing after data are available
% ============================================================

fprintf('\n============================================================\n');
fprintf('Starting SOLAR computation timing now\n');
fprintf('============================================================\n');

mem_snap = [];
mem_snap = add_mem_snapshot(mem_snap, 'before_SOLAR_start');

solar_total_tic = tic;

%% ============================================================
%  Center X and Y
% ============================================================

fprintf('\nCentering X and Y in place...\n');

center_tic = tic;

X_mean = mean(X, 1);
w_mean = mean(w_raw);

for jj = 1:p
    X(:, jj) = X(:, jj) - X_mean(jj);
    
    if mod(jj, 1000000) == 0
        fprintf('  Centered %d / %d columns\n', jj, p);
    end
end

w = w_raw - w_mean;

center_time_sec = toc(center_tic);

fprintf('Centering completed in %.2f seconds.\n', center_time_sec);

mem_snap = add_mem_snapshot(mem_snap, 'after_centering');

Xnorm2 = norm(X, 'fro')^2;
Xnorm  = sqrt(Xnorm2);

mem_snap = add_mem_snapshot(mem_snap, 'after_Xnorm');

%% ============================================================
%  SOLAR settings
% ============================================================

method_name = 'SOLAR_scalability_inmemory';

print_every = 2;

q_init = 10;
q_min  = 1;
q_max  = 15;

Num_iters = 500;

perform_gram_svd = 1;

T0 = 1.0;
T_min = 1e-2;
cool_pow = 0.7;

rho2 = 0.5^2;
g2   = 25.0;

c_kappa = 0.1;
zeta = 2/3;

kappa = c_kappa * ...
    (((sqrt(n) + sqrt(p))^2) / (n + p))^zeta;

q_max_eff  = min([q_max, n - 1, p]);
q_init_eff = min(q_init, q_max_eff);

fprintf('\nSOLAR settings:\n');
fprintf('method_name = %s\n', method_name);
fprintf('q_init = %d | q_min = %d | q_max = %d\n', q_init, q_min, q_max);
fprintf('q_init_eff = %d | q_max_eff = %d\n', q_init_eff, q_max_eff);
fprintf('Num_iters = %d | print_every = %d\n', Num_iters, print_every);
fprintf('perform_gram_svd = %d\n', perform_gram_svd);
fprintf('kappa = %.6f\n', kappa);

%% ============================================================
%  Precompute SVD using Gram matrix
% ============================================================

fprintf('\nPrecomputing leading SVD components...\n');

svd_tic = tic;

if perform_gram_svd == 0
    
    [U_full, S_full, V_full] = svds(X, q_max_eff);
    
else
    
    fprintf('Computing G = X * X''...\n');
    G = double(X * X');
    
    mem_snap = add_mem_snapshot(mem_snap, 'after_Gram_matrix');
    
    fprintf('Computing eigs(G)...\n');
    [U_full, S2] = eigs(G, q_max_eff);
    
    clear G
    
    s = sqrt(diag(S2));
    S_full = diag(s);
    
    clear S2
    
    fprintf('Computing V_full = X'' * U_full...\n');
    V_full = X' * U_full;
    
    for k = 1:q_max_eff
        V_full(:, k) = V_full(:, k) / s(k);
    end
end

svd_time_sec = toc(svd_tic);

fprintf('SVD/precomputation completed in %.2f seconds.\n', svd_time_sec);

mem_snap = add_mem_snapshot(mem_snap, 'after_SVD_precompute');

%% ============================================================
%  Constants
% ============================================================

fprintf('\nEstimating variance constants...\n');

const_tic = tic;

sigma2 = estimate_sigma2_resid_baseline(X, 2);
tau2   = var(w, 1);

prec_lik   = 1 / sigma2;
prec_prior = 1 / rho2;
prec_post  = prec_lik + prec_prior;

tau_inv2 = 1 / tau2;
g_inv2   = 1 / g2;
prec_b   = tau_inv2 + g_inv2;
beta_shrink = tau_inv2 / prec_b;

const_time_sec = toc(const_tic);

fprintf('sigma2_est = %.6f\n', sigma2);
fprintf('tau2_est   = %.6f\n', tau2);

mem_snap = add_mem_snapshot(mem_snap, 'after_constants');

%% ============================================================
%  Initialize
% ============================================================

q = q_init_eff;

H = U_full(:, 1:q);
V = V_full(:, 1:q);
D = (prec_lik / prec_post) * S_full(1:q, 1:q);
d = diag(D);

beta = beta_shrink * (H' * w);

obj_best  = -Inf;
q_best    = q;
H_best    = H;
V_best    = V;
D_best    = D;
beta_best = beta;

recon_best = Inf;

W_stop      = min(250, floor(0.5 * Num_iters));
epsJ_stop   = 1e-6;
epsAcc_stop = 0.01;
etaT_stop   = 0.10;

best_obj_hist = nan(Num_iters, 1);
best_q_hist   = nan(Num_iters, 1);
q_hist         = nan(Num_iters, 1);
q_accept_hist  = zeros(Num_iters, 1);

mem_snap = add_mem_snapshot(mem_snap, 'after_initialization');

%% ============================================================
%  Iterative MAP + trans-dimensional search
% ============================================================

fprintf('\nStarting SOLAR optimization...\n');

opt_tic = tic;
last_iter = Num_iters;

for iter = 1:Num_iters
    
    T = max(T_min, T0 * (iter^(-cool_pow)));
    
    Fv = X' * (H * D);
    [U1, ~, V1] = svd(Fv, 'econ');
    V = U1 * V1';
    
    Fh = (w * beta') / tau2 + (X * V * D) / sigma2;
    [U2, ~, V2] = svd(Fh, 'econ');
    H = U2 * V2';
    
    beta = beta_shrink * (H' * w);
    
    u = H' * w;
    
    if norm(u) > 1e-12
        
        a = u / norm(u);
        e1 = zeros(q, 1);
        e1(1) = 1;
        
        vv = a - e1;
        
        if norm(vv) > 1e-12
            vv = vv / norm(vv);
            Rw = eye(q) - 2 * (vv * vv');
        else
            Rw = eye(q);
        end
        
        H = H * Rw';
        V = V * Rw';
        beta = Rw * beta;
    end
    
    R = H' * X * V;
    [Ur, Sr, Vr] = svd(R, 'econ');
    
    H = H * Ur;
    V = V * Vr;
    
    d = (prec_lik / prec_post) * diag(Sr);
    D = diag(d);
    
    beta = beta_shrink * (H' * w);
    
    XV  = X * V;
    XVD = XV * D;
    
    residX2 = Xnorm2 + norm(D, 'fro')^2 - 2 * trace(H' * XVD);
    
    like_curr = ...
        -0.5 / sigma2 * residX2 ...
        -0.5 / tau2   * norm(w - H * beta)^2 ...
        -0.5 / rho2   * norm(d, 2)^2 ...
        -0.5 / g2     * norm(beta, 2)^2;
    
    dfq = q * (n + p - 2*q) + 2*q;
    pen_curr = kappa * 0.5 * log(n * p) * dfq;
    
    obj_curr = like_curr - pen_curr;
    
    if obj_curr > obj_best
        
        obj_best  = obj_curr;
        q_best    = q;
        H_best    = H;
        V_best    = V;
        D_best    = D;
        beta_best = beta;
        
        recon_best = sqrt(max(residX2, 0)) / Xnorm;
    end
    
    q_hist(iter) = q;
    
    if rand < 0.5
        
        if q == q_min
            q_prop = q + 1;
        elseif q == q_max_eff
            q_prop = q - 1;
        else
            q_prop = q + (2 * (rand < 0.5) - 1);
        end
        
        H_prop = U_full(:, 1:q_prop);
        V_prop = V_full(:, 1:q_prop);
        D_prop = (prec_lik / prec_post) * S_full(1:q_prop, 1:q_prop);
        d_prop = diag(D_prop);
        
        beta_prop = beta_shrink * (H_prop' * w);
        
        XVp  = X * V_prop;
        XVDp = XVp * D_prop;
        
        residXp2 = Xnorm2 + norm(D_prop, 'fro')^2 - ...
            2 * trace(H_prop' * XVDp);
        
        like_prop = ...
            -0.5 / sigma2 * residXp2 ...
            -0.5 / tau2   * norm(w - H_prop * beta_prop)^2 ...
            -0.5 / rho2   * norm(d_prop, 2)^2 ...
            -0.5 / g2     * norm(beta_prop, 2)^2;
        
        dfq_prop = q_prop * (n + p - 2*q_prop) + 2*q_prop;
        pen_prop = kappa * 0.5 * log(n * p) * dfq_prop;
        
        obj_prop = like_prop - pen_prop;
        dObj = obj_prop - obj_curr;
        
        if dObj >= 0 || rand < exp(dObj / T)
            
            q = q_prop;
            H = H_prop;
            V = V_prop;
            D = D_prop;
            d = d_prop;
            beta = beta_prop;
            
            q_accept_hist(iter) = 1;
            q_hist(iter) = q_prop;
        end
    end
    
    if mod(iter, print_every) == 0
        
        recon_iter = sqrt(max(residX2, 0)) / Xnorm;
        
        fprintf(['Iter %5d | q=%2d | recon=%.4f | obj=%.3e | ', ...
                 'best(q=%2d,recon=%.4f,obj=%.3e)\n'], ...
            iter, q, recon_iter, obj_curr, ...
            q_best, recon_best, obj_best);
        
        mem_snap = add_mem_snapshot(mem_snap, sprintf('iter_%05d', iter));
    end
    
    best_obj_hist(iter) = obj_best;
    best_q_hist(iter)   = q_best;
    
    if iter >= W_stop
        
        at_floor = (T <= (1 + etaT_stop) * T_min);
        
        q_win = best_q_hist(iter-W_stop+1:iter);
        q_stable = all(q_win == q_win(end));
        
        J_now  = best_obj_hist(iter);
        J_prev = best_obj_hist(iter-W_stop+1);
        
        rel_improve = abs(J_now - J_prev) / (1 + abs(J_now));
        J_stable = (rel_improve <= epsJ_stop);
        
        acc_rate = mean(q_accept_hist(iter-W_stop+1:iter));
        acc_small = (acc_rate <= epsAcc_stop);
        
        if at_floor && q_stable && J_stable && acc_small
            
            fprintf(['Early stop at iter=%d: T=%.2e, q_best=%d stable, ', ...
                     'rel_improve=%.2e, acc_rate=%.3f\n'], ...
                iter, T, q_best, rel_improve, acc_rate);
            
            last_iter = iter;
            break;
        end
    end
end

optimization_time_sec = toc(opt_tic);
solar_total_time_sec = toc(solar_total_tic);

q = q_best;

fprintf('\nSOLAR optimization complete.\n');
fprintf('Selected rank q_hat = %d\n', q);
fprintf('Optimization time = %.2f seconds\n', optimization_time_sec);
fprintf('Total SOLAR time after generation cleanup = %.2f seconds\n', ...
    solar_total_time_sec);

mem_snap = add_mem_snapshot(mem_snap, 'after_SOLAR_complete');

%% ============================================================
%  Memory summary
% ============================================================

if any(~isnan(mem_snap.matlab_mem_used_gb))
    matlab_peak_used_gb = max(mem_snap.matlab_mem_used_gb, [], 'omitnan');
else
    matlab_peak_used_gb = NaN;
end

%% ============================================================
%  Save computation summary
% ============================================================

summary_tbl = table;

summary_tbl.n = n;
summary_tbl.p = p;
summary_tbl.X_storage_GB_single = numel(X) * 4 / 1024^3;
summary_tbl.X_class = string(class(X));

summary_tbl.method = string(method_name);
summary_tbl.q_hat = q;
summary_tbl.Num_iters_requested = Num_iters;
summary_tbl.Num_iters_completed = last_iter;

summary_tbl.generation_time_sec = generation_time_sec;
summary_tbl.center_time_sec = center_time_sec;
summary_tbl.svd_time_sec = svd_time_sec;
summary_tbl.constants_time_sec = const_time_sec;
summary_tbl.optimization_time_sec = optimization_time_sec;
summary_tbl.total_SOLAR_time_after_generation_cleanup_sec = solar_total_time_sec;

summary_tbl.Peak_MATLAB_memory_footprint_GB = matlab_peak_used_gb;
summary_tbl.manual_peak_RAM_GB_to_fill = NaN;

out_file = fullfile(out_dir, sprintf('Comp_summary_n_%d_p_%d.csv', n, p));

writetable(summary_tbl, out_file);

fprintf('\nSaved computation summary:\n%s\n', out_file);

disp(summary_tbl);

fprintf('\nDone.\n');

%% ============================================================
%  Local helper: memory snapshot
% ============================================================

function mem_tbl = add_mem_snapshot(mem_tbl, label)
    
    matlab_mem_used_gb = NaN;
    
    try
        m = memory;
        matlab_mem_used_gb = m.MemUsedMATLAB / 1024^3;
    catch
    end
    
    new_row = table;
    new_row.label = string(label);
    new_row.matlab_mem_used_gb = matlab_mem_used_gb;
    
    if isempty(mem_tbl)
        mem_tbl = new_row;
    else
        mem_tbl = [mem_tbl; new_row];
    end
    
    if ~isnan(matlab_mem_used_gb)
        fprintf('  [Memory snapshot: %s] MATLAB used = %.3f GB\n', ...
            label, matlab_mem_used_gb);
    else
        fprintf('  [Memory snapshot: %s] MATLAB memory unavailable\n', ...
            label);
    end
end

%% ============================================================
%  Local function: singular values
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