%% ============================================================
%  SOLAR RAM/computation scalability run
%
%  Reads:
%      Output/X_n_<n>_p_<p>.mat
%      Output/Y_n_<n>_p_<p>.csv
%
%  Saves:
%      Output/Comp_summary_n_<n>_p_<p>.csv
%
%  Notes:
%  - X is kept in its stored precision, usually single.
%  - No full n x p double copy is intentionally created.
%  - Peak workspace memory is estimated from whos snapshots.
%  - True OS-level peak RAM can be filled manually if needed.
% ============================================================

clearvars %-except X
clc;
rng(123);

fprintf('\n============================================================\n');
fprintf('SOLAR scalability run\n');
fprintf('============================================================\n');

%% ============================================================
%  Settings
% ============================================================

n = 500;        % change: 500, 1000, 2000
p = 10^6;        

data_dir = 'Data';
out_dir  = 'Output';
supp_dir = 'SOLAR supp funs';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if ~exist(supp_dir, 'dir')
    error('SOLAR support-function folder not found: %s', supp_dir);
end

addpath(supp_dir);

X_mat_file = fullfile(data_dir, sprintf('X_n_%d_p_%d.mat', n, p));
Y_csv_file = fullfile(data_dir, sprintf('Y_n_%d_p_%d.csv', n, p));

if ~exist(X_mat_file, 'file')
    error('Cannot find X file: %s', X_mat_file);
end

if ~exist(Y_csv_file, 'file')
    error('Cannot find Y file: %s', Y_csv_file);
end

x_file_info = dir(X_mat_file);
x_file_size_gb = x_file_info.bytes / 1024^3;

fprintf('\nInput X file: %s\n', X_mat_file);
fprintf('X file size: %.3f GB\n', x_file_size_gb);

%% ============================================================
%  Load data
% ============================================================

fprintf('\nLoading X and Y...\n');

load_tic = tic;

tmp = load(X_mat_file);
var_names = fieldnames(tmp);

if numel(var_names) ~= 1
    disp(var_names);
    error('Expected exactly one variable in X .mat file.');
end

X = tmp.(var_names{1});
clear tmp

Y_tbl = readtable(Y_csv_file);
w_raw = Y_tbl{:, 1};

load_time_sec = toc(load_tic);

[n_check, p_check] = size(X);

if n_check ~= n || p_check ~= p
    error('Loaded X dimensions are %d x %d, expected %d x %d.', ...
        n_check, p_check, n, p);
end

if length(w_raw) ~= n
    error('Length of Y does not match n.');
end

fprintf('Loaded X: n = %d, p = %d, class = %s\n', n, p, class(X));
fprintf('Loaded Y length = %d\n', length(w_raw));
fprintf('Data loading time = %.2f seconds\n', load_time_sec);

mem_snap = [];
mem_snap = add_mem_snapshot(mem_snap, 'after_data_load');

%% ============================================================
%  Start SOLAR timing after data are available
% ============================================================

fprintf('\nStarting SOLAR computation timing now...\n');

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
    
    if mod(jj, 100000) == 0
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

method_name = 'SOLAR_scalability';

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
fprintf('Total SOLAR time after data load = %.2f seconds\n', solar_total_time_sec);

mem_snap = add_mem_snapshot(mem_snap, 'after_SOLAR_complete');

%% ============================================================
%  Compute memory summaries
% ============================================================

workspace_peak_gb = max(mem_snap.workspace_gb);
workspace_final_gb = mem_snap.workspace_gb(end);

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
summary_tbl.X_file_size_GB = x_file_size_gb;
summary_tbl.X_class = string(class(X));
summary_tbl.q_hat = q;
summary_tbl.Num_iters_requested = Num_iters;
summary_tbl.Num_iters_completed = last_iter;

summary_tbl.load_time_sec = load_time_sec;
summary_tbl.center_time_sec = center_time_sec;
summary_tbl.svd_time_sec = svd_time_sec;
summary_tbl.constants_time_sec = const_time_sec;
summary_tbl.optimization_time_sec = optimization_time_sec;
summary_tbl.total_SOLAR_time_after_load_sec = solar_total_time_sec;


summary_tbl.Peak_MATLAB_memory_footprint_GB = matlab_peak_used_gb;

summary_tbl.manual_peak_RAM_GB_to_fill = NaN;

out_file = fullfile(out_dir, sprintf('Comp_summary_n_%d_p_%d.csv', n, p));

writetable(summary_tbl, out_file);

fprintf('\nSaved computation summary:\n%s\n', out_file);

disp(summary_tbl);

%% ============================================================
%  Local helper: memory snapshot
% ============================================================

function mem_tbl = add_mem_snapshot(mem_tbl, label)

    s = whos;
    workspace_gb = sum([s.bytes]) / 1024^3;
    
    matlab_mem_used_gb = NaN;
    
    try
        m = memory;
        matlab_mem_used_gb = m.MemUsedMATLAB / 1024^3;
    catch
    end
    
    new_row = table;
    new_row.label = string(label);
    new_row.workspace_gb = workspace_gb;
    new_row.matlab_mem_used_gb = matlab_mem_used_gb;
    
    if isempty(mem_tbl)
        mem_tbl = new_row;
    else
        mem_tbl = [mem_tbl; new_row];
    end
    
    fprintf('  [Memory snapshot: %s] workspace = %.3f GB', ...
        label, workspace_gb);
    
    if ~isnan(matlab_mem_used_gb)
        fprintf(' | MATLAB used = %.3f GB', matlab_mem_used_gb);
    end
    
    fprintf('\n');
end