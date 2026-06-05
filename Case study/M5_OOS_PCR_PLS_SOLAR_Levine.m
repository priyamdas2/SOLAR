%% ============================================================
%  Real-data OOS validation: Levine residual acceleration
%  Methods: PCR, PLS, SOLAR
%
%  Same 80/20 split used for all three methods.
%
%  Outputs:
%      Output_OOS/
%        Output_from_M5_Levine_OOS_PCR.csv
%        Output_from_M5_Levine_OOS_PLS.csv
%        Output_from_M5_Levine_OOS_SOLAR.csv
%        Output_from_M5_Levine_OOS_combined_metrics.csv
%
%  Metrics saved:
%      RMSE, correlation, R2
% ============================================================

clearvars -except X
clc;
rng(123);

fprintf('\n============================================================\n');
fprintf('Real-data OOS validation: Levine residual acceleration\n');
fprintf('Methods: PCR, PLS, SOLAR\n');
fprintf('============================================================\n');

fprintf('\nMATLAB memory status at start:\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Directories and files
% ============================================================

data_dir = 'Age methylation data';
supp_dir = 'SOLAR supp funs';

out_dir = 'Output_OOS';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if ~exist(supp_dir, 'dir')
    error('SOLAR support-function folder not found: %s', supp_dir);
end

addpath(supp_dir);

X_mat_file = fullfile(data_dir, 'X_matrix_final.mat');

pheno_file = fullfile( ...
    data_dir, ...
    'Pheno_for_MATLAB_with_epi_ages.csv' ...
);

if ~exist(X_mat_file, 'file')
    error('X matrix file not found: %s', X_mat_file);
end

if ~exist(pheno_file, 'file')
    error('Phenotype file not found: %s', pheno_file);
end

%% ============================================================
%  Load phenotype
% ============================================================

fprintf('\n[Step 1] Reading phenotype file...\n');

pheno = readtable( ...
    pheno_file, ...
    'VariableNamingRule', 'preserve' ...
);

outcome_name = 'AgeAccelResidual_Levine';

if ~ismember(outcome_name, pheno.Properties.VariableNames)
    error('%s column not found in phenotype file.', outcome_name);
end

w_raw = pheno.(outcome_name);

fprintf('Outcome used: %s\n', outcome_name);
fprintf('Phenotype rows: %d\n', height(pheno));

%% ============================================================
%  Load methylation matrix X
% ============================================================

fprintf('\n[Step 2] Loading methylation matrix X...\n');

if exist('X', 'var')
    fprintf('Using X already available in memory.\n');
else
    fprintf('Loading X_matrix_final.mat. This may take time...\n');
    load_tic = tic;
    
    tmp = load(X_mat_file);
    var_names = fieldnames(tmp);
    
    if numel(var_names) ~= 1
        fprintf('\nVariables found in X_matrix_final.mat:\n');
        disp(var_names);
        error('Expected exactly one numeric matrix inside X_matrix_final.mat.');
    end
    
    X = tmp.(var_names{1});
    clear tmp
    
    fprintf('Finished loading X in %.2f seconds.\n', toc(load_tic));
end

[n_full, p] = size(X);

fprintf('\nLoaded methylation matrix:\n');
fprintf('n = %d samples\n', n_full);
fprintf('p = %d CpGs\n', p);
fprintf('Class of X: %s\n', class(X));

if height(pheno) ~= n_full
    error('Phenotype rows (%d) do not match X rows (%d).', height(pheno), n_full);
end

%% ============================================================
%  Remove missing outcomes
% ============================================================

fprintf('\n[Step 3] Filtering missing Levine residual outcomes...\n');

valid_idx = ~isnan(w_raw);

if any(~valid_idx)
    fprintf('Removing %d samples with missing Levine residual acceleration.\n', sum(~valid_idx));
    X = X(valid_idx, :);
    w_raw = w_raw(valid_idx);
    pheno = pheno(valid_idx, :);
else
    fprintf('No missing Levine residual outcomes detected.\n');
end

[n, p] = size(X);

fprintf('Analysis sample size after filtering: n = %d\n', n);
fprintf('CpGs retained: p = %d\n', p);

%% ============================================================
%  Fixed 80/20 train/test split
% ============================================================

fprintf('\n[Step 4] Creating fixed 80/20 train/test split...\n');

split_seed = 2026;
rng(split_seed);

train_frac = 0.80;

idx = randperm(n);
n_train = round(train_frac * n);

train_idx = idx(1:n_train);
test_idx  = idx(n_train+1:end);

w_train = w_raw(train_idx);
w_test  = w_raw(test_idx);

fprintf('Split seed = %d\n', split_seed);
fprintf('Training samples = %d\n', length(train_idx));
fprintf('Test samples     = %d\n', length(test_idx));

split_tbl = table;
split_tbl.sample_index = (1:n)';
split_tbl.is_training = false(n, 1);
split_tbl.is_test = false(n, 1);
split_tbl.is_training(train_idx) = true;
split_tbl.is_test(test_idx) = true;

writetable( ...
    split_tbl, ...
    fullfile(out_dir, 'Output_from_M5_Levine_OOS_train_test_split.csv') ...
);

fprintf('Saved train/test split file.\n');

%% ============================================================
%  Common settings
% ============================================================

q_grid = 5:5;

all_results = table;

%% ============================================================
%  Method 1: PCR
% ============================================================

fprintf('\n============================================================\n');
fprintf('[Method 1/3] PCR started\n');
fprintf('============================================================\n');

method_tic = tic;
method_name = 'PCR';

fprintf('Preparing training matrix for PCR...\n');

X_train = double(X(train_idx, :));
X_mean = mean(X_train, 1);
w_mean = mean(w_train);

fprintf('Centering PCR training matrix...\n');
X_train_c = X_train;
clear X_train

for jj = 1:p
    X_train_c(:, jj) = X_train_c(:, jj) - X_mean(jj);
    
    if mod(jj, 100000) == 0
        fprintf('  PCR centering progress: %d / %d CpGs\n', jj, p);
    end
end

w_train_c = w_train - w_mean;

fprintf('Computing truncated SVD for PCR with q_max = %d...\n', max(q_grid));
svd_tic = tic;

[U_pcr, S_pcr, V_pcr] = svds(X_train_c, max(q_grid));

fprintf('PCR SVD completed in %.2f seconds.\n', toc(svd_tic));

Scores_train = U_pcr * S_pcr;

fprintf('Selecting PCR rank by training BIC...\n');

bic_vals = nan(length(q_grid), 1);

for ii = 1:length(q_grid)
    
    q = q_grid(ii);
    
    Tq = Scores_train(:, 1:q);
    coef_q = Tq \ w_train_c;
    
    w_hat_q_centered = Tq * coef_q;
    w_hat_q = w_hat_q_centered + w_mean;
    
    rss_q = sum((w_train - w_hat_q).^2);
    
    bic_vals(ii) = n_train * log(rss_q / n_train) + q * log(n_train);
    
    fprintf('  PCR q = %2d | BIC = %.4f\n', q, bic_vals(ii));
end

[~, best_idx] = min(bic_vals);
q_hat_pcr = q_grid(best_idx);

fprintf('Selected PCR rank q_hat = %d\n', q_hat_pcr);

fprintf('Preparing test matrix for PCR prediction...\n');

X_test = double(X(test_idx, :));

fprintf('Centering PCR test matrix...\n');

for jj = 1:p
    X_test(:, jj) = X_test(:, jj) - X_mean(jj);
    
    if mod(jj, 100000) == 0
        fprintf('  PCR test-centering progress: %d / %d CpGs\n', jj, p);
    end
end

T_train_hat = Scores_train(:, 1:q_hat_pcr);
T_test_hat  = X_test * V_pcr(:, 1:q_hat_pcr);

beta_score = T_train_hat \ w_train_c;

w_hat_test_centered = T_test_hat * beta_score;
w_hat_test = w_hat_test_centered + w_mean;

[rmse_w, corr_w, R2_w] = compute_oos_metrics(w_test, w_hat_test);

runtime_sec = toc(method_tic);

result_pcr = table;
result_pcr.method = string(method_name);
result_pcr.selected_rank = q_hat_pcr;
result_pcr.rmse = rmse_w;
result_pcr.correlation = corr_w;
result_pcr.R2 = R2_w;
result_pcr.runtime_sec = runtime_sec;

writetable( ...
    result_pcr, ...
    fullfile(out_dir, 'Output_from_M5_Levine_OOS_PCR.csv') ...
);

fprintf('\nPCR completed.\n');
fprintf('PCR test RMSE = %.4f | corr = %.4f | R2 = %.4f | runtime = %.2f sec\n', ...
    rmse_w, corr_w, R2_w, runtime_sec);

all_results = [all_results; result_pcr];

clear X_train_c X_test X_mean w_train_c U_pcr S_pcr V_pcr Scores_train
clear T_train_hat T_test_hat beta_score w_hat_test_centered w_hat_test

fprintf('\nMemory after PCR cleanup:\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Method 2: PLS
% ============================================================

fprintf('\n============================================================\n');
fprintf('[Method 2/3] PLS started\n');
fprintf('============================================================\n');

method_tic = tic;
method_name = 'PLS';

Kfold = 5;

fprintf('PLS rank selection by %d-fold CV on training data.\n', Kfold);
fprintf('q_grid = ');
disp(q_grid);

cv = cvpartition(n_train, 'KFold', Kfold);

cv_rmse = nan(length(q_grid), 1);

for ii = 1:length(q_grid)
    
    q = q_grid(ii);
    fold_mse = nan(Kfold, 1);
    
    fprintf('\nPLS CV for q = %d\n', q);
    
    for kk = 1:Kfold
        
        fprintf('  Fold %d / %d: preparing data...\n', kk, Kfold);
        
        idx_tr_rel = training(cv, kk);
        idx_va_rel = test(cv, kk);
        
        idx_tr_global = train_idx(idx_tr_rel);
        idx_va_global = train_idx(idx_va_rel);
        
        X_cv_tr = double(X(idx_tr_global, :));
        w_cv_tr = w_raw(idx_tr_global);
        
        X_cv_mean = mean(X_cv_tr, 1);
        w_cv_mean = mean(w_cv_tr);
        
        fprintf('  Fold %d / %d: centering training fold...\n', kk, Kfold);
        
        for jj = 1:p
            X_cv_tr(:, jj) = X_cv_tr(:, jj) - X_cv_mean(jj);
        end
        
        w_cv_tr_c = w_cv_tr - w_cv_mean;
        
        fprintf('  Fold %d / %d: fitting PLS q = %d...\n', kk, Kfold, q);
        
        [~, ~, ~, ~, BETA_cv] = plsregress(X_cv_tr, w_cv_tr_c, q);
        
        clear X_cv_tr w_cv_tr_c
        
        X_cv_va = double(X(idx_va_global, :));
        w_cv_va = w_raw(idx_va_global);
        
        fprintf('  Fold %d / %d: centering validation fold and predicting...\n', kk, Kfold);
        
        for jj = 1:p
            X_cv_va(:, jj) = X_cv_va(:, jj) - X_cv_mean(jj);
        end
        
        w_cv_pred_c = [ones(length(idx_va_global), 1), X_cv_va] * BETA_cv;
        w_cv_pred = w_cv_pred_c + w_cv_mean;
        
        fold_mse(kk) = mean((w_cv_va - w_cv_pred).^2);
        
        fprintf('  Fold %d / %d: RMSE = %.4f\n', kk, Kfold, sqrt(fold_mse(kk)));
        
        clear X_cv_va X_cv_mean w_cv_va w_cv_pred_c w_cv_pred BETA_cv
    end
    
    cv_rmse(ii) = sqrt(mean(fold_mse, 'omitnan'));
    
    fprintf('PLS q = %d | CV RMSE = %.4f\n', q, cv_rmse(ii));
end

[~, best_idx] = min(cv_rmse);
q_hat_pls = q_grid(best_idx);

fprintf('\nSelected PLS rank q_hat = %d\n', q_hat_pls);

fprintf('Preparing full training matrix for final PLS...\n');

X_train = double(X(train_idx, :));
X_mean = mean(X_train, 1);
w_mean = mean(w_train);

fprintf('Centering final PLS training matrix...\n');

for jj = 1:p
    X_train(:, jj) = X_train(:, jj) - X_mean(jj);
    
    if mod(jj, 100000) == 0
        fprintf('  PLS final training-centering progress: %d / %d CpGs\n', jj, p);
    end
end

w_train_c = w_train - w_mean;

fprintf('Fitting final PLS model with q_hat = %d...\n', q_hat_pls);

[~, ~, ~, ~, BETA_hat] = plsregress(X_train, w_train_c, q_hat_pls);

clear X_train w_train_c

fprintf('Preparing test matrix for PLS prediction...\n');

X_test = double(X(test_idx, :));

fprintf('Centering PLS test matrix...\n');

for jj = 1:p
    X_test(:, jj) = X_test(:, jj) - X_mean(jj);
    
    if mod(jj, 100000) == 0
        fprintf('  PLS test-centering progress: %d / %d CpGs\n', jj, p);
    end
end

w_hat_test_centered = [ones(length(test_idx), 1), X_test] * BETA_hat;
w_hat_test = w_hat_test_centered + w_mean;

[rmse_w, corr_w, R2_w] = compute_oos_metrics(w_test, w_hat_test);

runtime_sec = toc(method_tic);

result_pls = table;
result_pls.method = string(method_name);
result_pls.selected_rank = q_hat_pls;
result_pls.rmse = rmse_w;
result_pls.correlation = corr_w;
result_pls.R2 = R2_w;
result_pls.runtime_sec = runtime_sec;

writetable( ...
    result_pls, ...
    fullfile(out_dir, 'Output_from_M5_Levine_OOS_PLS.csv') ...
);

fprintf('\nPLS completed.\n');
fprintf('PLS test RMSE = %.4f | corr = %.4f | R2 = %.4f | runtime = %.2f sec\n', ...
    rmse_w, corr_w, R2_w, runtime_sec);

all_results = [all_results; result_pls];

clear X_test X_mean BETA_hat w_hat_test_centered w_hat_test cv cv_rmse

fprintf('\nMemory after PLS cleanup:\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Method 3: SOLAR
% ============================================================

fprintf('\n============================================================\n');
fprintf('[Method 3/3] SOLAR started\n');
fprintf('============================================================\n');

method_tic = tic;
method_name = 'SOLAR';

fprintf('Preparing SOLAR training matrix...\n');

X_train = double(X(train_idx, :));
w_train = w_raw(train_idx);

[n_train_check, p_check] = size(X_train);

if n_train_check ~= n_train || p_check ~= p
    error('SOLAR training matrix dimensions are inconsistent.');
end

fprintf('SOLAR training n = %d | p = %d\n', n_train, p);

fprintf('Centering SOLAR training matrix and outcome...\n');

center_tic = tic;

X_mean = mean(X_train, 1);
w_mean = mean(w_train);

for jj = 1:p
    X_train(:, jj) = X_train(:, jj) - X_mean(jj);
    
    if mod(jj, 100000) == 0
        fprintf('  SOLAR centering progress: %d / %d CpGs\n', jj, p);
    end
end

w = w_train - w_mean;

center_time_sec = toc(center_tic);

fprintf('SOLAR centering completed in %.2f seconds.\n', center_time_sec);

Xnorm2 = norm(X_train, 'fro')^2;
Xnorm  = sqrt(Xnorm2);

%% ----- SOLAR settings -----

print_every = 2;

q_init = 5;
q_min  = q_grid(1);
q_max  = q_grid(end);

Num_iters = 50;

perform_gram_svd = 1;

T0 = 1.0;
T_min = 1e-2;
cool_pow = 0.7;

rho2 = 0.5^2;
g2   = 25.0;

c_kappa = 1;
zeta = 2/3;

kappa = c_kappa * ...
    (((sqrt(n_train) + sqrt(p))^2) / (n_train + p))^zeta;

q_max_eff  = min([q_max, n_train - 1, p]);
q_init_eff = min(q_init, q_max_eff);

fprintf('\nSOLAR tuning/settings:\n');
fprintf('q_init = %d | q_min = %d | q_max = %d\n', q_init, q_min, q_max);
fprintf('q_init_eff = %d | q_max_eff = %d\n', q_init_eff, q_max_eff);
fprintf('Num_iters = %d | print_every = %d\n', Num_iters, print_every);
fprintf('perform_gram_svd = %d\n', perform_gram_svd);
fprintf('kappa = %.6f\n', kappa);
fprintf('rho2 = %.6f | g2 = %.6f\n', rho2, g2);

%% ----- Precompute SVD(X_train) -----

fprintf('\nSOLAR precomputing leading SVD components...\n');

svd_tic = tic;

if perform_gram_svd == 0
    
    [U_full, S_full, V_full] = svds(X_train, q_max_eff);
    
else
    
    fprintf('Computing Gram matrix X_train * X_train''...\n');
    
    G = double(X_train * X_train');
    
    fprintf('Computing eigs of Gram matrix...\n');
    
    [U_full, S2] = eigs(G, q_max_eff);
    
    clear G
    
    s = sqrt(diag(S2));
    S_full = diag(s);
    
    clear S2
    
    fprintf('Computing right singular vectors V_full...\n');
    
    V_full = X_train' * U_full;
    
    for k = 1:q_max_eff
        V_full(:, k) = V_full(:, k) / s(k);
    end
end

svd_time_sec = toc(svd_tic);

fprintf('SOLAR SVD/precomputation completed in %.2f seconds.\n', svd_time_sec);

%% ----- Constants -----

fprintf('\nSOLAR estimating variance constants...\n');

sigma2 = estimate_sigma2_resid_baseline(X_train, 2);
tau2   = var(w, 1);

prec_lik   = 1 / sigma2;
prec_prior = 1 / rho2;
prec_post  = prec_lik + prec_prior;

tau_inv2 = 1 / tau2;
g_inv2   = 1 / g2;
prec_b   = tau_inv2 + g_inv2;
beta_shrink = tau_inv2 / prec_b;

fprintf('sigma2_est = %.6f\n', sigma2);
fprintf('tau2_est   = %.6f\n', tau2);

%% ----- Initialize -----

q = q_init_eff;

H = U_full(:, 1:q);
V = V_full(:, 1:q);
D = (prec_lik / prec_post) * S_full(1:q, 1:q);
d = diag(D);

beta = beta_shrink * (H' * w);

%% ----- Best-state tracking -----

obj_best  = -Inf;
q_best    = q;
H_best    = H;
V_best    = V;
D_best    = D;
beta_best = beta;

recon_best = Inf;

%% ----- Early stopping -----

W_stop      = min(30, floor(0.5 * Num_iters));
epsJ_stop   = 1e-6;
epsAcc_stop = 0.01;
etaT_stop   = 0.10;

best_obj_hist = nan(Num_iters, 1);
best_q_hist   = nan(Num_iters, 1);
q_hist         = nan(Num_iters, 1);
q_accept_hist  = zeros(Num_iters, 1);

%% ----- Iterative MAP + trans-dimensional search -----

fprintf('\nStarting SOLAR optimization...\n');

opt_tic = tic;
last_iter = Num_iters;

for iter = 1:Num_iters
    
    T = max(T_min, T0 * (iter^(-cool_pow)));
    
    Fv = X_train' * (H * D);
    [U1, ~, V1] = svd(Fv, 'econ');
    V = U1 * V1';
    
    Fh = (w * beta') / tau2 + (X_train * V * D) / sigma2;
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
    
    R = H' * X_train * V;
    [Ur, Sr, Vr] = svd(R, 'econ');
    
    H = H * Ur;
    V = V * Vr;
    
    d = (prec_lik / prec_post) * diag(Sr);
    D = diag(d);
    
    beta = beta_shrink * (H' * w);
    
    XV  = X_train * V;
    XVD = XV * D;
    
    residX2 = Xnorm2 + norm(D, 'fro')^2 - 2 * trace(H' * XVD);
    
    like_curr = ...
        -0.5 / sigma2 * residX2 ...
        -0.5 / tau2   * norm(w - H * beta)^2 ...
        -0.5 / rho2   * norm(d, 2)^2 ...
        -0.5 / g2     * norm(beta, 2)^2;
    
    dfq = q * (n_train + p - 2*q) + 2*q;
    pen_curr = kappa * 0.5 * log(n_train * p) * dfq;
    
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
        
        %%% Added this part only for OOS performance under fixed q
        if q_min == q_max_eff
            q_prop = q_min;
        end
        
        H_prop = U_full(:, 1:q_prop);
        V_prop = V_full(:, 1:q_prop);
        D_prop = (prec_lik / prec_post) * S_full(1:q_prop, 1:q_prop);
        d_prop = diag(D_prop);
        
        beta_prop = beta_shrink * (H_prop' * w);
        
        XVp  = X_train * V_prop;
        XVDp = XVp * D_prop;
        
        residXp2 = Xnorm2 + norm(D_prop, 'fro')^2 - ...
            2 * trace(H_prop' * XVDp);
        
        like_prop = ...
            -0.5 / sigma2 * residXp2 ...
            -0.5 / tau2   * norm(w - H_prop * beta_prop)^2 ...
            -0.5 / rho2   * norm(d_prop, 2)^2 ...
            -0.5 / g2     * norm(beta_prop, 2)^2;
        
        dfq_prop = q_prop * (n_train + p - 2*q_prop) + 2*q_prop;
        pen_prop = kappa * 0.5 * log(n_train * p) * dfq_prop;
        
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
        
        XV  = X_train * V;
        XVD = XV * D;
        
        residX2_iter = Xnorm2 + norm(D, 'fro')^2 - ...
            2 * trace(H' * XVD);
        
        recon_iter = sqrt(max(residX2_iter, 0)) / Xnorm;
        
        fprintf(['Iter %5d | q=%2d | recon=%.4f | obj=%.3e | ', ...
                 'best(q=%2d,recon=%.4f,obj=%.3e)\n'], ...
            iter, q, recon_iter, obj_curr, ...
            q_best, recon_best, obj_best);
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

q = q_best;
H = H_best;
V = V_best;
D = D_best;
beta = beta_best;
d = diag(D);

fprintf('\nSOLAR optimization complete.\n');
fprintf('Selected SOLAR rank q_hat = %d\n', q);
fprintf('SOLAR optimization time = %.2f seconds\n', optimization_time_sec);

clear H_best V_best D_best beta_best U_full S_full V_full
clear H D XV XVD Fv Fh R U1 V1 U2 V2 Ur Sr Vr

fprintf('\nPreparing SOLAR test matrix for OOS prediction...\n');

X_test = double(X(test_idx, :));

fprintf('Centering SOLAR test matrix using training means...\n');

for jj = 1:p
    X_test(:, jj) = X_test(:, jj) - X_mean(jj);
    
    if mod(jj, 100000) == 0
        fprintf('  SOLAR test-centering progress: %d / %d CpGs\n', jj, p);
    end
end

fprintf('Projecting test samples into SOLAR latent space...\n');

H_test = X_test * V;

for kk = 1:q
    if abs(d(kk)) > 1e-12
        H_test(:, kk) = H_test(:, kk) / d(kk);
    else
        H_test(:, kk) = 0;
    end
end

w_hat_test_centered = H_test * beta;
w_hat_test = w_hat_test_centered + w_mean;

[rmse_w, corr_w, R2_w] = compute_oos_metrics(w_test, w_hat_test);

runtime_sec = toc(method_tic);

result_solar = table;
result_solar.method = string(method_name);
result_solar.selected_rank = q;
result_solar.rmse = rmse_w;
result_solar.correlation = corr_w;
result_solar.R2 = R2_w;
result_solar.runtime_sec = runtime_sec;

writetable( ...
    result_solar, ...
    fullfile(out_dir, 'Output_from_M5_Levine_OOS_SOLAR.csv') ...
);

fprintf('\nSOLAR completed.\n');
fprintf('SOLAR test RMSE = %.4f | corr = %.4f | R2 = %.4f | runtime = %.2f sec\n', ...
    rmse_w, corr_w, R2_w, runtime_sec);

all_results = [all_results; result_solar];

clear X_train X_test X_mean H_test V beta d w_hat_test_centered w_hat_test

fprintf('\nMemory after SOLAR cleanup:\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Save combined OOS results
% ============================================================

fprintf('\n============================================================\n');
fprintf('Saving combined OOS validation results\n');
fprintf('============================================================\n');

writetable( ...
    all_results, ...
    fullfile(out_dir, 'Output_from_M5_Levine_OOS_combined_metrics.csv') ...
);

disp(all_results);

fprintf('\nSaved all OOS results in folder: %s\n', out_dir);
fprintf('\nDone.\n');

%% ============================================================
%  Local function: OOS metrics
% ============================================================

function [rmse_w, corr_w, R2_w] = compute_oos_metrics(y_true, y_pred)

    y_true = y_true(:);
    y_pred = y_pred(:);
    
    valid = ~isnan(y_true) & ~isnan(y_pred);
    
    y_true = y_true(valid);
    y_pred = y_pred(valid);
    
    rmse_w = sqrt(mean((y_true - y_pred).^2));
    
    if numel(y_true) > 2 && std(y_pred) > 0 && std(y_true) > 0
        corr_w = corr(y_pred, y_true, 'Rows', 'complete');
    else
        corr_w = NaN;
    end
    
    ss_res = sum((y_true - y_pred).^2);
    ss_tot = sum((y_true - mean(y_true)).^2);
    
    if ss_tot > 0
        R2_w = 1 - ss_res / ss_tot;
    else
        R2_w = NaN;
    end
end