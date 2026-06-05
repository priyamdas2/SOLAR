%% Run on > 400 GB RAM cluster (with 40 parallel threads)
%% ============================================================
%  SOLAR bootstrap rank stability: M1 Horvath residual acceleration
%
%  Output:
%      Output/Output_from_M1boot_onlyRanks_NumBoot_<Num_boot>.csv
% ============================================================

clearvars -except X pheno_FULL y_raw_FULL
clc;
rng(123);

fprintf('MATLAB mem (Windows only):\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Relative directories
% ============================================================

data_dir = 'Age methylation data';
supp_dir = 'SOLAR supp funs';
out_dir  = 'Output';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if ~exist(supp_dir, 'dir')
    error('SOLAR support-function folder not found: %s', supp_dir);
end

addpath(supp_dir);

X_mat_file = fullfile(data_dir, 'X_matrix_final.mat');
pheno_file = fullfile(data_dir, 'Pheno_for_MATLAB_with_epi_ages.csv');

if ~exist(X_mat_file, 'file')
    error('X matrix file not found: %s', X_mat_file);
end

if ~exist(pheno_file, 'file')
    error('Phenotype file not found: %s', pheno_file);
end

%% ============================================================
%  Load phenotype data
% ============================================================

if exist('pheno_FULL', 'var') && exist('y_raw_FULL', 'var')

    fprintf('\nUsing phenotype table and Horvath outcome already available in memory.\n');

    pheno = pheno_FULL;
    w_raw = y_raw_FULL;

else

    fprintf('\nReading phenotype file...\n');

    pheno_FULL = readtable(pheno_file, 'VariableNamingRule', 'preserve');

    if ~ismember('AgeAccelResidual_Horvath', pheno_FULL.Properties.VariableNames)
        error('AgeAccelResidual_Horvath column not found in phenotype file.');
    end

    y_raw_FULL = pheno_FULL.AgeAccelResidual_Horvath;

    pheno = pheno_FULL;
    w_raw = y_raw_FULL;
end

%% ============================================================
%  Load methylation matrix X
% ============================================================

if exist('X', 'var')

    fprintf('\nUsing X already available in memory.\n');

else

    fprintf('\nLoading X_matrix_final.mat. This may take time...\n');

    tic;

    tmp = load(X_mat_file);
    var_names = fieldnames(tmp);

    if numel(var_names) ~= 1
        fprintf('\nVariables found in X_matrix_final.mat:\n');
        disp(var_names);
        error('Expected exactly one numeric matrix inside X_matrix_final.mat.');
    end

    X = tmp.(var_names{1});
    clear tmp

    fprintf('Finished loading X in %.2f seconds.\n', toc);
end

%% ============================================================
%  Basic checks and filtering
% ============================================================

[n, p] = size(X);

fprintf('\nLoaded methylation matrix:\n');
fprintf('n = %d samples\n', n);
fprintf('p = %d CpGs\n', p);

fprintf('\nPhenotype table:\n');
fprintf('n_pheno = %d rows\n', height(pheno));

if height(pheno) ~= n
    error('Number of phenotype rows does not match number of X rows.');
end

if length(w_raw) ~= n
    error('Length of Horvath residual outcome does not match number of X rows.');
end

valid_idx = ~isnan(w_raw);

if any(~valid_idx)

    fprintf('\nRemoving %d samples with missing Horvath residual acceleration.\n', ...
        sum(~valid_idx));

    X = X(valid_idx, :);
    w_raw = w_raw(valid_idx);
    pheno = pheno(valid_idx, :);
end

[n, p] = size(X);

fprintf('\nAnalysis data after filtering:\n');
fprintf('n = %d samples\n', n);
fprintf('p = %d CpGs\n', p);

%% ============================================================
%  Center X and Horvath residual acceleration
% ============================================================

fprintf('\nCentering X and Horvath residual acceleration...\n');

center_tic = tic;

X_mean = mean(X, 1);
w_mean = mean(w_raw);

for jj = 1:p
    X(:, jj) = X(:, jj) - X_mean(jj);
end

w = w_raw - w_mean;

center_time_sec = toc(center_tic);

fprintf('Finished centering in %.2f seconds.\n', center_time_sec);

%% ============================================================
%  Bootstrap and SOLAR settings
% ============================================================

method_name = 'SOLAR_Horvath_bootstrap_rank';

Num_boot = 500;

print_every = 2;

q_init = 10;
q_min  = 1;
q_max  = 15;

Num_iters = 60;

perform_gram_svd = 1;

T0 = 1.0;
T_min = 1e-2;
cool_pow = 0.7;

rho2 = 0.5^2;
g2   = 25.0;

c_kappa = 1;
zeta = 2/3;

kappa = c_kappa * ...
    (((sqrt(n) + sqrt(p))^2) / (n + p))^zeta;

q_max_eff  = min([q_max, n - 1, p]);
q_init_eff = min(q_init, q_max_eff);

W_stop      = 40;
epsJ_stop   = 1e-6;
epsAcc_stop = 0.01;
etaT_stop   = 0.10;

fprintf('\nSOLAR bootstrap rank-stability settings:\n');
fprintf('method_name = %s\n', method_name);
fprintf('Num_boot = %d\n', Num_boot);
fprintf('q_init = %d | q_min = %d | q_max = %d\n', q_init, q_min, q_max);
fprintf('q_init_eff = %d | q_max_eff = %d\n', q_init_eff, q_max_eff);
fprintf('Num_iters = %d | W_stop = %d | print_every = %d\n', ...
    Num_iters, W_stop, print_every);
fprintf('perform_gram_svd = %d\n', perform_gram_svd);
fprintf('kappa = %.6f\n', kappa);
fprintf('rho2 = %.6f | g2 = %.6f\n', rho2, g2);

%% ============================================================
%  Start parallel pool
%  Prefer thread pool to reduce duplication of the large X matrix.
% ============================================================

poolobj = gcp('nocreate');

if isempty(poolobj)
    try
        parpool('threads');
        fprintf('\nStarted thread-based parallel pool.\n');
    catch
        warning(['Could not start thread-based pool. Starting default pool. ', ...
                 'This may duplicate X across workers and increase RAM use.']);
        parpool;
    end
else
    fprintf('\nUsing existing parallel pool with %d workers.\n', poolobj.NumWorkers);
end

%% ============================================================
%  Progress printing from parfor
% ============================================================

dq = parallel.pool.DataQueue;

afterEach(dq, @(msg) fprintf('%s\n', msg));

%% ============================================================
%  Bootstrap loop
% ============================================================

fprintf('\nStarting bootstrap rank-stability analysis...\n');

boot_tic = tic;

boot_rank = nan(Num_boot, 1);
boot_iter_completed = nan(Num_boot, 1);

boot_seeds = 123000 + (1:Num_boot)';

parfor bb = 1:Num_boot

    rng(boot_seeds(bb), 'twister');

    boot_idx = randsample(n, n, true);

    [q_hat_boot, last_iter_boot] = run_one_SOLAR_rank_bootstrap( ...
        X, ...
        w, ...
        boot_idx, ...
        q_init_eff, ...
        q_min, ...
        q_max_eff, ...
        Num_iters, ...
        print_every, ...
        perform_gram_svd, ...
        T0, ...
        T_min, ...
        cool_pow, ...
        rho2, ...
        g2, ...
        kappa, ...
        W_stop, ...
        epsJ_stop, ...
        epsAcc_stop, ...
        etaT_stop ...
    );

    boot_rank(bb) = q_hat_boot;
    boot_iter_completed(bb) = last_iter_boot;

    send(dq, sprintf( ...
        'Bootstrap iteration %d/%d completed; selected rank = %d.', ...
        bb, Num_boot, q_hat_boot));
end

boot_time_sec = toc(boot_tic);

fprintf('\nBootstrap rank-stability analysis complete.\n');
fprintf('Total bootstrap runtime = %.2f seconds\n', boot_time_sec);

%% ============================================================
%  Save one compact CSV output
% ============================================================

boot_tbl = table;

boot_tbl.boot_iter = (1:Num_boot)';
boot_tbl.selected_rank = boot_rank;
boot_tbl.Num_iters_completed = boot_iter_completed;

out_file = fullfile( ...
    out_dir, ...
    sprintf('Output_from_M1boot_onlyRanks_NumBoot_%d.csv', Num_boot) ...
);

writetable(boot_tbl, out_file);

fprintf('\nSaved bootstrap rank output:\n%s\n', out_file);

fprintf('\nBootstrap selected-rank frequencies:\n');
disp(tabulate(boot_rank));

try
    disp(memory);
catch
end

%% ============================================================
%  Local function: one bootstrap SOLAR rank fit
% ============================================================

function [q_best, last_iter] = run_one_SOLAR_rank_bootstrap( ...
    X_full, ...
    w_full, ...
    boot_idx, ...
    q_init_eff, ...
    q_min, ...
    q_max_eff, ...
    Num_iters, ...
    print_every, ...
    perform_gram_svd, ...
    T0, ...
    T_min, ...
    cool_pow, ...
    rho2, ...
    g2, ...
    kappa, ...
    W_stop, ...
    epsJ_stop, ...
    epsAcc_stop, ...
    etaT_stop)

    Xb = X_full(boot_idx, :);
    wb = w_full(boot_idx);

    [n_b, p_b] = size(Xb);

    Xnorm2 = norm(Xb, 'fro')^2;
    Xnorm  = sqrt(Xnorm2);

    %% --------------------------------------------------------
    %  Precompute SVD
    %% --------------------------------------------------------

    if perform_gram_svd == 0

        [U_full, S_full, V_full] = svds(Xb, q_max_eff);

    else

        G = double(Xb * Xb');

        [U_full, S2] = eigs(G, q_max_eff);

        s = sqrt(diag(S2));
        S_full = diag(s);

        V_full = Xb' * U_full;

        for kk = 1:q_max_eff
            if s(kk) > 0
                V_full(:, kk) = V_full(:, kk) / s(kk);
            end
        end
    end

    %% --------------------------------------------------------
    %  Constants
    %% --------------------------------------------------------

    sigma2 = estimate_sigma2_resid_baseline(Xb, 2);
    tau2   = var(wb, 1);

    prec_lik   = 1 / sigma2;
    prec_prior = 1 / rho2;
    prec_post  = prec_lik + prec_prior;

    tau_inv2 = 1 / tau2;
    g_inv2   = 1 / g2;
    prec_b   = tau_inv2 + g_inv2;
    beta_shrink = tau_inv2 / prec_b;

    %% --------------------------------------------------------
    %  Initialize
    %% --------------------------------------------------------

    q = q_init_eff;

    H = U_full(:, 1:q);
    V = V_full(:, 1:q);
    D = (prec_lik / prec_post) * S_full(1:q, 1:q);
    d = diag(D);

    beta = beta_shrink * (H' * wb);

    %% --------------------------------------------------------
    %  Best-state tracking
    %% --------------------------------------------------------

    obj_best = -Inf;
    q_best = q;
    recon_best = Inf;

    best_obj_hist = nan(Num_iters, 1);
    best_q_hist   = nan(Num_iters, 1);
    q_accept_hist = zeros(Num_iters, 1);

    last_iter = Num_iters;

    %% --------------------------------------------------------
    %  Iterative MAP + rank search
    %% --------------------------------------------------------

    for iter = 1:Num_iters

        T = max(T_min, T0 * (iter^(-cool_pow)));

        %% ----- V update -----

        Fv = Xb' * (H * D);
        [U1, ~, V1] = svd(Fv, 'econ');
        V = U1 * V1';

        %% ----- H update -----

        Fh = (wb * beta') / tau2 + (Xb * V * D) / sigma2;
        [U2, ~, V2] = svd(Fh, 'econ');
        H = U2 * V2';

        %% ----- beta update -----

        beta = beta_shrink * (H' * wb);

        %% ----- Canonical supervised orientation -----

        u = H' * wb;

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

        %% ----- Re-diagonalize -----

        R = H' * Xb * V;
        [Ur, Sr, Vr] = svd(R, 'econ');

        H = H * Ur;
        V = V * Vr;

        d = (prec_lik / prec_post) * diag(Sr);
        D = diag(d);

        beta = beta_shrink * (H' * wb);

        %% ----- Objective -----

        XV  = Xb * V;
        XVD = XV * D;

        residX2 = Xnorm2 + norm(D, 'fro')^2 - 2 * trace(H' * XVD);

        like_curr = ...
            -0.5 / sigma2 * residX2 ...
            -0.5 / tau2   * norm(wb - H * beta)^2 ...
            -0.5 / rho2   * norm(d, 2)^2 ...
            -0.5 / g2     * norm(beta, 2)^2;

        dfq = q * (n_b + p_b - 2*q) + 2*q;
        pen_curr = kappa * 0.5 * log(n_b * p_b) * dfq;

        obj_curr = like_curr - pen_curr;

        %% ----- Update best state -----

        if obj_curr > obj_best
            obj_best = obj_curr;
            q_best = q;
            recon_best = sqrt(max(residX2, 0)) / Xnorm;
        end

        %% ----- Propose rank move -----

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

            beta_prop = beta_shrink * (H_prop' * wb);

            XVp  = Xb * V_prop;
            XVDp = XVp * D_prop;

            residXp2 = Xnorm2 + norm(D_prop, 'fro')^2 - ...
                2 * trace(H_prop' * XVDp);

            like_prop = ...
                -0.5 / sigma2 * residXp2 ...
                -0.5 / tau2   * norm(wb - H_prop * beta_prop)^2 ...
                -0.5 / rho2   * norm(d_prop, 2)^2 ...
                -0.5 / g2     * norm(beta_prop, 2)^2;

            dfq_prop = q_prop * (n_b + p_b - 2*q_prop) + 2*q_prop;
            pen_prop = kappa * 0.5 * log(n_b * p_b) * dfq_prop;

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
            end
        end

        %% ----- Optional internal progress, suppressed in parfor -----
        % Keep print_every as an input for consistency, but avoid printing
        % every inner iteration from parallel workers.

        if mod(iter, print_every) == 0
            
            fprintf(['Iter %5d | q=%2d | ', ...
                'best(q=%2d)\n'], ...
                iter, q, q_best);
        end

        best_obj_hist(iter) = obj_best;
        best_q_hist(iter)   = q_best;

        %% ----- Early stopping -----

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
                last_iter = iter;
                break;
            end
        end
    end
end