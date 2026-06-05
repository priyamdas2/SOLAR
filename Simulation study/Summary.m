clear all;

%% ============================================================
%  Scenario inputs
% ============================================================

n = 500;          % 100, 1000
p = 1000;         % 1000, 10000, 100000

q_true = 5;       % 3, 5

target_SNR_X = 0.05;   % 0.01, 0.05, 0.2
target_SNR_w = 0.2;    % 0.05, 0.2, 0.5

rep_vec = 1:10;

method_names = {'PCR', 'PLS', 'SPC','SOLAR'}; % 'PCR', 'PLS', 'SPC', 'SOLAR'
%% ============================================================
%  Directories
% ============================================================

out_dir = 'Output';
summary_dir = 'Output_Summary';

if ~exist(summary_dir, 'dir')
    mkdir(summary_dir);
end

%% ============================================================
%  Filename strings
% ============================================================

snrX_str = strrep(sprintf('%.3f', target_SNR_X), '.', 'p');
snrw_str = strrep(sprintf('%.3f', target_SNR_w), '.', 'p');

%% ============================================================
%  Quantities to summarize
% ============================================================

metric_names = { ...
    'q_hat', ...
    'rank_correct', ...
    'rank_abs_error', ...
    'H_projector_error', ...
    'V_projector_error', ...
    'X_signal_recon_error', ...
    'X_signal_rmse', ...
    'supervised_signal_rmse', ...
    'supervised_signal_corr', ...
    'supervised_signal_rel_error', ...
    'rmse_w', ...
    'mae_w', ...
    'corr_w', ...
    'R2_w', ...
    'beta_l2_error', ...
    'runtime_sec'};


for mmm = 1:length(method_names)
    
    method_name = method_names{mmm};
    
    fprintf('\n=========================================\n');
    fprintf('Processing method: %s\n', method_name);
    fprintf('=========================================\n');
    
    
    %% ============================================================
    %  Read all replication-level outputs
    % ============================================================
    
    all_results = table;
    
    for rep = rep_vec
        
        fname = sprintf( ...
            ['Output_%s_snrX_%s_snrw_%s_qTrue_%d_' ...
            'n_%d_p_%d_rep_%02d.csv'], ...
            method_name, snrX_str, snrw_str, q_true, n, p, rep);
        
        fpath = fullfile(out_dir, fname);
        
        if ~exist(fpath, 'file')
            warning('File not found: %s. Skipping.', fpath);
            continue;
        end
        
        T = readtable(fpath);
        all_results = [all_results; T];
        
    end
    
    if isempty(all_results)
        error('No output files were found for the specified scenario.');
    end
    
    %% ============================================================
    %  Create one-row summary table: mean (SE)
    % ============================================================
    
    summary = table;
    summary.method = string(method_name);
    
    for mm = 1:length(metric_names)
        
        varname = metric_names{mm};
        
        if ~ismember(varname, all_results.Properties.VariableNames)
            warning('Metric %s not found. Skipping.', varname);
            continue;
        end
        
        x = all_results.(varname);
        x = x(~isnan(x));
        
        if isempty(x)
            summary.(varname) = string("NA");
            continue;
        end
        
        m = mean(x);
        se = std(x) / sqrt(length(x));
        
        summary.(varname) = string(format_mean_se(m, se, varname));
        
    end
    
    %% ============================================================
    %  Save summary CSV
    % ============================================================
    
    summary_fname = sprintf( ...
        ['Summary_%s_snrX_%s_snrw_%s_qTrue_%d_' ...
        'n_%d_p_%d_reps_%02d_to_%02d.csv'], ...
        method_name, snrX_str, snrw_str, q_true, n, p, ...
        min(rep_vec), max(rep_vec));
    
    writetable(summary, fullfile(summary_dir, summary_fname));
    
    fprintf('Saved summary: %s\n', summary_fname);
    disp(summary);
    
end

%% ============================================================
%  Display combined summary across available methods
%  Looks for PCR, PLS, SPC, and SOLAR summary files for this scenario
% ============================================================

method_list = {'PCR', 'PLS', 'SPC', 'SOLAR'};

combined_summary = table;

for mmeth = 1:length(method_list)
    
    meth = method_list{mmeth};
    
    summary_fname_meth = sprintf( ...
        ['Summary_%s_snrX_%s_snrw_%s_qTrue_%d_' ...
        'n_%d_p_%d_reps_%02d_to_%02d.csv'], ...
        meth, snrX_str, snrw_str, q_true, n, p, ...
        min(rep_vec), max(rep_vec));
    
    summary_path_meth = fullfile(summary_dir, summary_fname_meth);
    
    if exist(summary_path_meth, 'file')
        
        Tmeth = readtable(summary_path_meth);
        
        combined_summary = [combined_summary; Tmeth];
        
    else
        
        fprintf('Summary file not found for method %s. Skipping.\n', meth);
        
    end
    
end

if ~isempty(combined_summary)
    
    fprintf('\n========================================================\n');
    fprintf('Combined method summary for this scenario\n');
    fprintf('SNR_X=%s, SNR_w=%s, qTrue=%d, n=%d, p=%d\n', ...
        snrX_str, snrw_str, q_true, n, p);
    fprintf('========================================================\n');
    
    disp(combined_summary);
    
    %% ============================================================
    %  Save combined summary CSV
    % ============================================================
    
    combined_summary_fname = sprintf( ...
        ['Summary_ALL_snrX_%s_snrw_%s_qTrue_%d_' ...
        'n_%d_p_%d_reps_%02d_to_%02d.csv'], ...
        snrX_str, snrw_str, q_true, n, p, ...
        min(rep_vec), max(rep_vec));
    
    writetable(combined_summary, ...
        fullfile(summary_dir, combined_summary_fname));
    
    fprintf('Saved combined summary: %s\n', ...
        combined_summary_fname);
    
else
    
    fprintf('\nNo method summaries found for combined display.\n');
    
end

%% ============================================================
%  Helper function
% ============================================================

function out = format_mean_se(m, se, varname)

% Default: 3 decimals
nd = 3;

% Rank/component quantities
if ismember(varname, {'q_hat', 'rank_abs_error'})
    nd = 2;
end

% Binary/rate quantities
if ismember(varname, {'rank_correct'})
    nd = 3;
end

% Runtime
if strcmp(varname, 'runtime_sec')
    nd = 2;
end

% Prediction correlations/R2
if ismember(varname, {'supervised_signal_corr', 'corr_w', 'R2_w'})
    nd = 3;
end

% Errors
if contains(varname, 'error') || contains(varname, 'rmse') || contains(varname, 'mae')
    nd = 3;
end

% X signal rmse
if ismember(varname, {'X_signal_rmse'})
    nd = 4;
end

fmt = sprintf('%%.%df (%%.%df)', nd, nd);
out = sprintf(fmt, m, se);

end