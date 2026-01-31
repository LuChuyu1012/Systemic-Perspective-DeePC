clc; clear; close all;

sys_list = 1:6;
f_target = 2;
SIs      = 0:3;
nruns    = 20;

Fs = 100;

% ======= Spread-only convergence threshold (NO reference) =======
eps_spread = 0.1;

g_list = [0 1 2 3];

SNRdB_target = 30;
extra_suffix = sprintf('_snr%02ddB', round(SNRdB_target));

lambda_list = [0.001 0.01 0.1 1];

% Prefer y_seq_opt; fall back to others if missing
y_fields = {'y_seq_opt','y_seq','y','y_out'};

warning('off','MATLAB:print:InvalidGraphicsState');
warning('off','MATLAB:graphics:SceneNode');

nsys  = numel(sys_list);
ncase = numel(SIs);

% ======= store best results =======
best_g        = NaN(nsys, ncase);
best_lam      = NaN(nsys, ncase);

best_RMSEbar  = NaN(nsys, ncase);   % RMSE(mean traj vs ref) for selected g
best_Tspr     = NaN(nsys, ncase);   % spread-only convergence time for selected g

% ======= NEW: run-wise RMSE stats at selected g =======
best_RMSErun_mean = NaN(nsys, ncase);
best_RMSErun_med  = NaN(nsys, ncase);

fprintf(['==== Best g selection: minimize RMSE(mean traj vs reference). ', ...
         'Convergence time: spread-only tail (eps_spread=%.3g). ', ...
         'After selecting g: compute median RMSE across runs (vs reference). ====\n'], eps_spread);

for ss = 1:nsys
    sysnum = sys_list(ss);

    for cc = 1:ncase
        si = SIs(cc);

        RMSEbar_by_g = NaN(1, numel(g_list));
        Tspr_by_g    = NaN(1, numel(g_list));
        RMSErun_med_by_g  = NaN(1, numel(g_list));
        RMSErun_mean_by_g = NaN(1, numel(g_list));

        for gg = 1:numel(g_list)
            g = g_list(gg);

            fname_big = sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                                sysnum, f_target, si, g, extra_suffix);

            if ~isfile(fname_big)
                continue;
            end

            Y_all = cell(1, nruns);
            lens  = zeros(1, nruns);

            for i = 1:nruns
                varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
                try
                    S = load(fname_big, varname);
                    if ~isfield(S, varname), continue; end
                    rec = S.(varname);

                    y = [];
                    for fn = y_fields
                        if isfield(rec, fn{1}) && ~isempty(rec.(fn{1}))
                            y = rec.(fn{1});
                            break;
                        end
                    end
                    if isempty(y), continue; end

                    y = squeeze(y);
                    y = y(:,:);          % [T x p]

                    Y_all{i} = y;
                    lens(i)  = size(y,1);
                catch
                end
            end

            % ---- selection metric: RMSE(mean traj vs reference) ----
            RMSEbar_by_g(gg) = local_rmse_meantraj_vs_ref(Y_all, lens, Fs, f_target);

            % ---- convergence: spread-only (NO reference) ----
            Tspr_by_g(gg) = local_convtime_spread_only(Y_all, lens, Fs, eps_spread);

            % ---- NEW: run-wise RMSE stats (vs reference) at this g ----
            [rmse_mean, rmse_med] = local_rmse_runs_vs_ref(Y_all, lens, Fs, f_target);
            RMSErun_mean_by_g(gg) = rmse_mean;
            RMSErun_med_by_g(gg)  = rmse_med;
        end

        % ======= pick g: minimize RMSEbar =======
        finite_rmsebar = isfinite(RMSEbar_by_g);
        if ~any(finite_rmsebar)
            fprintf('sys=%d case=%d | no valid RMSEbar -> skip\n', sysnum, si);
            continue;
        end

        candidates = RMSEbar_by_g;
        candidates(~finite_rmsebar) = inf;

        minRMSEbar = min(candidates);
        idxs = find(candidates == minRMSEbar);

        % tie-breaker: smaller Tspr if finite; else smaller RMSErun_mean; else first
        if numel(idxs) > 1
            ts = Tspr_by_g(idxs);
            if any(isfinite(ts))
                ts(~isfinite(ts)) = inf;
                [~, k2] = min(ts);
                idx_min = idxs(k2);
                why = 'RMSEbar, tie->min Tspr';
            else
                mr = RMSErun_mean_by_g(idxs);
                mr(~isfinite(mr)) = inf;
                [~, k2] = min(mr);
                idx_min = idxs(k2);
                why = 'RMSEbar, tie->min mean(RMSErun)';
            end
        else
            idx_min = idxs(1);
            why = 'RMSEbar';
        end

        gstar = g_list(idx_min);

        best_g(ss, cc)           = gstar;
        best_lam(ss, cc)         = lambda_list(gstar + 1);

        best_RMSEbar(ss, cc)     = RMSEbar_by_g(idx_min);
        best_Tspr(ss, cc)        = Tspr_by_g(idx_min);

        best_RMSErun_mean(ss,cc) = RMSErun_mean_by_g(idx_min);
        best_RMSErun_med(ss,cc)  = RMSErun_med_by_g(idx_min);

        fprintf(['sys=%d case=%d | best g=%d (lambda=%.3g) | ', ...
                 'RMSEbar=%.4g | Tspr=%s | mean(RMSErun)=%.4g | med(RMSErun)=%.4g | pick=%s\n'], ...
            sysnum, si, gstar, best_lam(ss,cc), ...
            best_RMSEbar(ss,cc), num2str(best_Tspr(ss,cc)), ...
            best_RMSErun_mean(ss,cc), best_RMSErun_med(ss,cc), why);
    end
end

% ======= Output tables =======
rowNames = arrayfun(@(s) sprintf('Sys%d', s), sys_list, 'UniformOutput', false);
colNames = arrayfun(@(c) sprintf('case_%d', c), SIs, 'UniformOutput', false);

TRMSEbar = array2table(best_RMSEbar, 'RowNames', rowNames, 'VariableNames', colNames);
TTspr    = array2table(best_Tspr,    'RowNames', rowNames, 'VariableNames', colNames);
Tg       = array2table(best_g,       'RowNames', rowNames, 'VariableNames', colNames);
Tlam     = array2table(best_lam,     'RowNames', rowNames, 'VariableNames', colNames);

TRMSErun_mean = array2table(best_RMSErun_mean, 'RowNames', rowNames, 'VariableNames', colNames);
TRMSErun_med  = array2table(best_RMSErun_med,  'RowNames', rowNames, 'VariableNames', colNames);

fprintf('\n=== RMSEbar table: RMSE(mean trajectory vs reference) ===\n');
disp(TRMSEbar);

fprintf('=== Tspr table: spread-only convergence time (seconds), eps_spread=%.3g ===\n', eps_spread);
disp(TTspr);

fprintf('=== mean(RMSErun) table: mean across runs, RMSE(run vs reference) ===\n');
disp(TRMSErun_mean);

fprintf('=== med(RMSErun) table: median across runs, RMSE(run vs reference) ===\n');
disp(TRMSErun_med);

fprintf('=== Best g table ===\n');
disp(Tg);

fprintf('=== Best lambda table ===\n');
disp(Tlam);

% ========================= FUNCTIONS =========================

function rmse_bar = local_rmse_meantraj_vs_ref(Y_all, lens, Fs, f_target)
    % Mean trajectory across runs, then RMSE vs reference.
    % Multi-output: RMSE = sqrt( mean( sum( (ybar-r).^2, 2 ) ) ).

    rmse_bar = NaN;

    good = find(~cellfun(@isempty, Y_all) & lens > 0);
    if isempty(good), return; end

    Nmin = min(lens(good));
    p    = size(Y_all{good(1)}, 2);

    Nr = numel(good);
    Ystack = NaN(Nmin, p, Nr);

    for k = 1:Nr
        y = Y_all{good(k)};
        y = y(1:Nmin, :);
        Ystack(:,:,k) = y;
    end

    ybar = mean(Ystack, 3, 'omitnan');   % [Nmin x p]

    n = (0:Nmin-1).';
    r = sin(2*pi*f_target/Fs * n);
    if p > 1
        r = repmat(r, 1, p);
    end

    diff = ybar - r;

    if p == 1
        rmse_bar = sqrt(mean(diff.^2, 'omitnan'));
    else
        rmse_bar = sqrt(mean(sum(diff.^2, 2), 'omitnan'));
    end
end

function [rmse_mean, rmse_med] = local_rmse_runs_vs_ref(Y_all, lens, Fs, f_target)
    % RMSE per run vs reference, then take mean + median across runs.
    % Multi-output RMSE definition per run:
    %   rmse_i = sqrt( mean( sum( (y_i-r).^2, 2 ) ) )

    rmse_mean = NaN;
    rmse_med  = NaN;

    good = find(~cellfun(@isempty, Y_all) & lens > 0);
    if isempty(good), return; end

    Nmin = min(lens(good));
    p    = size(Y_all{good(1)}, 2);

    n = (0:Nmin-1).';
    r = sin(2*pi*f_target/Fs * n);
    if p > 1
        r = repmat(r, 1, p);
    end

    rmse_vec = NaN(1, numel(good));

    for k = 1:numel(good)
        y = Y_all{good(k)};
        y = y(1:Nmin, :);

        diff = y - r;

        if p == 1
            rmse_vec(k) = sqrt(mean(diff.^2, 'omitnan'));
        else
            rmse_vec(k) = sqrt(mean(sum(diff.^2, 2), 'omitnan'));
        end
    end

    rmse_vec = rmse_vec(isfinite(rmse_vec));
    if isempty(rmse_vec), return; end

    rmse_mean = mean(rmse_vec);
    rmse_med  = median(rmse_vec);
end

function Tspr = local_convtime_spread_only(Y_all, lens, Fs, eps_spread)
    % Spread-only convergence time (NO reference):
    % Delta(k) = max_j y_j(k) - min_j y_j(k) (per output dim), then take max across dims.
    % Converged at earliest k_spr such that Delta(k) <= eps_spread for all k >= k_spr.
    % Return Tspr = (k_spr-1)/Fs in seconds.

    Tspr = NaN;

    good = find(~cellfun(@isempty, Y_all) & lens>0);
    if isempty(good), return; end

    Nmin = min(lens(good));
    p    = size(Y_all{good(1)}, 2);

    Delta = zeros(Nmin,1);

    for k = 1:Nmin
        yk = NaN(numel(good), p);
        for idx = 1:numel(good)
            y = Y_all{good(idx)}(1:Nmin, :);
            yk(idx,:) = y(k,:);
        end

        range_dim = max(yk, [], 1) - min(yk, [], 1);  % 1 x p
        Delta(k)  = max(range_dim);                   % scalar
    end

    below   = (Delta <= eps_spread);
    tail_ok = flipud(cummin(flipud(double(below)))) == 1;
    k_spr   = find(tail_ok, 1, 'first');

    if ~isempty(k_spr)
        Tspr = (k_spr - 1)/Fs;
    end
end
