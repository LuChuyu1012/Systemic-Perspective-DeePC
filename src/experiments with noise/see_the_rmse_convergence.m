clc; clear; close all;

sys_list = 1:5;
f_target = 2;
SIs      = 0:3;
nruns    = 50;
Fs       = 100;

eps_track = 0.3;

g_list = [0 1 2 3 5];

suffix_list = {'_sigmae5e-02', '_sigmae1e-02', '_snr30dB'};
label_list  = {'sigmae=0.05', 'sigmae=0.01', 'snr30dB'};

lambda_map = containers.Map( ...
    {0,    1,    2,   3,   5}, ...
    {0.001,0.01,0.1, 1,   100} );

y_fields = {'y_seq_opt','y_seq','y','y_out'};

warning('off','MATLAB:print:InvalidGraphicsState');
warning('off','MATLAB:graphics:SceneNode');

nsys  = numel(sys_list);
ncase = numel(SIs);

for mm = 1:numel(suffix_list)

    extra_suffix = suffix_list{mm};
    label_now    = label_list{mm};

    best_g            = NaN(nsys, ncase);
    best_lam          = NaN(nsys, ncase);
    best_Ttrk         = NaN(nsys, ncase);
    best_RMSErun_mean = NaN(nsys, ncase);

    fprintf('\n============================================================\n');
    fprintf('Dataset: %s   (suffix = %s)\n', label_now, extra_suffix);
    fprintf('============================================================\n');
    fprintf('==== Best g selection: minimize mean(RMSErun) across %d runs. Then report Ttrk at selected g. ====\n', nruns);

    for ss = 1:nsys
        sysnum = sys_list(ss);

        for cc = 1:ncase
            si = SIs(cc);

            RMSErun_mean_by_g = NaN(1, numel(g_list));
            Ttrk_by_g         = NaN(1, numel(g_list));

            for gg = 1:numel(g_list)
                g = g_list(gg);

                data_dir  = fullfile(pwd, 'data');
                fname_sel = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                                    sysnum, f_target, si, g, extra_suffix));

                if ~isfile(fname_sel)
                    continue;
                end

                [Y_all, lens] = local_load_Yall(fname_sel, sysnum, f_target, si, nruns, y_fields);

                RMSErun_mean_by_g(gg) = local_rmse_runs_vs_ref_mean(Y_all, lens, Fs, f_target);
                Ttrk_by_g(gg)         = local_convtime_track_to_ref(Y_all, lens, Fs, f_target, eps_track);
            end

            finite_mean = isfinite(RMSErun_mean_by_g);
            if ~any(finite_mean)
                fprintf('sys=%d case=%d | no valid mean(RMSErun) across g -> skip\n', sysnum, si);
                continue;
            end

            cand = RMSErun_mean_by_g;
            cand(~finite_mean) = inf;

            minMean = min(cand);
            idxs = find(cand == minMean);
            idx_min = idxs(1);

            if numel(idxs) > 1
                why = 'mean(RMSErun), tie->first';
            else
                why = 'mean(RMSErun)';
            end

            gstar = g_list(idx_min);

            if isKey(lambda_map, gstar)
                lam = lambda_map(gstar);
            else
                lam = NaN;
            end

            best_g(ss, cc)            = gstar;
            best_lam(ss, cc)          = lam;
            best_RMSErun_mean(ss, cc) = RMSErun_mean_by_g(idx_min);
            best_Ttrk(ss, cc)         = Ttrk_by_g(idx_min);

            fprintf('sys=%d case=%d | best g=%d (lambda=%s) | mean(RMSErun)=%s | Ttrk=%s | pick=%s\n', ...
                sysnum, si, gstar, fmt4f(lam), ...
                fmt4f(best_RMSErun_mean(ss,cc)), fmt4f(best_Ttrk(ss,cc)), why);
        end
    end

    rowNames = arrayfun(@(s) sprintf('Sys%d', s), sys_list, 'UniformOutput', false);
    colNames = arrayfun(@(c) sprintf('case_%d', c), SIs, 'UniformOutput', false);

    TTtrk         = array2table(best_Ttrk,         'RowNames', rowNames, 'VariableNames', colNames);
    Tg            = array2table(best_g,            'RowNames', rowNames, 'VariableNames', colNames);
    TRMSErun_mean = array2table(best_RMSErun_mean, 'RowNames', rowNames, 'VariableNames', colNames);

    fprintf('\n=== %s : Ttrk table ===\n', label_now);
    disp(table_format_4f(TTtrk));

    fprintf('=== %s : mean(RMSErun) table ===\n', label_now);
    disp(table_format_4f(TRMSErun_mean));

    fprintf('=== %s : Best g table ===\n', label_now);
    disp(Tg);
end

function [Y_all, lens] = local_load_Yall(fname_big, sysnum, f_target, si, nruns, y_fields)
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
            y = y(:,:);

            Y_all{i} = y;
            lens(i)  = size(y,1);
        catch
        end
    end
end

function rmse_mean = local_rmse_runs_vs_ref_mean(Y_all, lens, Fs, f_target)
    rmse_mean = NaN;

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
end

function Ttrk = local_convtime_track_to_ref(Y_all, lens, Fs, f_target, eps_track)
    Ttrk = NaN;

    good = find(~cellfun(@isempty, Y_all) & lens>0);
    if isempty(good), return; end

    Nmin = min(lens(good));
    p    = size(Y_all{good(1)}, 2);

    n = (0:Nmin-1).';
    r = sin(2*pi*f_target/Fs * n);
    if p > 1
        r = repmat(r, 1, p);
    end

    e_max = zeros(Nmin,1);

    for k = 1:Nmin
        e_run = zeros(numel(good),1);
        for idx = 1:numel(good)
            y = Y_all{good(idx)}(1:Nmin, :);
            ek = abs(y(k,:) - r(k,:));
            e_run(idx) = max(ek);
        end
        e_max(k) = max(e_run);
    end

    below   = (e_max <= eps_track);
    tail_ok = flipud(cummin(flipud(double(below)))) == 1;
    k0      = find(tail_ok, 1, 'first');

    if ~isempty(k0)
        Ttrk = (k0 - 1)/Fs;
    end
end

function T4 = table_format_4f(T)
    A = T{:,:};
    S = strings(size(A));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            S(i,j) = fmt4f(A(i,j));
        end
    end

    T4 = array2table(S, 'RowNames', T.Properties.RowNames, 'VariableNames', T.Properties.VariableNames);
end

function s = fmt4f(x)
    if isnan(x)
        s = "NaN";
    elseif isinf(x)
        if x > 0
            s = "Inf";
        else
            s = "-Inf";
        end
    else
        s = sprintf('%.4f', x);
    end
end