clc; clear; close all;

sys_list = 1:5;
f_target = 2;
SIs      = 0:3;
nruns    = 50;
Fs       = 100;

eps_track = 0.3;

y_fields = {'y_seq_opt','y_seq','y','y_out'};

warning('off','MATLAB:print:InvalidGraphicsState');
warning('off','MATLAB:graphics:SceneNode');

nsys  = numel(sys_list);
ncase = numel(SIs);

best_Ttrk         = NaN(nsys, ncase);
best_RMSErun_mean = NaN(nsys, ncase);


for ss = 1:nsys
    sysnum = sys_list(ss);

    for cc = 1:ncase
        si = SIs(cc);

        data_dir  = fullfile(pwd, 'data');
        fname_sel = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d.mat', ...
                            sysnum, f_target, si));

        if ~isfile(fname_sel)
            fprintf('sys=%d case=%d | file not found -> skip\n', sysnum, si);
            continue;
        end

        [Y_all, lens] = local_load_Yall(fname_sel, sysnum, f_target, si, nruns, y_fields);

        best_RMSErun_mean(ss, cc) = local_rmse_runs_vs_ref_mean(Y_all, lens, Fs, f_target);
        best_Ttrk(ss, cc)         = local_convtime_track_to_ref(Y_all, lens, Fs, f_target, eps_track);

        fprintf('sys=%d case=%d | mean(RMSErun)=%s | Ttrk=%s\n', ...
            sysnum, si, ...
            fmt4f(best_RMSErun_mean(ss,cc)), fmt4f(best_Ttrk(ss,cc)));
    end
end

rowNames = arrayfun(@(s) sprintf('Sys%d', s), sys_list, 'UniformOutput', false);
colNames = arrayfun(@(c) sprintf('case_%d', c), SIs, 'UniformOutput', false);

TTtrk         = array2table(best_Ttrk,         'RowNames', rowNames, 'VariableNames', colNames);
TRMSErun_mean = array2table(best_RMSErun_mean, 'RowNames', rowNames, 'VariableNames', colNames);

fprintf('\n=== Ttrk table ===\n');
disp(table_format_4f(TTtrk));

fprintf('=== mean(RMSErun) table ===\n');
disp(table_format_4f(TRMSErun_mean));

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