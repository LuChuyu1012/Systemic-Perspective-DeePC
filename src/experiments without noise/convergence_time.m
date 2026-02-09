clc; clear; close all;

% ======= Select ONE system here =======
sysnum    = 3;          % <<< change to 1..6
f_target  = 2;
SIs       = 7:9;
nruns     = 20;
Fs        = 100;

consensus_tol_abs = 0.1;   % threshold for (max-min) across runs
win_periods_chk   = 0;

P       = max(1, round(Fs / max(f_target, eps)));
win_chk = max(1, round(win_periods_chk * P));

y_fields = {'y_seq_opt','y_seq','y','y_out'};

ConvTime = NaN(1, numel(SIs));

for sidx = 1:numel(SIs)
    si = SIs(sidx);

    fname_big = sprintf('s%02d_r%02d_case%02d.mat', sysnum, f_target, si);
    if ~isfile(fname_big)
        warning('error: %s (skip sys=%d case=%d)', fname_big, sysnum, si);
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
            y = y(:,:);                 % [T x ny]

            lens(i)  = size(y,1);
            Y_all{i} = y;

        catch ME
            warning('error: %s | %s', varname, ME.message);
        end
    end

    good = find(~cellfun(@isempty, Y_all) & lens>0);
    if isempty(good)
        warning('sys=%d case=%d invalid y', sysnum, si);
        continue;
    end

    Nmin = min(lens(good));
    ny   = size(Y_all{good(1)}, 2);

    % For each time k, compute range per output dim:
    %   range_j(k) = max_r y_r(k,j) - min_r y_r(k,j)
    % Then aggregate into scalar:
    %   range(k) = max_j range_j(k)
    range_t = zeros(Nmin,1);

    for k = 1:Nmin
        yk = NaN(numel(good), ny);
        for g = 1:numel(good)
            y = Y_all{good(g)}(1:Nmin, :);
            yk(g,:) = y(k,:);
        end
        range_dim = max(yk, [], 1) - min(yk, [], 1);   % [1 x ny]
        range_t(k) = max(range_dim);                   % scalar
    end

    if win_chk > 1
        rchk = movmax(range_t, [win_chk-1, 0], 'Endpoints','shrink');
    else
        rchk = range_t;
    end

    below = (rchk <= consensus_tol_abs);

    tail_ok = flipud(cummin(flipud(double(below)))) == 1;
    n0 = find(tail_ok, 1, 'first');

    if ~isempty(n0)
        ConvTime(1, sidx) = (n0 - 1)/Fs;
    else
        ConvTime(1, sidx) = NaN;
    end
end

fprintf('\n=== Consensus time by range (sys=%d, f=%g Hz, tol=%g) ===\n', ...
    sysnum, f_target, consensus_tol_abs);

for sidx = 1:numel(SIs)
    si = SIs(sidx);
    if isnan(ConvTime(sidx))
        fprintf('case=%d : NaN (no convergence)\n', si);
    else
        fprintf('case=%d : %.6f s\n', si, ConvTime(sidx));
    end
end

T = table(SIs(:), ConvTime(:), 'VariableNames', {'case','conv_time_s'});
disp(T);
