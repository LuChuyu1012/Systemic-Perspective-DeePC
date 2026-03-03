clc; clear; close all;

% ======= Configuration =======
sysnum   = 5;
f_target = 2;
SIs      = 0:3;
nruns    = 50;

Fs = 100;

% ======= Candidate g list (SELECT BEST g by mean of 50 RMSErun) =======
g_list = [0 1 2 3 5];

% ======= Noise tag (must match saved files) =======
suffix_list = {'_sigmae5e-02','_sigmae1e-02','_snr30dB'};

% ======= g -> lambda mapping (so g=5 prints 100) =======
lambda_map = containers.Map( ...
    {0,    1,    2,   3,   5}, ...
    {0.001,0.01,0.1, 1,   100} );

% ======= data folder (CURRENT PATH /data) =======
data_dir = fullfile(pwd, 'data');

plt_set.y_font_dim     = 16;
plt_set.x_font_dim     = 16;
plt_set.thick          = 1.0;
plt_set.plot_unit      = 'centimeters';
plt_set.fontname       = 'Times';
plt_set.fontsize       = 12;
plt_set.plot_dim_x     = 9;
plt_set.plot_dim_y     = 6;

case_hex_all = {'#007191', '#62c8d3', '#f47a00', '#d31f11'};
use_cases = SIs(:).';

use_rgb = zeros(numel(use_cases),3);
for k = 1:numel(use_cases)
    si = use_cases(k);
    use_rgb(k,:) = hex2rgb(case_hex_all{si+1});
end

num_cases = numel(SIs);

for mm = 1:numel(suffix_list)

    extra_suffix = suffix_list{mm};

    best_g         = NaN(1, num_cases);
    best_lam       = NaN(1, num_cases);
    best_meanRMSE  = NaN(1, num_cases);
    best_Nvalid    = zeros(1, num_cases);

    fprintf('==== Select best g per case (min mean of %d run-wise RMSE) ====\n', nruns);

    % ======= 1) select best g per case: minimize mean(rmse_runs) across runs =======
    for c = 1:num_cases
        si = SIs(c);

        mean_rmse_by_g = NaN(1, numel(g_list));
        Nvalid_by_g    = zeros(1, numel(g_list));

        for gg = 1:numel(g_list)
            g = g_list(gg);

            fname_big = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                                sysnum, f_target, si, g, extra_suffix));
            if ~isfile(fname_big)
                continue;
            end

            rmse_runs = NaN(1, nruns);

            for i = 1:nruns
                varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
                try
                    S = load(fname_big, varname);
                    if ~isfield(S, varname), continue; end
                    rec = S.(varname);

                    if isfield(rec,'rms_error') && ~isempty(rec.rms_error)
                        rmse_runs(i) = mean(rec.rms_error(:));
                    end
                catch
                end
            end

            v = rmse_runs(isfinite(rmse_runs));
            Nvalid_by_g(gg) = numel(v);

            if isempty(v)
                continue;
            end

            mean_rmse_by_g(gg) = mean(v);

            if Nvalid_by_g(gg) ~= nruns
                fprintf('case=%d g=%d | valid rms_error=%d/%d | %s\n', si, g, Nvalid_by_g(gg), nruns, extra_suffix);
            end
        end

        finite_mean = isfinite(mean_rmse_by_g);
        if ~any(finite_mean)
            fprintf('case=%d | no valid rms_error across g (skip) | %s\n', si, extra_suffix);
            continue;
        end

        cand = mean_rmse_by_g;
        cand(~finite_mean) = inf;

        minMean = min(cand);
        idxs = find(cand == minMean);
        idx_min = idxs(1);

        gstar = g_list(idx_min);
        best_g(c)        = gstar;
        best_meanRMSE(c) = mean_rmse_by_g(idx_min);
        best_Nvalid(c)   = Nvalid_by_g(idx_min);

        if isKey(lambda_map, gstar)
            lam = lambda_map(gstar);
        else
            lam = NaN;
        end
        best_lam(c) = lam;

        fprintf('case=%d | best g=%d (lambda=%.3g) | meanRMSE=%g | valid=%d/%d | %s\n', ...
            si, gstar, lam, best_meanRMSE(c), best_Nvalid(c), nruns, extra_suffix);
    end

    % ======= 2) mean y_seq_opt for each case using best g =======
    Ymean_all = cell(1, num_cases);
    N_all     = NaN(1, num_cases);
    n_used    = zeros(1, num_cases);

    for c = 1:num_cases
        si = SIs(c);
        g  = best_g(c);

        if ~isfinite(g)
            warning('case=%d has no best g, skip.', si);
            continue;
        end

        fname_best = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                             sysnum, f_target, si, g, extra_suffix));
        if ~isfile(fname_best)
            warning('File not found: %s (case=%d), skip.', fname_best, si);
            continue;
        end

        Ys   = cell(1, nruns);
        lens = NaN(1, nruns);

        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname_best, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec,'y_seq_opt') && ~isempty(rec.y_seq_opt)
                    y = rec.y_seq_opt(:);
                    Ys{i}   = y;
                    lens(i) = numel(y);
                end
            catch
            end
        end

        good = find(~cellfun(@isempty, Ys) & isfinite(lens) & lens>0);
        if isempty(good)
            warning('case=%d | g=%d: no y_seq_opt found, skip.', si, g);
            continue;
        end

        N = min(lens(good));
        Ymat = NaN(N, numel(good));
        for j = 1:numel(good)
            y = Ys{good(j)};
            Ymat(:,j) = y(1:N);
        end

        Ymean_all{c} = mean(Ymat, 2, 'omitnan');
        N_all(c)     = N;
        n_used(c)    = numel(good);

        if isKey(lambda_map, g)
            lam = lambda_map(g);
        else
            lam = NaN;
        end

        fprintf('case=%d | best g=%d (lambda=%.3g) | used %d/%d runs | N=%d | %s\n', ...
            si, g, lam, n_used(c), nruns, N, extra_suffix);
    end

    valid = find(~cellfun(@isempty, Ymean_all) & isfinite(N_all));
    if isempty(valid)
        error('No case produced mean y_seq_opt (check files/y_seq_opt field).');
    end

    Nmin = min(N_all(valid));
    t = (0:Nmin-1).';

    % ======= Reference =======
    r = sin(2*pi*f_target/Fs * t);

    % ======= Plot =======
    fig = figure('Color','w','Name','Mean output (all cases)');
    fig.Units = plt_set.plot_unit;
    fig.Position(3) = plt_set.plot_dim_x;
    fig.Position(4) = plt_set.plot_dim_y;
    fig.Renderer = 'opengl';

    ax = gca; hold(ax,'on');

    for c = 1:num_cases
        if isempty(Ymean_all{c}), continue; end
        y_plot = Ymean_all{c}(1:Nmin);
        plot(t, y_plot, '-', 'Color', use_rgb(c,:), 'LineWidth', plt_set.thick);
    end

    plot(t, r, 'k-', 'LineWidth', plt_set.thick);

    ax.FontName  = plt_set.fontname;
    ax.FontSize  = plt_set.fontsize;
    ax.LineWidth = 1.0;

    xlab = xlabel('t [samples]');
    ylab = ylabel(ax, '$\bar{y}(t)$', 'Interpreter','latex');
    set(xlab, 'FontName', plt_set.fontname, 'FontSize', plt_set.x_font_dim);
    set(ylab, 'FontName', plt_set.fontname, 'FontSize', plt_set.y_font_dim);

    grid(ax,'on');
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    ax.GridAlpha = 0.15;
    xlim([t(1) t(end)]);

    set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));
    fig.PaperPositionMode = 'auto';
    drawnow;

    % ======= Save =======
    % outdir = fullfile(pwd, 'figures');
    % if ~exist(outdir,'dir'); mkdir(outdir); end
    % 
    % suffix_clean = regexprep(extra_suffix, '^_+', '');
    % suffix_clean = regexprep(suffix_clean, '[^A-Za-z0-9\.\-]+', '_');
    % 
    % base = sprintf('MeanY_seqopt_AllCases_sys%02d_ft%02d_bestGbyMeanRMSE_%s', ...
    %                sysnum, f_target, suffix_clean);
    % 
    % savefig(fig, fullfile(outdir, [base, '.fig']));
    % exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');
    % 
    % fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
    % fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));
end

%% ========================= FUNCTIONS =========================
function rgb = hex2rgb(h)
    if isstring(h), h = char(h); end
    h = strtrim(h);
    if isempty(h), error('hex2rgb:Empty', 'Empty hex string.'); end
    if h(1) == '#', h = h(2:end); end
    if numel(h) ~= 6
        error('hex2rgb:InvalidHex', 'Hex color must be 6 characters, got: %s', h);
    end
    r = hex2dec(h(1:2));
    g = hex2dec(h(3:4));
    b = hex2dec(h(5:6));
    rgb = [r g b] / 255;
end