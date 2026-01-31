clc; clear; close all;

% ======= Configuration =======
sysnum    = 6;
f_target  = 2;
SIs       = 0:3;
nruns     = 20;

Fs = 100;

g_list = 0:3;

SNRdB_target = 20;
extra_suffix = sprintf('_snr%02ddB', round(SNRdB_target));

lambda_list = [0.001 0.01 0.1 1];

plt_set.y_font_dim     = 16;
plt_set.x_font_dim     = 16;
plt_set.thick          = 1.0;
plt_set.ref_thick      = 1.0;
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

num_cases   = numel(SIs);
best_g      = NaN(1, num_cases);
best_rmsebar = NaN(1, num_cases);

fprintf('==== Search best g per case (ONLY: min RMSEbar(mean(y_seq_opt) vs ref)) ====\n');

% ======= 1) best g per case =======
for c = 1:num_cases
    si = SIs(c);

    rmsebar_by_g = NaN(1, numel(g_list));

    for gg = 1:numel(g_list)
        g = g_list(gg);

        fname_big = sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                            sysnum, f_target, si, g, extra_suffix);

        if ~isfile(fname_big)
            continue;
        end

        Ysopt   = cell(1, nruns);
        lensOpt = NaN(1, nruns);

        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname_big, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec,'y_seq_opt') && ~isempty(rec.y_seq_opt)
                    yy = rec.y_seq_opt(:);
                    Ysopt{i}   = yy;
                    lensOpt(i) = numel(yy);
                end
            catch
            end
        end

        rmsebar_by_g(gg) = local_rmsebar_meantraj_vs_ref(Ysopt, lensOpt, Fs, f_target);
    end

    finite_bar = isfinite(rmsebar_by_g);
    if ~any(finite_bar)
        fprintf('case=%d | no valid RMSEbar found across g (skip)\n', si);
        continue;
    end

    cand = rmsebar_by_g;
    cand(~finite_bar) = inf;

    [minBar, idx_min] = min(cand);
    gstar = g_list(idx_min);

    best_g(c)       = gstar;
    best_rmsebar(c) = minBar;

    lam = lambda_list(gstar+1);
    fprintf('case=%d | best g=%d (lambda=%.3g) | RMSEbar=%.6g\n', si, gstar, lam, minBar);
end

% ======= 2) build mean trajectories using best g =======
Ymean_all = cell(1, num_cases);
N_all     = NaN(1, num_cases);
n_used    = zeros(1, num_cases);

for c = 1:num_cases
    si = SIs(c);
    g  = best_g(c);

    if ~isfinite(g)
        warning('case=%d has no best g, skip mean trajectory.', si);
        continue;
    end

    fname_best = sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                         sysnum, f_target, si, g, extra_suffix);

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

            if ~isfield(rec,'y_seq_opt') || isempty(rec.y_seq_opt)
                continue;
            end

            y = rec.y_seq_opt(:);
            Ys{i}   = y;
            lens(i) = numel(y);
        catch
        end
    end

    good = find(~cellfun(@isempty, Ys) & isfinite(lens) & lens>0);
    if isempty(good)
        warning('case=%d | best g=%d: no y_seq_opt found, skip.', si, g);
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

    fprintf('case=%d | best g=%d (lambda=%.3g) | used %d/%d runs | N=%d\n', ...
        si, g, lambda_list(g+1), n_used(c), nruns, N);
end

valid = find(~cellfun(@isempty, Ymean_all) & isfinite(N_all));
if isempty(valid)
    error('No case produced a mean trajectory (check files/variable names/y_seq_opt field).');
end

Nmin = min(N_all(valid));
t = (0:Nmin-1).';

% ======= Plot =======
fig = figure('Color','w','Name','Mean trajectory (all cases)');
fig.Units = plt_set.plot_unit;
fig.Position(3) = plt_set.plot_dim_x;
fig.Position(4) = plt_set.plot_dim_y;
fig.Renderer = 'painters';

ax = gca; hold(ax,'on');

for c = 1:num_cases
    if isempty(Ymean_all{c}), continue; end
    ybar = Ymean_all{c}(1:Nmin);
    plot(t, ybar, '-', 'Color', use_rgb(c,:), 'LineWidth', plt_set.thick);
end

r = sin(2*pi*f_target*(t/Fs));
plot(t, r, 'k-', 'LineWidth', plt_set.ref_thick);

ax.FontName  = plt_set.fontname;
ax.FontSize  = plt_set.fontsize;
ax.LineWidth = 1.0;

xlab = xlabel('t [samples]');
ylab = ylabel('y(t)');
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
outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

suffix_clean = regexprep(extra_suffix, '^_+', '');
suffix_clean = regexprep(suffix_clean, '[^A-Za-z0-9\.\-]+', '_');

base = sprintf('MeanTraj_AllCases_sys%02d_ft%02d_%s', sysnum, f_target, suffix_clean);

savefig(fig, fullfile(outdir, [base, '.fig']));
exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');

fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

%% ======= local functions =======
function rmse_bar = local_rmsebar_meantraj_vs_ref(Ysopt, lensOpt, Fs, f_target)
    % RMSE between mean y_seq_opt across runs and reference r(t)=sin(2*pi*f*t/Fs).

    rmse_bar = NaN;

    good = find(~cellfun(@isempty, Ysopt) & isfinite(lensOpt) & lensOpt>0);
    if isempty(good), return; end

    Nmin = min(lensOpt(good));

    Ymat = NaN(Nmin, numel(good));
    for j = 1:numel(good)
        y = Ysopt{good(j)};
        Ymat(:,j) = y(1:Nmin);
    end

    ybar = mean(Ymat, 2, 'omitnan');

    n = (0:Nmin-1).';
    r = sin(2*pi*f_target/Fs * n);

    diff = ybar - r;
    rmse_bar = sqrt(mean(diff.^2, 'omitnan'));
end

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
