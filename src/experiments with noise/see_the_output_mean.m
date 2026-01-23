clc; clear; close all;

% ======= Configuration (same as your original) =======
sysnum    = 6;
f_target  = 2;
SIs       = 0:3;
nruns     = 20;

Fs      = 100;     % Sampling frequency
eps_tol = 0.1;     % Convergence threshold ε for selecting g (max distance to reference)

g_list       = 0:3;
extra_suffix = '_se1e-02';

% g -> lambda_g (same as your original)
lambda_list = [0.001 0.01 0.1 1];

% Candidate y fields (for robust convergence-time reading, same as your original)
y_fields = {'y_seq_opt','y_seq','y','y_out'};

% ======= Teacher style (your plt_set logic) =======
plt_set.y_font_dim     = 16;
plt_set.x_font_dim     = 16;
plt_set.title_font_dim = 10.5; 
plt_set.thick          = 1.0;      % Thinner lines
plt_set.ref_thick      = 1.0;      % Thinner reference
plt_set.plot_unit      = 'centimeters';
plt_set.fontname       = 'Times';
plt_set.fontsize       = 12;       % Tick labels
plt_set.plot_dim_x     = 9;        % Single figure, single column
plt_set.plot_dim_y     = 6;

% ======= Case colors (consistent with previous bar charts: case=0..3) =======
case_hex_all = {'#007191', '#62c8d3', '#f47a00', '#d31f11'};
use_cases = SIs(:).';

use_rgb = zeros(numel(use_cases),3);
for k = 1:numel(use_cases)
    si = use_cases(k);
    use_rgb(k,:) = hex2rgb(case_hex_all{si+1});
end

%% ======= Stage 1: search best g for each case =======
num_cases = numel(SIs);
best_g    = NaN(1, num_cases);
best_meanRMS = NaN(1, num_cases);
best_convT   = NaN(1, num_cases);

fprintf('==== Search best g per case (priority: min ConvTime; if all NaN: min mean(RMSE)) ====\n');

for c = 1:num_cases
    si = SIs(c);

    meanR_by_g = inf(1, numel(g_list));
    conv_by_g  = NaN(1, numel(g_list));

    for gg = 1:numel(g_list)
        g = g_list(gg);

        fname_big = sprintf('sb%02d_r%02d_case%02d_g%d%s.mat', ...
                            sysnum, f_target, si, g, extra_suffix);

        if ~isfile(fname_big)
            continue;
        end

        RMS_tmp = NaN(nruns, 1);
        Y_all   = cell(1, nruns);
        lens    = zeros(1, nruns);

        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname_big, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec, 'rms_error') && ~isempty(rec.rms_error)
                    RMS_tmp(i) = rec.rms_error;
                end

                y = [];
                for fn = y_fields
                    if isfield(rec, fn{1}) && ~isempty(rec.(fn{1}))
                        y = rec.(fn{1}); break;
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

        v = RMS_tmp(~isnan(RMS_tmp));
        if ~isempty(v)
            meanR_by_g(gg) = mean(v);
        end

        conv_by_g(gg) = local_convtime_maxoverruns(Y_all, lens, Fs, f_target, eps_tol);
    end

    finite_conv = isfinite(conv_by_g);
    if any(finite_conv)
        conv_candidates = conv_by_g;
        conv_candidates(~finite_conv) = inf;

        minConv = min(conv_candidates);
        idxs = find(conv_candidates == minConv);

        if numel(idxs) > 1
            [~, k2] = min(meanR_by_g(idxs));
            idx_min = idxs(k2);
        else
            idx_min = idxs(1);
        end
        why = 'ConvTime';
    else
        [~, idx_min] = min(meanR_by_g);
        why = 'RMSE(fallback)';
    end

    gstar = g_list(idx_min);
    best_g(c)       = gstar;
    best_meanRMS(c) = meanR_by_g(idx_min);
    best_convT(c)   = conv_by_g(idx_min);

    lam = lambda_list(gstar+1);
    fprintf('case=%d | best g=%d (λg=%.3g) | ConvTime=%s | mean(RMSE)=%.4g | pick=%s\n', ...
        si, gstar, lam, num2str(best_convT(c)), best_meanRMS(c), why);
end

%% ======= Stage 2: read y_seq_opt under best g -> pointwise mean =======
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

    fname_best = sprintf('sb%02d_r%02d_case%02d_g%d%s.mat', ...
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

            if ~isfield(rec, 'y_seq_opt') || isempty(rec.y_seq_opt)
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

    fprintf('case=%d | best g=%d (λg=%.3g) | used %d/%d runs | N=%d\n', ...
        si, g, lambda_list(g+1), n_used(c), nruns, N);
end

valid = find(~cellfun(@isempty, Ymean_all) & isfinite(N_all));
if isempty(valid)
    error('No case produced a mean trajectory (check files/variable names/y_seq_opt field).');
end

Nmin = min(N_all(valid));

% ======= Unified x-axis in samples (NOT seconds) =======
t = (0:Nmin-1).';

%% ======= Plot: mean trajectories of all cases (single figure) =======
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

% ======= Add reference (solid black line) =======
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

%% ======= Save to ./figures (FIG + PNG); filename includes extra_suffix =======
outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

suffix_clean = regexprep(extra_suffix, '^_+', '');
suffix_clean = regexprep(suffix_clean, '[^A-Za-z0-9\.\-]+', '_');

base = sprintf('MeanTraj_AllCases_sys%02d_ft%02d_eps%.2g_%s', sysnum, f_target, eps_tol, suffix_clean);

savefig(fig, fullfile(outdir, [base, '.fig']));
exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');

fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

%% ======= Convergence time: "max over runs" =======
function Tconv = local_convtime_maxoverruns(Y_all, lens, Fs, f_target, eps_tol)
    Tconv = NaN;

    good = find(~cellfun(@isempty, Y_all) & lens>0);
    if isempty(good), return; end

    Nmin = min(lens(good));
    p    = size(Y_all{good(1)}, 2);

    n = (0:Nmin-1).';
    r = sin(2*pi*f_target/Fs * n);
    if p > 1
        r = repmat(r, 1, p);
    end

    dmax = zeros(Nmin,1);
    for k = 1:Nmin
        dm = 0;
        for idx = 1:numel(good)
            y = Y_all{good(idx)}(1:Nmin,:);
            diff = y(k,:) - r(k,:);
            if p == 1
                dij = abs(diff);
            else
                dij = norm(diff, 2);
            end
            dm = max(dm, dij);
        end
        dmax(k) = dm;
    end

    below   = (dmax <= eps_tol);
    tail_ok = flipud(cummin(flipud(double(below)))) == 1;
    n0      = find(tail_ok, 1, 'first');

    if ~isempty(n0)
        Tconv = (n0 - 1)/Fs;
    end
end

%% ======= Local: hex -> rgb =======
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
