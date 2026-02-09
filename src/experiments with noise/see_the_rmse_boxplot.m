clc; clear; close all;

% ======= Configuration =======
sysnum   = 6;
f_target = 2;
SIs      = 0:3;
nruns    = 20;

Fs     = 100;
g_list = 0:3;

SNRdB_target = 30;
extra_suffix = sprintf('_snr%02ddB', round(SNRdB_target));

lambda_list = [0.001 0.01 0.1 1];

num_cases     = numel(SIs);
best_g        = NaN(1, num_cases);
best_RMS_cols = cell(1, num_cases);
best_RMSEbar  = NaN(1, num_cases);

fprintf('==== Search best g per case (ONLY: minimum RMSEbar(mean(y_seq_opt) vs ref)) ====\n');

% ======= 1) pick best g =======
for c = 1:num_cases
    si = SIs(c);

    rmsebar_by_g = NaN(1, numel(g_list));
    rms_by_g     = cell(1, numel(g_list));

    for gg = 1:numel(g_list)
        g = g_list(gg);

        fname_big = sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                            sysnum, f_target, si, g, extra_suffix);

        RMS_tmp = NaN(nruns, 1);

        Ysopt   = cell(1, nruns);
        lensOpt = NaN(1, nruns);

        if ~isfile(fname_big)
            rms_by_g{gg}     = RMS_tmp;
            rmsebar_by_g(gg) = NaN;
            continue;
        end

        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname_big, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                % ----- RMSE data for boxplot -----
                if isfield(rec, 'rms_error') && ~isempty(rec.rms_error)
                    RMS_tmp(i) = mean(rec.rms_error(:));
                end

                % ----- y_seq_opt for RMSEbar selection -----
                if isfield(rec, 'y_seq_opt') && ~isempty(rec.y_seq_opt)
                    yy = rec.y_seq_opt(:);
                    Ysopt{i}   = yy;
                    lensOpt(i) = numel(yy);
                end
            catch
            end
        end

        rms_by_g{gg}     = RMS_tmp;
        rmsebar_by_g(gg) = local_rmsebar_meantraj_vs_ref(Ysopt, lensOpt, Fs, f_target);
    end

    finite_bar = isfinite(rmsebar_by_g);
    if ~any(finite_bar)
        fprintf('case=%d | no valid RMSEbar across g (skip)\n', si);
        continue;
    end

    cand = rmsebar_by_g;
    cand(~finite_bar) = inf;

    [minBar, idx_min] = min(cand);
    gstar = g_list(idx_min);

    best_g(c)        = gstar;
    best_RMS_cols{c} = rms_by_g{idx_min};
    best_RMSEbar(c)  = minBar;

    fprintf('case=%d | best g=%d (lambda=%.3g) | RMSEbar=%.6g\n', ...
        si, gstar, lambda_list(gstar+1), minBar);
end

% ======= 2) pack data for boxplot =======
Data   = NaN(nruns, num_cases);

labels_map = {'$\mathbf{WN}$', '$\mathbf{IBW}$', '$\mathbf{IBN}$', '$\mathbf{OB}$'};
labels = cell(1, num_cases);

for c = 1:num_cases
    if isempty(best_RMS_cols{c}), continue; end
    Data(:, c) = best_RMS_cols{c}(:);

    si = SIs(c);
    if si >= 0 && si <= 3
        labels{c} = labels_map{si+1};
    else
        labels{c} = sprintf('case=%d', si);
    end
end

% ======= Plot style (keep your original) =======
plt_set.y_font_dim     = 16;
plt_set.x_font_dim     = 11;
plt_set.thick          = 1.6;
plt_set.plot_unit      = 'centimeters';
plt_set.fontname       = 'Times';
plt_set.fontsize       = 12;
plt_set.plot_dim_x     = 8;
plt_set.plot_dim_y     = 6;

case_hex = {'#007191', '#62c8d3', '#f47a00', '#d31f11'};
use_rgb = zeros(num_cases,3);
for k = 1:num_cases
    use_rgb(k,:) = hex2rgb(case_hex{k});
end

fig = figure('Color','w','Name','RMSE boxplot (g picked by RMSEbar)');
fig.Units = plt_set.plot_unit;
fig.Position(3) = plt_set.plot_dim_x;
fig.Position(4) = plt_set.plot_dim_y;
set(fig,'Renderer','opengl');

boxplot(Data, 'Labels', labels, 'Symbol','k+');
ax = gca;

ax.FontName  = plt_set.fontname;
ax.FontSize  = plt_set.fontsize;
ax.LineWidth = 1.0;

xlabel('');

ylab = ylabel('$\bigl\{\mathrm{RMSE}_y^{j}\bigr\}_{j=1}^{20}$', 'Interpreter', 'latex');
set(ylab, 'FontName', plt_set.fontname, 'FontSize', plt_set.y_font_dim);

grid on;
ax.XGrid = 'off';
ax.YGrid = 'on';
ax.GridAlpha = 0.15;

thick = plt_set.thick;

hBox  = findobj(ax,'Tag','Box');
hMed  = findobj(ax,'Tag','Median');
hWhi  = findobj(ax,'Tag','Whisker');
hOut  = findobj(ax,'Tag','Outliers');
hCapU = findobj(ax,'Tag','Upper Adjacent Value');
hCapL = findobj(ax,'Tag','Lower Adjacent Value');

if ~isempty(hBox),  set(hBox,  'LineWidth', thick); end
if ~isempty(hMed),  set(hMed,  'LineWidth', thick); end
if ~isempty(hWhi),  set(hWhi,  'LineWidth', thick); end
if ~isempty(hOut),  set(hOut,  'LineWidth', 1.0);  end
if ~isempty(hCapU), set(hCapU, 'LineWidth', thick); end
if ~isempty(hCapL), set(hCapL, 'LineWidth', thick); end

for k = 1:num_cases
    kk = num_cases - k + 1;

    if numel(hBox) >= k, set(hBox(kk), 'Color', use_rgb(k,:)); end
    if numel(hMed) >= k, set(hMed(kk), 'Color', use_rgb(k,:)); end

    idxw = (2*(kk-1)+1):(2*(kk-1)+2);
    if numel(hWhi) >= max(idxw)
        set(hWhi(idxw), 'Color', use_rgb(k,:));
    end

    if ~isempty(hCapU)
        idxc = (2*(kk-1)+1):(2*(kk-1)+2);
        if numel(hCapU) >= max(idxc)
            set(hCapU(idxc), 'Color', use_rgb(k,:));
        end
    end
    if ~isempty(hCapL)
        idxc2 = (2*(kk-1)+1):(2*(kk-1)+2);
        if numel(hCapL) >= max(idxc2)
            set(hCapL(idxc2), 'Color', use_rgb(k,:));
        end
    end

    if numel(hOut) >= k
        set(hOut(kk), 'MarkerEdgeColor', use_rgb(k,:), 'MarkerFaceColor', use_rgb(k,:));
    end
end

set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));
fig.PaperPositionMode = 'auto';
drawnow;

% --- Make x tick labels bold (TeX interpreter needed for \textbf) ---
set(ax, 'TickLabelInterpreter', 'latex');

% ======= Save =======
outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

suffix_clean = regexprep(extra_suffix, '^_+', '');
suffix_clean = regexprep(suffix_clean, '[^A-Za-z0-9\.\-]+', '_');

base = sprintf('Boxplot_RMSE_bestGbyRMSEbar_sys%02d_ft%02d_%s', sysnum, f_target, suffix_clean);

savefig(fig, fullfile(outdir,[base,'.fig']));
exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');

fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

fprintf('==== Best-setting statistics (by case) ====\n');
for c = 1:num_cases
    v = Data(:, c); v = v(~isnan(v));
    if ~isempty(v) && isfinite(best_g(c))
        lam = lambda_list(best_g(c)+1);
        fprintf(['case=%d | g=%d (lambda=%.3g) | RMSEbar=%.6g | N=%d | ', ...
                 'mean=%.4g | median=%.4g | min=%.4g | max=%.4g\n'], ...
            SIs(c), best_g(c), lam, best_RMSEbar(c), numel(v), ...
            mean(v), median(v), min(v), max(v));
    else
        fprintf('case=%d | no valid RMSE data | RMSEbar=%s\n', ...
            SIs(c), num2str(best_RMSEbar(c)));
    end
end

%% ========================= FUNCTIONS =========================

function rmse_bar = local_rmsebar_meantraj_vs_ref(Ysopt, lensOpt, Fs, f_target)
    % RMSEbar = RMSE( mean(y_seq_opt across runs) - reference )

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
    if h(1) == '#', h = h(2:end); end
    if numel(h) ~= 6
        error('hex2rgb:InvalidHex', 'Hex color must be 6 characters, got: %s', h);
    end
    r = hex2dec(h(1:2));
    g = hex2dec(h(3:4));
    b = hex2dec(h(5:6));
    rgb = [r g b] / 255;
end

