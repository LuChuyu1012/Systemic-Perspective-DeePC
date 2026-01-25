clc; clear; close all;

sysnum   = 3;
f_target = 2;
SIs      = 0:3;
nruns    = 20;

Fs      = 100;
eps_tol = 0.2;

g_list = [0 1 2 3];

SNRdB_target = 20;
extra_suffix = sprintf('_snr%02ddB', round(SNRdB_target));

lambda_list = [0.001 0.01 0.1 1];

y_fields = {'y_seq_opt','y_seq','y','y_out'};

warning('off','MATLAB:print:InvalidGraphicsState');
warning('off','MATLAB:graphics:SceneNode');

num_cases     = numel(SIs);
best_g        = NaN(1, num_cases);
best_RMS_cols = cell(1, num_cases);
best_meanRMS  = NaN(1, num_cases);
best_convT    = NaN(1, num_cases);

SNRYN_out     = NaN(1, num_cases);
SNRUN_out     = NaN(1, num_cases);

fprintf('==== Search best g per case (priority: minimum convergence time; if all NaN: minimum mean(RMSE)) ====\n');

for c = 1:num_cases
    si = SIs(c);

    rms_by_g   = cell(1, numel(g_list));
    meanR_by_g = inf(1, numel(g_list));
    conv_by_g  = NaN(1, numel(g_list));

    for gg = 1:numel(g_list)
        g = g_list(gg);

        fname_big = sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                            sysnum, f_target, si, g, extra_suffix);

        RMS_tmp = NaN(nruns, 1);
        Y_all   = cell(1, nruns);
        lens    = zeros(1, nruns);

        if ~isfile(fname_big)
            rms_by_g{gg}   = RMS_tmp;
            meanR_by_g(gg) = inf;
            conv_by_g(gg)  = NaN;
            continue;
        end

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

        rms_by_g{gg} = RMS_tmp;

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
    best_g(c)        = gstar;
    best_RMS_cols{c} = rms_by_g{idx_min};
    best_meanRMS(c)  = meanR_by_g(idx_min);
    best_convT(c)    = conv_by_g(idx_min);

    lam = lambda_list(gstar+1);
    fprintf('case=%d | best g=%d (λg=%.3g) | ConvTime=%s | mean(RMSE)=%.4g | pick=%s\n', ...
        si, gstar, lam, num2str(best_convT(c)), best_meanRMS(c), why);

    fname_best = sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                         sysnum, f_target, si, gstar, extra_suffix);
    [SNRYN_out(c), SNRUN_out(c)] = local_extract_snr(fname_best, sysnum, f_target, si);
end

Data   = NaN(nruns, num_cases);
labels = cell(1, num_cases);

for c = 1:num_cases
    Data(:, c) = best_RMS_cols{c}(:);
    labels{c} = sprintf('case=%d', SIs(c));
end

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

fig = figure('Color','w','Name','RMSE boxplot (g picked by ConvTime)');
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
ylab = ylabel('RMSE_{y}');
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

outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

suffix_clean = regexprep(extra_suffix, '^_+', '');
suffix_clean = regexprep(suffix_clean, '[^A-Za-z0-9\.\-]+', '_');

base = sprintf('Boxplot_RMSE_bestGbyConv_sys%02d_ft%02d_eps%g_%s', sysnum, f_target, eps_tol, suffix_clean);
base = strrep(base, '.', 'p');

savefig(fig, fullfile(outdir,[base,'.fig']));
exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');

fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

fprintf('==== Best-setting statistics (by case) ====\n');
for c = 1:num_cases
    v = Data(:, c); v = v(~isnan(v));
    if ~isempty(v) && isfinite(best_g(c))
        lam = lambda_list(best_g(c)+1);
        fprintf('case=%d | g=%d (λg=%.3g) | ConvTime=%s | N=%d | mean=%.4g | median=%.4g | min=%.4g | max=%.4g | SNRYN=%s | SNRUN=%s\n', ...
            SIs(c), best_g(c), lam, num2str(best_convT(c)), numel(v), mean(v), median(v), min(v), max(v), ...
            num2str(SNRYN_out(c)), num2str(SNRUN_out(c)));
    else
        fprintf('case=%d | no valid RMSE data | ConvTime=%s | SNRYN=%s | SNRUN=%s\n', ...
            SIs(c), num2str(best_convT(c)), num2str(SNRYN_out(c)), num2str(SNRUN_out(c)));
    end
end

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

function [snryn, snrun] = local_extract_snr(fname_best, sysnum, f_target, si)
    snryn = NaN; snrun = NaN;
    if ~isfile(fname_best), return; end

    info = whos('-file', fname_best);
    pat  = sprintf('^s%02d_r%02d_case%02d_run\\d+$', sysnum, f_target, si);
    names_all = {info.name};
    keep = ~cellfun(@isempty, regexp(names_all, pat, 'once'));
    names = names_all(keep);

    for k = 1:numel(names)
        S = load(fname_best, names{k});
        rec = S.(names{k});
        [a, b] = local_get_snr_pair(rec);
        if ~isnan(a) || ~isnan(b)
            snryn = a; snrun = b;
            return;
        end
    end

    Sfull = load(fname_best);
    [snryn, snrun] = local_get_snr_pair(Sfull);
end

function [snryn, snrun] = local_get_snr_pair(rec)
    snryn = NaN; snrun = NaN;

    cand_y = {'SNRYN','SNR_YN','snr_yn','snrYN','SNRyN','SNR_Y','SNR_Ynoise','SNR_Y'};
    cand_u = {'SNRUN','SNR_UN','snr_un','snrUN','SNRuN','SNR_U','SNR_Unoise','SNR_U'};

    try
        for k = 1:numel(cand_y)
            nm = cand_y{k};
            if isfield(rec, nm) && ~isempty(rec.(nm)) && isscalar(rec.(nm))
                snryn = double(rec.(nm));
                break;
            end
        end
        for k = 1:numel(cand_u)
            nm = cand_u{k};
            if isfield(rec, nm) && ~isempty(rec.(nm)) && isscalar(rec.(nm))
                snrun = double(rec.(nm));
                break;
            end
        end
    catch
    end
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
