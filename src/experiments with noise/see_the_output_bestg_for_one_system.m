clc; clear; close all;

sysnum       = 5;
f_target     = 2;
SIs          = 2:3;
nruns        = 20;
Fs           = 100;
g_list       = 0:3;

SNRdB_target = 20;
extra_suffix = sprintf('_snr%02ddB', round(SNRdB_target));

lambda_list = [0.001 0.01 0.1 1];

plt_set.y_font_dim     = 16;
plt_set.x_font_dim     = 16;
plt_set.thick          = 1.0;
plt_set.ref_thick      = 0.9;
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

outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

suffix_clean = regexprep(extra_suffix, '^_+', '');
suffix_clean = regexprep(suffix_clean, '[^A-Za-z0-9\.\-]+', '_');

num_cases = numel(SIs);
best_g    = NaN(1, num_cases);
best_bar  = NaN(1, num_cases);

fprintf('==== Pick best g per case (ONLY min RMSEbar(mean(y_seq_opt) vs ref)) ====\n');

% =====================================================================
% 1) For EACH case: pick its own best g by RMSEbar
% =====================================================================
for ci = 1:num_cases
    si = SIs(ci);

    bar_by_g = NaN(1, numel(g_list));

    for gg = 1:numel(g_list)
        g = g_list(gg);

        fname_big = sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                            sysnum, f_target, si, g, extra_suffix);

        if ~isfile(fname_big)
            continue;
        end

        Ysopt   = cell(1, nruns);
        lensOpt = NaN(1, nruns);

        for runi = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, runi);
            try
                S = load(fname_big, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec, 'y_seq_opt') && ~isempty(rec.y_seq_opt)
                    yy = rec.y_seq_opt(:);
                    Ysopt{runi}   = yy;
                    lensOpt(runi) = numel(yy);
                end
            catch
            end
        end

        bar_by_g(gg) = local_rmsebar_meantraj_vs_ref(Ysopt, lensOpt, Fs, f_target);
    end

    finite_bar = isfinite(bar_by_g);
    if ~any(finite_bar)
        fprintf('case=%d | no valid RMSEbar -> skip saving\n', si);
        continue;
    end

    cand = bar_by_g;
    cand(~finite_bar) = inf;

    [minBar, idx_min] = min(cand);

    best_g(ci)   = g_list(idx_min);
    best_bar(ci) = minBar;

    lam = lambda_list(best_g(ci)+1);
    fprintf('case=%d | best g=%d (λg=%.3g) | RMSEbar=%s\n', ...
        si, best_g(ci), lam, num2str(best_bar(ci)));
end


for ci = 1:num_cases
    si = SIs(ci);
    g  = best_g(ci);

    if ~isfinite(g)
        warning('Skip: case=%d has no best g.', si);
        continue;
    end

    fname_best = sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                         sysnum, f_target, si, g, extra_suffix);
    if ~isfile(fname_best)
        warning('File not found: %s (skip case=%d)', fname_best, si);
        continue;
    end

    allY = {};
    allY_len = [];

    for runi = 1:nruns
        varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, runi);
        try
            S = load(fname_best, varname);
            if ~isfield(S, varname), continue; end
            rec = S.(varname);

            if ~isfield(rec, 'y_seq_opt') || isempty(rec.y_seq_opt), continue; end
            y = rec.y_seq_opt(:);

            allY{end+1} = y;
            allY_len(end+1) = numel(y);
        catch
        end
    end

    if isempty(allY_len)
        warning('case=%d | best g=%d has no y_seq_opt.', si, g);
        continue;
    end

    N = min(allY_len);

    t = (0:N-1).';
    r = sin(2*pi*f_target*(t/Fs));

    lam = lambda_list(g+1);

    fig = figure('Color','w','Name',sprintf('sys%02d case%02d',sysnum,si));
    fig.Units = plt_set.plot_unit;
    fig.Position(3) = plt_set.plot_dim_x;
    fig.Position(4) = plt_set.plot_dim_y;
    fig.Renderer = 'opengl';

    ax = gca; hold(ax,'on');

    plot(t, r, 'k-', 'LineWidth', plt_set.ref_thick);

    C = use_rgb(ci,:);
    first = true;
    for j = 1:numel(allY)
        y = allY{j};
        y = y(1:N);

        if first
            plot(t, y, '-', 'Color', C, 'LineWidth', plt_set.thick);
            first = false;
        else
            plot(t, y, '-', 'Color', C, 'LineWidth', plt_set.thick, 'HandleVisibility','off');
        end
    end

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

    % ======= SAVE: ONLY this best-g figure for this case =======
    base = sprintf('TrajOverlay_sys%02d_ft%02d_case%02d_bestg%d_lam%.3g_%s', ...
                   sysnum, f_target, si, g, lam, suffix_clean);
    base = strrep(base, '.', 'p');

    savefig(fig, fullfile(outdir,[base,'.fig']));
    exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');

    fprintf('Saved (case=%d best g=%d): %s\n', si, g, base);
end

% ========================= FUNCTIONS =========================
function rmse_bar = local_rmsebar_meantraj_vs_ref(Ysopt, lensOpt, Fs, f_target)
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
