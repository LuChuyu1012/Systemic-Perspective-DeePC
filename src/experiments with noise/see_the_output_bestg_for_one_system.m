clc; clear; close all;

sysnum       = 5;
f_target     = 2;
SIs          = 2:3;
nruns        = 20;
Fs           = 100;
eps_tol      = 0.1;
g_list       = 0:3;
extra_suffix = '_se1e-02';

lambda_list = [0.001 0.01 0.1 1];

y_fields = {'y_seq_opt','y_seq','y','y_out'};

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
eps_tag = strrep(sprintf('%.2g', eps_tol), '.', 'p');

num_cases = numel(SIs);
best_g    = NaN(1, num_cases);
best_mean = NaN(1, num_cases);
best_conv = NaN(1, num_cases);

fprintf('==== Search best g per case (min ConvTime, eps=%g; fallback min mean RMSE) ====\n', eps_tol);

for ci = 1:num_cases
    si = SIs(ci);

    mean_by_g = inf(1, numel(g_list));
    conv_by_g = NaN(1, numel(g_list));

    for gg = 1:numel(g_list)
        g = g_list(gg);

        fname_big = sprintf('sb%02d_r%02d_case%02d_g%d%s.mat', ...
                            sysnum, f_target, si, g, extra_suffix);

        if ~isfile(fname_big)
            continue;
        end

        RMS_tmp = NaN(nruns,1);
        Y_all   = cell(1, nruns);
        lens    = zeros(1, nruns);

        for runi = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, runi);
            try
                S = load(fname_big, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec, 'rms_error') && ~isempty(rec.rms_error)
                    RMS_tmp(runi) = rec.rms_error;
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
                Y_all{runi} = y;
                lens(runi)  = size(y,1);
            catch
            end
        end

        v = RMS_tmp(~isnan(RMS_tmp));
        if ~isempty(v)
            mean_by_g(gg) = mean(v);
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
            [~, k2] = min(mean_by_g(idxs));
            idx_min = idxs(k2);
        else
            idx_min = idxs(1);
        end

        best_g(ci)    = g_list(idx_min);
        best_mean(ci) = mean_by_g(idx_min);
        best_conv(ci) = conv_by_g(idx_min);

        lam = lambda_list(best_g(ci)+1);
        fprintf('case=%d | best g=%d (λg=%.3g) | ConvTime=%s | mean(RMSE)=%.4g\n', ...
            si, best_g(ci), lam, num2str(best_conv(ci)), best_mean(ci));
    else
        [mval, idx_min] = min(mean_by_g);
        if isfinite(mval)
            best_g(ci)    = g_list(idx_min);
            best_mean(ci) = mval;
            best_conv(ci) = NaN;

            lam = lambda_list(best_g(ci)+1);
            fprintf('case=%d | ConvTime all NaN -> fallback | best g=%d (λg=%.3g) | mean(RMSE)=%.4g\n', ...
                si, best_g(ci), lam, best_mean(ci));
        else
            fprintf('case=%d | no valid data\n', si);
        end
    end
end

for ci = 1:num_cases
    si = SIs(ci);
    g  = best_g(ci);

    if isnan(g)
        warning('Skip: case=%d has no best g.', si);
        continue;
    end

    fname_best = sprintf('sb%02d_r%02d_case%02d_g%d%s.mat', ...
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

    % ======= CHANGED: samples axis (NOT seconds) =======
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

    % ======= CHANGED: xlabel text =======
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

    base = sprintf('TrajOverlay_sys%02d_ft%02d_case%02d_g%d_lam%.3g_eps%s_%s', ...
                   sysnum, f_target, si, g, lam, eps_tag, suffix_clean);
    base = strrep(base, '.', 'p');

    ws = warning; warning('off','all');
    cleanupObj = onCleanup(@() warning(ws));

    savefig(fig, fullfile(outdir,[base,'.fig']));
    exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');

    fprintf('Saved FIG: %s\n', fullfile(outdir,[base,'.fig']));
    fprintf('Saved PNG: %s\n', fullfile(outdir,[base,'.png']));
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
