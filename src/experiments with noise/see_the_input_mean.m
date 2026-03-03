clc; clear; close all;

sysnum   = 5;
f_target = 2;
SIs      = 0:3;
nruns    = 50;

Fs = 100;

g_list = [0 1 2 3 5];

suffix_list = {'_sigmae5e-02','_sigmae1e-02','_snr30dB'};

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
    best_g       = NaN(1, num_cases);

    fprintf('==== Best g selection per case: minimize mean(rms_error) across %d runs. ====\n', nruns);

    for c = 1:num_cases
        si = SIs(c);

        mean_rms_by_g = NaN(1, numel(g_list));
        Nvalid_by_g   = zeros(1, numel(g_list));

        for gg = 1:numel(g_list)
            g = g_list(gg);

            fname_sel = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                                sysnum, f_target, si, g, extra_suffix));
            if ~isfile(fname_sel)
                continue;
            end

            rms_runs = NaN(1, nruns);

            for i = 1:nruns
                varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
                try
                    S = load(fname_sel, varname);
                    if ~isfield(S, varname), continue; end
                    rec = S.(varname);

                    if isfield(rec,'rms_error') && ~isempty(rec.rms_error)
                        rms_runs(i) = mean(rec.rms_error(:));
                    end
                catch
                end
            end

            v = rms_runs(isfinite(rms_runs));
            Nvalid_by_g(gg) = numel(v);

            if ~isempty(v)
                mean_rms_by_g(gg) = mean(v);
            end

            if Nvalid_by_g(gg) ~= nruns
                fprintf('case=%d g=%d | valid rms_error=%d/%d | %s\n', si, g, Nvalid_by_g(gg), nruns, extra_suffix);
            end
        end

        finite_mean = isfinite(mean_rms_by_g);
        if ~any(finite_mean)
            fprintf('case=%d | no valid mean(rms_error) across g (skip) | %s\n', si, extra_suffix);
            continue;
        end

        cand = mean_rms_by_g;
        cand(~finite_mean) = inf;

        minMean = min(cand);
        idxs = find(cand == minMean);
        idx_min = idxs(1);

        gstar = g_list(idx_min);
        best_g(c) = gstar;

        fprintf('case=%d | best g=%d | mean(rms_error)=%g | valid=%d/%d | %s\n', ...
            si, gstar, mean_rms_by_g(idx_min), Nvalid_by_g(idx_min), nruns, extra_suffix);
    end

    Umean_all = cell(1, num_cases);
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

        Us   = cell(1, nruns);
        lens = NaN(1, nruns);

        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname_best, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec,'u_seq_opt') && ~isempty(rec.u_seq_opt)
                    u = rec.u_seq_opt(:);
                    Us{i}   = u;
                    lens(i) = numel(u);
                end
            catch
            end
        end

        good = find(~cellfun(@isempty, Us) & isfinite(lens) & lens>0);
        if isempty(good)
            warning('case=%d | g=%d: no u_seq_opt found, skip.', si, g);
            continue;
        end

        N = min(lens(good));
        Umat = NaN(N, numel(good));
        for j = 1:numel(good)
            u = Us{good(j)};
            Umat(:,j) = u(1:N);
        end

        Umean_all{c} = mean(Umat, 2, 'omitnan');
        N_all(c)     = N;
        n_used(c)    = numel(good);

        fprintf('case=%d | best g=%d | used %d/%d runs | N=%d | %s\n', ...
            si, g, n_used(c), nruns, N, extra_suffix);
    end

    valid = find(~cellfun(@isempty, Umean_all) & isfinite(N_all));
    if isempty(valid)
        error('No case produced mean u_seq_opt (check files/u_seq_opt field).');
    end

    Nmin = min(N_all(valid));
    t = (0:Nmin-1).';

    fig = figure('Color','w','Name','Mean input (all cases)');
    fig.Units = plt_set.plot_unit;
    fig.Position(3) = plt_set.plot_dim_x;
    fig.Position(4) = plt_set.plot_dim_y;
    fig.Renderer = 'opengl';

    ax = gca; hold(ax,'on');

    for c = 1:num_cases
        if isempty(Umean_all{c}), continue; end
        u_plot = Umean_all{c}(1:Nmin);
        plot(t, u_plot, '-', 'Color', use_rgb(c,:), 'LineWidth', plt_set.thick);
    end

    ax.FontName  = plt_set.fontname;
    ax.FontSize  = plt_set.fontsize;
    ax.LineWidth = 1.0;

    xlab = xlabel('t [samples]');
    ylab = ylabel(ax, '$\bar{u}(t)$', 'Interpreter','latex');
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
%     outdir = fullfile(pwd, 'figures');
%     if ~exist(outdir,'dir'); mkdir(outdir); end
% 
%     suffix_clean = regexprep(extra_suffix, '^_+', '');
%     suffix_clean = regexprep(suffix_clean, '[^A-Za-z0-9\.\-]+', '_');
% 
%     base = sprintf('MeanU_seqopt_AllCases_sys%02d_ft%02d_bestGbyMeanRMSErun_%s', ...
%                    sysnum, f_target, suffix_clean);
% 
%     savefig(fig, fullfile(outdir, [base, '.fig']));
%     exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');
% 
%     fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
%     fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));
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