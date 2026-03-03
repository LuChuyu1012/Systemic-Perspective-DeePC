clc; clear; close all;

sysnum   = 5;
f_target = 2;
SIs      = 0:3;
nruns    = 50;
Fs       = 100;

conf = 0.95;

suffix_list = {'_sigmae5e-02','_sigmae1e-02','_snr30dB'};

g_list = [0 1 2 3 5];

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

outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

data_dir = fullfile(pwd, 'data');

for mm = 1:numel(suffix_list)

    extra_suffix = suffix_list{mm};

    for kcase = 1:numel(use_cases)
        si = use_cases(kcase);

        rmse_mean_by_g = NaN(1, numel(g_list));
        n_used_by_g    = zeros(1, numel(g_list));

        for gg = 1:numel(g_list)
            g = g_list(gg);

            fname_big = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
                sysnum, f_target, si, g, extra_suffix));

            if ~isfile(fname_big)
                continue;
            end

            rmse_vec = NaN(1, nruns);

            for i = 1:nruns
                varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
                try
                    S = load(fname_big, varname);
                    if ~isfield(S, varname), continue; end
                    rec = S.(varname);

                    if ~isfield(rec, 'y_seq_opt') || isempty(rec.y_seq_opt)
                        continue;
                    end

                    y = rec.y_seq_opt(:);
                    N = numel(y);
                    if N <= 0, continue; end

                    t = (0:N-1).';
                    r = sin(2*pi*f_target*(t/Fs));

                    diff = y - r;
                    rmse_vec(i) = sqrt(mean(diff.^2, 'omitnan'));
                catch
                end
            end

            rmse_vec = rmse_vec(isfinite(rmse_vec));
            n_used_by_g(gg) = numel(rmse_vec);

            if ~isempty(rmse_vec)
                rmse_mean_by_g(gg) = mean(rmse_vec);
            end
        end

        finite_g = isfinite(rmse_mean_by_g);
        if ~any(finite_g)
            fprintf('Skip (no valid g found): case=%d | %s\n', si, extra_suffix);
            continue;
        end

        cand = rmse_mean_by_g;
        cand(~finite_g) = inf;

        minMean = min(cand);
        idxs = find(cand == minMean);

        idx_best = idxs(1);
        why = 'min mean RMSE';

        g_fixed = g_list(idx_best);

        fprintf('case=%d | best g=%d | meanRMSE=%g | used=%d/%d | pick=%s | %s\n', ...
            si, g_fixed, rmse_mean_by_g(idx_best), n_used_by_g(idx_best), nruns, why, extra_suffix);

        fname_big = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d_g%d%s.mat', ...
            sysnum, f_target, si, g_fixed, extra_suffix));

        if ~isfile(fname_big)
            fprintf('Skip (file not found after selection): %s\n', fname_big);
            continue;
        end

        Ys   = cell(1, nruns);
        Us   = cell(1, nruns);
        ly   = NaN(1, nruns);
        lu   = NaN(1, nruns);

        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname_big, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec, 'y_seq_opt') && ~isempty(rec.y_seq_opt)
                    y = rec.y_seq_opt(:);
                    Ys{i} = y;
                    ly(i) = numel(y);
                end

                if isfield(rec, 'u_seq_opt') && ~isempty(rec.u_seq_opt)
                    u = rec.u_seq_opt(:);
                    Us{i} = u;
                    lu(i) = numel(u);
                end
            catch
            end
        end

        goodY = find(~cellfun(@isempty, Ys) & isfinite(ly) & ly > 0);
        if isempty(goodY)
            fprintf('Skip (no valid y_seq_opt): case=%d | %s\n', si, extra_suffix);
            continue;
        end
        Ny    = min(ly(goodY));
        nYeff = numel(goodY);

        Ymat = NaN(Ny, nYeff);
        for j = 1:nYeff
            y = Ys{goodY(j)};
            Ymat(:,j) = y(1:Ny);
        end

        [Data_y_mean, Data_y_var, CI_y] = Mean_Variance_and_CI(Ymat, nYeff, conf, 2);

        goodU = find(~cellfun(@isempty, Us) & isfinite(lu) & lu > 0);
        hasU  = ~isempty(goodU);

        if hasU
            Nu    = min(lu(goodU));
            nUeff = numel(goodU);

            Umat = NaN(Nu, nUeff);
            for j = 1:nUeff
                u = Us{goodU(j)};
                Umat(:,j) = u(1:Nu);
            end

            [Data_u_mean, Data_u_var, CI_u] = Mean_Variance_and_CI(Umat, nUeff, conf, 2);
        end

        t_y = (0:Ny-1).';
        if hasU
            t_u = (0:Nu-1).';
        end

        r_y = sin(2*pi*f_target*(t_y/Fs));

        fprintf('case=%d | y used %d/%d runs | Ny=%d | g=%d | %s\n', si, nYeff, nruns, Ny, g_fixed, extra_suffix);
        if hasU
            fprintf('case=%d | u used %d/%d runs | Nu=%d | g=%d | %s\n', si, nUeff, nruns, Nu, g_fixed, extra_suffix);
        else
            fprintf('case=%d | u not found (no u_seq_opt in file) | %s\n', si, extra_suffix);
        end

        figy = figure('Color','w','Visible','on', ...
            'Name', sprintf('Output mean and CI95 case=%d %s', si, extra_suffix));
        figy.Units = plt_set.plot_unit;
        figy.Position(3) = plt_set.plot_dim_x;
        figy.Position(4) = plt_set.plot_dim_y;
        figy.Renderer = 'opengl';

        axy = axes(figy); hold(axy,'on');

        col = use_rgb(kcase,:);
        col_fill = to_white(col, 0.75);
        y_lo = Data_y_mean - CI_y;
        y_hi = Data_y_mean + CI_y;

        patch(axy, [t_y; flipud(t_y)], [y_lo; flipud(y_hi)], col_fill, ...
            'EdgeColor', 'none', 'FaceAlpha', 1.0);

        plot(axy, t_y, Data_y_mean, 'LineWidth', plt_set.thick, ...
            'Color', col);

        plot(axy, t_y, r_y, 'k', 'LineWidth', plt_set.thick);

        axy.FontName  = plt_set.fontname;
        axy.FontSize  = plt_set.fontsize;
        axy.LineWidth = 1.0;

        xlab = xlabel(axy, 't [samples]');
        ylab = ylabel(axy, 'y(t)');
        set(xlab, 'FontName', plt_set.fontname, 'FontSize', plt_set.x_font_dim);
        set(ylab, 'FontName', plt_set.fontname, 'FontSize', plt_set.y_font_dim);

        grid(axy,'on');
        axy.XGrid = 'off';
        axy.YGrid = 'on';
        axy.GridAlpha = 0.15;

        xlim(axy, [t_y(1) t_y(end)]);
        ylim(axy, [-1 2]);

        set(axy,'LooseInset', max(get(axy,'TightInset'), 0.02));
        figy.PaperPositionMode = 'auto';
        drawnow;

%         basey = sprintf('MeanTraj_CI95_y_sys%02d_ft%02d_case%02d_g%d_%s', ...
%             sysnum, f_target, si, g_fixed, regexprep(extra_suffix,'^_',''));
%         savefig(figy, fullfile(outdir, [basey, '.fig']));
%         exportgraphics(figy, fullfile(outdir, [basey, '.png']), ...
%             'Resolution', 300, 'BackgroundColor', 'white');
% 
%         fprintf('Y FIG saved to: %s\n', fullfile(outdir,[basey,'.fig']));
%         fprintf('Y PNG saved to: %s\n', fullfile(outdir,[basey,'.png']));

    end
end

function [Data_mean, Data_var, Data_CI] = Mean_Variance_and_CI(data, n_data, conf, direction)
    alpha = 1 - conf;
    pUp   = 1 - alpha/2;

    Data_mean = (1/n_data) * sum(data, direction);
    Data_var  = var(data, 0, direction);

    Data_SEM = ((1/n_data) * Data_var).^(1/2);

    Data_CI = Data_SEM * tinv(pUp, n_data - 1);
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

function Cw = to_white(C, a)
    Cw = (1-a)*C + a*[1 1 1];
end