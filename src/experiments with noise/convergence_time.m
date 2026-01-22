clc; clear; close all;

sys_list  = 1:6;
f_target  = 2;
nruns     = 20;

Fs        = 100;
eps_pick  = 0.1;
eps_eval  = 0.3;

g_list       = [0 1 2 3];
lambda_list  = [0.001 0.01 0.1 1];
extra_suffix = '_se5e-03';

SIs_all  = 0:3;
SIs_plot = 0:2;

y_fields = {'y_seq_opt','y_seq','y','y_out'};

nsys      = numel(sys_list);
ncase_all = numel(SIs_all);

BestG     = NaN(nsys, ncase_all);
BestLam   = NaN(nsys, ncase_all);
BestMeanR = NaN(nsys, ncase_all);
ConvTime  = NaN(nsys, ncase_all);

fprintf('==== Pick g by ConvTime(eps=%.3g); then evaluate ConvTime(eps=%.3g) ====\n', eps_pick, eps_eval);

for ss = 1:nsys
    sysnum = sys_list(ss);

    for cc = 1:ncase_all
        si = SIs_all(cc);

        mean_by_g = inf(1, numel(g_list));
        conv_by_g = NaN(1, numel(g_list));

        for gg = 1:numel(g_list)
            g = g_list(gg);

            fname = sprintf('sb%02d_r%02d_case%02d_g%d%s.mat', sysnum, f_target, si, g, extra_suffix);
            if ~isfile(fname), continue; end

            RMS   = NaN(nruns,1);
            Y_all = cell(1, nruns);
            lens  = zeros(1, nruns);

            for i = 1:nruns
                varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
                try
                    S = load(fname, varname);
                    if ~isfield(S, varname), continue; end
                    rec = S.(varname);

                    if isfield(rec,'rms_error') && ~isempty(rec.rms_error)
                        RMS(i) = rec.rms_error;
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

            v = RMS(~isnan(RMS));
            if ~isempty(v)
                mean_by_g(gg) = mean(v);
            end

            conv_by_g(gg) = local_convtime_maxoverruns(Y_all, lens, Fs, f_target, eps_pick);
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
        else
            [mval, idx_min] = min(mean_by_g);
            if ~isfinite(mval)
                fprintf('sys=%d case=%d | no valid data, skip\n', sysnum, si);
                continue;
            end
        end

        gstar = g_list(idx_min);
        BestG(ss,cc)     = gstar;
        BestLam(ss,cc)   = lambda_list(gstar+1);
        BestMeanR(ss,cc) = mean_by_g(idx_min);

        fname_best = sprintf('sb%02d_r%02d_case%02d_g%d%s.mat', sysnum, f_target, si, gstar, extra_suffix);
        if isfile(fname_best)
            Y_all = cell(1, nruns);
            lens  = zeros(1, nruns);

            for i = 1:nruns
                varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
                try
                    S = load(fname_best, varname);
                    if ~isfield(S, varname), continue; end
                    rec = S.(varname);

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

            ConvTime(ss,cc) = local_convtime_maxoverruns(Y_all, lens, Fs, f_target, eps_eval);
        end

        fprintf('sys=%d case=%d | best g=%d (λg=%.3g) | meanRMSE=%.4g | ConvPick(ε=%.1f)=%s | ConvFinal(ε=%.1f)=%s\n', ...
            sysnum, si, gstar, BestLam(ss,cc), BestMeanR(ss,cc), eps_pick, num2str(conv_by_g(idx_min)), eps_eval, num2str(ConvTime(ss,cc)));
    end
end

idx_plot = ismember(SIs_all, SIs_plot);
ConvPlot = ConvTime(:, idx_plot);

labels_all = arrayfun(@(x) sprintf('case=%d', x), SIs_all, 'UniformOutput', false);

plt_set.y_font_dim     = 11;
plt_set.x_font_dim     = 10;
plt_set.thick          = 1.6;
plt_set.plot_unit      = 'centimeters';
plt_set.fontname       = 'Times';
plt_set.fontsize       = 10;
plt_set.plot_dim_x     = 9;
plt_set.plot_dim_y     = 6;

case_hex = ["#007191"; "#62c8d3"; "#f47a00"; "#d31f11"];
use_rgb = zeros(4,3);
for j = 1:4
    use_rgb(j,:) = hex2rgb(case_hex(j));
end

fig = figure('Color','w','Name',sprintf('ConvTime pick=%.1f eval=%.1f',eps_pick,eps_eval));
fig.Units = plt_set.plot_unit;
fig.Position(3) = plt_set.plot_dim_x;
fig.Position(4) = plt_set.plot_dim_y;

fig.Renderer = 'opengl';

b = bar(ConvPlot, 'grouped'); hold on;

for k = 1:numel(b)
    b(k).FaceColor = use_rgb(k,:);
    b(k).EdgeColor = 'none';
end

ax = gca;
ax.FontName  = plt_set.fontname;
ax.FontSize  = plt_set.fontsize;
ax.LineWidth = 1.0;
ax.XTick = 1:nsys;
ax.XTickLabel = arrayfun(@(s) sprintf('Sys=%d', s), sys_list, 'UniformOutput', false);

xlabel('');
ylab = ylabel('c [s]');
set(ylab, 'FontName', plt_set.fontname, 'FontSize', plt_set.y_font_dim);

grid on;
ax.XGrid = 'off';
ax.YGrid = 'on';
ax.GridAlpha = 0.15;

hleg = gobjects(4,1);
for j = 1:4
    hleg(j) = line(NaN, NaN, ...
        'LineStyle','none', ...
        'Marker','s', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', use_rgb(j,:), ...
        'MarkerEdgeColor', use_rgb(j,:));
end

lgd = legend(hleg, labels_all, 'Location','northoutside', 'Orientation','horizontal', 'NumColumns', 2);
set(lgd, 'FontName', plt_set.fontname, 'FontSize', 10, 'Box','off');

set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));
fig.PaperPositionMode = 'auto';
drawnow;
hold off;

outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

suffix_clean = regexprep(extra_suffix, '^_+', '');
suffix_clean = regexprep(suffix_clean, '[^A-Za-z0-9\.\-]+', '_');

base = sprintf('ConvTime_Pick%g_Eval%g_ft%02d_%s', eps_pick, eps_eval, f_target, suffix_clean);
base = strrep(base, '.', 'p');

ws = warning; warning('off','all');
cleanupObj = onCleanup(@() warning(ws));

drawnow;

savefig(fig, fullfile(outdir,[base,'.fig']));
exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');

fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

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
    if h(1) == '#', h = h(2:end); end
    if numel(h) ~= 6
        error('hex2rgb:InvalidHex', 'Hex color must be 6 characters, got: %s', h);
    end
    r = hex2dec(h(1:2));
    g = hex2dec(h(3:4));
    b = hex2dec(h(5:6));
    rgb = [r g b] / 255;
end
