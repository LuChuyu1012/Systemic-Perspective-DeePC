clc; clear; close all;

% ===================== Settings =====================
sys_list = 1:5;
f_target = 2;
SIs      = 0:3;
nruns    = 50;

Fs = 100;

y_fields = {'y_seq_opt','y_seq','y','y_out'};
data_dir = fullfile(pwd, 'data');

nsys  = numel(sys_list);
ncase = numel(SIs);

% ===================== Plot style =====================
plt_set.plot_unit      = 'centimeters';
plt_set.fontname       = 'Times';

plt_set.fontsize       = 8;
plt_set.x_font_dim     = 8;
plt_set.y_font_dim     = 10;

plt_set.plot_dim_x     = 10;
plt_set.plot_dim_y     = 6;

case_hex = {'#007191', '#62c8d3', '#f47a00', '#d31f11'};
use_rgb = zeros(ncase,3);
for k = 1:ncase
    use_rgb(k,:) = hex2rgb(case_hex{k});
end

gap_sys   = 0.55;
width_box = 0.32;

RMS = NaN(nsys, ncase, nruns);

fprintf('==== Load noise-free RMS error per (system, case) ====\n');

% ===================== Main loops =====================
for ss = 1:nsys
    sysnum = sys_list(ss);

    for cc = 1:ncase
        si = SIs(cc);

        fname = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d.mat', ...
                            sysnum, f_target, si));
        if ~isfile(fname)
            fprintf('System %d case=%d | file not found -> skip\n', sysnum, si);
            continue;
        end

        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec, 'rms_error') && ~isempty(rec.rms_error)
                    RMS(ss, cc, i) = mean(rec.rms_error(:));
                end
            catch
            end
        end
    end
end

% ===================== build one long vector + positions grouped by system =====================
X = [];
pos = [];

for ss = 1:nsys
    base = (ss-1)*(ncase + gap_sys);

    for cc = 1:ncase
        p = base + cc;

        v = squeeze(RMS(ss, cc, :));
        v = v(~isnan(v));
        if isempty(v), continue; end

        X   = [X; v];
        pos = [pos; repmat(p, numel(v), 1)];
    end
end

pos_unique = unique(pos);
pos_unique = sort(pos_unique, 'ascend');

fig = figure('Color','w','Name','RMSE boxplots grouped by system');
fig.Units = plt_set.plot_unit;
fig.Position(3) = plt_set.plot_dim_x;
fig.Position(4) = plt_set.plot_dim_y;
fig.Renderer = 'opengl';

boxplot(X, pos, 'Positions', pos_unique, 'Widths', width_box, 'Symbol','k+');

ax = gca;
ax.FontName  = plt_set.fontname;
ax.FontSize  = plt_set.fontsize;
ax.LineWidth = 1.0;

xlabel('');

ylab = ylabel('$\mathrm{RMSE}_y^i$', 'Interpreter', 'latex');
set(ylab, 'FontName', plt_set.fontname, 'FontSize', plt_set.y_font_dim);

grid on;
ax.XGrid = 'off';
ax.YGrid = 'on';
ax.GridAlpha = 0.15;

x_min = min(pos_unique) - 0.75;
x_max = max(pos_unique) + 0.75;
xlim(ax, [x_min, x_max]);

line_w_box = 1.0;
line_w_mid = 1.0;
line_w_whi = 0.9;
line_w_out = 0.9;

get_case_from_pos = @(p) round(p - (ncase+gap_sys)*floor((p-1)/(ncase+gap_sys)));

hLines = findobj(ax, 'Type', 'Line');

for h = reshape(hLines,1,[])
    xd = get(h, 'XData');
    if isempty(xd) || all(isnan(xd)), continue; end

    xmid = mean(xd(~isnan(xd)));
    [~, idxp] = min(abs(pos_unique - xmid));
    p = pos_unique(idxp);

    cc = get_case_from_pos(p);
    cc = max(1, min(ncase, cc));
    col = use_rgb(cc,:);

    tg = get(h, 'Tag');
    switch tg
        case 'Box'
            set(h, 'Color', col, 'LineWidth', line_w_box);
        case 'Median'
            set(h, 'Color', col, 'LineWidth', line_w_mid);
        case 'Whisker'
            set(h, 'Color', col, 'LineWidth', line_w_whi);
        case 'Upper Adjacent Value'
            set(h, 'Color', col, 'LineWidth', line_w_whi);
        case 'Lower Adjacent Value'
            set(h, 'Color', col, 'LineWidth', line_w_whi);
        case 'Outliers'
            set(h, 'Color', col, 'LineWidth', line_w_out, ...
                   'MarkerEdgeColor', col, 'MarkerFaceColor', col);
        otherwise
    end
end

% ===================== X ticks: only System labels at group centers (bottom) =====================
centers = NaN(nsys,1);
for ss = 1:nsys
    base = (ss-1)*(ncase + gap_sys);
    centers(ss) = base + (ncase+1)/2;
end

set(ax, 'XTick', centers);
set(ax, 'TickLabelInterpreter', 'latex');
xticklabels(arrayfun(@(k) sprintf('$\\mathbf{System\\ %d}$', k), sys_list, 'UniformOutput', false));
ax.XTickLabelRotation = 0;

% ===================== legend on top =====================
case_names = {'WN','IBW','IBN','OB'};
hLeg = gobjects(1,ncase);
for c = 1:ncase
    hLeg(c) = line(ax, [NaN NaN], [NaN NaN], ...
        'Color', use_rgb(c,:), 'LineWidth', 6);
end
lgd = legend(hLeg, case_names, 'Location','northoutside');
set(lgd, 'Orientation','horizontal', 'Box','off');
lgd.FontName = plt_set.fontname;
lgd.FontSize = plt_set.fontsize;

set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));
fig.PaperPositionMode = 'auto';
drawnow;

% ===================== Save =====================
outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

    base = sprintf('Boxplot_RMSE_allSys_ft%02d', f_target);

    savefig(fig, fullfile(outdir, [base, '.fig']));
    exportgraphics(fig, fullfile(outdir, [base, '.png']), 'Resolution', 300, 'BackgroundColor', 'white');

    fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
    fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

% ===================== Print stats =====================
fprintf('==== RMS error stats (by system, case) ====\n');
for ss = 1:nsys
    for cc = 1:ncase
        v = squeeze(RMS(ss, cc, :)); 
        v = v(~isnan(v));
        if ~isempty(v)
            fprintf('System %d case=%d | N=%d | mean=%.4f | min=%.4f | max=%.4f\n', ...
                sys_list(ss), SIs(cc), numel(v), mean(v), min(v), max(v));
        else
            fprintf('System %d case=%d | no valid RMS data\n', sys_list(ss), SIs(cc));
        end
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