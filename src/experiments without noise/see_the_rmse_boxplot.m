clc; clear; close all;

sysnum   = 6;
f_target = 2;
SIs      = 0:3;
nruns    = 20;

RMS     = NaN(numel(SIs), nruns);
COND    = NaN(numel(SIs), nruns);
COVRMS  = NaN(numel(SIs), nruns);

for sidx = 1:numel(SIs)
    si = SIs(sidx);
    fname_big = sprintf('s%02d_r%02d_case%02d.mat', sysnum, f_target, si);

    if ~isfile(fname_big)
        warning('Big file not found: %s (skip case=%d)', fname_big, si);
        continue;
    end

    for i = 1:nruns
        varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
        try
            S = load(fname_big, varname);
            rec = S.(varname);

            if isfield(rec, 'rms_error')
                RMS(sidx, i) = rec.rms_error;
            end
            if isfield(rec, 'cond_num')
                COND(sidx, i) = rec.cond_num;
            end
            if isfield(rec, 'coverage') && ~isempty(rec.coverage)
                v = rec.coverage(:);
                COVRMS(sidx, i) = sqrt(mean(v.^2));
            end
        catch ME
            warning('Read failed: %s | %s', varname, ME.message);
        end
    end
end

plt_set.thick          = 1.6;
plt_set.plot_unit      = 'centimeters';
plt_set.fontname       = 'Times';

plt_set.fontsize       = 10.5;
plt_set.x_font_dim     = 11;
plt_set.y_font_dim     = 16;

plt_set.plot_dim_x     = 8;
plt_set.plot_dim_y     = 6;

case_hex = {'#007191', '#62c8d3', '#f47a00', '#d31f11'};
ncase = numel(SIs);
use_rgb = zeros(ncase,3);
for k = 1:ncase
    use_rgb(k,:) = hex2rgb(case_hex{k});
end

labels = arrayfun(@(x) sprintf('case=%d', x), SIs, 'UniformOutput', false);

fig = figure('Color','w','Name','RMSE (boxplot) by case');
fig.Units = plt_set.plot_unit;
fig.Position(3) = plt_set.plot_dim_x;
fig.Position(4) = plt_set.plot_dim_y;

fig.Renderer = 'opengl';

boxplot(RMS.', 'Labels', labels, 'Symbol','k+');
ax = gca;

ax.FontName  = plt_set.fontname;
ax.FontSize  = plt_set.fontsize;
ax.LineWidth = 1.0;

xlabel('');

% ======= Y label (UPDATED) =======
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

set(hBox, 'LineWidth', thick);
set(hMed, 'LineWidth', thick);
set(hWhi, 'LineWidth', thick);
set(hOut, 'LineWidth', 1.2);
if ~isempty(hCapU), set(hCapU,'LineWidth', thick); end
if ~isempty(hCapL), set(hCapL,'LineWidth', thick); end

for k = 1:ncase
    kk = ncase - k + 1;

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

% ======= X tick labels (UPDATED) =======
set(ax, 'TickLabelInterpreter', 'latex');
xticklabels({'$\mathbf{WN}$', '$\mathbf{IBW}$', '$\mathbf{IBN}$', '$\mathbf{OB}$'});

set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));
fig.PaperPositionMode = 'auto';
drawnow;

outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

base = sprintf('Boxplot_RMSE_sys%02d_ft%02d_sidebyside', sysnum, f_target);

set(fig, 'Visible', 'on');
savefig(fig, fullfile(outdir, [base, '.fig']));

exportgraphics(fig, fullfile(outdir, [base, '.png']), 'Resolution', 300, 'BackgroundColor', 'white');

fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

fprintf('==== RMS error ====\n');
for sidx = 1:numel(SIs)
    v = RMS(sidx, :); v = v(~isnan(v));
    if ~isempty(v)
        fprintf('case=%d | N=%d | mean=%.4f | median=%.4f | min=%.4f | max=%.4f\n', ...
            SIs(sidx), numel(v), mean(v), median(v), min(v), max(v));
    else
        fprintf('case=%d | no valid RMS data\n', SIs(sidx));
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
