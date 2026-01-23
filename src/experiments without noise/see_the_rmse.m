clc; clear; close all;

sys_list = 1:6;
f_target = 2;
SIs      = 0:3;
nruns    = 20;

medRMSE = NaN(numel(sys_list), numel(SIs));

for k = 1:numel(sys_list)
    sysnum = sys_list(k);

    for sidx = 1:numel(SIs)
        si = SIs(sidx);
        fname_big = sprintf('s%02d_r%02d_case%02d.mat', sysnum, f_target, si);

        if ~isfile(fname_big)
            warning('File not found: %s (sys=%d case=%d)', fname_big, sysnum, si);
            continue;
        end

        v = NaN(1, nruns);
        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname_big, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                if isfield(rec,'rms_error') && ~isempty(rec.rms_error)
                    v(i) = rec.rms_error;
                end
            catch
            end
        end

        v = v(~isnan(v));
        if ~isempty(v)
            medRMSE(k, sidx) = median(v);
        else
            warning('No valid rms_error: sys=%d case=%d', sysnum, si);
        end
    end
end

plt_set.y_font_dim     = 11;
plt_set.x_font_dim     = 10;
plt_set.title_font_dim = 10.5;
plt_set.thick          = 1.6;
plt_set.plot_unit      = 'centimeters';
plt_set.fontname       = 'Times';
plt_set.fontsize       = 10;

plt_set.plot_dim_x     = 9;
plt_set.plot_dim_y     = 6;

case_colors = ["#007191"; "#62c8d3"; "#f47a00"; "#d31f11"];

ncase = numel(SIs);
if ncase <= numel(case_colors)
    use_colors = case_colors(1:ncase);
else
    use_colors = repmat(case_colors, ceil(ncase/numel(case_colors)), 1);
    use_colors = use_colors(1:ncase);
end

fig = figure('Color','w','Name','Median RMSE (sys x case)');
fig.Units = plt_set.plot_unit;
fig.Position(3) = plt_set.plot_dim_x;
fig.Position(4) = plt_set.plot_dim_y;

fig.Renderer = 'opengl';

b = bar(medRMSE, 'grouped');
grid on;

for j = 1:ncase
    b(j).FaceColor = hex2rgb(use_colors(j));
    b(j).EdgeColor = 'none';
end

ax = gca;
ax.FontName  = plt_set.fontname;
ax.FontSize  = plt_set.fontsize;
ax.LineWidth = 1.0;
ax.XTick = 1:numel(sys_list);
ax.XTickLabel = arrayfun(@(s) sprintf('Sys=%d', s), sys_list, 'UniformOutput', false);

xlabel('');
ylab = ylabel('med(RMSE_{y})');
set(ylab, 'FontName', plt_set.fontname, 'FontSize', plt_set.y_font_dim);

lgd = legend(arrayfun(@(x) sprintf('case=%d', x), SIs, 'UniformOutput', false), ...
             'Location','northoutside', ...
             'Orientation','horizontal', ...
             'NumColumns', 2);
set(lgd, 'FontName', plt_set.fontname, 'FontSize', 10, 'Box','off');

try
    lgd.ItemTokenSize = [12 8];
catch
end

set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));
fig.PaperPositionMode = 'auto';
drawnow;

outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

base = sprintf('MedianRMSE_GroupedBar_ft%02d', f_target);

savefig(fig, fullfile(outdir, [base, '.fig']));
exportgraphics(fig, fullfile(outdir, [base, '.png']), 'Resolution',300, 'BackgroundColor','white');

fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

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
