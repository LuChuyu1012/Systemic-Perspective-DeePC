clc; clear; close all;

sysnum   = 6;
f_target = 2;
SIs      = 0:3;
nruns    = 20;
Fs       = 100;

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

Umean_all = cell(1, numel(use_cases));
N_all     = NaN(1, numel(use_cases));
n_used    = zeros(1, numel(use_cases));

for kcase = 1:numel(use_cases)
    si = use_cases(kcase);

    fname_big = sprintf('sf%02d_r%02d_case%02d.mat', sysnum, f_target, si);
    if ~isfile(fname_big)
        fprintf('Skip (file not found): %s\n', fname_big);
        continue;
    end

    Us   = cell(1, nruns);
    lens = NaN(1, nruns);

    for i = 1:nruns
        varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
        try
            S = load(fname_big, varname);
            if ~isfield(S, varname), continue; end
            rec = S.(varname);

            if ~isfield(rec, 'u_seq_opt') || isempty(rec.u_seq_opt), continue; end
            u = rec.u_seq_opt(:);

            Us{i}   = u;
            lens(i) = numel(u);
        catch
        end
    end

    good = find(~cellfun(@isempty, Us) & isfinite(lens) & lens > 0);
    if isempty(good)
        fprintf('Skip (no valid u_seq_opt): case=%d\n', si);
        continue;
    end

    N = min(lens(good));

    Umat = NaN(N, numel(good));
    for j = 1:numel(good)
        u = Us{good(j)};
        Umat(:,j) = u(1:N);
    end

    Umean_all{kcase} = mean(Umat, 2, 'omitnan');
    N_all(kcase)     = N;
    n_used(kcase)    = numel(good);

    fprintf('case=%d | used %d/%d runs | N=%d\n', si, n_used(kcase), nruns, N);
end

valid = find(~cellfun(@isempty, Umean_all) & isfinite(N_all));
if isempty(valid)
    error('No valid data read for any case.');
end

Nmin = min(N_all(valid));
t = (0:Nmin-1).';

fig = figure('Color','w','Name','Mean input (all cases)','Visible','on');
fig.Units = plt_set.plot_unit;
fig.Position(3) = plt_set.plot_dim_x;
fig.Position(4) = plt_set.plot_dim_y;

fig.Renderer = 'opengl';

ax = axes(fig); hold(ax,'on');

for kcase = 1:numel(use_cases)
    if isempty(Umean_all{kcase}), continue; end
    ubar = Umean_all{kcase}(1:Nmin);
    plot(ax, t, ubar, '-', 'Color', use_rgb(kcase,:), 'LineWidth', plt_set.thick);
end

ax.FontName  = plt_set.fontname;
ax.FontSize  = plt_set.fontsize;
ax.LineWidth = 1.0;

xlab = xlabel(ax, 't [samples]');
ylab = ylabel(ax, 'u(t)');
set(xlab, 'FontName', plt_set.fontname, 'FontSize', plt_set.x_font_dim);
set(ylab, 'FontName', plt_set.fontname, 'FontSize', plt_set.y_font_dim);

grid(ax,'on');
ax.XGrid = 'off';
ax.YGrid = 'on';
ax.GridAlpha = 0.15;

xlim(ax, [t(1) t(end)]);

set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));
fig.PaperPositionMode = 'auto';

drawnow;

outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

base = sprintf('MeanTrajU_AllCases_sys%02d_ft%02d', sysnum, f_target);

savefig(fig, fullfile(outdir, [base, '.fig']));
exportgraphics(fig, fullfile(outdir, [base, '.png']), ...
    'Resolution', 300, 'BackgroundColor', 'white');

fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

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
