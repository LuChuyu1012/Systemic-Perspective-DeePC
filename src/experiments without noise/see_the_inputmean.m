clc; clear; close all;

sysnum   = 5;
f_target = 2;
SIs      = 0:3;
nruns    = 50;
Fs       = 100;

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

Umean_all = cell(1, num_cases);
N_all     = NaN(1, num_cases);
n_used    = zeros(1, num_cases);

for c = 1:num_cases
    si = SIs(c);

    fname_best = fullfile(data_dir, sprintf('s%02d_r%02d_case%02d.mat', ...
                         sysnum, f_target, si));
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
        warning('case=%d: no u_seq_opt found, skip.', si);
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

    fprintf('case=%d | used %d/%d runs | N=%d\n', si, n_used(c), nruns, N);
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
%     base = sprintf('MeanTrajU_AllCases_sys%02d_ft%02d', ...
%                    sysnum, f_target);
% 
%     savefig(fig, fullfile(outdir, [base, '.fig']));
%     exportgraphics(fig, fullfile(outdir,[base,'.png']), 'Resolution',300, 'BackgroundColor','white');
% 
%     fprintf('FIG saved to: %s\n', fullfile(outdir,[base,'.fig']));
%     fprintf('PNG saved to: %s\n', fullfile(outdir,[base,'.png']));

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