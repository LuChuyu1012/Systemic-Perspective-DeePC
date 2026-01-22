clc; clear; close all;

sys_list  = 1:1;
f_target  = 2;
SIs       = 4:11;
nruns     = 20;
Fs        = 100;

consensus_tol_abs = 0.04;
win_periods_chk   = 0;

P       = max(1, round(Fs / max(f_target, eps)));
win_chk = max(1, round(win_periods_chk * P));

y_fields = {'y_seq_opt','y_seq','y','y_out'};

ConvTime = NaN(numel(sys_list), numel(SIs));

for kk = 1:numel(sys_list)
    sysnum = sys_list(kk);

    for sidx = 1:numel(SIs)
        si = SIs(sidx);

        fname_big = sprintf('sf%02d_r%02d_case%02d.mat', sysnum, f_target, si);
        if ~isfile(fname_big)
            warning('error: %s (skip sys=%d case=%d)', fname_big, sysnum, si);
            continue;
        end

        Y_all = cell(1, nruns);
        lens  = zeros(1, nruns);

        for i = 1:nruns
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            try
                S = load(fname_big, varname);
                if ~isfield(S, varname), continue; end
                rec = S.(varname);

                y = [];
                for fn = y_fields
                    if isfield(rec, fn{1}) && ~isempty(rec.(fn{1}))
                        y = rec.(fn{1});
                        break;
                    end
                end
                if isempty(y), continue; end

                y = squeeze(y);
                y = y(:,:);

                lens(i)  = size(y,1);
                Y_all{i} = y;

            catch ME
                warning('error: %s | %s', varname, ME.message);
            end
        end

        good = find(~cellfun(@isempty, Y_all) & lens>0);
        if isempty(good)
            warning('sys=%d case=%d invalid y', sysnum, si);
            continue;
        end

        Nmin = min(lens(good));

        n = (0:Nmin-1).';
        yref = sin(2*pi*f_target/Fs * n);

        ny = size(Y_all{good(1)}, 2);
        if ny > 1
            yref = repmat(yref, 1, ny);
        end

        D = NaN(Nmin, numel(good));
        for g = 1:numel(good)
            y = Y_all{good(g)}(1:Nmin, :);

            diff = y - yref;
            if size(diff,2) == 1
                D(:,g) = abs(diff);
            else
                D(:,g) = sqrt(sum(diff.^2, 2));
            end
        end

        dmax = max(D, [], 2);

        if win_chk > 1
            dchk = movmax(dmax, [win_chk-1, 0], 'Endpoints','shrink');
        else
            dchk = dmax;
        end

        below = (dchk <= consensus_tol_abs);

        tail_ok = flipud(cummin(flipud(double(below)))) == 1;
        n0 = find(tail_ok, 1, 'first');

        if ~isempty(n0)
            ConvTime(kk, sidx) = (n0 - 1)/Fs;
        else
            ConvTime(kk, sidx) = NaN;
        end
    end
end

labels_case = arrayfun(@(x) sprintf('case=%d', x), SIs, 'UniformOutput', false);

plt_set.y_font_dim     = 11;
plt_set.x_font_dim     = 10;
plt_set.thick          = 1.6;
plt_set.plot_unit      = 'centimeters';
plt_set.fontname       = 'Times';
plt_set.fontsize       = 10;

plt_set.plot_dim_x     = 9;
plt_set.plot_dim_y     = 6;

case_hex = ["#007191"; "#62c8d3"; "#f47a00"; "#d31f11"];

ncase = numel(SIs);
if ncase <= numel(case_hex)
    use_hex = case_hex(1:ncase);
else
    use_hex = repmat(case_hex, ceil(ncase/numel(case_hex)), 1);
    use_hex = use_hex(1:ncase);
end

use_rgb = zeros(ncase,3);
for j = 1:ncase
    use_rgb(j,:) = hex2rgb(use_hex(j));
end

fig = figure('Color','w','Name','Consensus-to-reference time (sys x case)');
fig.Units = plt_set.plot_unit;
fig.Position(3) = plt_set.plot_dim_x;
fig.Position(4) = plt_set.plot_dim_y;

fig.Renderer = 'opengl';

b = bar(ConvTime, 'grouped'); hold on;

for j = 1:ncase
    b(j).FaceColor = use_rgb(j,:);
    b(j).EdgeColor = 'none';
end

ax = gca;
ax.FontName = plt_set.fontname;
ax.FontSize = plt_set.fontsize;
ax.LineWidth = 1.0;
ax.XTick = 1:numel(sys_list);
ax.XTickLabel = arrayfun(@(s) sprintf('Sys=%d', s), sys_list, 'UniformOutput', false);

ylab = ylabel('c [s]');
set(ylab, 'FontName', plt_set.fontname, 'FontSize', plt_set.y_font_dim);

grid on;
ax.XGrid = 'off';
ax.YGrid = 'on';
ax.GridAlpha = 0.15;

hleg = gobjects(ncase,1);
for j = 1:ncase
    hleg(j) = line(NaN, NaN, ...
        'LineStyle','none', ...
        'Marker','s', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', use_rgb(j,:), ...
        'MarkerEdgeColor', use_rgb(j,:));
end

lgd = legend(hleg, labels_case, ...
    'Location','northoutside', 'Orientation','horizontal', 'NumColumns',2);
set(lgd, 'FontName', plt_set.fontname, 'FontSize', 10, 'Box','off');

set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));
fig.PaperPositionMode = 'auto';
drawnow;

outdir = fullfile(pwd, 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

base = sprintf('ConvTime_GroupedBar_ft%02d_tol%g', f_target, consensus_tol_abs);
base = strrep(base, '.', 'p');

savefig(fig, fullfile(outdir, [base, '.fig']));

exportgraphics(fig, fullfile(outdir, [base, '.png']), ...
    'Resolution',300, 'BackgroundColor','white');

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
