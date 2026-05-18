% Plot the axial magnetic field profile B0(z) at r = 0.
% Default: use B0_load.txt plus 1parameter.dat, i.e. the B0 field interpolated
% onto the main simulation grid. Set use_dense_output = true to plot the
% high-resolution auxiliary field from background_b0 on its fixed [-1, 1] m grid.
clear; clc; close all;

use_dense_output = false;
i_switch_print = 0;
resolution = '-r450';

script_path = mfilename('fullpath');
if isempty(script_path)
    script_dir = pwd;
else
    script_dir = fileparts(script_path);
end

candidate_dirs = unique({pwd, script_dir, fileparts(script_dir)}, 'stable');
data_dir = '';

for k = 1:numel(candidate_dirs)
    if ~use_dense_output && isfile(fullfile(candidate_dirs{k}, 'B0_load.txt')) && ...
            isfile(fullfile(candidate_dirs{k}, '1parameter.dat'))
        data_dir = candidate_dirs{k};
        break;
    end
end

if isempty(data_dir)
    for k = 1:numel(candidate_dirs)
        if use_dense_output && isfile(fullfile(candidate_dirs{k}, '2for_plot_B0.txt')) && ...
                isfile(fullfile(candidate_dirs{k}, '2for_plot_z_B0.txt'))
            data_dir = candidate_dirs{k};
            break;
        end
    end
end

if isempty(data_dir)
    if use_dense_output
        error('Cannot find 2for_plot_B0.txt and 2for_plot_z_B0.txt.');
    else
        error('Cannot find B0_load.txt and 1parameter.dat.');
    end
end

if use_dense_output
    z = load(fullfile(data_dir, '2for_plot_z_B0.txt'));
    z = z(:).';

    raw = load(fullfile(data_dir, '2for_plot_B0.txt'));
    mz = numel(z);
    if size(raw, 1) < 2 * mz
        error('2for_plot_B0.txt has %d rows, but z file indicates mz=%d.', size(raw, 1), mz);
    end

    br = raw(1:mz, :);
    bz = raw(mz + (1:mz), :);
    br_axis = br(:, 1).';
    bz_axis = bz(:, 1).';
    input_name = '2for\_plot\_B0.txt';
else
    para = load(fullfile(data_dir, '1parameter.dat'));
    nr = round(para(1));
    nz = round(para(2));
    zs = para(8);
    zl = para(10);
    zd = zs + zl;

    raw = load(fullfile(data_dir, 'B0_load.txt'));
    if size(raw, 1) < 2 * nz || size(raw, 2) < nr
        error('B0_load.txt shape is %dx%d, but nr=%d and nz=%d.', ...
            size(raw, 1), size(raw, 2), nr, nz);
    end

    z = zs + (0:nz-1) * (zd - zs) / (nz - 1);
    br = raw(1:nz, 1:nr);
    bz = raw(nz + (1:nz), 1:nr);
    br_axis = br(:, 1).';
    bz_axis = bz(:, 1).';
    input_name = 'B0\_load.txt';
end

b0_axis = sqrt(br_axis.^2 + bz_axis.^2);

figure('Color', 'w');
plot(z, b0_axis, 'k-', 'LineWidth', 2.2);
hold on;
plot(z, bz_axis, 'r--', 'LineWidth', 1.6);
grid on;
xlim([min(z), max(z)]);
xlabel('z (m)');
ylabel('B_0 (T)');
legend('|B_0| at r=0', 'B_z at r=0', 'Location', 'best');
title(sprintf('Axial magnetic field: %s', input_name));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

fprintf('Loaded B0 data from %s\n', data_dir);
fprintf('max |B0| on axis = %.6g T\n', max(b0_axis));

if i_switch_print > 0.1
    print(fullfile(data_dir, 'B0_axis.png'), '-dpng', resolution);
end
