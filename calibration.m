clc;
clear;
close all;

%% Step 1: Read Excel file
data = readtable('calibration3.xlsx');

imageNames = data.filename;
concentration = data.concentration;

n = length(imageNames);

%% Step 2: Read reference image (first image, no dye)
refPath = fullfile('calibration_photos_3', imageNames{1});
refImg = im2double(imread(refPath));

% Crop: remove top 1000 pixels and bottom 200 pixels
refImg = refImg(1001:end-200, :, :);

% Avoid division by zero
epsilon = 1e-6;
refImg = refImg + epsilon;

%% Step 3: Initialize storage
R_mean = zeros(n,1);
G_mean = zeros(n,1);
B_mean = zeros(n,1);

R_std = zeros(n,1);
G_std = zeros(n,1);
B_std = zeros(n,1);

%% Step 4: Loop through images
for i = 1:n
    
    imgPath = fullfile('calibration_photos_3', imageNames{i});
    img = im2double(imread(imgPath));
    
    % Apply same cropping
    img = img(1001:end-200, :, :);
    
    % Ensure same size
    if ~isequal(size(img), size(refImg))
        error('Image sizes do not match after cropping!');
    end
    
    % Pixel-wise calibration
    calibratedImg = img ./ refImg;
    
    % Extract channels
    R = calibratedImg(:,:,1);
    G = calibratedImg(:,:,2);
    B = calibratedImg(:,:,3);
    
    % Mean values
    R_mean(i) = mean(R(:));
    G_mean(i) = mean(G(:));
    B_mean(i) = mean(B(:));
    
    % Spatial uncertainty (pixel variation)
    R_std(i) = std(R(:));
    G_std(i) = std(G(:));
    B_std(i) = std(B(:));
end

%% Step 4.5: Group repeated measurements
[unique_conc, ~, idx_group] = unique(concentration);
n_groups = length(unique_conc);

R_group_mean = zeros(n_groups,1);
G_group_mean = zeros(n_groups,1);
B_group_mean = zeros(n_groups,1);

R_group_std = zeros(n_groups,1);
G_group_std = zeros(n_groups,1);
B_group_std = zeros(n_groups,1);

N_repeats = zeros(n_groups,1);

for j = 1:n_groups
    group_idx = (idx_group == j);
    N_repeats(j) = sum(group_idx);
    
    % Mean across repeated images
    R_group_mean(j) = mean(R_mean(group_idx));
    G_group_mean(j) = mean(G_mean(group_idx));
    B_group_mean(j) = mean(B_mean(group_idx));
    
    % Spatial uncertainty only
    R_group_std(j) = sqrt(mean(R_std(group_idx).^2));
    G_group_std(j) = sqrt(mean(G_std(group_idx).^2));
    B_group_std(j) = sqrt(mean(B_std(group_idx).^2));
end

%% Step 5: Write calibration CSV
% Main calibration table
calibrationTable = table( ...
    unique_conc, ...
    R_group_mean, G_group_mean, B_group_mean, ...
    R_group_std,  G_group_std,  B_group_std, ...
    N_repeats, ...
    'VariableNames', { ...
        'concentration', ...
        'R', 'G', 'B', ...
        'R_std', 'G_std', 'B_std', ...
        'n_repeats'});

writetable(calibrationTable, 'rgb_concentration_calibration.csv');

disp('Saved rgb_concentration_calibration.csv');

%% Step 6: Log-scale plots
figure;

subplot(3,2,1);
plot(concentration, R_mean, 'r-o', 'LineWidth', 1.5);
xlabel('Dye Concentration');
ylabel('Red Channel');
title('Red vs Concentration');
grid on;

subplot(3,2,3);
plot(concentration, G_mean, 'g-o', 'LineWidth', 1.5);
xlabel('Dye Concentration');
ylabel('Green Channel');
title('Green vs Concentration');
grid on;

subplot(3,2,5);
plot(concentration, B_mean, 'b-o', 'LineWidth', 1.5);
xlabel('Dye Concentration');
ylabel('Blue Channel');
title('Blue vs Concentration');
grid on;

sgtitle('RGB Values vs Dye Concentration');

epsilon = 1e-10;
R_log = R_mean + epsilon;
G_log = G_mean + epsilon;
B_log = B_mean + epsilon;

subplot(3,2,2);
semilogy(concentration, R_log, 'r-o', 'LineWidth', 1.5);
xlabel('Dye Concentration');
ylabel('Red Channel (log)');
title('Red vs Concentration (Log Scale)');
grid on;

subplot(3,2,4);
semilogy(concentration, G_log, 'g-o', 'LineWidth', 1.5);
xlabel('Dye Concentration');
ylabel('Green Channel (log)');
title('Green vs Concentration (Log Scale)');
grid on;

%% Step 7: Weighted linear regression on Blue (log scale)
fit_idx = unique_conc >= 0 & unique_conc <= 0.8e-3;

x = unique_conc(fit_idx);
y = log(B_group_mean(fit_idx));

sigma = B_group_std(fit_idx) ./ B_group_mean(fit_idx);

valid = isfinite(x) & isfinite(y) & isfinite(sigma) & sigma > 0;
x = x(valid);
y = y(valid);
sigma = sigma(valid);

w = 1 ./ (sigma.^2);

X = [x ones(size(x))];
W = diag(w);

beta = (X' * W * X) \ (X' * W * y);

m = beta(1);
b = beta(2);

cov_beta = inv(X' * W * X);
m_err = sqrt(cov_beta(1,1));
b_err = sqrt(cov_beta(2,2));

x_line = linspace(min(x), max(x), 100);
y_line = m*x_line + b;

y_pred = X * beta;
SS_res = sum(w .* (y - y_pred).^2);
SS_tot = sum(w .* (y - mean(y)).^2);
R2 = 1 - SS_res / SS_tot;

subplot(3,2,6);

errorbar(unique_conc, B_group_mean, B_group_std, ...
    'bo', 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
hold on;

semilogy(x_line, exp(y_line), 'k--', 'LineWidth', 2);

xlabel('Dye Concentration');
ylabel('Blue Channel (log)');
title('Blue vs Concentration (Log Scale)');
grid on;

eqn_text = sprintf(['ln(B) = (%.3e ± %.1e)x + (%.3e ± %.1e)\n' ...
                    'R^2 = %.4f'], ...
                    m, m_err, b, b_err, R2);

x_pos = min(unique_conc) + 0.4*(max(unique_conc)-min(unique_conc));
y_pos = max(B_group_mean) * 0.4;

text(x_pos, y_pos, eqn_text, ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k');

annotation('textbox', [0.15, 0.01, 0.7, 0.05], ...
    'String', 'Error bars: Top = mean + standard deviation, Bottom = mean - standard deviation', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 13);

hold off;

%% Step 8: Save blue-fit model to CSV
blueFitTable = table(m, m_err, b, b_err, R2, ...
    'VariableNames', {'slope_m', 'slope_unc', 'intercept_b', 'intercept_unc', 'R_squared'});

writetable(blueFitTable, 'blue_log_fit_parameters.csv');

disp('Saved blue_log_fit_parameters.csv');