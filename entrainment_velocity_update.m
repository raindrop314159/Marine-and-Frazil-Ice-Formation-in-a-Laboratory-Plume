clear; clc; close all;

%% folders to analyse

datasets = {
    'plume',             'No bubbles, 15.8 Brix',       0,     15.8, false;
    'plume_46V_bubbles', 'Bubble flux = 4.42 ml/s',     4.42,  15.8, true;
    'plume_63V_bubbles', 'Bubble flux = 6.50 ml/s',     6.50,  15.8, true;
    'plume_73V_bubbles', 'Bubble flux = 8.32 ml/s',     8.32,  15.8, true;
};

imagePattern = 'DSC_*.JPG';

dt = 2;
pumpRate_ml_s = 2;
tankInitialVolume_ml = 8000;

%% constants

g = 9.81;
rho_fresh = 1000;
rho_air = 1.2;

Q_water_m3_s = pumpRate_ml_s * 1e-6;

%% analysis settings

mean_end_time = 50;
maxFrames = 15;

Bdiff_threshold_low = 0.35;
Bdiff_threshold_high = 3;

strip_x1 = 3100;
strip_x2 = 3800;

tank_depth_cm = 22;
cm_per_pixel = 0.0111111;
cm2_per_pixel2 = cm_per_pixel^2;

crop_top = 1000;
crop_bottom = 300;

smooth_window = 5;

gaussian_sigma = 2;
min_object_size_px = 80;
closing_disk_size = 4;

hBottom_min_px = 150;

%% start parallel pool

if isempty(gcp('nocreate'))
    parpool;
end

%% storage

nDatasets = size(datasets,1);

emptyResult = struct( ...
    'folder', [], ...
    'label', [], ...
    'bubbleFlux_ml_s', [], ...
    'salinity_brix', [], ...
    'hasBubbles', [], ...
    'rho_ambient', [], ...
    'g_reduced', [], ...
    'B_freshwater', [], ...
    'B_bubbles', [], ...
    'B_total', [], ...
    'meanEndTime_s', [], ...
    'meanEntrainmentRate_ml_s', [], ...
    'peakEntrainmentRate_ml_s', [], ...
    'meanUe_eff_cm_s', [], ...
    'peakUe_eff_cm_s', [], ...
    'time_s', [], ...
    'pumpedVolume_ml', [], ...
    'thresholdedVolume_ml', [], ...
    'entrainedVolume_ml', [], ...
    'entrained_smoothed', [], ...
    'time_rate_s', [], ...
    'entrainmentRate_ml_s', [], ...
    'ue_eff_cm_s', [], ...
    'plumeArea_px', [], ...
    'plumeArea_cm2', [], ...
    'plumeHeight_px', [], ...
    'plumeHeight_cm', [], ...
    'hTop_px', [], ...
    'hTop_cm', [], ...
    'hBottom_px', [], ...
    'hBottom_cm', [], ...
    'interfaceArea_cm2', [], ...
    'stripBdiffPixels_px', [], ...
    'stripTotalPixels_px', [], ...
    'stripBdiffPercentage', [], ...
    'strip_x1', [], ...
    'strip_x2', [] ...
);

allResults = repmat(emptyResult, nDatasets, 1);

%% process folders in parallel

parfor d = 1:nDatasets

    imageFolder = datasets{d,1};
    labelName = datasets{d,2};
    bubbleFlux_ml_s = datasets{d,3};
    salinity_brix = datasets{d,4};
    hasBubbles = datasets{d,5};

    %% dataset-specific buoyancy flux

    salinity_ppt = salinity_brix * 10;
    rho_ambient_local = 1000 + 0.75 * salinity_ppt;

    g_reduced_local = g * (rho_ambient_local - rho_fresh) / rho_ambient_local;
    B_freshwater_local = g_reduced_local * Q_water_m3_s;

    Q_bubble_m3_s = bubbleFlux_ml_s * 1e-6;
    B_bubbles = g * (rho_ambient_local - rho_air) / rho_ambient_local * Q_bubble_m3_s;

    B_total = B_freshwater_local + B_bubbles;

    fprintf('\nProcessing folder: %s\n', imageFolder);
    fprintf('Salinity/Brix = %.1f\n', salinity_brix);
    fprintf('Ambient density = %.2f kg/m^3\n', rho_ambient_local);
    fprintf('Reduced gravity = %.4f m/s^2\n', g_reduced_local);
    fprintf('Freshwater buoyancy flux = %.4e m^4/s^3\n', B_freshwater_local);
    fprintf('Bubble flux = %.2f ml/s\n', bubbleFlux_ml_s);
    fprintf('Bubble buoyancy flux = %.4e m^4/s^3\n', B_bubbles);
    fprintf('Total buoyancy flux = %.4e m^4/s^3\n', B_total);

    %% output folders

    overlayFolder = fullfile(imageFolder, 'threshold_overlays');
    if ~exist(overlayFolder, 'dir')
        mkdir(overlayFolder);
    end

    diagnosticFolder = fullfile(imageFolder, 'diagnostic_deltaB');
    if ~exist(diagnosticFolder, 'dir')
        mkdir(diagnosticFolder);
    end

    %% load files

    files = dir(fullfile(imageFolder, imagePattern));
    [~, idx] = sort({files.name});
    files = files(idx);

    nFrames = min(length(files), maxFrames);
    files = files(1:nFrames);

    if nFrames == 0
        error('No images found in folder: %s', imageFolder);
    end

    fprintf('Using first %d images from %s\n', nFrames, imageFolder);

    %% reference image

    I0 = im2double(imread(fullfile(imageFolder, files(1).name)));

    if crop_bottom > 0
        I0 = I0(crop_top+1:end-crop_bottom, :, :);
    else
        I0 = I0(crop_top+1:end, :, :);
    end

    B0 = I0(:,:,3);
    waterMask = true(size(B0));

    %% strip limits

    imageWidth = size(B0, 2);

    local_strip_x1 = strip_x1;
    local_strip_x2 = strip_x2;

    if local_strip_x2 > imageWidth
        local_strip_x2 = imageWidth;
    end

    if local_strip_x1 < 1
        local_strip_x1 = 1;
    end

    if local_strip_x1 >= local_strip_x2
        error('Invalid strip range in %s.', imageFolder);
    end

    stripMask = false(size(B0));
    stripMask(:, local_strip_x1:local_strip_x2) = true;

    %% effective tank depth

    tankArea_px = nnz(waterMask);
    tankArea_cm2 = tankArea_px * cm2_per_pixel2;
    effective_depth_cm = tankInitialVolume_ml / tankArea_cm2;

    %% allocate variables

    time_s = (0:nFrames-1)' * dt;
    pumpedVolume_ml = time_s * pumpRate_ml_s;

    plumeArea_px = zeros(nFrames,1);
    plumeArea_cm2 = zeros(nFrames,1);

    plumeHeight_px = zeros(nFrames,1);
    plumeHeight_cm = zeros(nFrames,1);

    hTop_px = zeros(nFrames,1);
    hTop_cm = zeros(nFrames,1);

    hBottom_px = zeros(nFrames,1);
    hBottom_cm = zeros(nFrames,1);

    interfaceArea_cm2 = zeros(nFrames,1);

    thresholdedVolume_ml = zeros(nFrames,1);
    entrainedVolume_ml = zeros(nFrames,1);

    stripBdiffPixels_px = zeros(nFrames,1);
    stripTotalPixels_px = zeros(nFrames,1);
    stripBdiffPercentage = zeros(nFrames,1);

    %% process images

    for i = 1:nFrames

        I = im2double(imread(fullfile(imageFolder, files(i).name)));

        if crop_bottom > 0
            I = I(crop_top+1:end-crop_bottom, :, :);
        else
            I = I(crop_top+1:end, :, :);
        end

        B = I(:,:,3);

        %% background-subtracted blue channel

        deltaB = B0 - B;

        deltaB(isnan(deltaB)) = 0;
        deltaB(isinf(deltaB)) = 0;
        deltaB(deltaB < 0) = 0;

        deltaB_smooth = imgaussfilt(deltaB, gaussian_sigma);
        deltaB_enhanced = adapthisteq(mat2gray(deltaB_smooth));

        %% threshold plume mask

        dyeMask = ...
            (deltaB_enhanced > Bdiff_threshold_low) & ...
            (deltaB_enhanced < Bdiff_threshold_high);

        dyeMask = dyeMask & waterMask;

        dyeMask = bwareaopen(dyeMask, min_object_size_px);
        dyeMask = imclose(dyeMask, strel('disk', closing_disk_size));
        dyeMask = imfill(dyeMask, 'holes');

        %% diagram-based plume interface area

        [imgH, imgW] = size(dyeMask);

        % 1. Estimate total plume height H from middle of image
        mid_x = round(imgW/2);
        mid_band_width = 100;

        mid_x1 = max(1, mid_x - mid_band_width);
        mid_x2 = min(imgW, mid_x + mid_band_width);

        middleMask = dyeMask(:, mid_x1:mid_x2);
        rows_middle = find(any(middleMask, 2));

        if isempty(rows_middle)
            plumeHeight_px(i) = NaN;
        else
            plume_top_row = min(rows_middle);
            plume_bottom_row = max(rows_middle);
            plumeHeight_px(i) = plume_bottom_row - plume_top_row + 1;
        end

        plumeHeight_cm(i) = plumeHeight_px(i) * cm_per_pixel;

        % 2. Estimate h_top from red pixels in selected side strip
        stripMaskNow = dyeMask(:, local_strip_x1:local_strip_x2);
        rows_strip = find(any(stripMaskNow, 2));

        if isempty(rows_strip)
            hTop_px(i) = 0;
        else
            strip_top_row = min(rows_strip);
            strip_bottom_row = max(rows_strip);
            hTop_px(i) = strip_bottom_row - strip_top_row + 1;
        end

        hTop_cm(i) = hTop_px(i) * cm_per_pixel;

        % 3. Calculate h_bottom = H - h_top
        hBottom_px(i) = plumeHeight_px(i) - hTop_px(i);

        if isnan(hBottom_px(i))
            hBottom_px(i) = NaN;
        else
            hBottom_px(i) = max(hBottom_px(i), hBottom_min_px);
        end

        hBottom_cm(i) = hBottom_px(i) * cm_per_pixel;

        % 4. Interface area: two plume sides times tank thickness
        interfaceArea_cm2(i) = 2 * hBottom_cm(i) * tank_depth_cm;

        % Boundary only for overlay display
        boundaryMask = bwperim(dyeMask);

        %% strip percentage

        stripRangeMask = ...
            (deltaB_enhanced > Bdiff_threshold_low) & ...
            (deltaB_enhanced < Bdiff_threshold_high) & ...
            stripMask;

        stripBdiffPixels_px(i) = nnz(stripRangeMask);
        stripTotalPixels_px(i) = nnz(stripMask);

        stripBdiffPercentage(i) = ...
            100 * stripBdiffPixels_px(i) / stripTotalPixels_px(i);

        %% cyan overlay

        overlayImg = I;

        redLayer = overlayImg(:,:,1);
        greenLayer = overlayImg(:,:,2);
        blueLayer = overlayImg(:,:,3);

        alpha = 0.45;

        redLayer(dyeMask)   = (1-alpha)*redLayer(dyeMask);
        greenLayer(dyeMask) = (1-alpha)*greenLayer(dyeMask) + alpha;
        blueLayer(dyeMask)  = (1-alpha)*blueLayer(dyeMask)  + alpha;

        overlayImg(:,:,1) = redLayer;
        overlayImg(:,:,2) = greenLayer;
        overlayImg(:,:,3) = blueLayer;

        overlayImg(:,:,1) = overlayImg(:,:,1) .* ~boundaryMask;
        overlayImg(:,:,2) = overlayImg(:,:,2) + boundaryMask;
        overlayImg(:,:,3) = overlayImg(:,:,3) + boundaryMask;

        %% strip boundary in green

        overlayImg(:, local_strip_x1, 1) = 0;
        overlayImg(:, local_strip_x1, 2) = 1;
        overlayImg(:, local_strip_x1, 3) = 0;

        overlayImg(:, local_strip_x2, 1) = 0;
        overlayImg(:, local_strip_x2, 2) = 1;
        overlayImg(:, local_strip_x2, 3) = 0;

        overlayImg = min(overlayImg, 1);

        [~, baseName, ~] = fileparts(files(i).name);
        outName = fullfile(overlayFolder, [baseName '_threshold_overlay.png']);
        imwrite(overlayImg, outName);

        %% diagnostic image

        diagnosticName = fullfile(diagnosticFolder, [baseName '_deltaB.png']);
        imwrite(deltaB_enhanced, diagnosticName);

        %% area and volume

        plumeArea_px(i) = nnz(dyeMask);
        plumeArea_cm2(i) = plumeArea_px(i) * cm2_per_pixel2;

        thresholdedVolume_ml(i) = plumeArea_cm2(i) * effective_depth_cm;

        maxTankVolume_ml = tankInitialVolume_ml + pumpedVolume_ml(i);
        thresholdedVolume_ml(i) = min(thresholdedVolume_ml(i), maxTankVolume_ml);

        entrainedVolume_ml(i) = thresholdedVolume_ml(i) - pumpedVolume_ml(i);
        entrainedVolume_ml(i) = max(entrainedVolume_ml(i), 0);

        fprintf(['Dataset %d/%d | Frame %3d/%3d  %s  time = %.1f s  ' ...
                 'thresholded = %.2f ml  pumped = %.2f ml  entrained = %.2f ml  ' ...
                 'H = %.2f cm  h_top = %.2f cm  h_bottom = %.2f cm  A_p = %.2f cm^2\n'], ...
            d, nDatasets, i, nFrames, files(i).name, time_s(i), ...
            thresholdedVolume_ml(i), pumpedVolume_ml(i), entrainedVolume_ml(i), ...
            plumeHeight_cm(i), hTop_cm(i), hBottom_cm(i), interfaceArea_cm2(i));

    end

    %% entrainment rate

    entrained_smoothed = movmean(entrainedVolume_ml, smooth_window);

    entrainmentRate_ml_s = diff(entrained_smoothed) / dt;
    time_rate_s = time_s(1:end-1) + dt/2;

    entrainmentRate_ml_s(entrainmentRate_ml_s < 0) = 0;

    %% effective entrainment velocity

    interfaceArea_rate_cm2 = ...
        0.5 * (interfaceArea_cm2(1:end-1) + interfaceArea_cm2(2:end));

    ue_eff_cm_s = entrainmentRate_ml_s ./ interfaceArea_rate_cm2;

    ue_eff_cm_s(interfaceArea_rate_cm2 <= 0) = NaN;
    ue_eff_cm_s(isinf(ue_eff_cm_s)) = NaN;

    %% mean values up to 50 s

    validMean = time_rate_s <= mean_end_time;

    meanEntrainmentRate_ml_s = ...
        mean(entrainmentRate_ml_s(validMean), 'omitnan');

    peakEntrainmentRate_ml_s = max(entrainmentRate_ml_s);

    meanUe_eff_cm_s = ...
        mean(ue_eff_cm_s(validMean), 'omitnan');

    peakUe_eff_cm_s = ...
        max(ue_eff_cm_s, [], 'omitnan');

    %% result struct

    result = emptyResult;

    result.folder = imageFolder;
    result.label = labelName;

    result.bubbleFlux_ml_s = bubbleFlux_ml_s;
    result.salinity_brix = salinity_brix;
    result.hasBubbles = hasBubbles;

    result.rho_ambient = rho_ambient_local;
    result.g_reduced = g_reduced_local;
    result.B_freshwater = B_freshwater_local;
    result.B_bubbles = B_bubbles;
    result.B_total = B_total;

    result.meanEndTime_s = mean_end_time;
    result.meanEntrainmentRate_ml_s = meanEntrainmentRate_ml_s;
    result.peakEntrainmentRate_ml_s = peakEntrainmentRate_ml_s;

    result.meanUe_eff_cm_s = meanUe_eff_cm_s;
    result.peakUe_eff_cm_s = peakUe_eff_cm_s;

    result.time_s = time_s;
    result.pumpedVolume_ml = pumpedVolume_ml;
    result.thresholdedVolume_ml = thresholdedVolume_ml;
    result.entrainedVolume_ml = entrainedVolume_ml;
    result.entrained_smoothed = entrained_smoothed;
    result.time_rate_s = time_rate_s;
    result.entrainmentRate_ml_s = entrainmentRate_ml_s;
    result.ue_eff_cm_s = ue_eff_cm_s;

    result.plumeArea_px = plumeArea_px;
    result.plumeArea_cm2 = plumeArea_cm2;

    result.plumeHeight_px = plumeHeight_px;
    result.plumeHeight_cm = plumeHeight_cm;

    result.hTop_px = hTop_px;
    result.hTop_cm = hTop_cm;

    result.hBottom_px = hBottom_px;
    result.hBottom_cm = hBottom_cm;

    result.interfaceArea_cm2 = interfaceArea_cm2;

    result.stripBdiffPixels_px = stripBdiffPixels_px;
    result.stripTotalPixels_px = stripTotalPixels_px;
    result.stripBdiffPercentage = stripBdiffPercentage;

    result.strip_x1 = local_strip_x1;
    result.strip_x2 = local_strip_x2;

    allResults(d) = result;

    %% save csv files

    results = table(time_s, pumpedVolume_ml, thresholdedVolume_ml, ...
        entrainedVolume_ml, entrained_smoothed, plumeArea_px, plumeArea_cm2, ...
        plumeHeight_px, plumeHeight_cm, hTop_px, hTop_cm, ...
        hBottom_px, hBottom_cm, interfaceArea_cm2, ...
        stripBdiffPixels_px, stripTotalPixels_px, stripBdiffPercentage);

    writetable(results, fullfile(imageFolder, 'entrained_volume_results.csv'));

    rate_results = table(time_rate_s, entrainmentRate_ml_s, ue_eff_cm_s);
    writetable(rate_results, fullfile(imageFolder, 'entrainment_rate_results.csv'));

    strip_results = table(time_s, stripBdiffPixels_px, stripTotalPixels_px, ...
        stripBdiffPercentage);

    writetable(strip_results, fullfile(imageFolder, 'strip_bdiff_percentage_results.csv'));

    buoyancy_results = table( ...
        salinity_brix, ...
        bubbleFlux_ml_s, ...
        rho_ambient_local, ...
        g_reduced_local, ...
        B_freshwater_local, ...
        B_bubbles, ...
        B_total, ...
        mean_end_time, ...
        meanEntrainmentRate_ml_s, ...
        peakEntrainmentRate_ml_s, ...
        meanUe_eff_cm_s, ...
        peakUe_eff_cm_s);

    writetable(buoyancy_results, fullfile(imageFolder, 'buoyancy_flux_results.csv'));

end

%% plot all datasets together

figure('Units','centimeters', ...
       'Position',[2 2 29.7 10.5], ...
       'PaperUnits','centimeters', ...
       'PaperSize',[29.7 21.0], ...
       'PaperPosition',[1.5 5.5 26.7 9.5]);

tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

%% plot 1: thresholded volume

nexttile;
hold on;

for d = 1:length(allResults)
    plot(allResults(d).time_s, allResults(d).thresholdedVolume_ml, ...
        'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
end

plot(allResults(1).time_s, allResults(1).pumpedVolume_ml, ...
    'k--', 'LineWidth', 1.2);

xlabel('Time (s)');
ylabel('Volume (ml)');
title('(a) Volume passing threshold');

legend([string({allResults.label}), "Pumped volume"], 'Location', 'best');

grid on; box on;

%% plot 2: entrained volume

nexttile;
hold on;

for d = 1:length(allResults)
    plot(allResults(d).time_s, allResults(d).entrainedVolume_ml, ...
        'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
end

xlabel('Time (s)');
ylabel('Entrained volume (ml)');
title('(b) Entrained volume');

legend(string({allResults.label}), 'Location', 'best');

grid on; box on;

%% plot 3: effective entrainment velocity

nexttile;
hold on;

for d = 1:length(allResults)
    plot(allResults(d).time_rate_s, allResults(d).ue_eff_cm_s, ...
        'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
end

xlabel('Time (s)');
ylabel('u_{e,eff} (cm/s)');
title('(c) Effective entrainment velocity');

legend(string({allResults.label}), 'Location', 'best');

grid on; box on;

sgtitle('Plume Entrainment Analysis', 'FontWeight', 'bold');

%% strip Bdiff percentage plot

figure('Units','centimeters', ...
       'Position',[3 3 18 12], ...
       'PaperUnits','centimeters', ...
       'PaperSize',[21.0 29.7], ...
       'PaperPosition',[2 8 17 12]);

hold on;

for d = 1:length(allResults)
    plot(allResults(d).time_s, allResults(d).stripBdiffPercentage, ...
        'o-', 'LineWidth', 1.5, 'MarkerSize', 5);
end

xlabel('Time (s)');
ylabel('Pixels in Bdiff range within strip (%)');

title(sprintf('Strip filling: x = %d to %d, %.2f < Bdiff < %.2f', ...
    strip_x1, strip_x2, Bdiff_threshold_low, Bdiff_threshold_high));

legend(string({allResults.label}), 'Location', 'best');

grid on; box on;

%% buoyancy flux plot using effective entrainment velocity

B_total_all = [allResults.B_total];
meanUe_all = [allResults.meanUe_eff_cm_s];
peakUe_all = [allResults.peakUe_eff_cm_s];
hasBubbles_all = [allResults.hasBubbles];

noBubbleIdx = ~hasBubbles_all;
bubbleIdx = hasBubbles_all;

figure('Units','centimeters', ...
       'Position',[4 4 20 13], ...
       'PaperUnits','centimeters', ...
       'PaperSize',[21.0 29.7], ...
       'PaperPosition',[1 8 19 13]);

hold on;

plot(B_total_all(noBubbleIdx), meanUe_all(noBubbleIdx), ...
    'bo-', 'LineWidth', 1.7, 'MarkerSize', 8);

plot(B_total_all(noBubbleIdx), peakUe_all(noBubbleIdx), ...
    'bs--', 'LineWidth', 1.7, 'MarkerSize', 8);

plot(B_total_all(bubbleIdx), meanUe_all(bubbleIdx), ...
    'ro-', 'LineWidth', 1.7, 'MarkerSize', 8);

plot(B_total_all(bubbleIdx), peakUe_all(bubbleIdx), ...
    'rs--', 'LineWidth', 1.7, 'MarkerSize', 8);

%% scaling curves: u_e = k B_total^(1/3)

B_fit = linspace(min(B_total_all(noBubbleIdx)), max(B_total_all), 300);

k_mean_noBubble = mean(meanUe_all(noBubbleIdx) ./ ...
    (B_total_all(noBubbleIdx).^(1/3)), 'omitnan');

ue_fit_mean = k_mean_noBubble * B_fit.^(1/3);

plot(B_fit, ue_fit_mean, 'k-', 'LineWidth', 2);

k_peak_noBubble = mean(peakUe_all(noBubbleIdx) ./ ...
    (B_total_all(noBubbleIdx).^(1/3)), 'omitnan');

ue_fit_peak = k_peak_noBubble * B_fit.^(1/3);

plot(B_fit, ue_fit_peak, 'k--', 'LineWidth', 2);

for d = 1:length(allResults)
    text(B_total_all(d), peakUe_all(d), "  " + string(allResults(d).label), ...
        'FontSize', 8);
end

xlabel('Total buoyancy flux, B_{total} (m^4 s^{-3})');
ylabel('Effective entrainment velocity, u_{e,eff} (cm/s)');

title(sprintf('Effective entrainment velocity as a function of total buoyancy flux; mean over t \\leq %.0f s', mean_end_time));

legend('Mean u_{e,eff}, no bubbles', ...
       'Peak u_{e,eff}, no bubbles', ...
       'Mean u_{e,eff}, bubbles', ...
       'Peak u_{e,eff}, bubbles', ...
       'u_{e,eff} = kB_{total}^{1/3}, fitted to no-bubble mean', ...
       'u_{e,eff} = kB_{total}^{1/3}, fitted to no-bubble peak', ...
       'Location', 'best');

grid on; box on;

%% save combined buoyancy summary

datasetLabel = string({allResults.label})';
bubbleFlux_ml_s_all = [allResults.bubbleFlux_ml_s]';
salinity_brix_all = [allResults.salinity_brix]';
hasBubbles_all_col = [allResults.hasBubbles]';
rho_ambient_all = [allResults.rho_ambient]';
g_reduced_all = [allResults.g_reduced]';
B_freshwater_all = [allResults.B_freshwater]';
B_bubbles_all = [allResults.B_bubbles]';
B_total_all_col = [allResults.B_total]';
meanEndTime_s_all = [allResults.meanEndTime_s]';
meanE_all_col = [allResults.meanEntrainmentRate_ml_s]';
peakE_all_col = [allResults.peakEntrainmentRate_ml_s]';
meanUe_all_col = [allResults.meanUe_eff_cm_s]';
peakUe_all_col = [allResults.peakUe_eff_cm_s]';

combined_buoyancy_summary = table( ...
    datasetLabel, ...
    salinity_brix_all, ...
    hasBubbles_all_col, ...
    bubbleFlux_ml_s_all, ...
    rho_ambient_all, ...
    g_reduced_all, ...
    B_freshwater_all, ...
    B_bubbles_all, ...
    B_total_all_col, ...
    meanEndTime_s_all, ...
    meanE_all_col, ...
    peakE_all_col, ...
    meanUe_all_col, ...
    peakUe_all_col);

writetable(combined_buoyancy_summary, 'combined_buoyancy_flux_summary.csv');