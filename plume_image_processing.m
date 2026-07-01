clear; clc; close all;

%% folders to analyse

datasets = {
    'plume',          'No bubbles';
    'plume_46V_bubbles', 'Bubbles Flux = 4.42ml/s (4.6V)';
    'plume_63V_bubbles', 'Bubbles Flux = 6.50ml/s (6.3V)'; 
    'plume_73V_bubbles', 'Bubble Flux = 8.32ml/s (7.3V)';
};

imagePattern = 'DSC_*.JPG';

dt = 2;
pumpRate_ml_s = 2;
tankInitialVolume_ml = 8000;

%% Bdiff threshold range

Bdiff_threshold_low = 0.3;
Bdiff_threshold_high = 4;

%% strip region after y-cropping

strip_x1 = 3100;
strip_x2 = 3800;

tank_depth_cm = 22;
cm_per_pixel = 0.0111111;
cm2_per_pixel2 = cm_per_pixel^2;

crop_top = 780;
crop_bottom = 300;

smooth_window = 5;

%% image processing parameters

gaussian_sigma = 2;
min_object_size_px = 80;
closing_disk_size = 4;

%% storage for all datasets

allResults = struct();

%% process each folder

for d = 1:size(datasets,1)

    imageFolder = datasets{d,1};
    labelName = datasets{d,2};

    fprintf('\nProcessing folder: %s\n', imageFolder);

    overlayFolder = fullfile(imageFolder, 'threshold_overlays');
    if ~exist(overlayFolder, 'dir')
        mkdir(overlayFolder);
    end

    diagnosticFolder = fullfile(imageFolder, 'diagnostic_deltaB');
    if ~exist(diagnosticFolder, 'dir')
        mkdir(diagnosticFolder);
    end

    files = dir(fullfile(imageFolder, imagePattern));
    [~, idx] = sort({files.name});
    files = files(idx);

    nFrames = length(files);

    if nFrames == 0
        error('No images found in folder: %s', imageFolder);
    end

    %% reference image

    I0 = im2double(imread(fullfile(imageFolder, files(1).name)));

    if crop_bottom > 0
        I0 = I0(crop_top+1:end-crop_bottom, :, :);
    else
        I0 = I0(crop_top+1:end, :, :);
    end

    B0 = I0(:,:,3);

    waterMask = true(size(B0));

    %% check strip limits

    imageWidth = size(B0, 2);

    if strip_x2 > imageWidth
        warning('strip_x2 exceeds image width. Changing strip_x2 from %d to %d.', ...
            strip_x2, imageWidth);
        strip_x2 = imageWidth;
    end

    if strip_x1 < 1
        strip_x1 = 1;
    end

    if strip_x1 >= strip_x2
        error('Invalid strip range. Check strip_x1 and strip_x2.');
    end

    stripMask = false(size(B0));
    stripMask(:, strip_x1:strip_x2) = true;

    %% effective tank depth

    tankArea_px = nnz(waterMask);
    tankArea_cm2 = tankArea_px * cm2_per_pixel2;

    effective_depth_cm = tankInitialVolume_ml / tankArea_cm2;

    %% allocate variables

    time_s = (0:nFrames-1)' * dt;
    pumpedVolume_ml = time_s * pumpRate_ml_s;

    plumeArea_px = zeros(nFrames,1);
    plumeArea_cm2 = zeros(nFrames,1);
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

        %% background-subtracted blue-channel detection

        deltaB = B0 - B;

        deltaB(isnan(deltaB)) = 0;
        deltaB(isinf(deltaB)) = 0;
        deltaB(deltaB < 0) = 0;

        deltaB_smooth = imgaussfilt(deltaB, gaussian_sigma);

        deltaB_enhanced = adapthisteq(mat2gray(deltaB_smooth));

        %% dual-threshold detection for full plume

        dyeMask = ...
            (deltaB_enhanced > Bdiff_threshold_low) & ...
            (deltaB_enhanced < Bdiff_threshold_high);

        dyeMask = dyeMask & waterMask;

        dyeMask = bwareaopen(dyeMask, min_object_size_px);
        dyeMask = imclose(dyeMask, strel('disk', closing_disk_size));
        dyeMask = imfill(dyeMask, 'holes');

        %% strip percentage calculation

        stripRangeMask = ...
            (deltaB_enhanced > Bdiff_threshold_low) & ...
            (deltaB_enhanced < Bdiff_threshold_high) & ...
            stripMask;

        stripBdiffPixels_px(i) = nnz(stripRangeMask);
        stripTotalPixels_px(i) = nnz(stripMask);

        stripBdiffPercentage(i) = ...
            100 * stripBdiffPixels_px(i) / stripTotalPixels_px(i);

        %% save overlay image

        overlayImg = I;
        
        redLayer = overlayImg(:,:,1);
        greenLayer = overlayImg(:,:,2);
        blueLayer = overlayImg(:,:,3);
        
        alpha = 0.45;
        
        % Cyan overlay: RGB = [0, 1, 1]
        redLayer(dyeMask)   = (1-alpha)*redLayer(dyeMask);
        greenLayer(dyeMask) = (1-alpha)*greenLayer(dyeMask) + alpha;
        blueLayer(dyeMask)  = (1-alpha)*blueLayer(dyeMask)  + alpha;
        
        overlayImg(:,:,1) = redLayer;
        overlayImg(:,:,2) = greenLayer;
        overlayImg(:,:,3) = blueLayer;
        
        % Cyan boundary
        boundaryMask = bwperim(dyeMask);
        
        overlayImg(:,:,1) = overlayImg(:,:,1) .* ~boundaryMask;
        overlayImg(:,:,2) = overlayImg(:,:,2) + boundaryMask;
        overlayImg(:,:,3) = overlayImg(:,:,3) + boundaryMask;
        
        %% draw strip boundary in green
        
        overlayImg(:, strip_x1, 1) = 0;
        overlayImg(:, strip_x1, 2) = 1;
        overlayImg(:, strip_x1, 3) = 0;
        
        overlayImg(:, strip_x2, 1) = 0;
        overlayImg(:, strip_x2, 2) = 1;
        overlayImg(:, strip_x2, 3) = 0;
        
        overlayImg = min(overlayImg, 1);
        
        [~, baseName, ~] = fileparts(files(i).name);
        outName = fullfile(overlayFolder, [baseName '_threshold_overlay.png']);
        
        imwrite(overlayImg, outName);

        %% save diagnostic image

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

        fprintf(['Frame %3d/%3d  %s  time = %.1f s  thresholded = %.2f ml  ' ...
                 'pumped = %.2f ml  entrained = %.2f ml  strip = %.2f %%\n'], ...
            i, nFrames, files(i).name, time_s(i), thresholdedVolume_ml(i), ...
            pumpedVolume_ml(i), entrainedVolume_ml(i), stripBdiffPercentage(i));

    end

    %% entrainment rate

    entrained_smoothed = movmean(entrainedVolume_ml, smooth_window);

    entrainmentRate_ml_s = diff(entrained_smoothed) / dt;
    time_rate_s = time_s(1:end-1) + dt/2;

    entrainmentRate_ml_s(entrainmentRate_ml_s < 0) = 0;

    %% store results

    allResults(d).folder = imageFolder;
    allResults(d).label = labelName;
    allResults(d).time_s = time_s;
    allResults(d).pumpedVolume_ml = pumpedVolume_ml;
    allResults(d).thresholdedVolume_ml = thresholdedVolume_ml;
    allResults(d).entrainedVolume_ml = entrainedVolume_ml;
    allResults(d).entrained_smoothed = entrained_smoothed;
    allResults(d).time_rate_s = time_rate_s;
    allResults(d).entrainmentRate_ml_s = entrainmentRate_ml_s;
    allResults(d).plumeArea_px = plumeArea_px;
    allResults(d).plumeArea_cm2 = plumeArea_cm2;

    allResults(d).stripBdiffPixels_px = stripBdiffPixels_px;
    allResults(d).stripTotalPixels_px = stripTotalPixels_px;
    allResults(d).stripBdiffPercentage = stripBdiffPercentage;

    %% save individual csv files

    results = table(time_s, pumpedVolume_ml, thresholdedVolume_ml, ...
        entrainedVolume_ml, entrained_smoothed, plumeArea_px, plumeArea_cm2, ...
        stripBdiffPixels_px, stripTotalPixels_px, stripBdiffPercentage);

    writetable(results, fullfile(imageFolder, 'entrained_volume_results.csv'));

    rate_results = table(time_rate_s, entrainmentRate_ml_s);
    writetable(rate_results, fullfile(imageFolder, 'entrainment_rate_results.csv'));

    strip_results = table(time_s, stripBdiffPixels_px, stripTotalPixels_px, ...
        stripBdiffPercentage);

    writetable(strip_results, fullfile(imageFolder, 'strip_bdiff_percentage_results.csv'));

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
ylabel('Entrained Volume (ml)');
title('(b) Entrained Volume');

legend(string({allResults.label}), 'Location', 'best');

grid on; box on;

%% plot 3: entrainment rate

nexttile;
hold on;

for d = 1:length(allResults)
    plot(allResults(d).time_rate_s, allResults(d).entrainmentRate_ml_s, ...
        'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
end

xlabel('Time (s)');
ylabel('Entrainment Rate (ml/s)');
title('(c) Entrainment Rate');

legend(string({allResults.label}), 'Location', 'best');

grid on; box on;

sgtitle('Plume Entrainment Analysis', 'FontWeight', 'bold');

%% new plot: strip Bdiff percentage

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

title(sprintf('Strip Filling: x = %d to %d, %.2f < Bdiff < %.2f', ...
    strip_x1, strip_x2, Bdiff_threshold_low, Bdiff_threshold_high));

legend(string({allResults.label}), 'Location', 'best');

grid on; box on;
