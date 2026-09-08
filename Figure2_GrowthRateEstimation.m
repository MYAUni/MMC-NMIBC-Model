% Source-specific estimation of the total tumor-growth rate.
%
% Model:
%   dT/dt = r*T*(1 - T/k)
%
%
% Statistical analysis:
%   * Each untreated/control series is fitted separately.
%   * The first retained positive observation initializes the solution.
%   * r minimizes unweighted squared residuals on the log10 cell-count scale.
%   * r is optimized on the real line; the data determine its sign and value.
%   * Recorded zeros are omitted only because a logistic solution initialized
%     at T=0 remains zero and cannot represent later positive observations.
%   * RMSE and R^2 are reported only for series with at least three retained
%     time points. Two-point series still provide an r estimate.
%   * The primary fit fixes k at the exact largest converted literature value.
%   * Fixed-k sensitivity uses the exact prespecified k interval.
%   * Joint r-k fitting is a diagnostic. k is restricted only to the
%     prespecified interval; allowed endpoint solutions are retained.
%   * Two-point series are not jointly fitted because, after fixing the first
%     observation as the initial condition, one remaining observation cannot
%     identify both r and k.
%
% Uncertainty:
% Source-reported measurement bounds are plotted without assuming a common
% measurement-error distribution. Confidence intervals for fitted parameters
% are not calculated.

clear; clc; close all;

%========================================================================
% DISPLAY SETTINGS 
% ========================================================================

set(groot, 'defaultAxesFontSize', 15);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.08);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.05);
set(groot, 'defaultAxesLineWidth', 1.35);
set(groot, 'defaultLineLineWidth', 2.5);
set(groot, 'defaultTextFontSize', 15);
set(groot, 'defaultLegendFontSize', 13);
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultLegendFontName', 'Arial');
set(groot, 'defaultLegendInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');

% Figure 1 font sizes.
Fig1Font.tick = 14.8;
Fig1Font.panelTitle = 15.6;
Fig1Font.panelLetter = 20;
Fig1Font.axisLabel = 20;
Fig1Font.legend = 20.0;

C = getPlotColors();


% Panel order.
plotOrder = [ ...
    "DS158_MRI_MBT2"; ...
    "DS163_UMUC3_noMMC"; ...
    "DS164_UMUC3_PBS"; ...
    "DS160_HRUS_MB49"; ...
    "DS162_Mouse01_vehicle"; ...
    "DS162_Mouse02_vehicle"; ...
    "DS162_Mouse03_vehicle"; ...
    "DS162_Mouse04_vehicle"];

%% ========================================================================
% ANALYSIS SETTINGS
% ========================================================================

cellsPerMM3 = 1e6;

kMain = 288737497606.131;          % cells
kMin  = 9.00e7;                    % cells
kMax  = 288737497606.131;          % cells

rLiterature = [6.25e-3, 0.5];      % day^-1; comparison/simulation interval

% Grid sizes for fixed-k sensitivity and joint r-k profiling.
nKSensitivity = 9;
nKProfile     = 201;

kSensitivity = logspace(log10(kMin), log10(kMax), nKSensitivity);
kSensitivity([1 end]) = [kMin kMax];

assert(kMain == kMax, 'The primary k must equal the upper k diagnostic bound.');
assert(kMin == 9.00e7, 'The lower k diagnostic bound has changed.');

%% ========================================================================
% DATA ENTRY
% ========================================================================

S = struct( ...
    'id', {}, 'shortLabel', {}, 'group', {}, ...
    'timeDays', {}, 'cells', {}, 'lowerCells', {}, 'upperCells', {}, ...
    'useForFit', {}, 'notes', {});

S(end+1) = makeSeries( ...
    "DS158_MRI_MBT2", "MBT-2 MRI", "untreated_control", ...
    [10 14 17 24], ...
    [7.20e4 3.24e5 4.95e5 9.54e5], ...
    [5.70e4 1.77e5 4.35e5 8.85e5], ...
    [8.70e4 4.71e5 5.55e5 1.02e6], ...
    true, "Converted MRI tumor-area measurements.");

S(end+1) = makeSeries( ...
    "DS163_UMUC3_noMMC", "UMUC3 untreated", "untreated_control", ...
    [0 30], ...
    [6240 1.31e7], ...
    [2456 5.16e6], ...
    [10000 2.10e7], ...
    true, "Converted bioluminescence measurements; two time points.");

S(end+1) = makeSeries( ...
    "DS164_UMUC3_PBS", "UMUC3 PBS", "untreated_control", ...
    [4 21], ...
    [8.36e2 1.68e5], ...
    [4.48e2 1.04e5], ...
    [1.22e3 2.31e5], ...
    true, "Converted bioluminescence measurements; two time points.");

S(end+1) = makeSeries( ...
    "DS160_HRUS_MB49", "MB49 HRUS", "untreated_control", ...
    [0 11 32], ...
    [3.00e3 0.4*cellsPerMM3 308.6*cellsPerMM3], ...
    [], [], ...
    true, "Reported inoculum followed by converted HRUS tumor volumes.");

S(end+1) = makeSeries( ...
    "DS162_Mouse01_vehicle", "MB49 vehicle 1", "vehicle_control", ...
    0:6, ...
    [0 28720 47981 510020 3513702 10118921 18899380], ...
    [], [], true, "Individual vehicle-control mouse.");

S(end+1) = makeSeries( ...
    "DS162_Mouse02_vehicle", "MB49 vehicle 2", "vehicle_control", ...
    0:6, ...
    [5012 28718 58773 412388 2244575 8273819 15008719], ...
    [], [], true, "Individual vehicle-control mouse.");

S(end+1) = makeSeries( ...
    "DS162_Mouse03_vehicle", "MB49 vehicle 3", "vehicle_control", ...
    0:6, ...
    [0 13091 49840 482910 1872512 6998162 13298610], ...
    [], [], true, "Individual vehicle-control mouse.");

S(end+1) = makeSeries( ...
    "DS162_Mouse04_vehicle", "MB49 vehicle 4", "vehicle_control", ...
    0:6, ...
    [0 29102 54972 719261 2569162 7397218 16009271], ...
    [], [], true, "Individual vehicle-control mouse.");

% Converted treated series retained in the audit but excluded from baseline
% untreated/control growth estimation because it is MMC-treated.
S(end+1) = makeSeries( ...
    "DS169_UMUC3_MMC_excluded", "UMUC3 MMC-treated", "treated_excluded", ...
    [4 21], ...
    [6.44e2 1.55e3], ...
    [4.96e2 9.24e2], ...
    [7.92e2 2.18e3], ...
    false, "Excluded from baseline growth estimation because it is treated.");

seriesIDVector = vertcat(S.id);
if numel(unique(seriesIDVector)) ~= numel(seriesIDVector)
    error('Series identifiers must be unique.');
end

%% ========================================================================
% CLEANING AUDIT AND PRIMARY FIXED-k FITS
% ========================================================================

nSeries = numel(S);
Clean = cell(nSeries, 1);

fitRows = struct([]);
auditRows = struct([]);

for i = 1:nSeries
    [cleanSeries, audit] = prepareSeries(S(i));
    Clean{i} = cleanSeries;
    auditRows = appendStruct(auditRows, audit);

    if ~S(i).useForFit
        continue;
    end

    fit = fitFixedK(Clean{i}.timeDays, Clean{i}.cells, kMain);
    goodnessAssessed = Clean{i}.nPositive >= 3;

    if goodnessAssessed
        rmseReported = fit.RMSElog10;
        r2Reported = fit.R2log10;
    else
        rmseReported = NaN;
        r2Reported = NaN;
    end

    fitRows = appendStruct(fitRows, struct( ...
        'SeriesID', S(i).id, ...
        'Label', S(i).shortLabel, ...
        'Group', string(S(i).group), ...
        'k_cells', kMain, ...
        'r_per_day', fit.r, ...
        'SSE_log10', fit.SSElog10, ...
        'RMSE_log10_raw', fit.RMSElog10, ...
        'RMSE_log10_reported', rmseReported, ...
        'R2_log10_raw', fit.R2log10, ...
        'R2_log10_reported', r2Reported, ...
        'N_raw_points', Clean{i}.nRaw, ...
        'N_positive_points', Clean{i}.nPositive, ...
        'N_recorded_zeros_omitted', Clean{i}.nZeroOmitted, ...
        'Initial_time_day', Clean{i}.timeDays(1), ...
        'Initial_burden_cells', Clean{i}.cells(1), ...
        'Goodness_of_fit_assessed', goodnessAssessed, ...
        'Optimizer_exitflag', fit.exitflag, ...
        'Maximum_start_solution_difference', fit.startSpread));
end

AuditTable = struct2table(auditRows);
PrimaryFitTable = struct2table(fitRows);

if height(AuditTable) ~= nSeries
    error('The cleaning audit does not contain one row per entered series.');
end
if height(PrimaryFitTable) ~= sum([S.useForFit])
    error('The primary-fit table does not contain one row per fitted series.');
end
if any(~isfinite(PrimaryFitTable.r_per_day)) || ...
        any(~isfinite(PrimaryFitTable.SSE_log10))
    error('The primary-fit table contains a non-finite fitted result.');
end

writetable(AuditTable, 'GrowthRateEstimation_data_and_cleaning_audit.csv');
writetable(PrimaryFitTable, 'GrowthRateEstimation_primary_fixed_k_fits.csv');

disp('=== Data and cleaning audit ===');
disp(AuditTable);
disp('=== Primary fixed-k fits ===');
disp(PrimaryFitTable);

%% ========================================================================
% DESCRIPTIVE RANGES
% ========================================================================

isUntreated = PrimaryFitTable.Group == "untreated_control";
isVehicle   = PrimaryFitTable.Group == "vehicle_control";

if any(~(isUntreated | isVehicle))
    error('An unexpected group is present in the primary-fit table.');
end
if ~any(isUntreated) || ~any(isVehicle)
    error('Both untreated/control and vehicle-control estimates are required.');
end

RangeTable = table( ...
    ["untreated_control"; "vehicle_control"; "all_fitted_estimates"], ...
    [sum(isUntreated); sum(isVehicle); height(PrimaryFitTable)], ...
    [min(PrimaryFitTable.r_per_day(isUntreated)); ...
     min(PrimaryFitTable.r_per_day(isVehicle)); ...
     min(PrimaryFitTable.r_per_day)], ...
    [median(PrimaryFitTable.r_per_day(isUntreated)); ...
     median(PrimaryFitTable.r_per_day(isVehicle)); ...
     median(PrimaryFitTable.r_per_day)], ...
    [max(PrimaryFitTable.r_per_day(isUntreated)); ...
     max(PrimaryFitTable.r_per_day(isVehicle)); ...
     max(PrimaryFitTable.r_per_day)], ...
    'VariableNames', {'EstimateGroup','N_estimates','r_min_day','r_median_day','r_max_day'});

writetable(RangeTable, 'GrowthRateEstimation_descriptive_ranges.csv');

LiteratureTable = table( ...
    rLiterature(1), rLiterature(2), ...
    'VariableNames', {'r_min_day','r_max_day'});
writetable(LiteratureTable, 'GrowthRateEstimation_literature_interval.csv');

%% ========================================================================
% FIXED-k SENSITIVITY
% ========================================================================

sensitivityRows = struct([]);

for i = 1:nSeries
    if ~S(i).useForFit
        continue;
    end

    for j = 1:numel(kSensitivity)
        fit = fitFixedK(Clean{i}.timeDays, Clean{i}.cells, kSensitivity(j));

        sensitivityRows = appendStruct(sensitivityRows, struct( ...
            'SeriesID', S(i).id, ...
            'Label', S(i).shortLabel, ...
            'Group', string(S(i).group), ...
            'k_cells', kSensitivity(j), ...
            'r_per_day', fit.r, ...
            'SSE_log10', fit.SSElog10, ...
            'RMSE_log10_raw', fit.RMSElog10, ...
            'RMSE_log10_reported', ...
                conditionalValue(Clean{i}.nPositive >= 3, fit.RMSElog10, NaN), ...
            'N_positive_points', Clean{i}.nPositive, ...
            'Goodness_of_fit_assessed', Clean{i}.nPositive >= 3, ...
            'Optimizer_exitflag', fit.exitflag, ...
            'Maximum_start_solution_difference', fit.startSpread));
    end
end

SensitivityTable = struct2table(sensitivityRows);
expectedSensitivityRows = sum([S.useForFit]) * numel(kSensitivity);
if height(SensitivityTable) ~= expectedSensitivityRows
    error('The fixed-k sensitivity table is incomplete.');
end
if any(~isfinite(SensitivityTable.r_per_day)) || ...
        any(~isfinite(SensitivityTable.SSE_log10))
    error('The fixed-k sensitivity table contains a non-finite fitted result.');
end
writetable(SensitivityTable, 'GrowthRateEstimation_fixed_k_sensitivity.csv');

%% ========================================================================
% JOINT r-k DIAGNOSTIC
% ========================================================================

jointRows = struct([]);

for i = 1:nSeries
    if ~S(i).useForFit
        continue;
    end

    if Clean{i}.nPositive < 3
        jointRows = appendStruct(jointRows, struct( ...
            'SeriesID', S(i).id, ...
            'Label', S(i).shortLabel, ...
            'Group', string(S(i).group), ...
            'N_positive_points', Clean{i}.nPositive, ...
            'r_fixed_k_per_day', lookupPrimaryR(PrimaryFitTable, S(i).id), ...
            'r_joint_per_day', NaN, ...
            'k_joint_cells', NaN, ...
            'SSE_joint_log10', NaN, ...
            'k_at_lower_bound', false, ...
            'k_at_upper_bound', false, ...
            'Diagnostic_status', "not_estimable_from_two_points"));
        continue;
    end

    joint = fitJointRKProfile( ...
        Clean{i}.timeDays, Clean{i}.cells, kMin, kMax, nKProfile);

    jointRows = appendStruct(jointRows, struct( ...
        'SeriesID', S(i).id, ...
        'Label', S(i).shortLabel, ...
        'Group', string(S(i).group), ...
        'N_positive_points', Clean{i}.nPositive, ...
        'r_fixed_k_per_day', lookupPrimaryR(PrimaryFitTable, S(i).id), ...
        'r_joint_per_day', joint.r, ...
        'k_joint_cells', joint.k, ...
        'SSE_joint_log10', joint.SSElog10, ...
        'k_at_lower_bound', joint.atLowerBound, ...
        'k_at_upper_bound', joint.atUpperBound, ...
        'Diagnostic_status', "estimated"));
end

JointTable = struct2table(jointRows);
if height(JointTable) ~= sum([S.useForFit])
    error('The joint-fit diagnostic table is incomplete.');
end
estimatedJointRows = JointTable.Diagnostic_status == "estimated";
if any(~isfinite(JointTable.r_joint_per_day(estimatedJointRows))) || ...
        any(~isfinite(JointTable.k_joint_cells(estimatedJointRows))) || ...
        any(~isfinite(JointTable.SSE_joint_log10(estimatedJointRows)))
    error('An estimated joint-fit row contains a non-finite result.');
end
writetable(JointTable, 'GrowthRateEstimation_joint_r_k_diagnostic.csv');

disp('=== Joint r-k diagnostic ===');
disp(JointTable);

%% ========================================================================
% MACHINE-READABLE ANALYSIS SETTINGS
% ========================================================================

SettingsTable = table( ...
    ["model"; "objective"; "initial_condition"; "r_domain"; ...
     "primary_k_cells"; "k_diagnostic_min_cells"; "k_diagnostic_max_cells"; ...
     "literature_r_min_day"; "literature_r_max_day"; ...
     "two_point_primary_fit"; "two_point_goodness_of_fit"; ...
     "two_point_joint_r_k"; "uncertainty_intervals"; ...
     "n_k_sensitivity_values"; "n_k_profile_values"], ...
    ["dT/dt = r*T*(1-T/k)"; ...
     "unweighted sum of squared log10 cell-count residuals"; ...
     "first retained positive observation"; ...
     "unbounded real line"; ...
     string(sprintf('%.15g', kMain)); ...
     string(sprintf('%.15g', kMin)); ...
     string(sprintf('%.15g', kMax)); ...
     string(sprintf('%.15g', rLiterature(1))); ...
     string(sprintf('%.15g', rLiterature(2))); ...
     "included"; "not assessed"; "not estimable"; ...
     "not calculated without a common measurement-error model"; ...
     string(nKSensitivity); string(nKProfile)], ...
    'VariableNames', {'Setting','Value'});

writetable(SettingsTable, 'GrowthRateEstimation_analysis_settings.csv');

%% ========================================================================
% MAIN FIGURE: ALL SOURCE-SPECIFIC FITS
% ========================================================================

fitIndices = find([S.useForFit]);
PlotTable = reorderTableByIDs(PrimaryFitTable, plotOrder);

if ~isempty(PlotTable)

    Fig1 = figure('Color','w', ...
        'Name','All source total tumor-burden fits');
    set(Fig1, 'Units','inches', 'Position',[1 1 14.6 9.6]);

    TL = tiledlayout(Fig1, 2, 4, ...
        'TileSpacing','compact', 'Padding','loose');
    try
        TL.Units = 'normalized';
        TL.Position = [0.070 0.205 0.905 0.705];
    catch
    end

    panelLetters = char('A' + (0:25));

    for p = 1:height(PlotTable)
        idx = find(seriesIDVector == PlotTable.SeriesID(p), 1);
        if isempty(idx)
            error('No source definition was found for %s.', ...
                PlotTable.SeriesID(p));
        end

        t = Clean{idx}.timeDays;
        T = Clean{idx}.cells;
        L = Clean{idx}.lowerCells;
        U = Clean{idx}.upperCells;

        tf = linspace(min(t), max(t), 300);
        centralFit = logisticSolution( ...
            tf, t(1), T(1), PlotTable.r_per_day(p), kMain);

        ax = nexttile(TL);
        hold(ax, 'on');

        plot(ax, tf, centralFit, '-', ...
            'Color',C.fit, ...
            'LineWidth',2.7);

        if Clean{idx}.hasBounds
            errorbar(ax, t, T, T-L, U-T, 'o', ...
                'Color',C.data, ...
                'MarkerEdgeColor',C.data, ...
                'MarkerFaceColor','w', ...
                'LineWidth',1.55, ...
                'MarkerSize',6.8, ...
                'CapSize',5);
        else
            plot(ax, t, T, 'o', ...
                'Color',C.data, ...
                'MarkerEdgeColor',C.data, ...
                'MarkerFaceColor','w', ...
                'LineWidth',1.55, ...
                'MarkerSize',6.8);
        end

        set(ax, ...
            'YScale','log', ...
            'FontSize',Fig1Font.tick, ...
            'FontWeight','bold', ...
            'LineWidth',1.20, ...
            'TickDir','out', ...
            'Box','on');

        grid(ax, 'on');
        ax.GridAlpha = 0.14;
        ax.XMinorGrid = 'off';
        ax.YMinorGrid = 'off';
        ax.XMinorTick = 'off';
        ax.YMinorTick = 'off';

        title(ax, S(idx).shortLabel, ...
            'Interpreter','tex', ...
            'FontSize',Fig1Font.panelTitle, ...
            'FontWeight','bold');

        text(ax, 0.03, 0.92, panelLetters(p), ...
            'Units','normalized', ...
            'FontSize',Fig1Font.panelLetter, ...
            'FontWeight','bold', ...
            'BackgroundColor','w', ...
            'Margin',1.2);
    end

    xlabel(TL, 'Time (days)', ...
        'FontSize',Fig1Font.axisLabel, ...
        'FontWeight','bold');
    ylabel(TL, 'Tumor burden (cells)', ...
        'FontSize',Fig1Font.axisLabel, ...
        'FontWeight','bold');

    addMainFitLegend(Fig1, C, Fig1Font.legend);

    exportJournalFigure(Fig1, ...
        '01_Main_Figure_Logistic_Fits', 14.6, 9.6);
end

%% ========================================================================
% SUPPLEMENTARY FIGURE S1: SOURCE ESTIMATES AND LITERATURE COMPARISON
% ========================================================================

FigS1 = figure('Color','w', ...
    'Name','Supplementary Figure S1: source estimates and literature comparison');
set(FigS1, 'Units','inches', 'Position',[1 1 13.2 8.8]);

TLS1 = tiledlayout(FigS1, 5, 2, ...
    'TileSpacing','compact', 'Padding','loose');
try
    TLS1.Units = 'normalized';
    TLS1.Position = [0.180 0.105 0.790 0.820];
catch
end

SourcePlotTable = reorderTableByIDs(PrimaryFitTable, plotOrder);

axS1a = nexttile(TLS1, 1, [3 1]);
plotEstimatePanelExactStyle( ...
    axS1a, SourcePlotTable(SourcePlotTable.Group == "untreated_control",:), ...
    'Untreated and control estimates', C);
addPanelLetter(axS1a, 'A', 15);

axS1b = nexttile(TLS1, 2, [3 1]);
plotEstimatePanelExactStyle( ...
    axS1b, SourcePlotTable(SourcePlotTable.Group == "vehicle_control",:), ...
    'MB49 vehicle-control mice', C);
addPanelLetter(axS1b, 'B', 15);

axS1c = nexttile(TLS1, 7, [2 2]);
plotRangeSummaryExactStyle(axS1c, RangeTable, rLiterature, C);
title(axS1c, 'Literature comparison', ...
    'FontSize',13, 'FontWeight','bold');
addPanelLetter(axS1c, 'C', 15);

exportJournalFigure(FigS1, ...
    '02_Supplementary_Figure_S1_source_estimates_and_literature_comparison', ...
    13.2, 8.8);

%% ========================================================================
% SUPPLEMENTARY FIGURE S2: FIXED-k SENSITIVITY
% ========================================================================

availableIDs = unique(SensitivityTable.SeriesID, 'stable');
uniqueIDs = plotOrder(ismember(plotOrder, availableIDs));
lineColors = lines(numel(uniqueIDs));

FigS2 = figure('Color','w', ...
    'Name','Supplementary Figure S2: fixed-k sensitivity');
set(FigS2, 'Units','inches', 'Position',[1 1 13.2 7.2]);

TLS2 = tiledlayout(FigS2, 1, 2, ...
    'TileSpacing','compact', 'Padding','loose');
try
    TLS2.Units = 'normalized';
    TLS2.Position = [0.085 0.300 0.885 0.590];
catch
end

axS2a = nexttile(TLS2, 1);
plotFixedKSensitivityPanelExactStyle( ...
    axS2a, SensitivityTable, uniqueIDs, lineColors, ...
    kSensitivity, kMain, C);
title(axS2a, 'Growth-rate estimates', ...
    'FontSize',14, 'FontWeight','bold');
addPanelLetter(axS2a, 'A');

axS2b = nexttile(TLS2, 2);
plotFixedKRMSEPanelExactStyle( ...
    axS2b, SensitivityTable, uniqueIDs, lineColors, ...
    kSensitivity, kMain, C);
title(axS2b, 'Fit error', ...
    'FontSize',14, 'FontWeight','bold');
addPanelLetter(axS2b, 'B');

addFixedKSensitivityLegendExactStyle( ...
    FigS2, uniqueIDs, lineColors, C);

exportJournalFigure(FigS2, ...
    '03_Supplementary_Figure_S2_fixed_k_sensitivity', ...
    13.2, 7.2);

%% ========================================================================
% SUPPLEMENTARY FIGURE S3: JOINT r-k DIAGNOSTIC
% ========================================================================

jointPlotted = JointTable(JointTable.Diagnostic_status == "estimated", :);
if isempty(jointPlotted)
    error('No joint r-k estimates are available for Supplementary Figure S3.');
end

jointPlotted = reorderTableByIDs(jointPlotted, plotOrder);

FigS3 = figure('Color','w', ...
    'Name','Supplementary Figure S3: fixed-k and joint r-k estimates');
set(FigS3, 'Units','inches', 'Position',[1 1 8.8 5.6]);

axS3 = axes(FigS3, 'Position',[0.115 0.160 0.820 0.760]);
plotFixedKVsJointRKExactStyle(axS3, jointPlotted, C);

exportJournalFigure(FigS3, ...
    '04_Supplementary_Figure_S3_joint_r_k_diagnostic', ...
    8.8, 5.6);

fprintf('\nCompleted without finite r bounds or fit-rejection thresholds.\n');
fprintf('Primary k: %.15g cells\n', kMain);
fprintf('Diagnostic k interval: %.15g to %.15g cells\n', kMin, kMax);
fprintf('No inferential uncertainty interval was calculated.\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function C = getPlotColors()

    C.fit = [0.0000 0.4470 0.7410];
    C.data = [0.8500 0.3250 0.0980];
    C.slow = [0.0000 0.4470 0.7410];
    C.rapid = [0.8500 0.3250 0.0980];
    C.interval = [0.25 0.25 0.25];
    C.range = [0.0000 0.4470 0.7410];
    C.median = [0.8500 0.3250 0.0980];
    C.literature = [0.95 0.86 0.50];
    C.literatureEdge = [0.60 0.42 0.05];
    C.ref = [0.15 0.15 0.15];
    C.diagonal = [0.15 0.15 0.15];
end

function s = makeSeries(id, shortLabel, group, timeDays, cells, ...
    lowerCells, upperCells, useForFit, notes)

    s = struct( ...
        'id', string(id), ...
        'shortLabel', string(shortLabel), ...
        'group', string(group), ...
        'timeDays', double(timeDays(:)), ...
        'cells', double(cells(:)), ...
        'lowerCells', double(lowerCells(:)), ...
        'upperCells', double(upperCells(:)), ...
        'useForFit', logical(useForFit), ...
        'notes', string(notes));
end

function [clean, audit] = prepareSeries(s)

    t = s.timeDays(:);
    y = s.cells(:);

    if numel(t) ~= numel(y)
        error('%s: time and cell-count vectors have different lengths.', s.id);
    end
    if any(~isfinite(t)) || any(~isfinite(y))
        error('%s: missing or non-finite time/cell values require explicit review.', s.id);
    end
    if any(y < 0)
        error('%s: negative cell counts are not permitted.', s.id);
    end
    if any(diff(t) <= 0)
        error('%s: time points must be strictly increasing.', s.id);
    end

    zeroMask = y == 0;
    firstPositiveIndex = find(y > 0, 1, 'first');
    if isempty(firstPositiveIndex)
        if s.useForFit
            error('%s: no positive observations remain for fitting.', s.id);
        end
    elseif any(zeroMask(firstPositiveIndex:end))
        error(['%s: a recorded zero occurs after the first positive ' ...
            'observation and requires explicit review.'], s.id);
    end
    keep = ~zeroMask;

    if s.useForFit && sum(keep) < 2
        error('%s: fewer than two positive observations remain.', s.id);
    end

    clean.timeDays = t(keep);
    clean.cells = y(keep);
    clean.nRaw = numel(y);
    clean.nPositive = sum(keep);
    clean.nZeroOmitted = sum(zeroMask);

    hasLower = ~isempty(s.lowerCells);
    hasUpper = ~isempty(s.upperCells);
    if xor(hasLower, hasUpper)
        error('%s: both lower and upper values must be supplied together.', s.id);
    end

    clean.hasBounds = hasLower && hasUpper;
    if clean.hasBounds
        if numel(s.lowerCells) ~= numel(y) || numel(s.upperCells) ~= numel(y)
            error('%s: bound vectors do not match the observations.', s.id);
        end
        if any(s.lowerCells <= 0) || any(s.upperCells <= 0)
            error('%s: source-reported bounds must be positive for log plotting.', s.id);
        end
        if any(s.lowerCells > y) || any(s.upperCells < y)
            error('%s: central values fall outside source-reported bounds.', s.id);
        end
        clean.lowerCells = s.lowerCells(keep);
        clean.upperCells = s.upperCells(keep);
    else
        clean.lowerCells = [];
        clean.upperCells = [];
    end

    audit = struct( ...
        'SeriesID', s.id, ...
        'Label', s.shortLabel, ...
        'Group', s.group, ...
        'Used_for_baseline_fit', s.useForFit, ...
        'N_raw_points', clean.nRaw, ...
        'N_positive_points', clean.nPositive, ...
        'N_recorded_zeros_omitted', clean.nZeroOmitted, ...
        'Has_source_reported_bounds', clean.hasBounds, ...
        'Notes', s.notes);
end

function fit = fitFixedK(t, y, k)
% Estimate r without a finite parameter bound.

    t = t(:);
    y = y(:);

    if any(y <= 0)
        error('fitFixedK requires positive observations.');
    end
    if y(1) >= k
        error(['The initial burden must be below k for the increasing ' ...
            'logistic solution used here.']);
    end

    objective = @(r) sum(( ...
        logisticLog10(t, t(1), y(1), r, k) - log10(y)).^2);

    elapsed = t(end) - t(1);
    if elapsed <= 0
        error('At least two distinct time points are required.');
    end

    endpointSlope = log(y(end)/y(1)) / elapsed;
    starts = unique([endpointSlope, 0, -abs(endpointSlope), abs(endpointSlope)]);

    options = optimset( ...
        'Display', 'off', ...
        'TolX', 1e-12, ...
        'TolFun', 1e-12, ...
        'MaxIter', 2e4, ...
        'MaxFunEvals', 5e4);

    candidateR = nan(numel(starts),1);
    candidateSSE = inf(numel(starts),1);
    candidateFlag = zeros(numel(starts),1);

    for q = 1:numel(starts)
        [candidateR(q), candidateSSE(q), candidateFlag(q)] = ...
            fminsearch(objective, starts(q), options);
    end

    validSolution = isfinite(candidateR) & isfinite(candidateSSE) & ...
        candidateFlag > 0;
    if ~any(validSolution)
        error('Fixed-k optimization failed for every starting value.');
    end

    candidateSSE(~validSolution) = Inf;
    [~, bestIndex] = min(candidateSSE);
    bestR = candidateR(bestIndex);

    successfulSolutions = candidateR(validSolution);
    startSpread = max(successfulSolutions) - min(successfulSolutions);

    predictedLog10 = logisticLog10(t, t(1), y(1), bestR, k);
    observedLog10 = log10(y);
    residuals = predictedLog10 - observedLog10;
    bestSSE = sum(residuals.^2);

    sst = sum((observedLog10 - mean(observedLog10)).^2);
    if sst > 0
        r2 = 1 - bestSSE/sst;
    else
        r2 = NaN;
    end

    fit = struct( ...
        'r', bestR, ...
        'SSElog10', bestSSE, ...
        'RMSElog10', sqrt(mean(residuals.^2)), ...
        'R2log10', r2, ...
        'exitflag', candidateFlag(bestIndex), ...
        'startSpread', startSpread);
end

function joint = fitJointRKProfile(t, y, kLower, kUpper, nGrid)
% Profile over the prespecified k interval; estimate unbounded r at each k.
% Include both k bounds among the candidate solutions.
    if ~(isfinite(kLower) && isfinite(kUpper) && ...
            kLower > 0 && kUpper > kLower)
        error('The joint-fit k interval must contain two finite positive bounds.');
    end
    if ~(isscalar(nGrid) && isfinite(nGrid) && nGrid >= 3 && ...
            nGrid == floor(nGrid))
        error('The joint-fit profile grid size must be an integer of at least 3.');
    end

    logKGrid = linspace(log(kLower), log(kUpper), nGrid);
    profileSSE = nan(size(logKGrid));

    for j = 1:numel(logKGrid)
        fit = fitFixedK(t, y, exp(logKGrid(j)));
        profileSSE(j) = fit.SSElog10;
    end

    [~, gridBest] = min(profileSSE);

    if gridBest == 1
        refinedLogK = logKGrid(1);
    elseif gridBest == numel(logKGrid)
        refinedLogK = logKGrid(end);
    else
        profileObjective = @(logK) ...
            fixedKObjectiveAtLogK(logK, t, y);
        [refinedLogK, ~, refinementExitflag] = fminbnd(profileObjective, ...
            logKGrid(gridBest-1), logKGrid(gridBest+1), ...
            optimset('Display','off','TolX',1e-11,'MaxIter',1e3));
        if refinementExitflag <= 0 || ~isfinite(refinedLogK)
            error('Joint r-k profile refinement did not converge.');
        end
    end

    candidateLogK = [log(kLower), refinedLogK, log(kUpper)];
    candidateSSE = nan(size(candidateLogK));
    candidateR = nan(size(candidateLogK));

    for j = 1:numel(candidateLogK)
        fit = fitFixedK(t, y, exp(candidateLogK(j)));
        candidateSSE(j) = fit.SSElog10;
        candidateR(j) = fit.r;
    end

    validCandidate = isfinite(candidateSSE) & isfinite(candidateR);
    if ~any(validCandidate)
        error('Joint r-k fitting produced no finite candidate solution.');
    end
    candidateSSE(~validCandidate) = Inf;
    [bestSSE, bestIndex] = min(candidateSSE);
    bestK = exp(candidateLogK(bestIndex));

    relativeTolerance = 1e-8; % Relative tolerance for identifying solutions at a k bound.
    joint = struct( ...
        'r', candidateR(bestIndex), ...
        'k', bestK, ...
        'SSElog10', bestSSE, ...
        'atLowerBound', abs(bestK-kLower)/kLower <= relativeTolerance, ...
        'atUpperBound', abs(bestK-kUpper)/kUpper <= relativeTolerance);
end

function value = fixedKObjectiveAtLogK(logK, t, y)
    fit = fitFixedK(t, y, exp(logK));
    value = fit.SSElog10;
end

function y = logisticSolution(t, t0, y0, r, k)
    log10y = logisticLog10(t, t0, y0, r, k);
    y = 10.^log10y;
end

function log10y = logisticLog10(t, t0, y0, r, k)
% Numerically stable log10 of the analytical logistic solution.

    if y0 <= 0 || y0 >= k
        error('The analytical form requires 0 < initial burden < k.');
    end

    a = log(k/y0 - 1) - r.*(t(:)-t0);
    logDenominator = max(a,0) + log1p(exp(-abs(a)));
    log10y = (log(k) - logDenominator) / log(10);
end

function out = appendStruct(in, row)
    if isempty(in)
        out = row;
    else
        out = in;
        out(end+1) = row;
    end
end

function value = lookupPrimaryR(T, seriesID)
    row = T(T.SeriesID == string(seriesID), :);
    if height(row) ~= 1
        error('Expected one primary-fit row for %s.', seriesID);
    end
    value = row.r_per_day(1);
end

function value = conditionalValue(condition, valueIfTrue, valueIfFalse)
    if condition
        value = valueIfTrue;
    else
        value = valueIfFalse;
    end
end

function Tordered = reorderTableByIDs(T, idOrder)

    Tordered = T([], :);
    if isempty(T)
        return;
    end

    used = false(height(T), 1);
    for i = 1:numel(idOrder)
        rows = T.SeriesID == string(idOrder(i));
        if any(rows)
            Tordered = [Tordered; T(rows,:)]; %#ok<AGROW>
            used = used | rows;
        end
    end
    Tordered = [Tordered; T(~used,:)];
end

function labels = makeCompactSourceLabels(seriesIDs)

    seriesIDs = string(seriesIDs);
    labels = cell(numel(seriesIDs), 1);

    for i = 1:numel(seriesIDs)
        switch seriesIDs(i)
            case "DS158_MRI_MBT2"
                labels{i} = 'MBT-2 MRI';
            case "DS160_HRUS_MB49"
                labels{i} = 'MB49 HRUS';
            case "DS162_Mouse01_vehicle"
                labels{i} = 'MB49 vehicle 1';
            case "DS162_Mouse02_vehicle"
                labels{i} = 'MB49 vehicle 2';
            case "DS162_Mouse03_vehicle"
                labels{i} = 'MB49 vehicle 3';
            case "DS162_Mouse04_vehicle"
                labels{i} = 'MB49 vehicle 4';
            case "DS163_UMUC3_noMMC"
                labels{i} = 'UMUC3 untreated';
            case "DS164_UMUC3_PBS"
                labels{i} = 'UMUC3 PBS';
            otherwise
                labels{i} = char(seriesIDs(i));
        end
    end
end

function addMainFitLegend(figHandle, C, legendFontSize)

    axLeg = axes(figHandle, ...
        'Position',[0.260 0.035 0.480 0.060], ...
        'Visible','off');
    hold(axLeg, 'on');

    hFit = plot(axLeg, NaN, NaN, '-', ...
        'Color',C.fit, ...
        'LineWidth',2.2);
    hData = plot(axLeg, NaN, NaN, 'o', ...
        'Color',C.data, ...
        'MarkerEdgeColor',C.data, ...
        'MarkerFaceColor','w', ...
        'LineWidth',1.35, ...
        'MarkerSize',5.8);

    legend(axLeg, [hFit hData], ...
        {'logistic fit', 'converted data'}, ...
        'Orientation','horizontal', ...
        'Location','south', ...
        'FontSize',legendFontSize, ...
        'FontWeight','bold', ...
        'Box','on');
end

function plotEstimatePanelExactStyle(ax, T, panelTitle, C)

    hold(ax, 'on');
    if isempty(T)
        text(ax, 0.5, 0.5, 'No estimates available', ...
            'Units','normalized', 'HorizontalAlignment','center');
        axis(ax, 'off');
        return;
    end

    y = 1:height(T);
    for i = 1:height(T)
        if T.Group(i) == "vehicle_control"
            rowColor = C.rapid;
        else
            rowColor = C.slow;
        end

        rowMarker = 'o';
        if T.N_positive_points(i) == 2
            rowMarker = 's';
        end

        plot(ax, T.r_per_day(i), y(i), rowMarker, ...
            'Color',rowColor, ...
            'MarkerEdgeColor',rowColor, ...
            'MarkerFaceColor','w', ...
            'MarkerSize',9.0, ...
            'LineWidth',1.65);
    end

    yticks(ax, y);
    yticklabels(ax, makeCompactSourceLabels(T.SeriesID));
    ylim(ax, [0.5, height(T)+0.5]);

    set(ax, ...
        'TickLabelInterpreter','tex', ...
        'FontSize',11.5, ...
        'FontWeight','normal', ...
        'LineWidth',1.15, ...
        'TickDir','out', ...
        'XScale','linear', ...
        'YDir','reverse', ...
        'Box','on');

    xlabel(ax, 'Estimated growth rate, r (day^{-1})', ...
        'FontSize',12.0, 'FontWeight','normal');
    title(ax, panelTitle, ...
        'FontSize',13.0, 'FontWeight','bold');

    grid(ax, 'on');
    ax.GridAlpha = 0.16;
    ax.YMinorGrid = 'off';
    ax.MinorGridAlpha = 0.06;

    if any(T.Group == "vehicle_control")
        xlim(ax, [1.0 2.5]);
        xticks(ax, 1.0:0.5:2.5);
    else
        xlim(ax, [0.1 0.5]);
        xticks(ax, 0.1:0.1:0.5);
    end
    ax.XMinorTick = 'on';
    ax.XMinorGrid = 'on';
end

function plotRangeSummaryExactStyle(ax, R, literature, C)

    hold(ax, 'on');
    if isempty(R)
        text(ax, 0.5, 0.5, 'No range summary available', ...
            'Units','normalized', 'HorizontalAlignment','center');
        axis(ax, 'off');
        return;
    end

    R = R(R.N_estimates > 0, :);
    y = 1:height(R);

    litMin = literature(1);
    litMax = literature(2);
    validLit = all(isfinite(literature)) && ...
        litMin > 0 && litMax > litMin;

    hLit = [];
    hRange = [];
    hMedian = [];

    if validLit
        hLit = patch(ax, [litMin litMax litMax litMin], ...
            [0.5 0.5 height(R)+0.5 height(R)+0.5], ...
            C.literature, ...
            'EdgeColor',C.literatureEdge, ...
            'LineWidth',0.95, ...
            'FaceAlpha',0.46);
    end

    for i = 1:height(R)
        if R.EstimateGroup(i) == "untreated_control"
            rowColor = C.slow;
        elseif R.EstimateGroup(i) == "vehicle_control"
            rowColor = C.rapid;
        else
            rowColor = C.interval;
        end

        hR = plot(ax, [R.r_min_day(i), R.r_max_day(i)], ...
            [y(i), y(i)], '-', ...
            'Color',rowColor, 'LineWidth',3.2);
        hM = plot(ax, R.r_median_day(i), y(i), 'o', ...
            'Color',rowColor, ...
            'MarkerEdgeColor',rowColor, ...
            'MarkerFaceColor','w', ...
            'MarkerSize',9.0, ...
            'LineWidth',1.65);

        if isempty(hRange)
            hRange = hR;
            hMedian = hM;
        end
    end

    labels = cell(height(R),1);
    for i = 1:height(R)
        switch R.EstimateGroup(i)
            case "untreated_control"
                labels{i} = 'Untreated/control estimates';
            case "vehicle_control"
                labels{i} = 'MB49 vehicle-control mice';
            otherwise
                labels{i} = 'All fitted estimates';
        end
    end

    yticks(ax, y);
    yticklabels(ax, labels);
    ylim(ax, [0.5, height(R)+0.5]);

    set(ax, ...
        'TickLabelInterpreter','tex', ...
        'FontSize',11.5, ...
        'FontWeight','normal', ...
        'LineWidth',1.10, ...
        'TickDir','out', ...
        'XScale','log', ...
        'YDir','reverse', ...
        'XAxisLocation','bottom', ...
        'Box','off');

    xlabel(ax, 'Estimated growth rate, r (day^{-1})', ...
        'FontSize',12.0, ...
        'FontWeight','normal');

    vals = [R.r_min_day; R.r_max_day; R.r_median_day; ...
        litMin; litMax];
    setGrowthRateXAxis(ax, vals);

    grid(ax, 'on');
    ax.GridAlpha = 0.16;
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'off';
    ax.MinorGridAlpha = 0.05;
    ax.XAxisLocation = 'bottom';

    addUntickedTopAndRightBorders(ax);

    if validLit
        yl = ylim(ax);
        yLabel = yl(2) - 0.04*(yl(2)-yl(1));

        text(ax, litMin*1.08, yLabel, sprintf('%.5g', litMin), ...
            'FontSize',10.0, ...
            'FontWeight','bold', ...
            'Color',C.literatureEdge, ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom', ...
            'BackgroundColor','w', ...
            'Margin',0.8);

        text(ax, litMax/1.08, yLabel, sprintf('%.3g', litMax), ...
            'FontSize',10.0, ...
            'FontWeight','bold', ...
            'Color',C.literatureEdge, ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','bottom', ...
            'BackgroundColor','w', ...
            'Margin',0.8);

        legend(ax, [hLit hRange hMedian], ...
            {'literature-supported interval', ...
             'minimum to maximum', ...
             'median'}, ...
            'Location','northwest', ...
            'FontSize',9.5, ...
            'FontWeight','normal', ...
            'Box','on');
    end
end

function addUntickedTopAndRightBorders(ax)
    xl = xlim(ax);
    yl = ylim(ax);

    line(ax, xl, [yl(1) yl(1)], ...
        'Color', ax.XColor, ...
        'LineWidth', ax.LineWidth, ...
        'HandleVisibility', 'off', ...
        'Clipping', 'off');

    line(ax, [xl(2) xl(2)], yl, ...
        'Color', ax.YColor, ...
        'LineWidth', ax.LineWidth, ...
        'HandleVisibility', 'off', ...
        'Clipping', 'off');
end

function plotFixedKSensitivityPanelExactStyle( ...
    ax, FitTable, uniqueIDs, lineColors, kGrid, kRef, C)

    hold(ax, 'on');
    for i = 1:numel(uniqueIDs)
        rows = FitTable.SeriesID == uniqueIDs(i) & ...
            isfinite(FitTable.r_per_day);
        if ~any(rows)
            continue;
        end

        T = sortrows(FitTable(rows,:), 'k_cells');
        plot(ax, T.k_cells, T.r_per_day, '-o', ...
            'Color',lineColors(i,:), ...
            'MarkerEdgeColor',lineColors(i,:), ...
            'MarkerFaceColor','w', ...
            'LineWidth',1.3, ...
            'MarkerSize',5.5);
    end

    set(ax, 'XScale','log', ...
        'FontSize',12.5, ...
        'LineWidth',1.15, ...
        'TickDir','out', ...
        'Box','on');
    setCarryingCapacityXAxis(ax, kGrid);

    ylNow = ylim(ax);
    plot(ax, [kRef kRef], ylNow, '--', ...
        'Color',C.ref, 'LineWidth',1.25);
    ylim(ax, ylNow);

    xlabel(ax, 'Fixed carrying capacity, k (cells)', ...
        'FontSize',13.5, 'FontWeight','bold');
    ylabel(ax, 'Estimated growth rate, r (day^{-1})', ...
        'FontSize',13.5, 'FontWeight','bold');

    grid(ax, 'on');
    ax.GridColor = [0.45 0.45 0.45];
    ax.GridAlpha = 0.18;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    ax.Layer = 'top';
end

function plotFixedKRMSEPanelExactStyle( ...
    ax, FitTable, uniqueIDs, lineColors, kGrid, kRef, C)

    hold(ax, 'on');
    for i = 1:numel(uniqueIDs)
        rows = FitTable.SeriesID == uniqueIDs(i) & ...
            isfinite(FitTable.RMSE_log10_reported);
        if ~any(rows)
            continue;
        end

        T = sortrows(FitTable(rows,:), 'k_cells');
        plot(ax, T.k_cells, T.RMSE_log10_reported, '-o', ...
            'Color',lineColors(i,:), ...
            'MarkerEdgeColor',lineColors(i,:), ...
            'MarkerFaceColor','w', ...
            'LineWidth',1.3, ...
            'MarkerSize',5.5);
    end

    set(ax, 'XScale','log', ...
        'FontSize',12.5, ...
        'LineWidth',1.15, ...
        'TickDir','out', ...
        'Box','on');
    setCarryingCapacityXAxis(ax, kGrid);

    ylNow = ylim(ax);
    plot(ax, [kRef kRef], ylNow, '--', ...
        'Color',C.ref, 'LineWidth',1.25);
    ylim(ax, ylNow);

    xlabel(ax, 'Fixed carrying capacity, k (cells)', ...
        'FontSize',13.5, 'FontWeight','bold');
    ylabel(ax, 'RMSE on log_{10} cell-count scale', ...
        'FontSize',13.5, 'FontWeight','bold');

    grid(ax, 'on');
    ax.GridColor = [0.45 0.45 0.45];
    ax.GridAlpha = 0.18;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    ax.Layer = 'top';
end

function addFixedKSensitivityLegendExactStyle( ...
    figHandle, uniqueIDs, lineColors, C)

    axLeg = axes(figHandle, ...
        'Position',[0.100 0.035 0.800 0.180], ...
        'Visible','off');
    hold(axLeg, 'on');

    h = gobjects(numel(uniqueIDs)+1, 1);
    labels = makeCompactSourceLabels(uniqueIDs);

    for i = 1:numel(uniqueIDs)
        h(i) = plot(axLeg, NaN, NaN, '-o', ...
            'Color',lineColors(i,:), ...
            'MarkerEdgeColor',lineColors(i,:), ...
            'MarkerFaceColor','w', ...
            'LineWidth',1.3, ...
            'MarkerSize',5.5);
    end

    h(end) = plot(axLeg, NaN, NaN, '--', ...
        'Color',C.ref, 'LineWidth',1.25);
    labels{end+1} = 'k used for main fits';

    try
        lgd = legend(axLeg, h, labels, ...
            'Orientation','horizontal', ...
            'NumColumns',4, ...
            'Location','south', ...
            'FontSize',10.8, ...
            'Box','on');
    catch
        lgd = legend(axLeg, h, labels, ...
            'Orientation','horizontal', ...
            'Location','south', ...
            'FontSize',10.8, ...
            'Box','on');
    end

    try
        lgd.Units = 'normalized';
        lgd.Position = [0.285 0.035 0.430 0.135];
    catch
    end
end

function plotFixedKVsJointRKExactStyle(ax, JointTable, C)

    hold(ax, 'on');
    xVals = JointTable.r_fixed_k_per_day;
    yVals = JointTable.r_joint_per_day;

    valid = isfinite(xVals) & isfinite(yVals);
    xVals = xVals(valid);
    yVals = yVals(valid);
    labels = cellstr(JointTable.Label(valid));

    if isempty(xVals)
        title(ax, 'No joint r-k fits available');
        box(ax, 'on');
        return;
    end

    mn = min([xVals; yVals]);
    mx = max([xVals; yVals]);
    span = max(mx-mn, 0.1);
    lo = max(0, mn-0.10*span);
    hi = mx+0.10*span;

    hPoints = scatter(ax, xVals, yVals, 80, ...
        'MarkerEdgeColor',C.fit, ...
        'MarkerFaceColor',C.fit, ...
        'LineWidth',1.0, ...
        'DisplayName','source-specific estimates');
    hDiagonal = plot(ax, [lo hi], [lo hi], '--', ...
        'Color',C.diagonal, ...
        'LineWidth',1.3, ...
        'DisplayName','equal estimates');

   
    labelOffsets = [ ...
         0.012  0.012; ...
         0.012 -0.018; ...
        -0.012  0.016; ...
         0.012  0.017; ...
         0.012  0.018; ...
        -0.012 -0.020];
    for i = 1:numel(xVals)
        j = min(i, size(labelOffsets,1));
        if labelOffsets(j,1) < 0
            hAlign = 'right';
        else
            hAlign = 'left';
        end
        text(ax, ...
            xVals(i) + labelOffsets(j,1)*span, ...
            yVals(i) + labelOffsets(j,2)*span, ...
            labels{i}, ...
            'FontSize',8.5, ...
            'HorizontalAlignment',hAlign, ...
            'VerticalAlignment','middle', ...
            'Clipping','on');
    end

    xlim(ax, [lo hi]);
    ylim(ax, [lo hi]);
    set(ax, ...
        'FontSize',12.5, ...
        'LineWidth',1.15, ...
        'TickDir','out', ...
        'Box','on');

    xlabel(ax, 'r estimated with fixed k (day^{-1})', ...
        'FontSize',13.5, ...
        'FontWeight','bold');
    ylabel(ax, 'r estimated jointly with k (day^{-1})', ...
        'FontSize',13.5, ...
        'FontWeight','bold');

    grid(ax, 'on');
    ax.GridAlpha = 0.16;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';

    legend(ax, [hPoints hDiagonal], ...
        {'source-specific estimates', 'equal estimates'}, ...
        'Location','southeast', ...
        'FontSize',11.5, ...
        'Box','on');
end

function addPanelLetter(ax, panelLetter, fontSize)

    if nargin < 3
        fontSize = 17;
    end

    text(ax, -0.11, 1.05, panelLetter, ...
        'Units', 'normalized', ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'Clipping', 'off');
end

function setForestGrowthRateXAxis(ax, vals, isRapidPanel)

    vals = vals(isfinite(vals) & vals > 0);
    if isempty(vals)
        return;
    end

    if isRapidPanel
        tickStep = 0.5;
        lo = max(0, floor(min(vals)/tickStep)*tickStep);
        hi = ceil(max(vals)*1.08/tickStep)*tickStep;
        hi = max(hi, lo + 3*tickStep);
    else
        tickStep = 0.1;
        lo = max(0, floor(min(vals)/tickStep)*tickStep);
        hi = ceil(max(vals)*1.10/tickStep)*tickStep;
        hi = max(hi, lo + 4*tickStep);
    end

    if lo == 0
        lo = tickStep;
    end

    ticks = lo:tickStep:hi;
    xlim(ax, [lo hi]);
    xticks(ax, ticks);
    xticklabels(ax, arrayfun(@(x) sprintf('%.1f', x), ...
        ticks, 'UniformOutput',false));
    ax.XMinorTick = 'on';
    ax.XMinorGrid = 'on';
end

function setGrowthRateXAxis(ax, vals)

    vals = vals(isfinite(vals) & vals > 0);
    if isempty(vals)
        return;
    end

    lo = max(min(vals)/1.45, 1e-5);
    hi = max(vals)*1.45;
    if hi <= lo
        hi = lo*10;
    end

    xlim(ax, [lo hi]);
    candidateTicks = [1e-5 2e-5 5e-5 ...
                      1e-4 2e-4 5e-4 ...
                      1e-3 2e-3 5e-3 ...
                      1e-2 2e-2 5e-2 ...
                      1e-1 2e-1 5e-1 ...
                      1 2 5 10];
    ticks = candidateTicks( ...
        candidateTicks >= lo & candidateTicks <= hi);
    if numel(ticks) < 3
        ticks = logspace(log10(lo), log10(hi), 4);
    end

    xticks(ax, ticks);
    xticklabels(ax, arrayfun(@(x) sprintf('%.3g', x), ...
        ticks, 'UniformOutput',false));
    ax.XMinorTick = 'on';
    ax.XMinorGrid = 'on';
end

function setCarryingCapacityXAxis(ax, kGrid)

    vals = kGrid(isfinite(kGrid) & kGrid > 0);
    if isempty(vals)
        return;
    end

    ticks = [min(vals), 1e9, 1e10, 1e11, max(vals)];
    ticks = unique(ticks(ticks >= min(vals) & ticks <= max(vals)), ...
        'stable');

    xlim(ax, [min(vals)*0.85, max(vals)*1.15]);
    xticks(ax, ticks);
    xticklabels(ax, arrayfun(@formatKTickLabel, ...
        ticks, 'UniformOutput',false));
    ax.TickLabelInterpreter = 'tex';
    ax.XMinorGrid = 'on';
    ax.XMinorTick = 'on';
    ax.GridAlpha = 0.20;
    ax.MinorGridAlpha = 0.08;
end

function label = formatKTickLabel(kValue)

    if kValue <= 0 || ~isfinite(kValue)
        label = '';
        return;
    end

    exponent = floor(log10(kValue));
    mantissa = kValue/10^exponent;
    if abs(mantissa - 1) < 1e-10
        label = sprintf('10^{%d}', exponent);
    else
        label = sprintf('%.2g\\times10^{%d}', ...
            mantissa, exponent);
    end
end

function exportJournalFigure( ...
    figHandle, baseName, widthInches, heightInches)

    if nargin < 3
        widthInches = 11;
    end
    if nargin < 4
        heightInches = 6.5;
    end

    set(figHandle, 'Color','w');
    set(figHandle, 'Units','inches');
    set(figHandle, 'Position', ...
        [1 1 widthInches heightInches]);
    set(figHandle, 'InvertHardcopy','off');

    axAll = findall(figHandle, 'Type','axes');
    for ii = 1:numel(axAll)
        try
            disableDefaultInteractivity(axAll(ii));
        catch
        end
        try
            axAll(ii).Toolbar.Visible = 'off';
        catch
        end
    end

    drawnow;
    try
        exportgraphics(figHandle, ...
            [baseName '.pdf'], ...
            'ContentType','vector');
        exportgraphics(figHandle, ...
            [baseName '.tif'], ...
            'Resolution',600);
        exportgraphics(figHandle, ...
            [baseName '.png'], ...
            'Resolution',1200, ...
            'BackgroundColor','white');
    catch
        warning('exportgraphics failed; using saveas for %s.', ...
            baseName);
        saveas(figHandle, [baseName '.png']);
        saveas(figHandle, [baseName '.fig']);
    end
end
