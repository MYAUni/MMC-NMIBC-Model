clear; clc; close all;
% Figure 9: injury restoration, LI response and cancer growth.

outDir = fullfile(fileparts(mfilename('fullpath')), 'Fig9_selected_panels');
if ~exist(outDir, 'dir'), mkdir(outDir); end
dpi = 600;
blue = [0.00 0.35 0.70];
red = [0.75 0.00 0.00];
bandColor = [0.65 0.82 0.95];
options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5, ...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off');

%% LI data and delayed-logistic derivative fits
% LI(t) = A*du/dt; u is a normalized logistic profile, A has units % days.
tLI_h = [1 4 8 16 24 72 110 240]';
tLI = tLI_h / 24;
LI = [0 0 1.9 31.2 21.9 3.5 2.1 0.2]';
LI_sd = [0 0 0.8 8 8 1 0.1 0.1]';
fitMask = tLI_h <= 110; % The 240 h point is excluded from fitting and these plots.
area0 = 0.1; % Fixed LI-integral offset (% days); preserves the supplied fit.
tauBounds = [4 8] / 24;
sigmaStarts = [1 2 5 10 15 20 30 50 80 120 200];
AStarts = [1 5 10 20 30 50 80 120 200 400 800];
tauStarts_h = [4.2 5 6 7 7.8];
LI_targets = [LI, max(LI - LI_sd, 0), LI + LI_sd];
LI_fits = cell(1, 3);
fitNames = {'Mean', 'Mean minus SD', 'Mean plus SD'};

for j = 1:3
    fprintf('Fitting LI: %s...\n', fitNames{j});
    LI_fits{j} = fitLI(tLI(fitMask), LI_targets(fitMask,j), ...
        area0, tauBounds, sigmaStarts, AStarts, tauStarts_h, options);
end

tLI_smooth = linspace(0, 120/24, 1500)';
LI_curves = zeros(numel(tLI_smooth), 3);
u_curves = zeros(numel(tLI_smooth), 3);
du_curves = zeros(numel(tLI_smooth), 3);
for j = 1:3
    f = LI_fits{j};
    [u_curves(:,j), du_curves(:,j), LI_curves(:,j)] = delayedLIProfile( ...
        f.sigma, f.A, f.tau, area0, tLI_smooth);
    fprintf('\nLI %s: sigma = %.6f /day, A = %.6f %% days, tau = %.6f h\n', ...
        fitNames{j}, f.sigma, f.A, 24*f.tau);
    fprintf('SSE = %.6f, RMSE = %.6f percentage points, R^2 = %.6f\n', ...
        f.SSE, f.RMSE, f.R2);
    tPeak = f.tau + log((f.A-area0)/area0)/f.sigma;
    fprintf('Fitted peak = %.6f %% at %.6f h\n', f.sigma*f.A/4, 24*tPeak);
    disp(table(tLI_h(fitMask), LI_targets(fitMask,j), f.fitted, f.residuals, ...
        'VariableNames', {'Time_h','LI_percent','Fitted_percent','Residual'}));
end

% These envelopes join fits to mean +/- SD; they are not confidence intervals.
LI_lo = min(LI_curves(:,2:3), [], 2);
LI_hi = max(LI_curves(:,2:3), [], 2);
LI_ylim = [0, max(55, 1.10*max([LI_curves(:); LI+LI_sd]))];
schematicScale = 1; % a.u. per (% day); one common scale for all profiles.
profileScales = cellfun(@(f) f.A, LI_fits)*schematicScale;
restoration_curves = bsxfun(@times, u_curves, profileScales);
restoration_lo = min(restoration_curves(:,2:3), [], 2);
restoration_hi = max(restoration_curves(:,2:3), [], 2);

% Injury restoration with envelope
[fig, ax] = panel('Fig9_LI_integrated_envelope', ...
    {'Logistic replacement profile','with variability envelope'}, ...
    'Time after injury (h)', 'Schematic cell-count increase (a.u.)');
hBand = fill(ax, [24*tLI_smooth; flipud(24*tLI_smooth)], ...
    [restoration_lo; flipud(restoration_hi)], bandColor, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.50);
hMean = plot(ax, 24*tLI_smooth, restoration_curves(:,1), '-', ...
    'Color', blue, 'LineWidth', 2.2);
panelLegend(ax, [hBand hMean], ...
    {'Fits to mean LI \pm SD','Profile from mean LI'}, 'southeast');
injuryAxes(ax, [0 1.10*max(restoration_curves(:))]);
savePanel(fig, outDir, dpi);

% LI response with envelope
[fig, ax] = panel('Fig9_LI_fit_envelope', ...
    {'Proliferative response','with variability envelope'}, ...
    'Time after injury (h)', 'Labeling index (%)');
plotLIEnvelope(ax, tLI_h(fitMask), LI(fitMask), LI_sd(fitMask), ...
    24*tLI_smooth, LI_curves(:,1), LI_lo, LI_hi, blue, bandColor);
injuryAxes(ax, LI_ylim);
savePanel(fig, outDir, dpi);

% Cancer cell-count data and logistic fit
% Counts were supplied after area-to-cell conversion; fit log10(count) residuals.
% Four-point phenomenological fit with K > the largest observed count.
tCancer = [10 14 17 24]';
N = [7.20e4 3.24e5 4.95e5 9.54e5]';
Kmin = 1.001*max(N);
rStarts = [0.03 0.05 0.08 0.10 0.15 0.20 0.30 0.50 0.80];
KMultipliers = [1.02 1.05 1.10 1.25 1.50 2.00 3.00 5.00 10.00];
tInfStarts = [8 10 12 14 16 17 18 20 22 24];
model = @(q,t) logisticN(exp(q(1)), Kmin+exp(q(2)), q(3), t);
objective = @(q) sum((log10(N)-log10(model(q,tCancer))).^2);
bestSSE = Inf;
bestQ = [];
fprintf('\nFitting cancer cell counts...\n');
for r0 = rStarts
    for km = KMultipliers
        for t0 = tInfStarts
            K0 = max(km*max(N), Kmin+1);
            q0 = [log(r0), log(K0-Kmin), t0];
            q = fminsearch(objective, q0, options);
            sse = objective(q);
            if isfinite(sse) && sse < bestSSE
                bestSSE = sse;
                bestQ = q;
            end
        end
    end
end
assert(~isempty(bestQ), 'No finite cancer fit was found.');
cancer.r = exp(bestQ(1));
cancer.K = Kmin+exp(bestQ(2));
cancer.tInf = bestQ(3);
cancer.fitted = model(bestQ, tCancer);
cancer.residuals = N-cancer.fitted;
cancer.SSE_linear = sum(cancer.residuals.^2);
cancer.RMSE_linear = sqrt(mean(cancer.residuals.^2));
cancer.R2_linear = 1-cancer.SSE_linear/sum((N-mean(N)).^2);
cancer.SSE_log10 = bestSSE;
cancer.RMSE_log10 = sqrt(bestSSE/numel(N));
tCancer_smooth = linspace(0, 30, 1500)';
N_smooth = logisticN(cancer.r, cancer.K, cancer.tInf, tCancer_smooth);
N_rate = cancer.r*N_smooth.*(1-N_smooth/cancer.K);

fprintf('\nCancer: r = %.6f /day, K = %.6e cells, t_inf = %.6f days\n', ...
    cancer.r, cancer.K, cancer.tInf);
fprintf('K/2 = %.6e cells; peak dN/dt = %.6e cells/day\n', ...
    cancer.K/2, cancer.r*cancer.K/4);
fprintf('Linear SSE = %.6e, RMSE = %.6e cells, R^2 = %.6f\n', ...
    cancer.SSE_linear, cancer.RMSE_linear, cancer.R2_linear);
fprintf('Log10 SSE = %.6e, RMSE = %.6e\n', cancer.SSE_log10, cancer.RMSE_log10);
disp(table(tCancer, N, cancer.fitted, cancer.residuals, ...
    'VariableNames', {'Time_days','Data_cells','Fitted_cells','Residual_cells'}));

% Cancer growth on a logarithmic axis
[fig, ax] = panel('Fig9_Cancer_logscale', {'Cancer-growth logistic fit','Log scale'}, ...
    'Time (days)', 'Tumor cell count');
hFit = plot(ax, tCancer_smooth, N_smooth, '-', 'Color', red, 'LineWidth', 2.3);
hData = plot(ax, tCancer, N, 'ko', ...
    'MarkerFaceColor', 'w', 'MarkerSize', 7, 'LineWidth', 1.4);
set(ax, 'YScale', 'log');
panelLegend(ax, [hFit hData], ...
    {'Fitted logistic profile','Converted cell-count data'}, 'southeast');
xlim(ax, [0 30]);
savePanel(fig, outDir, dpi);

save(fullfile(outDir, 'Fig9_fit_results.mat'), 'tLI_h', 'LI', 'LI_sd', ...
    'fitMask', 'area0', 'tauBounds', 'LI_fits', 'tLI_smooth', 'LI_curves', ...
    'u_curves', 'du_curves', 'schematicScale', 'profileScales', 'restoration_curves', ...
    'tCancer', 'N', 'cancer', 'tCancer_smooth', 'N_smooth', 'N_rate');
fprintf('\nThree separate figures saved as PNG (%d dpi), PDF and FIG in:\n%s\n', dpi, outDir);

%% Functions
function fit = fitLI(t, y, area0, bounds, sigmaStarts, AStarts, tauStarts_h, options)
    model = @(q) liModel(q, t, area0, bounds);
    objective = @(q) sum((y-model(q)).^2);
    bestSSE = Inf;
    bestQ = [];
    for sigma0 = sigmaStarts
        for A0 = AStarts
            for tau0_h = tauStarts_h
                A0 = max(A0, area0+1e-6);
                p0 = (tau0_h/24-bounds(1))/diff(bounds);
                p0 = min(max(p0, 1e-6), 1-1e-6);
                q0 = [log(sigma0), log(A0-area0), log(p0/(1-p0))];
                q = fminsearch(objective, q0, options);
                sse = objective(q);
                if isfinite(sse) && sse < bestSSE
                    bestSSE = sse;
                    bestQ = q;
                end
            end
        end
    end
    assert(~isempty(bestQ), 'No finite LI fit was found.');
    fit.sigma = exp(bestQ(1));
    fit.A = area0+exp(bestQ(2));
    fit.tau = bounds(1)+diff(bounds)/(1+exp(-bestQ(3)));
    fit.fitted = model(bestQ);
    fit.residuals = y-fit.fitted;
    fit.SSE = bestSSE;
    fit.RMSE = sqrt(mean(fit.residuals.^2));
    sst = sum((y-mean(y)).^2);
    fit.R2 = NaN;
    if sst > 0, fit.R2 = 1-bestSSE/sst; end
end

function y = liModel(q, t, area0, bounds)
    tau = bounds(1)+diff(bounds)/(1+exp(-q(3)));
    [~, ~, y] = delayedLIProfile(exp(q(1)), area0+exp(q(2)), tau, area0, t);
end

function [u, du, LI] = delayedLIProfile(sigma, A, tau, area0, t)
    z = area0*ones(size(t));
    active = t > tau;
    z(active) = A./(1+((A-area0)/area0).*exp(-sigma.*(t(active)-tau)));
    LI = zeros(size(t));
    LI(active) = sigma.*z(active).*(1-z(active)/A);
    u = z/A;
    du = LI/A;
end

function N = logisticN(r, K, tInf, t)
    N = K./(1+exp(-r.*(t-tInf)));
end

function [fig, ax] = panel(name, heading, xLabel, yLabel)
    fig = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 12.5 9.5]);
    ax = axes('Parent', fig);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold', ...
        'LineWidth', 1.1, 'GridAlpha', 0.22);
    xlabel(ax, xLabel, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, yLabel, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
    %title(ax, heading, 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
end

function panelLegend(ax, handles, labels, location)
    legend(ax, handles, labels, 'Location', location, 'Box', 'off', ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'AutoUpdate', 'off');
end

function injuryAxes(ax, limits)
    xlim(ax, [0 120]); xticks(ax, [0 24 48 72 96 120]); ylim(ax, limits);
end

function plotLIEnvelope(ax, t, y, sd, ts, curve, lo, hi, color, bandColor)
    hBand = fill(ax, [ts; flipud(ts)], [lo; flipud(hi)], bandColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.50);
    hData = errorbar(ax, t, y, sd, 'ko', 'MarkerFaceColor', 'w', ...
        'MarkerSize', 6.5, 'LineWidth', 1.2, 'CapSize', 6, 'LineStyle', 'none');
    hFit = plot(ax, ts, curve, '-', 'Color', color, 'LineWidth', 2.2);
    panelLegend(ax, [hData hBand hFit], ...
        {'Measured LI \pm SD','Fits to mean \pm SD','Fitted logistic derivative'}, 'northeast');
end

function savePanel(fig, outDir, dpi)
    drawnow;
    base = fullfile(outDir, get(fig, 'Name'));
    savefig(fig, [base '.fig']);
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(fig, [base '.png'], 'Resolution', dpi, 'BackgroundColor', 'white');
end
