% -------------------------------------------------------------------------
% PURPOSE
%   Post-process the saved paper-mode Sobol GSA results without rerunning
%   the ODE model and without changing the system, bounds, sampling rules,
%   pulse schedule, output definition, or Sobol estimators.
%
% WORKFLOW
%   1. Prompts the user to select the completed paper-mode
%      Sobol_GSA_Day360_complete_results*.mat file.
%   2. Computes a nested-prefix convergence diagnostic at
%      N = 128, 256, 512 and 1024 from the already saved YA, YB and YAB.
%   3. Regenerates publication Figure 2 with nonoverlapping panel labels
%      A and B placed in the figure margin (not inside either axes).
%   4. Regenerates the two supplementary GSA figures with legends
%      outside the data region, and creates one convergence figure.
%
%
% OUTPUTS
%   A timestamped folder is created beside the selected MAT file containing:
%     Fig02_GSA_Day360_Sobol_S1_ST.pdf/.png
%     FigS_GSA_Output_and_Combined_Indices.pdf/.png
%     FigS_GSA_r_Sampling_Scale_Robustness.pdf/.png
%     FigS_GSA_Nested_Sample_Convergence.pdf/.png
%     Sobol_nested_sample_convergence.csv
%     Sobol_nested_sample_convergence_summary.csv
%     GSA_postprocessing_notes.txt
%
% -------------------------------------------------------------------------

clearvars
clc
close all

fprintf('\n============================================================\n')
fprintf('GSA POST-PROCESSING FROM SAVED PAPER-MODE RESULTS\n')
fprintf('No ODE simulations will be run.\n')
fprintf('============================================================\n\n')

[matName,matFolder] = uigetfile( ...
    {'*.mat','MATLAB result files (*.mat)'}, ...
    'Select the completed paper-mode Sobol_GSA_Day360_complete_results MAT file');

if isequal(matName,0)
    error('No MAT file was selected. Post-processing was cancelled.')
end

matPath = fullfile(matFolder,matName);
fprintf('Selected saved result:\n%s\n\n',matPath)

saved = load(matPath,'primary','rLog','factors','settings');
requiredTopLevel = {'primary','rLog','factors','settings'};
for j = 1:numel(requiredTopLevel)
    if ~isfield(saved,requiredTopLevel{j})
        error('Selected MAT file lacks required variable "%s".', ...
            requiredTopLevel{j})
    end
end

primary = saved.primary;
rLog = saved.rLog;
factors = saved.factors;
sourceSettings = saved.settings;

cfg = makePostprocessSettings(matFolder,sourceSettings);
validateSavedPaperResults(primary,rLog,factors,sourceSettings,cfg)

if ~exist(cfg.outputDir,'dir')
    mkdir(cfg.outputDir)
end

fprintf('Output folder:\n%s\n\n',cfg.outputDir)

%% Nested-prefix convergence diagnostics from saved model evaluations
fprintf('Computing nested-prefix convergence diagnostics from saved outputs ...\n')
[convergenceTable,convergenceSummary,convergenceSobol] = ...
    buildConvergenceDiagnostics(primary,rLog,factors,cfg.sampleSizes);

writetable(convergenceTable,fullfile(cfg.outputDir, ...
    'Sobol_nested_sample_convergence.csv'))
writetable(convergenceSummary,fullfile(cfg.outputDir, ...
    'Sobol_nested_sample_convergence_summary.csv'))

% Regenerate  GSA figures from the saved final estimates
fprintf('Regenerating Figure 2 and supplementary GSA figures ...\n')
makeSectionFigure2(primary.sobol,factors,cfg)
makeDiagnosticSobolFigure(primary,factors,cfg)

rhoPrimaryVsRlog = localSpearman(primary.sobol.ST,rLog.sobol.ST);
makeRScaleRobustnessFigure(primary.sobol,rLog.sobol, ...
    factors,cfg,rhoPrimaryVsRlog)
makeConvergenceFigure(convergenceSobol,factors,cfg)

writePostprocessingNotes(matPath,primary,rLog,cfg,rhoPrimaryVsRlog)

fprintf('\n============================================================\n')
fprintf('POST-PROCESSING COMPLETED SUCCESSFULLY\n')
fprintf('No ODE simulations were rerun.\n')
fprintf('Outputs saved in:\n%s\n',cfg.outputDir)
fprintf('============================================================\n')

%% ========================================================================
% LOCAL FUNCTIONS: SETTINGS AND VALIDATION
% ========================================================================

function cfg = makePostprocessSettings(matFolder,sourceSettings)
    cfg.sampleSizes = [128 256 512 1024];
    cfg.requiredPaperBaseN = 1024;
    cfg.topN = 15;
    cfg.convergenceTopN = 6;
    cfg.exportResolution = 1200;
    cfg.figureFont = 'Arial';
    cfg.mainFigureWidth = 13.4;
    cfg.mainFigureHeight = 6.2;
    cfg.detectionLimit = pi*(1/2)^2*(3*0.01)*1e6;

    if isfield(sourceSettings,'topN') && isfinite(sourceSettings.topN)
        cfg.topN = sourceSettings.topN;
    end
    if isfield(sourceSettings,'figureFont') && ...
            (ischar(sourceSettings.figureFont) || ...
             isstring(sourceSettings.figureFont))
        cfg.figureFont = char(sourceSettings.figureFont);
    end
    if isfield(sourceSettings,'detectionLimit') && ...
            isfinite(sourceSettings.detectionLimit) && ...
            sourceSettings.detectionLimit>0
        cfg.detectionLimit = sourceSettings.detectionLimit;
    end

    runStamp = datestr(now,'yyyymmdd_HHMMSS');
    cfg.outputDir = fullfile(matFolder, ...
        ['GSA_Postprocessed_Publication_' runStamp]);
end

function validateSavedPaperResults(primary,rLog,factors,sourceSettings,cfg)
    if isfield(sourceSettings,'runMode')
        if ~strcmpi(char(string(sourceSettings.runMode)),'paper')
            error(['The selected MAT file is not marked as a paper-mode ' ...
                'result. Select the completed paper-mode MAT file.'])
        end
    end

    validateScenario(primary,'primary',factors)
    validateScenario(rLog,'rLog',factors)

    nPrimary = numel(primary.YA);
    nRlog = numel(rLog.YA);
    if nPrimary ~= nRlog
        error('Primary and r-log scenarios have different base sample sizes.')
    end
    if nPrimary < cfg.requiredPaperBaseN
        error(['The selected result has N=%d. Publication post-processing ' ...
            'requires the completed paper-mode result with N>=1024.'],nPrimary)
    end
    if any(cfg.sampleSizes>nPrimary)
        error('A requested convergence size exceeds the saved base sample size.')
    end

    % Recompute the final point estimates once from the saved outputs and
    % require exact estimator agreement (within floating-point tolerance).
    checkPrimary = computeSobolSubset(primary.YA,primary.YB, ...
        primary.YAB,factors);
    checkRlog = computeSobolSubset(rLog.YA,rLog.YB,rLog.YAB,factors);
    requireEstimatorAgreement(checkPrimary,primary.sobol,'primary')
    requireEstimatorAgreement(checkRlog,rLog.sobol,'r-log robustness')
end

function validateScenario(scenario,scenarioName,factors)
    required = {'YA','YB','YAB','sobol'};
    for j = 1:numel(required)
        if ~isfield(scenario,required{j})
            error('Scenario %s lacks required field "%s".', ...
                scenarioName,required{j})
        end
    end
    n = numel(scenario.YA);
    D = numel(factors);
    if numel(scenario.YB)~=n || size(scenario.YAB,1)~=n || ...
            size(scenario.YAB,2)~=D
        error('Scenario %s has incompatible YA, YB or YAB dimensions.', ...
            scenarioName)
    end
    if any(~isfinite(scenario.YA(:))) || any(~isfinite(scenario.YB(:))) || ...
            any(~isfinite(scenario.YAB(:)))
        error('Scenario %s contains nonfinite saved model outputs.',scenarioName)
    end

    requiredSobol = {'S1','ST','S1ciLow','S1ciHigh', ...
        'STciLow','STciHigh','rankByS1','rankByST'};
    for j = 1:numel(requiredSobol)
        if ~isfield(scenario.sobol,requiredSobol{j})
            error('Scenario %s lacks saved Sobol field "%s".', ...
                scenarioName,requiredSobol{j})
        end
    end
end

function requireEstimatorAgreement(recomputed,savedSobol,scenarioName)
    scale = max(1,max(abs([savedSobol.S1(:);savedSobol.ST(:)])));
    tol = 5e-12*scale;
    maxDifference = max(abs([recomputed.S1(:)-savedSobol.S1(:); ...
        recomputed.ST(:)-savedSobol.ST(:)]));
    if ~(isfinite(maxDifference) && maxDifference<=tol)
        error(['Saved and recomputed Sobol point estimates disagree for ' ...
            '%s (maximum difference %.4g).'],scenarioName,maxDifference)
    end
end

%% ========================================================================
% LOCAL FUNCTIONS: NESTED-SAMPLE CONVERGENCE
% ========================================================================

function [longTable,summaryTable,store] = ...
        buildConvergenceDiagnostics(primary,rLog,factors,sampleSizes)
    scenarios = {primary,rLog};
    scenarioNames = ["Primary_uniform_r","Robustness_log_uniform_r"];
    nScenarios = numel(scenarios);
    nSizes = numel(sampleSizes);
    D = numel(factors);

    store = cell(nScenarios,nSizes);
    nLong = nScenarios*nSizes*D;
    Scenario = strings(nLong,1);
    BaseSampleSize = zeros(nLong,1);
    Factor = strings(nLong,1);
    S1_Raw = nan(nLong,1);
    ST_Raw = nan(nLong,1);
    S1_Rank = nan(nLong,1);
    ST_Rank = nan(nLong,1);

    nSummary = nScenarios*nSizes;
    SummaryScenario = strings(nSummary,1);
    SummaryN = zeros(nSummary,1);
    OutputVariance = nan(nSummary,1);
    SumS1 = nan(nSummary,1);
    SumST = nan(nSummary,1);
    ST_RankSpearmanToFinal = nan(nSummary,1);
    Top5_ST_OverlapWithFinal = nan(nSummary,1);
    MaximumAbsoluteSTDifferenceFromFinal = nan(nSummary,1);

    longRow = 0;
    summaryRow = 0;
    factorNames = string({factors.name}).';

    for s = 1:nScenarios
        sc = scenarios{s};
        for q = 1:nSizes
            N = sampleSizes(q);
            store{s,q} = computeSobolSubset(sc.YA(1:N),sc.YB(1:N), ...
                sc.YAB(1:N,:),factors);
        end
        finalSobol = store{s,end};

        for q = 1:nSizes
            N = sampleSizes(q);
            sobol = store{s,q};
            rankS1 = inverseRank(sobol.rankByS1,D);
            rankST = inverseRank(sobol.rankByST,D);

            rows = longRow+(1:D);
            Scenario(rows) = scenarioNames(s);
            BaseSampleSize(rows) = N;
            Factor(rows) = factorNames;
            S1_Raw(rows) = sobol.S1;
            ST_Raw(rows) = sobol.ST;
            S1_Rank(rows) = rankS1;
            ST_Rank(rows) = rankST;
            longRow = longRow+D;

            summaryRow = summaryRow+1;
            SummaryScenario(summaryRow) = scenarioNames(s);
            SummaryN(summaryRow) = N;
            OutputVariance(summaryRow) = sobol.outputVariance;
            SumS1(summaryRow) = sobol.sumS1;
            SumST(summaryRow) = sobol.sumST;
            ST_RankSpearmanToFinal(summaryRow) = ...
                localSpearman(sobol.ST,finalSobol.ST);
            topCount = min(5,D);
            Top5_ST_OverlapWithFinal(summaryRow) = numel(intersect( ...
                sobol.rankByST(1:topCount), ...
                finalSobol.rankByST(1:topCount)));
            MaximumAbsoluteSTDifferenceFromFinal(summaryRow) = ...
                max(abs(sobol.ST-finalSobol.ST));
        end
    end

    longTable = table(Scenario,BaseSampleSize,Factor,S1_Raw,ST_Raw, ...
        S1_Rank,ST_Rank);
    summaryTable = table(SummaryScenario,SummaryN,OutputVariance, ...
        SumS1,SumST,ST_RankSpearmanToFinal, ...
        Top5_ST_OverlapWithFinal, ...
        MaximumAbsoluteSTDifferenceFromFinal, ...
        'VariableNames',{'Scenario','BaseSampleSize','OutputVariance', ...
        'Sum_S1','Sum_ST','ST_Rank_Spearman_to_Final_N', ...
        'Top5_ST_Overlap_with_Final_N', ...
        'Maximum_Absolute_ST_Difference_from_Final_N'});
end

function sobol = computeSobolSubset(YA,YB,YAB,factors)
    YA = YA(:);
    YB = YB(:);
    outputVariance = var([YA;YB],1);
    if ~(isfinite(outputVariance) && outputVariance>0)
        error('A nested subset has zero or nonfinite output variance.')
    end

    D = size(YAB,2);
    if D~=numel(factors)
        error('Saved hybrid-output width does not match the factor list.')
    end
    S1 = nan(D,1);
    ST = nan(D,1);
    for j = 1:D
        % Exact estimators used in the GSA code:
        % Saltelli first-order and Jansen total-order.
        S1(j) = mean(YB.*(YAB(:,j)-YA))/outputVariance;
        ST(j) = mean((YA-YAB(:,j)).^2)/(2*outputVariance);
    end

    [~,rankByS1] = sort(S1,'descend');
    [~,rankByST] = sort(ST,'descend');
    sobol = struct('S1',S1,'ST',ST, ...
        'outputVariance',outputVariance, ...
        'sumS1',sum(S1),'sumST',sum(ST), ...
        'rankByS1',rankByS1,'rankByST',rankByST);
end

function ranks = inverseRank(order,D)
    ranks = nan(D,1);
    ranks(order) = 1:D;
end

%% ========================================================================
% LOCAL FUNCTIONS: PUBLICATION FIGURE 2
% ========================================================================

function makeSectionFigure2(sobol,factors,cfg)
    fig = figure('Color','w','Name', ...
        'Figure 2: day-360 global sensitivity analysis');
    set(fig,'Units','inches','Position',[1 1 ...
        cfg.mainFigureWidth cfg.mainFigureHeight]);
    tl = tiledlayout(fig,1,2,'TileSpacing','loose','Padding','compact');

    topN = min(cfg.topN,numel(factors));
    xCommon = 1.12*max([sobol.S1ciHigh(:);sobol.STciHigh(:); ...
        sobol.S1(:);sobol.ST(:);0.05]);
    xCommon = min(1,max(0.05,xCommon));

    axA = nexttile(tl,1);
    idxS1 = flip(sobol.rankByS1(1:topN));
    plotSingleIndexPanel(axA,sobol.S1,sobol.S1ciLow, ...
        sobol.S1ciHigh,idxS1,factors, ...
        'First-order Sobol index, S_1','Independent effects', ...
        [0.82 0.22 0.13],xCommon,cfg)

    axB = nexttile(tl,2);
    idxST = flip(sobol.rankByST(1:topN));
    plotSingleIndexPanel(axB,sobol.ST,sobol.STciLow, ...
        sobol.STciHigh,idxST,factors, ...
        'Total-order Sobol index, S_T', ...
        'Total effects including interactions', ...
        [0.22 0.43 0.70],xCommon,cfg)


    drawnow
    placeFigure2PanelLabels(fig,[axA axB],{'A','B'})

    exportPublicationFigure(fig,cfg,'Fig02_GSA_Day360_Sobol_S1_ST')
end

function plotSingleIndexPanel(ax,values,ciLow,ciHigh,idx,factors, ...
        xLabel,panelTitle,faceColor,xCommon,cfg)
    y = 1:numel(idx);

    % Keep raw estimates in exported CSVs. Truncate only the visual bars
    % and visual interval endpoints at zero.
    vals = max(values(idx),0);
    lo = max(ciLow(idx),0);
    hi = max(ciHigh(idx),0);

    hold(ax,'on')
    barh(ax,y,vals,0.68,'FaceColor',faceColor,'EdgeColor','none', ...
        'FaceAlpha',0.90)
    for i = 1:numel(idx)
        drawHorizontalInterval(ax,y(i),lo(i),hi(i))
    end
    yticks(ax,y)
    yticklabels(ax,{factors(idx).label})
    set(ax,'TickLabelInterpreter','tex')
    xlabel(ax,xLabel,'FontWeight','bold','Interpreter','tex')
    ylabel(ax,'Model input','FontWeight','bold')
    title(ax,panelTitle,'FontWeight','bold')
    xlim(ax,[0 xCommon])
    ylim(ax,[0.35 numel(idx)+0.65])
    formatAxes(ax,cfg)
end

function placeFigure2PanelLabels(fig,axesHandles,letters)
    for j = 1:numel(axesHandles)
        pos = axesHandles(j).Position;
        x = max(0.006,pos(1)-0.052);
        y = min(0.955,pos(2)+pos(4)+0.010);
        annotation(fig,'textbox',[x y 0.035 0.040], ...
            'String',letters{j},'EdgeColor','none', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontName','Arial','FontSize',23,'FontWeight','bold', ...
            'Margin',0,'FitBoxToText','off');
    end
end

%% ========================================================================
% LOCAL FUNCTIONS: SUPPLEMENTARY FIGURES 
% ========================================================================

function makeDiagnosticSobolFigure(primary,factors,cfg)
    Ydist = [primary.YA;primary.YB];
    fig = figure('Color','w','Name','Supplementary GSA diagnostics');
    set(fig,'Units','inches','Position',[1 1 13.6 6.3]);
    tl = tiledlayout(fig,1,2,'TileSpacing','loose','Padding','compact');

    axA = nexttile(tl,1);
    histogram(axA,Ydist,32,'FaceColor',[0.35 0.35 0.35], ...
        'EdgeColor','none','FaceAlpha',0.88);
    hold(axA,'on')
    xRef = log10(cfg.detectionLimit+1);
    yLimits = ylim(axA);
    plot(axA,[xRef xRef],yLimits,'--','Color',[0.1 0.1 0.1], ...
        'LineWidth',1.35)
    ylim(axA,yLimits)
    text(axA,xRef+0.08,yLimits(1)+0.82*(yLimits(2)-yLimits(1)), ...
        '1-mm detection threshold','FontWeight','bold', ...
        'FontSize',12,'BackgroundColor','w','Margin',1, ...
        'Clipping','on')
    xlabel(axA,'log_{10}(T_s(360)+T_u(360)+1)', ...
        'FontWeight','bold','Interpreter','tex')
    ylabel(axA,'Number of simulations','FontWeight','bold')
    title(axA,'Day-360 GSA output','FontWeight','bold')
    formatAxes(axA,cfg)

    axB = nexttile(tl,2);
    plotCombinedSobolPanel(axB,primary.sobol,factors,cfg)
    title(axB,'Ranked sensitivity indices','FontWeight','bold')


    drawnow
    placeFigure2PanelLabels(fig,[axA axB],{'A','B'})

    exportPublicationFigure(fig,cfg, ...
        'FigS_GSA_Output_and_Combined_Indices')
end

function plotCombinedSobolPanel(ax,sobol,factors,cfg)
    topN = min(cfg.topN,numel(factors));
    idx = flip(sobol.rankByST(1:topN));
    y = 1:topN;
    ST = max(sobol.ST(idx),0);
    S1 = max(sobol.S1(idx),0);
    STlo = max(sobol.STciLow(idx),0);
    SThi = max(sobol.STciHigh(idx),0);
    S1lo = max(sobol.S1ciLow(idx),0);
    S1hi = max(sobol.S1ciHigh(idx),0);

    hold(ax,'on')
    hST = barh(ax,y,ST,0.78,'FaceColor',[0.22 0.43 0.70], ...
        'EdgeColor','none','FaceAlpha',0.88);
    hS1 = barh(ax,y,S1,0.38,'FaceColor',[0.82 0.22 0.13], ...
        'EdgeColor','none','FaceAlpha',0.92);
    for i = 1:topN
        drawHorizontalInterval(ax,y(i)+0.18,STlo(i),SThi(i))
        drawHorizontalInterval(ax,y(i)-0.18,S1lo(i),S1hi(i))
    end
    yticks(ax,y)
    yticklabels(ax,{factors(idx).label})
    set(ax,'TickLabelInterpreter','tex')
    xlabel(ax,'Sobol index','FontWeight','bold','Interpreter','tex')
    ylabel(ax,'Model input','FontWeight','bold')
    xMax = max([ST(:);S1(:);SThi(:);S1hi(:);0.05]);
    xlim(ax,[0 min(1,1.12*xMax)])
    ylim(ax,[0.35 topN+0.65])
    legend(ax,[hST hS1],{'Total order, S_T','First order, S_1'}, ...
        'Location','southoutside','Orientation','horizontal', ...
        'Box','off','Interpreter','tex')
    formatAxes(ax,cfg)
end

function makeRScaleRobustnessFigure(primarySobol,rLogSobol, ...
        factors,cfg,rho)
    D = numel(factors);
    topN = min(cfg.topN,D);
    [~,rankByMax] = sort(max([primarySobol.ST rLogSobol.ST],[],2), ...
        'descend');
    idx = flip(rankByMax(1:topN));
    y = 1:topN;
    values = [max(primarySobol.ST(idx),0) max(rLogSobol.ST(idx),0)];

   fig = figure('Color','w','Name', ...
    'Robustness to computational sampling measure for r');
    set(fig,'Units','inches','Position',[1 1 9.2 7.1]);
    ax = axes(fig);
    barHandles = barh(ax,y,values,0.75,'grouped','EdgeColor','none');
    barHandles(1).FaceColor = [0.10 0.45 0.74];
    barHandles(2).FaceColor = [0.85 0.33 0.10];
    yticks(ax,y)
    yticklabels(ax,{factors(idx).label})
    set(ax,'TickLabelInterpreter','tex')
    xlabel(ax,'Total-order Sobol index, S_T','FontWeight','bold', ...
        'Interpreter','tex')
    ylabel(ax,'Model input','FontWeight','bold')
  title(ax,{ ...
    'Robustness to the computational sampling measure for r', ...
    sprintf('(rank \\rho = %.3f)',rho)}, ...
    'FontWeight','bold','Interpreter','tex')
    legend(ax,{'Primary: uniform r','Robustness: log-uniform r'}, ...
        'Location','southoutside','Orientation','horizontal','Box','off')
    xMax = max(values(:));
    xlim(ax,[0 min(1,max(0.05,1.12*xMax))])
    formatAxes(ax,cfg)
    exportPublicationFigure(fig,cfg, ...
        'FigS_GSA_r_Sampling_Scale_Robustness')
end

%% ========================================================================
% LOCAL FUNCTIONS: CONVERGENCE FIGURE
% ========================================================================

function makeConvergenceFigure(store,factors,cfg)
    sampleSizes = cfg.sampleSizes;
    finalPrimary = store{1,end};
    topN = min(cfg.convergenceTopN,numel(factors));
    topIdx = finalPrimary.rankByST(1:topN);

    fig = figure('Color','w','Name','Nested-sample convergence diagnostic');
    set(fig,'Units','inches','Position',[1 1 13.2 5.9]);
    tl = tiledlayout(fig,1,2,'TileSpacing','loose','Padding','compact');

    axA = nexttile(tl,1);
    hold(axA,'on')
    colors = lines(topN);
    lineHandles = gobjects(topN,1);
    for j = 1:topN
        values = nan(size(sampleSizes));
        for q = 1:numel(sampleSizes)
            values(q) = max(store{1,q}.ST(topIdx(j)),0);
        end
        lineHandles(j) = semilogx(axA,sampleSizes,values,'-o', ...
            'Color',colors(j,:),'LineWidth',1.6,'MarkerSize',5.5, ...
            'MarkerFaceColor',colors(j,:));
    end
    xticks(axA,sampleSizes)
    xlim(axA,[0.9*sampleSizes(1) 1.1*sampleSizes(end)])
    xlabel(axA,'Nested base sample size, N','FontWeight','bold')
    ylabel(axA,'Total-order Sobol index, S_T', ...
        'FontWeight','bold','Interpreter','tex')
    title(axA,'Primary total-order index stability','FontWeight','bold')
    factorLabels = {factors(topIdx).label};
    leg = legend(axA,lineHandles,factorLabels,'Location','southoutside', ...
        'Orientation','horizontal','Box','off','Interpreter','tex');
    try
        leg.NumColumns = 3;
    catch
    end
    formatAxes(axA,cfg)

    axB = nexttile(tl,2);
    hold(axB,'on')
    rhoPrimary = nan(size(sampleSizes));
    rhoRlog = nan(size(sampleSizes));
    for q = 1:numel(sampleSizes)
        rhoPrimary(q) = localSpearman(store{1,q}.ST,store{1,end}.ST);
        rhoRlog(q) = localSpearman(store{2,q}.ST,store{2,end}.ST);
    end
    plot(axB,sampleSizes,rhoPrimary,'-o','LineWidth',1.8, ...
        'MarkerSize',6,'Color',[0.10 0.45 0.74], ...
        'MarkerFaceColor',[0.10 0.45 0.74])
    plot(axB,sampleSizes,rhoRlog,'-s','LineWidth',1.8, ...
        'MarkerSize',6,'Color',[0.85 0.33 0.10], ...
        'MarkerFaceColor',[0.85 0.33 0.10])
    xticks(axB,sampleSizes)
    xlim(axB,[0.9*sampleSizes(1) 1.1*sampleSizes(end)])
    ylim(axB,[0 1.03])
    xlabel(axB,'Nested base sample size, N','FontWeight','bold')
    ylabel(axB,'S_T rank correlation with N=1024', ...
        'FontWeight','bold','Interpreter','tex')
    title(axB,'Rank agreement with the final design','FontWeight','bold')
    legend(axB,{'Primary: uniform r','Robustness: log-uniform r'}, ...
        'Location','southoutside','Orientation','horizontal','Box','off')
    formatAxes(axB,cfg)

      drawnow
    placeFigure2PanelLabels(fig,[axA axB],{'A','B'})

    exportPublicationFigure(fig,cfg, ...
        'FigS_GSA_Nested_Sample_Convergence')
end

%% ========================================================================
% LOCAL FUNCTIONS: COMMON PLOTTING, RANKS AND NOTES
% ========================================================================

function drawHorizontalInterval(ax,y,lo,hi)
    if ~(isfinite(lo) && isfinite(hi))
        return
    end
    if hi<lo
        tmp = lo;
        lo = hi;
        hi = tmp;
    end
    lineColor = [0.05 0.05 0.05];
    plot(ax,[lo hi],[y y],'-','Color',lineColor,'LineWidth',1.1)
    plot(ax,[lo lo],[y-0.11 y+0.11],'-','Color',lineColor, ...
        'LineWidth',1.1)
    plot(ax,[hi hi],[y-0.11 y+0.11],'-','Color',lineColor, ...
        'LineWidth',1.1)
end

function formatAxes(ax,cfg)
    set(ax,'FontName',cfg.figureFont,'FontSize',15, ...
        'LineWidth',1.25,'TickDir','out','Box','off')
    grid(ax,'off')
end

function exportPublicationFigure(fig,cfg,fileStem)
    pngFile = fullfile(cfg.outputDir,[fileStem '.png']);
    pdfFile = fullfile(cfg.outputDir,[fileStem '.pdf']);
    if exist('exportgraphics','file')==2
        exportgraphics(fig,pngFile,'Resolution',cfg.exportResolution, ...
            'BackgroundColor','white')
        exportgraphics(fig,pdfFile,'ContentType','vector', ...
            'BackgroundColor','white')
    else
        print(fig,pngFile,'-dpng',sprintf('-r%d',cfg.exportResolution))
        print(fig,pdfFile,'-dpdf','-painters')
    end
end

function rho = localSpearman(x,y)
    x = x(:);
    y = y(:);
    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);
    if numel(x)<2
        rho = NaN;
        return
    end
    rx = localTiedRank(x);
    ry = localTiedRank(y);
    if std(rx)==0 || std(ry)==0
        rho = NaN;
        return
    end
    C = corrcoef(rx,ry);
    rho = C(1,2);
end

function ranks = localTiedRank(x)
    [sortedX,order] = sort(x);
    n = numel(x);
    sortedRanks = nan(n,1);
    first = 1;
    while first<=n
        last = first;
        while last<n && sortedX(last+1)==sortedX(first)
            last = last+1;
        end
        sortedRanks(first:last) = mean(first:last);
        first = last+1;
    end
    ranks = nan(n,1);
    ranks(order) = sortedRanks;
end

function writePostprocessingNotes(matPath,primary,rLog,cfg,rho)
    notesFile = fullfile(cfg.outputDir,'GSA_postprocessing_notes.txt');
    fid = fopen(notesFile,'w');
    if fid<0
        warning('Could not create the post-processing notes file.')
        return
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid,'Selected saved MAT file: %s\n',matPath);
    fprintf(fid,'No ODE simulations were rerun.\n');
    fprintf(fid,['No model equation, parameter bound, sampling rule, ' ...
        'pulse schedule, output definition or Sobol estimator was changed.\n']);
    fprintf(fid,'Nested convergence sizes: ');
    fprintf(fid,'%d ',cfg.sampleSizes);
    fprintf(fid,'\n');
    fprintf(fid,['Convergence interpretation: nested prefixes of one ' ...
        'scrambled Sobol design; numerical stability diagnostic only, ' ...
        'not independent-scramble uncertainty quantification.\n']);
    fprintf(fid,['Saved error bars: paired-row percentile-bootstrap ' ...
        'resampling intervals from the original saved design; not ' ...
        'biological, clinical, patient-level or inferential CIs.\n']);
    fprintf(fid,['Display rule: raw finite-sample Sobol estimates remain ' ...
        'in CSV output; negative estimates are shown as zero in figures.\n']);
    fprintf(fid,'Primary-versus-r-log ST rank Spearman rho: %.10g\n',rho);
    fprintf(fid,'Primary sum(S1): %.10g\n',primary.sobol.sumS1);
    fprintf(fid,'Primary sum(ST): %.10g\n',primary.sobol.sumST);
    fprintf(fid,'r-log robustness sum(S1): %.10g\n',rLog.sobol.sumS1);
    fprintf(fid,'r-log robustness sum(ST): %.10g\n',rLog.sobol.sumST);
       fprintf(fid,['Panel letters A and B appear in main figure and the ' ...
        'supplementary diagnostic and convergence figures, and are ' ...
        'positioned in the outer figure margin.\n']);
end
