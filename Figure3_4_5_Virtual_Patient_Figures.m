% FIGURE-ONLY SCRIPT
%
% Required input to run this code:
%   VP_MMC_continuous_independent_p8_dual_rAlpha_results.mat
% generated in PAPER mode by
%   Figure3_4_5_Virtual_Patient_Analyses.m
%
% Outputs (vector PDF; any rasterized content is exported at resolution of 1200 dpi):
%   FigVP01_Cumulative_Ever_Below_Detection.pdf
%   FigVP02_Continuous_Total_and_Ts_Burden_Refined.pdf
%   FigVP03_rAlpha_Through_0p5_2D.pdf
%
% Figure definitions:
%   * Figure 1 begins at day 30. Baseline day 0 is not counted.
%   * At month j, a patient is counted if the tumor has been below the
%     1-mm-equivalent detection threshold at least once during months 1:j.
%   * Baseline group membership never changes during follow-up.
%   * Figure 2 uses the stored continuous-ODE trajectories and contains no
%     one-cell absorbing rule or reporting-stage reset.
%   * Figure 3 uses only the prespecified r,alpha <= 0.5 day^-1 result.
%     Filled contours provide standard piecewise-linear visualization of
%     the simulated grid; no new parameter values or ODE solutions are
%     created. Both panels use one common color scale and one shared
%     colorbar. The two dashed lines mark differentiation-range boundaries,
%     not outcome thresholds. 

clearvars
close all
clc

%% Locate and load the completed PAPER results
resultName = 'VP_MMC_continuous_independent_p8_dual_rAlpha_results.mat';
candidateFiles = { ...
    fullfile(pwd,resultName), ...
    fullfile(pwd,'VP_MMC_Continuous_IndependentP8_Dual_rAlpha_Outputs', ...
        resultName), ...
    fullfile(pwd,'VP_MMC_continuous_independent_p8_dual_rAlpha_results(1).mat')};

matFile = '';
for j = 1:numel(candidateFiles)
    if exist(candidateFiles{j},'file')==2
        matFile = candidateFiles{j};
        break
    end
end
if isempty(matFile)
    [selectedFile,selectedPath] = uigetfile('*.mat', ...
        'Select the completed PAPER result MAT file');
    if isequal(selectedFile,0)
        error('No result MAT file was selected.')
    end
    matFile = fullfile(selectedPath,selectedFile);
end

loaded = load(matFile,'settings','sizeAnalysis','p8Analysis', ...
    'persistenceAnalysis','rAlphaPrimary');
requiredVariables = {'settings','sizeAnalysis','p8Analysis', ...
    'persistenceAnalysis','rAlphaPrimary'};
for j = 1:numel(requiredVariables)
    if ~isfield(loaded,requiredVariables{j})
        error('The selected MAT file does not contain %s.', ...
            requiredVariables{j})
    end
end

settings = loaded.settings;
sizeAnalysis = loaded.sizeAnalysis;
p8Analysis = loaded.p8Analysis;
persistenceAnalysis = loaded.persistenceAnalysis;
rAlphaPrimary = loaded.rAlphaPrimary;

validatePaperResults(settings,sizeAnalysis,p8Analysis, ...
    persistenceAnalysis,rAlphaPrimary)

[matFolder,~,~] = fileparts(matFile);
outputDir = fullfile(matFolder,'VP_MMC_Publication_Figures_Only');
if exist(outputDir,'dir')~=7
    mkdir(outputDir)
end

fprintf('\nCreating three figures from the completed PAPER results ...\n')

% Figure 1: cumulative post-baseline attainment of below-detection burden
fig1 = makeCumulativeBelowDetectionFigure( ...
    sizeAnalysis,p8Analysis,settings);
installFigureAspectLock(fig1,15.8/7.4)
exportPublicationPDF(fig1,fullfile(outputDir, ...
    'FigVP01_Cumulative_Ever_Below_Detection.pdf'))

% Figure 2: temporal continuous total and Ts burdens
fig2 = makeTemporalPersistenceFigure(persistenceAnalysis,settings);
installFigureAspectLock(fig2,15.8/7.5)
exportPublicationPDF(fig2,fullfile(outputDir, ...
    'FigVP02_Continuous_Total_and_Ts_Burden_Refined.pdf'))

% Figure 3: two-dimensional r-alpha maps through 0.5 day^-1
fig3 = makeRAlpha2DFigure(rAlphaPrimary,settings);
installFigureAspectLock(fig3,10.5/4.75)
exportPublicationPDF(fig3,fullfile(outputDir, ...
    'FigVP03_rAlpha_Through_0p5_2D.pdf'))

fprintf('\nAll three figure-only outputs were created successfully.\n')
fprintf('Output directory: %s\n',outputDir)

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function validatePaperResults(settings,sizeAnalysis,p8Analysis, ...
        persistenceAnalysis,rAlphaPrimary)
    if ~isfield(settings,'runMode') || ~strcmpi(settings.runMode,'paper')
        error(['The selected MAT file is not a PAPER-mode result. ' ...
            'No figure was created.'])
    end
    if settings.nSizePerGroup~=1000 || settings.nP8PerGroup~=1000 || ...
            settings.nPersistence~=1000 || settings.nBackgrounds~=50
        error('The selected result does not contain the agreed PAPER cohorts.')
    end
    if numel(sizeAnalysis.monthlyDetectableStatus)~=2 || ...
            numel(p8Analysis.monthlyDetectableStatus)~=3
        error('Stored monthly detectability matrices are incomplete.')
    end
    for g = 1:2
        status = sizeAnalysis.monthlyDetectableStatus{g};
        if ~isequal(size(status),[1000 12])
            error('Initial-size monthly status matrix %d is not 1000 x 12.',g)
        end
    end
    for g = 1:3
        status = p8Analysis.monthlyDetectableStatus{g};
        if ~isequal(size(status),[1000 12])
            error('Initial-Ts-proportion status matrix %d is not 1000 x 12.',g)
        end
    end
    if ~isequal(persistenceAnalysis.visitTimes,0:30:360)
        error('The temporal analysis does not contain baseline plus 12 months.')
    end
    if abs(rAlphaPrimary.rGrid(1)-0.00625)>1e-12 || ...
            abs(rAlphaPrimary.rGrid(end)-0.5)>1e-12 || ...
            abs(rAlphaPrimary.alphaGrid(1)-1/280)>1e-12 || ...
            abs(rAlphaPrimary.alphaGrid(end)-0.5)>1e-12
        error('Figure 3 requires the agreed r-alpha domain through 0.5 day^-1.')
    end
    if any(~isfinite(rAlphaPrimary.meanTotalTumorCells(:))) || ...
            any(~isfinite(rAlphaPrimary.meanTsCells(:)))
        error('Figure-3 response matrices contain nonfinite values.')
    end
end

function fig = makeCumulativeBelowDetectionFigure( ...
        sizeAnalysis,p8Analysis,settings)
    days = settings.monthDays:settings.monthDays:settings.tEnd;
    sizeColors = [0.0000 0.4470 0.7410;0.8500 0.3250 0.0980];
    p8Colors = [0.00 0.52 0.32;0.49 0.18 0.56;0.93 0.55 0.02];

    cumulativeSize = zeros(2,numel(days));
    for g = 1:2
        belowAtVisit = ~sizeAnalysis.monthlyDetectableStatus{g};
        cumulativeSize(g,:) = mean(cumsum(belowAtVisit,2)>0,1);
    end
    cumulativeP8 = zeros(3,numel(days));
    for g = 1:3
        belowAtVisit = ~p8Analysis.monthlyDetectableStatus{g};
        cumulativeP8(g,:) = mean(cumsum(belowAtVisit,2)>0,1);
    end
    if any(diff(cumulativeSize,1,2)<-1e-14,'all') || ...
            any(diff(cumulativeP8,1,2)<-1e-14,'all')
        error('A cumulative-ever-below curve decreased unexpectedly.')
    end

    fig = figure('Color','w','Name','Cumulative below-detection burden', ...
        'Units','inches','Position',[0.3 0.3 15.8 7.4], ...
        'Renderer','painters');
    axA = axes(fig,'Units','normalized','Position',[0.075 0.15 0.395 0.74]);
    axB = axes(fig,'Units','normalized','Position',[0.565 0.15 0.395 0.74]);
    hold(axA,'on'); hold(axB,'on')

    for g = 1:2
        stairs(axA,days,cumulativeSize(g,:),'-','LineWidth',2.7, ...
            'Color',sizeColors(g,:));
    end
    title(axA,'Initial tumor size')
    ylabel(axA,{'Cumulative proportion reaching', ...
        'below-detection burden'})
    legend(axA,{sprintf('<3 cm (n=%d)', ...
        sizeAnalysis.resultsTable.N_Small(1)), ...
        sprintf('>=3 cm (n=%d)', ...
        sizeAnalysis.resultsTable.N_Large(1))}, ...
        'Location','southeast','Box','off','Interpreter','none', ...
        'FontName','Arial','FontSize',12.5,'FontWeight','bold')

    thirdNames = {'Lower third','Middle third','Upper third'};
    p8Legend = cell(3,1);
    for g = 1:3
        stairs(axB,days,cumulativeP8(g,:),'-','LineWidth',2.5, ...
            'Color',p8Colors(g,:));
        p8Legend{g} = sprintf( ...
            '%s of initial T_s fraction (p_8) (n=%d)', ...
            thirdNames{g},p8Analysis.groupSizes(g));
    end
    title(axB,'Initial tumor stem-cell fraction (p_8)','Interpreter','tex')
    legend(axB,p8Legend,'Location','southeast','Box','off', ...
        'Interpreter','tex','FontName','Arial','FontSize',11.5, ...
        'FontWeight','bold')

    styleMonthlyAxes(axA,settings)
    styleMonthlyAxes(axB,settings)
    linkaxes([axA axB],'xy')
    addSharedXLabel(fig,'Follow-up time [days]')
    addPanelLetter(fig,[0.010 0.900 0.035 0.05],'A')
    addPanelLetter(fig,[0.500 0.900 0.035 0.05],'B')
end

function styleMonthlyAxes(ax,settings)
    xlim(ax,[settings.monthDays settings.tEnd])
    xticks(ax,[30 60 120 180 240 300 360])
    ylim(ax,[0.50 0.85])
    yticks(ax,0.50:0.05:0.85)
    style2DAxes(ax)
end

function fig = makeTemporalPersistenceFigure(analysis,settings)
    fig = figure('Color','w','Name','Continuous temporal tumor burden', ...
        'Units','inches','Position',[0.3 0.3 15.8 7.5], ...
        'Renderer','painters');

    axA = axes(fig,'Units','normalized','Position',[0.075 0.15 0.395 0.73]);
    axB = axes(fig,'Units','normalized','Position',[0.565 0.15 0.395 0.73]);
    hold(axA,'on'); hold(axB,'on')

    colors = [0.8500 0.3250 0.0980;0.0000 0.4470 0.7410];
    labels = cell(1,2);
    lineHandles = gobjects(2,1);
    for g = 1:2
        summary = analysis.groupSummary(g);
        labels{g} = sprintf('%s (n=%d)', ...
            conciseOutcomeLabel(summary.label),summary.n);
        lineHandles(g) = plotMedianIQR(axA,analysis.visitTimes, ...
            summary.medianTotal,summary.q25Total,summary.q75Total, ...
            colors(g,:));
        plotMedianIQR(axB,analysis.visitTimes,summary.medianTs, ...
            summary.q25Ts,summary.q75Ts,colors(g,:));
    end

    thresholdY = log10(settings.detectionLimit+1);
    yline(axA,thresholdY,'k:','LineWidth',1.7,'HandleVisibility','off');
    text(axA,settings.tEnd-5,thresholdY+0.08,'1-mm detection threshold', ...
        'HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'FontName','Arial','FontSize',14,'FontWeight','bold', ...
        'Color',[0.12 0.12 0.12]);

    title(axA,'Total tumor burden')
    legend(axA,lineHandles,labels,'Location','northwest','Box','off', ...
        'Interpreter','none','FontName','Arial','FontSize',15, ...
        'FontWeight','bold')

    title(axB,'Tumor stem-cell burden')

    styleTemporalAxes(axA,settings)
    styleTemporalAxes(axB,settings)
    sharedYMaximum = sharedTemporalYMaximum(analysis,thresholdY);
    ylim(axA,[0 sharedYMaximum])
    ylim(axB,[0 sharedYMaximum])
    yticks(axA,0:2:sharedYMaximum)
    yticks(axB,0:2:sharedYMaximum)
    linkaxes([axA axB],'xy')
    addSharedXLabel(fig,'Follow-up time [days]')

    ylabel(axA,{'Modeled tumor burden', ...
        'log_{10}(cell count + 1) [cells]'}, ...
        'Interpreter','tex','FontName','Arial','FontSize',15, ...
        'FontWeight','bold')
    addPanelLetter(fig,[0.010 0.890 0.035 0.05],'A')
    addPanelLetter(fig,[0.500 0.890 0.035 0.05],'B')
end

function sharedMaximum = sharedTemporalYMaximum(analysis,thresholdY)
    upperValues = thresholdY;
    for g = 1:numel(analysis.groupSummary)
        summary = analysis.groupSummary(g);
        upperValues = [upperValues, ... %#ok<AGROW>
            log10(1+reshape(summary.q75Total,1,[])), ...
            log10(1+reshape(summary.q75Ts,1,[]))]; %#ok<AGROW>
    end
    sharedMaximum = max(2,2*ceil(max(upperValues)/2));
end

function label = conciseOutcomeLabel(fullLabel)
    if contains(string(fullLabel),'Detectable tumor')
        label = 'Detectable at day 360';
    else
        label = 'Below detection at day 360';
    end
end

function lineHandle = plotMedianIQR(ax,x,medianValues,q25Values,q75Values,color)
    x = reshape(x,1,[]);
    transformedMedian = log10(1+reshape(medianValues,1,[]));
    transformedLower = log10(1+reshape(q25Values,1,[]));
    transformedUpper = log10(1+reshape(q75Values,1,[]));
    fill(ax,[x fliplr(x)],[transformedLower fliplr(transformedUpper)], ...
        color,'FaceAlpha',0.16,'EdgeColor','none','HandleVisibility','off');
    lineHandle = plot(ax,x,transformedMedian,'-o','Color',color, ...
        'LineWidth',2.5,'MarkerSize',5.2,'MarkerFaceColor',color);
end

function styleTemporalAxes(ax,settings)
    xlim(ax,[0 settings.tEnd])
    xticks(ax,0:60:settings.tEnd)
    style2DAxes(ax)
end

function fig = makeRAlpha2DFigure(result,settings)
    r = reshape(result.rGrid,1,[]);
    alpha = reshape(result.alphaGrid,1,[]);
    logAlpha = log10(alpha);
    totalCells = result.meanTotalTumorCells;
    tsCells = result.meanTsCells;
    sharedMaximum = max([totalCells(:);tsCells(:)]);
    if sharedMaximum<=0; sharedMaximum=1; end

    if ~isequal(size(totalCells),[numel(alpha) numel(r)]) || ...
            ~isequal(size(tsCells),[numel(alpha) numel(r)])
        error('The r-alpha matrices do not match the stored grids.')
    end

    fig = figure('Color','w','Name','Two-dimensional r-alpha response', ...
        'Units','inches','Position',[0.25 0.25 10.5 4.75], ...
        'Renderer','painters');
    axA = axes(fig,'Units','normalized','Position',[0.090 0.14 0.355 0.75]);
    axB = axes(fig,'Units','normalized','Position',[0.505 0.14 0.355 0.75]);
    hold(axA,'on'); hold(axB,'on')

    plotResponseMap(axA,r,logAlpha,totalCells,sharedMaximum, ...
        'Mean total tumor-cell count',settings,true)
    plotResponseMap(axB,r,logAlpha,tsCells,sharedMaximum, ...
        'Mean T_s-compartment count',settings,false)

    % One colorbar is used because both panels encode the same dimensional
    % quantity (mean modeled cell count) on exactly the same color limits.
    cb = colorbar(axB,'Location','eastoutside');
    cb.Label.String = 'Mean modeled cell count at day 360 [cells]';
    styleColorbar(cb)
    cb.Position = [0.890 0.14 0.018 0.75];
    axA.Position = [0.090 0.14 0.355 0.75];
    axB.Position = [0.505 0.14 0.355 0.75];

    colormap(fig,vividMap(256))
    linkaxes([axA axB],'xy')
    yticks(axB,[])
addSharedXLabel(fig,'Tumor growth rate r [day^{-1}]',13.5,0.015)

    ylabel(axA,'Differentiation-like transition rate \alpha [day^{-1}]', ...
        'Interpreter','tex','FontName','Arial','FontSize',13, ...
        'FontWeight','bold')
addPanelLetter(fig,[0.040 0.910 0.035 0.05],'A')
addPanelLetter(fig,[0.455 0.910 0.035 0.05],'B')
    drawnow
    updateDomainTextBackgrounds(axA)
    setappdata(fig,'DomainTextAxes',axA)
end

function plotResponseMap(ax,r,logAlpha,zCells,sharedMaximum,titleText, ...
        settings,showDomainLabels)
    % contourf linearly interpolates between the actual simulated grid
    % coordinates. It smooths display only and does not create new outcomes.
    contourf(ax,r,logAlpha,zCells,40,'LineColor','none');
    caxis(ax,[0 sharedMaximum])

    boundary28 = log10(settings.normalAlphaLower);
    boundary4 = log10(settings.normalAlphaUpper);
    yline(ax,boundary28,'k--','LineWidth',2.3,'HandleVisibility','off');
    yline(ax,boundary4,'k--','LineWidth',2.3,'HandleVisibility','off');
    xlim(ax,[r(1) r(end)])
    ylim(ax,[logAlpha(1) logAlpha(end)])

    if showDomainLabels
        xText = r(end)-0.018*(r(end)-r(1));
        yCancer = mean([logAlpha(1),boundary28]);
        yNormal = mean([boundary28,boundary4]);
        yHigh = mean([boundary4,logAlpha(end)]);
        addDomainText(ax,xText,yCancer, ...
            {'Hypothetical cancer-associated', ...
            'loss-of-differentiation range (>28 days)'});
        addDomainText(ax,xText,yNormal, ...
            {'Normal-tissue reference','4-28 days'});
        addDomainText(ax,xText,yHigh, ...
            {'Hypothetical faster range','<4 days'});
    end

    title(ax,titleText,'Interpreter','tex')
    xticks(ax,[r(1) 0.1 0.2 0.3 0.4 0.5])
    xticklabels(ax,{'6.25\times10^{-3}','0.1','0.2','0.3','0.4','0.5'})
    alphaTicks = [1/280 1/28 1/4 0.5];
    yticks(ax,log10(alphaTicks))
    yticklabels(ax,{'1/280','1/28','1/4','0.5'})
    style2DAxes(ax)
    ax.XTickLabelRotation = 0;
    ax.FontSize = 12;
    ax.Title.FontSize = 15;
end

function addDomainText(ax,xValue,yValue,labelLines)
    text(ax,xValue,yValue,labelLines, ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'FontName','Arial','FontSize',14,'FontWeight','bold', ...
        'Color',[0.05 0.05 0.05],'BackgroundColor','none', ...
        'Margin',2,'Clipping','on','Tag','DomainRangeLabel');
end

function updateDomainTextBackgrounds(ax)
    if ~isgraphics(ax,'axes')
        return
    end
    delete(findall(ax,'Type','Patch','Tag','DomainRangeBackground'))
    labelTexts = findall(ax,'Type','Text','Tag','DomainRangeLabel');
    xPadding = 0.006*diff(xlim(ax));
    yPadding = 0.010*diff(ylim(ax));
    for j = 1:numel(labelTexts)
        labelExtent = labelTexts(j).Extent;
        xBox = [labelExtent(1)-xPadding, ...
            labelExtent(1)+labelExtent(3)+xPadding, ...
            labelExtent(1)+labelExtent(3)+xPadding, ...
            labelExtent(1)-xPadding];
        yBox = [labelExtent(2)-yPadding, ...
            labelExtent(2)-yPadding, ...
            labelExtent(2)+labelExtent(4)+yPadding, ...
            labelExtent(2)+labelExtent(4)+yPadding];
        patch(ax,xBox,yBox,[1 1 1], ...
            'FaceAlpha',0.55,'EdgeColor','none', ...
            'HandleVisibility','off','Clipping','on', ...
            'HitTest','off','PickableParts','none', ...
            'Tag','DomainRangeBackground');
    end
    for j = 1:numel(labelTexts)
        uistack(labelTexts(j),'top')
    end
end

function installFigureAspectLock(fig,aspectRatio)
    if ~isgraphics(fig,'figure')
        return
    end
    fig.WindowState = 'normal';
    fig.Units = 'pixels';
    positionValue = fig.Position;
    if positionValue(3)/positionValue(4) < aspectRatio
        positionValue(4) = positionValue(3)/aspectRatio;
    else
        positionValue(3) = positionValue(4)*aspectRatio;
    end
    fig.Position = positionValue;
    axesHandles = findall(fig,'Type','axes');
    for j = 1:numel(axesHandles)
        axesHandles(j).Toolbar = [];
    end
    setappdata(fig,'LockedAspectRatio',aspectRatio)
    setappdata(fig,'PreviousLockedSize',positionValue(3:4))
    setappdata(fig,'AspectRatioUpdateBusy',false)
    drawnow
    refreshFigureDependentLayout(fig)
    fig.SizeChangedFcn = @maintainFigureAspectRatio;
end

function maintainFigureAspectRatio(fig,~)
    if ~isgraphics(fig,'figure')
        return
    end
    isBusy = getappdata(fig,'AspectRatioUpdateBusy');
    if ~isempty(isBusy) && isBusy
        return
    end
    setappdata(fig,'AspectRatioUpdateBusy',true)
    cleanupObject = onCleanup(@()setappdata( ...
        fig,'AspectRatioUpdateBusy',false)); %#ok<NASGU>

    aspectRatio = getappdata(fig,'LockedAspectRatio');
    previousSize = getappdata(fig,'PreviousLockedSize');
    positionValue = fig.Position;
    currentSize = positionValue(3:4);
    if isempty(previousSize)
        previousSize = currentSize;
    end

    relativeWidthChange = abs(currentSize(1)-previousSize(1))/ ...
        max(previousSize(1),1);
    relativeHeightChange = abs(currentSize(2)-previousSize(2))/ ...
        max(previousSize(2),1);
    if relativeWidthChange >= relativeHeightChange
        positionValue(4) = positionValue(3)/aspectRatio;
    else
        positionValue(3) = positionValue(4)*aspectRatio;
    end

    if any(abs(fig.Position(3:4)-positionValue(3:4))>0.5)
        fig.Position = positionValue;
    end
    setappdata(fig,'PreviousLockedSize',positionValue(3:4))
    drawnow limitrate nocallbacks
    refreshFigureDependentLayout(fig)
end

function refreshFigureDependentLayout(fig)
    if isappdata(fig,'DomainTextAxes')
        domainAxes = getappdata(fig,'DomainTextAxes');
        if isgraphics(domainAxes,'axes')
            updateDomainTextBackgrounds(domainAxes)
        end
    end
end

function style2DAxes(ax)
    grid(ax,'on')
    ax.GridAlpha = 0.12;
    ax.MinorGridAlpha = 0.06;
    ax.FontName = 'Arial';
    ax.FontSize = 13;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.1;
    ax.TickDir = 'out';
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    box(ax,'off')
    ax.XLabel.FontSize = 14;
    ax.YLabel.FontSize = 14;
    ax.Title.FontSize = 15;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontWeight = 'bold';
    ax.Title.FontWeight = 'bold';
end

function styleColorbar(cb)
    cb.FontName = 'Arial';
    cb.FontSize = 12.5;
    cb.FontWeight = 'bold';
    cb.Label.FontName = 'Arial';
    cb.Label.FontSize = 13.5;
    cb.Label.FontWeight = 'bold';
end

function map = vividMap(n)
    if exist('turbo','file')==2 || exist('turbo','builtin')==5
        map = turbo(n);
    else
        map = jet(n);
    end
end

function addFigureTitle(fig,titleText,yPosition,fontSize)
    annotation(fig,'textbox',[0.08 yPosition-0.025 0.84 0.045], ...
        'String',titleText,'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Arial','FontSize',fontSize,'FontWeight','bold', ...
        'Interpreter','tex');
end

function addSharedXLabel(fig,labelText,fontSize,yPosition)
    if nargin<3
        fontSize = 15;
    end
    if nargin<4
        yPosition = 0.035;
    end

    annotation(fig,'textbox',[0.35 yPosition 0.30 0.045], ...
        'String',labelText,'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Arial','FontSize',fontSize,'FontWeight','bold');
end

function addPanelLetter(fig,positionValue,letter)
    annotation(fig,'textbox',positionValue,'String',letter, ...
        'EdgeColor','none','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontName','Arial', ...
        'FontSize',21,'FontWeight','bold','Color',[0.05 0.05 0.05]);
end

function exportPublicationPDF(fig,outputFile)
    if exist(outputFile,'file')==2
        delete(outputFile)
    end
    exportgraphics(fig,outputFile,'ContentType','vector', ...
        'BackgroundColor','white','Resolution',1200);
    fprintf('Saved: %s\n',outputFile)
end
