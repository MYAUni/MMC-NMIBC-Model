% Generate Figures 7 and 8 from saved paper-mode analysis results.
% Figure 7 draws the mechanism network and saved transition counts.
% Figure 8 computes the selected patient's continuous parameter path
% (41 ODE simulations) and checks its endpoints before saving Figure 8.
%
% Place this script beside the file MMC_mechanism_failure_final_paper_results.
% Run the script to save figures in the results folder's
% publication_figures subfolder.

clearvars
clc
close all

%% USER SETTINGS -- figure export only
pngDPI = 1200;
maximumPngDimensionPixels = 4800;
keepFiguresOpen = true;
patientFigureViewAngles = [125 24]; % [azimuth elevation], degrees

scriptPath = [mfilename('fullpath') '.m'];
scriptFolder = fileparts(scriptPath);
inputFolder = locatePaperResultFolder(scriptFolder);
outputFolder = fullfile(inputFolder,'publication_figures');
if exist(outputFolder,'dir')~=7
    mkdir(outputFolder)
end

resultPath = locateOneFile(inputFolder, ...
    'MMC_mechanism_failure_results.mat', ...
    'MMC_mechanism_failure_results*.mat');

fprintf('\n===============================================================\n')
fprintf('MMC MECHANISM-OF-FAILURE FIGURES\n')
fprintf('===============================================================\n')
fprintf('Loading completed results: %s\n',resultPath)
loaded = load(resultPath);
requiredTopLevel = {'Results','patients','perturbationPlans','spec', ...
    'fixed','mechanisms','settings'};
for q = 1:numel(requiredTopLevel)
    if ~isfield(loaded,requiredTopLevel{q})
        error('Saved results are missing top-level variable %s.', ...
            requiredTopLevel{q})
    end
end

Results = loaded.Results;
patients = loaded.patients;
plans = loaded.perturbationPlans;
spec = loaded.spec;
fixed = loaded.fixed;
mechanisms = loaded.mechanisms;
settings = loaded.settings;

validateLoadedPaperResults( ...
    Results,patients,plans,mechanisms,spec,fixed,settings)
simulation = Results.Simulation;
summaryTable = Results.MechanismSummary;
selectedCase = Results.SelectedCases;

settings.outputDir = outputFolder;
settings.figureDPI = pngDPI;
settings.keepFiguresOpen = keepFiguresOpen;
settings.makeFigures = true;
settings.makePatientFigures = true;

odeOpt = odeset('AbsTol',settings.AbsTol,'RelTol',settings.RelTol, ...
    'MaxStep',settings.MaxStep,'NonNegative',1:7);

% FIGURE 7- draw the network and plot the saved transition counts
fprintf('\nDrawing Figure 9 from the saved numerical results...\n')
[fig9,figure9Audit] = makeMechanismFigure(summaryTable,settings);
figure9Stem = 'Fig09_Mechanism_Network_and_Transition_Counts';
figure9FigPath = fullfile(outputFolder,[figure9Stem '.fig']);
savefig(fig9,figure9FigPath)
validateEditableFigureFile(figure9FigPath,2)
figure9Export = exportMechanismFigure( ...
    fig9,outputFolder,figure9Stem,pngDPI,maximumPngDimensionPixels);
if ~keepFiguresOpen
    close(fig9)
end

% FIGURE 8-selected VP 127/M4 path only
fprintf('\nRebuilding the selected paper patient figure...\n')
patientIndex = selectedCase.VirtualPatient(1);
mechanismNumber = selectedCase.MechanismNumber(1);
sampleIndex = selectedCase.FailureSample(1);
fprintf('Selected paper case: VP %d, %s, draw %d (%d/%d failures).\n', ...
    patientIndex,mechanisms(mechanismNumber).id,sampleIndex, ...
    selectedCase.FailureDraws(1),selectedCase.ValidDraws(1))

patientScreen = buildPatientMechanismScreen( ...
    patientIndex,simulation,mechanisms,settings);
[Surface,pathTable,endpointTable,surfaceAudit] = computeCaseSurface( ...
    patients{patientIndex},plans(mechanismNumber),sampleIndex,spec, ...
    fixed,settings,odeOpt,simulation.ReferenceFinalTumor(patientIndex), ...
    simulation.PerturbedFinalTumor( ...
    patientIndex,mechanismNumber,sampleIndex));

[fig10,figure10Audit] = makePatientCaseFigure( ...
    mechanismNumber,patientScreen,Surface,patients{patientIndex}, ...
    patientFigureViewAngles,settings);
figure10Stem = sprintf('Fig10_Selected_Patient_VP%04d_%s', ...
    patientIndex,mechanisms(mechanismNumber).id);
writetable(patientScreen,fullfile(outputFolder, ...
    [figure10Stem '_patient_mechanism_screen.csv']))
writetable(pathTable,fullfile(outputFolder, ...
    [figure10Stem '_continuous_parameter_path.csv']))
writetable(endpointTable,fullfile(outputFolder, ...
    [figure10Stem '_reference_and_failure_parameters.csv']))
save(fullfile(outputFolder,[figure10Stem '_surface_data.mat']), ...
    'Surface','-v7.3')
figure10FigPath = fullfile(outputFolder,[figure10Stem '.fig']);
savefig(fig10,figure10FigPath)
validateEditableFigureFile(figure10FigPath,2)
figure10Export = exportMemorySafeFigure( ...
    fig10,outputFolder,figure10Stem,pngDPI, ...
    maximumPngDimensionPixels,'image');
if ~keepFiguresOpen
    close(fig10)
end

CaseRank = ones(height(surfaceAudit)+height(figure10Audit),1);
CaseID = repmat(string(figure10Stem),numel(CaseRank),1);
selectedAudit = [surfaceAudit;figure10Audit];
selectedAudit = addvars(selectedAudit,CaseRank,CaseID,'Before',1);
writetable(selectedAudit,fullfile(outputFolder, ...
    'Selected_case_figure_and_surface_audit.csv'))

Figure = [repmat("Figure 9",height(figure9Audit),1); ...
    repmat("Figure 10",height(selectedAudit),1)];
contentAudit = [figure9Audit;selectedAudit(:,{'Check','Passed','Detail'})];
contentAudit = addvars(contentAudit,Figure,'Before',1);
writetable(contentAudit,fullfile(outputFolder, ...
    'Figure_content_audit.csv'))

exportAudit = [figure9Export;figure10Export];
writetable(exportAudit,fullfile(outputFolder,'Figure_export_audit.csv'))

fprintf('\n===============================================================\n')
fprintf('PAPER FIGURES COMPLETED WITHOUT PRIMARY-ANALYSIS RECOMPUTATION\n')
fprintf('===============================================================\n')
fprintf('Figure 9: %s\n',fullfile(outputFolder,[figure9Stem '.png']))
fprintf('Figure 10: %s\n',fullfile(outputFolder,[figure10Stem '.png']))
fprintf('Editable .fig and PDF versions are in the same folder.\n')
fprintf('Primary simulations rerun: 0\nBootstrap replicates rerun: 0\n')
if ispc
    winopen(outputFolder)
end

% ========================================================================
% INPUT LOCATION AND PAPER-RESULT VALIDATION
% ========================================================================

function folder = locatePaperResultFolder(scriptFolder)
candidate = fullfile(scriptFolder, ...
    'MMC_mechanism_failure_final_paper_results');
if exist(fullfile(candidate,'MMC_mechanism_failure_results.mat'),'file')==2
    folder = candidate;
elseif exist(fullfile(scriptFolder, ...
        'MMC_mechanism_failure_results.mat'),'file')==2
    folder = scriptFolder;
else
    exact = dir(fullfile(scriptFolder, ...
        'MMC_mechanism_failure_results*.mat'));
    if numel(exact)==1
        folder = scriptFolder;
    else
        error(['Could not identify the paper-results folder. Put this ' ...
            'script in that folder or immediately beside it.'])
    end
end
end

function path = locateOneFile(folder,exactName,wildcard)
exactPath = fullfile(folder,exactName);
if exist(exactPath,'file')==2
    path = exactPath;
    return
end
matches = dir(fullfile(folder,wildcard));
matches = matches(~[matches.isdir]);
if numel(matches)~=1
    error('Expected exactly one file matching %s in %s; found %d.', ...
        wildcard,folder,numel(matches))
end
path = fullfile(matches(1).folder,matches(1).name);
end

function validateLoadedPaperResults( ...
        Results,patients,plans,mechanisms,spec,fixed,settings)
requiredResults = {'Simulation','MechanismSummary','PatientMechanismSummary', ...
    'LongTable','PatientOverlap','Bootstrap','CohortStability', ...
    'SolverAudit','CandidateCases','SelectedCases'};
for q = 1:numel(requiredResults)
    if ~isfield(Results,requiredResults{q})
        error('Saved Results structure is missing %s.',requiredResults{q})
    end
end
if ~strcmp(settings.runMode,'paper') || settings.nPatients~=500 || ...
        settings.nPerturbations~=20
    error('This figure-only code requires the completed 500 x 20 paper run.')
end
if ~strcmp(settings.codeVersion,'MMC_MF_FINAL_CORRECTED_2026-08-23_B') || ...
        settings.surfaceLevels~=41 || settings.surfaceTimePoints~=181 || ...
        settings.nSelectedPatients~=1 || ...
        abs(settings.detectionLimit-23561.9449019234)>1e-9
    error('Saved paper analysis or selected-surface settings changed.')
end
if numel(patients)~=500 || numel(plans)~=5 || numel(mechanisms)~=5
    error('Saved cohort or mechanism dimensions are incorrect.')
end
expectedPlanInputs = {{'p1','p2'}, ...
    {'alpha','theta1','theta2'}, ...
    {'beta1','beta2','p3','p4'}, ...
    {'gamma','p5'}, {'eta'}};
for g = 1:5
    if ~isequal(plans(g).inputNames,expectedPlanInputs{g}) || ...
            size(plans(g).values,1)~=20
        error('Saved perturbation plan M%d does not match the specified inputs or draw count.',g)
    end
end
gammaSpec = getSpec(spec,'gamma');
p5Spec = getSpec(spec,'p5');
if ~strcmp(gammaSpec.measure,'log') || ...
        ~strcmp(p5Spec.measure,'linear') || ...
        fixed.M0~=0 || fixed.Dm0~=1 || fixed.E0~=1 || fixed.R0~=1 || ...
        fixed.p6~=1.44e-5 || fixed.p7~=3.110210655
    error('Selected-path input mapping or fixed values changed.')
end
simulation = Results.Simulation;
if ~all(simulation.ReferenceOK) || ~all(simulation.PerturbedOK(:))
    error('The loaded paper run contains an incomplete ODE simulation.')
end
if height(Results.LongTable)~=50000 || ...
        height(Results.PatientMechanismSummary)~=2500 || ...
        height(Results.PatientOverlap)~=500
    error('A saved paper-mode table has an incorrect number of rows.')
end
S = Results.MechanismSummary;
if ~isequal(S.MechanismID,string({'M1','M2','M3','M4','M5'}).') || ...
        any(S.EligibleReferencePatients~=322) || ...
        any(S.ValidPerturbations~=6440) || ...
        ~isequal(S.TransitionEvents,[33;57;21;142;0]) || ...
        ~isequal(S.AffectedPatients,[5;16;7;19;0])
    error('Saved paper mechanism summaries do not match the reconciled run.')
end
Selected = Results.SelectedCases;
if height(Selected)~=1 || Selected.VirtualPatient(1)~=127 || ...
        Selected.MechanismNumber(1)~=4 || ...
        Selected.FailureSample(1)~=16 || ...
        Selected.FailureDraws(1)~=17 || Selected.ValidDraws(1)~=20
    error('Saved primary patient case does not match VP 127/M4/draw 16.')
end
if abs(Selected.ReferenceTumorCells(1)-3462.70272649384)>1e-8 || ...
        abs(Selected.FailureTumorCells(1)-526977373.158697)>1e-4
    error('Saved selected-case endpoints changed.')
end
if settings.AbsTol~=1e-9 || settings.RelTol~=1e-6 || ...
        settings.MaxStep~=1 || settings.oneCellLimit~=1 || ...
        ~isequal(settings.treatmentTimes, ...
        [0 7 14 21 28 35 70 98 126 154 182 210]) || ...
        settings.doseDuration~=2/24
    error('Saved model or solver settings do not match the expected values.')
end
fprintf(['Loaded paper results validated: 500 references, 50,000 ' ...
    'perturbations, 322 eligible patients, 253 transitions.\n'])
end

% ========================================================================
% FIGURE 7-MECHANISM NETWORK AND TRANSITION COUNTS
% ========================================================================

function [fig,Audit] = makeMechanismFigure(S,settings)
% Draw the mechanism network and plot the saved transition counts.
fig = figure('Color','w','Units','pixels','Position',[40 40 1400 1200], ...
    'Renderer','painters','Visible','off','NumberTitle','off', ...
    'Name','Mechanism network and control-to-failure transitions');
annotation(fig,'textbox',[0.025 0.940 0.035 0.045], ...
    'String','A','EdgeColor','none','FontName',settings.figureFont, ...
    'FontSize',20,'FontWeight','bold','Color',[0.1294 0.1294 0.1294]);
annotation(fig,'textbox',[0.025 0.335 0.035 0.045], ...
    'String','B','EdgeColor','none','FontName',settings.figureFont, ...
    'FontSize',20,'FontWeight','bold','Color',[0.1294 0.1294 0.1294]);

networkAxis = axes('Parent',fig,'Units','normalized', ...
    'Position',[0.065 0.405 0.870 0.545]);
drawMechanismNetwork(networkAxis,settings.figureFont);

barAxis = axes('Parent',fig,'Units','normalized', ...
    'Position',[0.260 0.085 0.480 0.245], ...
    'FontName',settings.figureFont,'FontSize',11,'FontWeight','bold', ...
    'LineWidth',1.2,'Box','on','TickDir','in', ...
    'XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off');
hold(barAxis,'on')
counts = double(S.TransitionEvents(:));
barObject = bar(barAxis,1:numel(counts),counts,0.62, ...
    'FaceColor',[0.06 0.24 0.48],'EdgeColor',[0.02 0.10 0.24], ...
    'LineWidth',1.3);
barAxis.XLim = [0.45 5.55];
barAxis.XTick = 1:5;
barAxis.XTickLabel = cellstr(S.MechanismID);
xlabel(barAxis,'Mechanism','FontSize',11,'FontWeight','bold');
ylabel(barAxis,'Control-to-failure transitions', ...
    'FontSize',11,'FontWeight','bold');
[upperLimit,ticks] = niceCountAxis(max(counts));
barAxis.YLim = [0 upperLimit];
barAxis.YTick = ticks;

labelOffset = 0.025*upperLimit;
for g = 1:5
    if counts(g)==1
        eventWord = 'event';
    else
        eventWord = 'events';
    end
    label = sprintf('%d %s\n%d/%d VPs',counts(g),eventWord, ...
        S.AffectedPatients(g),S.EligibleReferencePatients(g));
    text(barAxis,g,counts(g)+labelOffset,label, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontName',settings.figureFont,'FontSize',9.5, ...
        'FontWeight','bold','Color',[0 0 0], ...
        'Interpreter','none','Clipping','off');
end

axesObjects = [networkAxis;barAxis];
enableEditableFigureControls(fig,axesObjects)
set(fig,'Visible','on')
drawnow

% Fit the 1400-by-1200 design to the available window without stretching.
canvas = double(fig.Position(3:4));
canvasScale = min([1 canvas./[1400 1200]]);
fig.Position(3:4) = canvasScale*[1400 1200];
scaleMechanismGraphics(fig,canvasScale)
drawnow

allStrings = collectFigureStrings(fig);
Check = ["Two editable axes created"; ...
    "Exactly five mechanism bars created"; ...
    "Paper transition counts loaded"; ...
    "Paper affected-patient counts loaded"; ...
    "No automatic data labels present"; ...
    "M1--M5 order retained"];
Passed = [numel(axesObjects)==2;numel(barObject.YData)==5; ...
    isequal(double(barObject.YData(:)),[33;57;21;142;0]); ...
    isequal(double(S.AffectedPatients(:)),[5;16;7;19;0]); ...
    ~any(startsWith(lower(allStrings),'data')); ...
    isequal(string(barAxis.XTickLabel(:)),string({'M1';'M2';'M3';'M4';'M5'}))];
Detail = [string(numel(axesObjects))+" axes"; ...
    string(numel(barObject.YData))+" bars"; ...
    "33, 57, 21, 142, 0";"5, 16, 7, 19, 0 of 322"; ...
    "legend strings checked";"M1, M2, M3, M4, M5"];
Audit = table(Check,Passed,Detail);
if ~all(Passed)
    error('Figure 9 content audit failed: %s.', ...
        char(strjoin(Check(~Passed),'; ')))
end
end

function Audit = exportMechanismFigure(fig,outputFolder,stem,dpi,maximumPixels)
% Export Figure 9 with the same proportions as its editable canvas.
pngPath = fullfile(outputFolder,[stem '.png']);
pdfPath = fullfile(outputFolder,[stem '.pdf']);
canvas = double(fig.Position(3:4));
paperSize = (maximumPixels/dpi)*canvas/max(canvas);

exportFigure = copyobj(fig,groot);
cleanup = onCleanup(@()close(exportFigure)); %#ok<NASGU>
set(exportFigure,'Visible','off','MenuBar','none','ToolBar','none');
axesObjects = findall(exportFigure,'Type','axes');
for q = 1:numel(axesObjects)
    axesObjects(q).Toolbar.Visible = 'off';
end

% Scale fonts, line widths and markers with the export canvas.
pixelsPerInch = get(groot,'ScreenPixelsPerInch');
scale = paperSize(1)*pixelsPerInch/canvas(1);
scaleMechanismGraphics(exportFigure,scale)
set(exportFigure,'Units','inches', ...
    'Position',[0.5 0.5 paperSize], ...
    'PaperUnits','inches','PaperPositionMode','manual', ...
    'PaperPosition',[0 0 paperSize],'PaperSize',paperSize, ...
    'Color','w','InvertHardcopy','off');
drawnow

deleteIfFileExists(pngPath)
print(exportFigure,pngPath,'-dpng',sprintf('-r%d',dpi),'-painters');
[valid,detail,width,height,blackFraction,whiteCorners] = ...
    validatePngAcrossCanvas(pngPath,maximumPixels);
expectedPixels = round(dpi*paperSize);
if any(abs([width height]-expectedPixels)>8)
    valid = false;
    detail = detail+"; exported dimensions differ from the canvas";
end
if ~valid
    error('Saved PNG validation failed for %s: %s.',stem,char(detail))
end
fprintf('Validated PNG: %s (%s).\n',pngPath,char(detail))

deleteIfFileExists(pdfPath)
try
    print(exportFigure,pdfPath,'-dpdf','-painters');
catch ME
    warning('PDF export failed for %s; .fig and validated PNG remain: %s', ...
        stem,ME.message)
end

Figure = string(stem);
PngPath = string(pngPath);
WidthPixels = width;
HeightPixels = height;
RequestedDPI = dpi;
MaximumDimensionPixels = maximumPixels;
BlackPixelPercent = 100*blackFraction;
WhiteCorners = whiteCorners;
Passed = valid;
Detail = detail;
Audit = table(Figure,PngPath,WidthPixels,HeightPixels,RequestedDPI, ...
    MaximumDimensionPixels,BlackPixelPercent,WhiteCorners,Passed,Detail);
end

function scaleMechanismGraphics(fig,scale)
% Keep text, lines and node sizes proportional when the canvas changes size.
fontObjects = findall(fig,'-property','FontSize');
fontSizes = zeros(numel(fontObjects),1);
for q = 1:numel(fontObjects)
    if isprop(fontObjects(q),'FontUnits')
        fontObjects(q).FontUnits = 'points';
    end
    fontSizes(q) = fontObjects(q).FontSize;
end
% Capture all sizes first: changing axes fonts can also change their labels.
for q = 1:numel(fontObjects)
    fontObjects(q).FontSize = scale*fontSizes(q);
end

properties = {'LineWidth','MarkerSize','Margin'};
for p = 1:numel(properties)
    objects = findall(fig,'-property',properties{p});
    values = cell(numel(objects),1);
    for q = 1:numel(objects)
        values{q} = get(objects(q),properties{p});
    end
    for q = 1:numel(objects)
        set(objects(q),properties{p},scale*values{q});
    end
end
nodes = findall(fig,'Type','scatter');
for q = 1:numel(nodes)
    nodes(q).SizeData = scale^2*nodes(q).SizeData;
end
end

function drawMechanismNetwork(ax,fontName)
% Node positions in diagram coordinates.
MMC = [0.75 3.45];
M1 = [5.00 6.10];
M2 = [3.20 4.65];
M3 = [3.55 3.35];
M4 = [6.45 3.35];
M5 = [5.00 1.05];
day360 = [9.55 3.40];
primaryColor = [0.43 0.67 0.88];
modulatoryColor = [0.84 0.12 0.10];

hold(ax,'on')
set(ax,'XLim',[0 10.4],'YLim',[0 7.15], ...
    'DataAspectRatio',[1 1 1],'Visible','off','SortMethod','childorder');

% Each control point sets the curvature between the two named nodes.
drawNetworkArrow(ax,MMC,M1,[2.133 5.965],primaryColor,2.5,'-');
drawNetworkArrow(ax,M1,day360,[7.896 5.7965],primaryColor,2.5,'-');
drawNetworkArrow(ax,MMC,M3,[2.150 3.400],primaryColor,2.5,'-');
drawNetworkArrow(ax,M3,M4,[5.000 3.350],primaryColor,2.5,'-');
drawNetworkArrow(ax,M4,day360,[8.000 3.375],primaryColor,2.5,'-');
drawNetworkArrow(ax,MMC,M5,[2.155 0.975],primaryColor,2.5,'-');
drawNetworkArrow(ax,M5,M4,[6.231 1.881],modulatoryColor,2.2,'--');
drawNetworkArrow(ax,M2,M1,[3.984 5.519],modulatoryColor,2.1,'--');
drawNetworkArrow(ax,M2,M3,[3.310 3.9825],modulatoryColor,2.1,'--');
drawNetworkArrow(ax,M2,M4,[4.565 3.350],modulatoryColor,2.1,'--');

drawNetworkNode(ax,MMC,[0.10 0.45 0.85]);
drawNetworkNode(ax,M1,[0.10 0.45 0.85]);
drawNetworkNode(ax,M2,[0.78 0.25 0.20]);
drawNetworkNode(ax,M3,[0.20 0.65 0.35]);
drawNetworkNode(ax,M4,[0.20 0.65 0.35]);
drawNetworkNode(ax,M5,[0.60 0.25 0.75]);
drawNetworkNode(ax,day360,[0.98 0.72 0.08]);

drawNetworkLabel(ax,[0.43 3.45],'MMC',fontName,12,'right','middle');
drawNetworkLabel(ax,[5.00 7.08],'M1',fontName,12,'center','top');
drawNetworkLabel(ax,[3.52 4.65],'M2',fontName,12,'left','middle');
drawNetworkLabel(ax,[3.55 2.96],'M3',fontName,12,'center','top');
drawNetworkLabel(ax,[6.45 4.44],'M4',fontName,12,'center','top');
drawNetworkLabel(ax,[5.00 0.72],'M5',fontName,12,'center','top');
drawNetworkLabel(ax,[9.98 3.40],'Day-360',fontName,12,'left','middle');
drawNetworkLabel(ax,[5.00 6.70],'Direct MMC effects', ...
    fontName,11.5,'center','top');
drawNetworkLabel(ax,[3.98 4.65],'Tumor heterogeneity', ...
    fontName,11.5,'left','middle');
drawNetworkLabel(ax,[3.55 2.55],'ICD/DC maturation', ...
    fontName,11.5,'center','top');
drawNetworkLabel(ax,[6.45 4.06],'Antitumor immune response', ...
    fontName,11.5,'center','top');
drawNetworkLabel(ax,[5.00 0.33],'Treg-mediated suppression', ...
    fontName,11.5,'center','top');

plot(ax,[6.65 7.25],[6.95 6.95],'-', ...
    'Color',primaryColor,'LineWidth',2.5);
plot(ax,[6.65 7.25],[6.65 6.65],'--', ...
    'Color',modulatoryColor,'LineWidth',2.2);
drawNetworkLabel(ax,[7.37 6.95],'primary model effect', ...
    fontName,9.5,'left','middle');
drawNetworkLabel(ax,[7.37 6.65],'modulatory or suppressive effect', ...
    fontName,9.5,'left','middle');
end

function drawNetworkArrow(ax,startNode,endNode,controlPoint,color,width,style)
% Quadratic Bezier curve with clearance around each node.
nodeClearance = 0.27;
startDirection = controlPoint-startNode;
endDirection = endNode-controlPoint;
startPoint = startNode+nodeClearance*startDirection/norm(startDirection);
endPoint = endNode-nodeClearance*endDirection/norm(endDirection);
t = linspace(0,1,220).';
curve = (1-t).^2*startPoint+2*(1-t).*t*controlPoint+t.^2*endPoint;
plot(ax,curve(:,1),curve(:,2),'Color',color, ...
    'LineWidth',width,'LineStyle',style);

direction = curve(end,:)-curve(end-6,:);
direction = direction/norm(direction);
normal = [-direction(2) direction(1)];
headLength = 0.22;
headHalfWidth = 0.09;
base = endPoint-headLength*direction;
vertices = [endPoint; ...
    base+headHalfWidth*normal;base-headHalfWidth*normal];
patch(ax,vertices(:,1),vertices(:,2),color,'EdgeColor',color);
end

function drawNetworkNode(ax,position,color)
scatter(ax,position(1),position(2),360,color,'filled', ...
    'MarkerEdgeColor',[0.08 0.08 0.08],'LineWidth',1.1);
end

function drawNetworkLabel(ax,position,label,fontName,fontSize,horizontal,vertical)
text(ax,position(1),position(2),label,'FontName',fontName, ...
    'FontSize',fontSize,'FontWeight','bold','Color',[0 0 0], ...
    'HorizontalAlignment',horizontal,'VerticalAlignment',vertical, ...
    'Interpreter','none','Clipping','off');
end

function [upperLimit,tickValues] = niceCountAxis(maximumCount)
if maximumCount<=1
    upperLimit = 1.30;
    tickValues = [0 0.5 1];
    return
end
rawUpper = 1.24*maximumCount;
magnitude = 10^floor(log10(rawUpper));
normalized = rawUpper/magnitude;
if normalized<=1
    niceNormalized = 1;
elseif normalized<=2
    niceNormalized = 2;
elseif normalized<=5
    niceNormalized = 5;
else
    niceNormalized = 10;
end
upperLimit = niceNormalized*magnitude;
tickStep = upperLimit/5;
tickValues = 0:tickStep:upperLimit;
end

function values = collectFigureStrings(fig)
objects = findall(fig,'-property','String');
values = strings(0,1);
for q = 1:numel(objects)
    try
        current = string(objects(q).String);
        values = [values;current(:)]; %#ok<AGROW>
    catch
    end
end
values = values(strlength(values)>0);
end

% ========================================================================
% SELECTED PATIENT SCREEN AND CONTINUOUS PATH
% ========================================================================

function T = buildPatientMechanismScreen( ...
        patientIndex,simulation,mechanisms,settings)
G = numel(mechanisms);
Sdraw = settings.nPerturbations;
MechanismNumber = (1:G).';
MechanismID = string({mechanisms.id}).';
Mechanism = string({mechanisms.title}).';
FailureDraws = zeros(G,1);
ValidDraws = zeros(G,1);
FailurePercent = nan(G,1);
for g = 1:G
    ok = reshape(simulation.PerturbedOK(patientIndex,g,:),Sdraw,1);
    tumor = reshape(simulation.PerturbedFinalTumor( ...
        patientIndex,g,:),Sdraw,1);
    valid = simulation.ReferenceControlled(patientIndex) & ok & ...
        isfinite(tumor);
    failure = valid & tumor>=settings.detectionLimit;
    ValidDraws(g) = sum(valid);
    FailureDraws(g) = sum(failure);
    if ValidDraws(g)>0
        FailurePercent(g) = 100*FailureDraws(g)/ValidDraws(g);
    end
end
T = table(MechanismNumber,MechanismID,Mechanism,FailureDraws, ...
    ValidDraws,FailurePercent);
if ~isequal(T.FailureDraws,[0;4;4;17;0]) || any(T.ValidDraws~=20)
    error('VP 127 mechanism screen does not reconcile with saved results.')
end
end

function [Surface,PathTable,EndpointTable,AuditTable] = ...
        computeCaseSurface(patient,plan,targetSample,spec,fixed, ...
        settings,odeOpt,expectedReferenceTumor,expectedFailureTumor)
pathPosition = linspace(0,1,settings.surfaceLevels);
timeGrid = linspace(0,settings.tEnd,settings.surfaceTimePoints);
tumor = nan(numel(pathPosition),numel(timeGrid));
finalTumor = nan(numel(pathPosition),1);
TsAbsorbed = false(numel(pathPosition),1);
TotalTumorAbsorbed = false(numel(pathPosition),1);
TsAbsorptionTime = nan(numel(pathPosition),1);
TotalTumorAbsorptionTime = nan(numel(pathPosition),1);
Completed = false(numel(pathPosition),1);

for q = 1:numel(pathPosition)
    patientUse = interpolatePatientToDraw( ...
        patient,plan,targetSample,pathPosition(q),spec);
    y0 = [fixed.M0;patientUse.Ts0;patientUse.Tu0;patientUse.Di0; ...
        fixed.Dm0;fixed.E0;fixed.R0];
    result = simulateScheduleWithAbsorption( ...
        y0,patientUse,settings,odeOpt,true);
    Completed(q) = result.ok;
    if ~result.ok
        continue
    end
    [uniqueTime,uniqueIndex] = unique(result.time,'stable');
    uniqueState = result.state(uniqueIndex,:);
    stateOnGrid = interp1(uniqueTime,uniqueState,timeGrid,'linear');
    if result.TsAbsorbed
        stateOnGrid(timeGrid>=result.TsAbsorptionTime,2) = 0;
    end
    if result.totalAbsorbed
        stateOnGrid(timeGrid>=result.totalAbsorptionTime,2:3) = 0;
    end
    tumor(q,:) = stateOnGrid(:,2).'+stateOnGrid(:,3).';
    finalTumor(q) = result.finalState(2)+result.finalState(3);
    TsAbsorbed(q) = result.TsAbsorbed;
    TotalTumorAbsorbed(q) = result.totalAbsorbed;
    TsAbsorptionTime(q) = result.TsAbsorptionTime;
    TotalTumorAbsorptionTime(q) = result.totalAbsorptionTime;
end

if ~all(Completed)
    error('Selected-case path failed at indices %s.', ...
        mat2str(find(~Completed).'))
end
if any(~isfinite(tumor(:))) || any(tumor(:)<0)
    error('Selected-case surface contains invalid tumor values.')
end
surfaceEndpointDifference = max(abs(tumor(:,end)-finalTumor));
surfaceTolerance = 1e-7*max(1,max(abs(finalTumor)));
referenceDifference = abs(finalTumor(1)-expectedReferenceTumor);
failureDifference = abs(finalTumor(end)-expectedFailureTumor);
referenceTolerance = 1e-7*max(1,abs(expectedReferenceTumor));
failureTolerance = 1e-7*max(1,abs(expectedFailureTumor));
if surfaceEndpointDifference>surfaceTolerance || ...
        referenceDifference>referenceTolerance || ...
        failureDifference>failureTolerance
    error(['Selected-case surface does not reproduce the saved paper ' ...
        'endpoints.'])
end
if finalTumor(1)>=settings.detectionLimit || ...
        finalTumor(end)<settings.detectionLimit
    error('Selected path is not a control-to-failure path.')
end

Surface = struct('PathPosition',pathPosition,'TimeGrid',timeGrid, ...
    'Tumor',tumor,'FinalTumor',finalTumor,'TsAbsorbed',TsAbsorbed, ...
    'TotalTumorAbsorbed',TotalTumorAbsorbed, ...
    'TsAbsorptionTime',TsAbsorptionTime, ...
    'TotalTumorAbsorptionTime',TotalTumorAbsorptionTime, ...
    'TargetSample',targetSample,'InputNames',{plan.inputNames});

Parameter = string(plan.inputNames(:));
ReferenceValue = nan(numel(plan.inputNames),1);
FailureValue = nan(numel(plan.inputNames),1);
Direction = strings(numel(plan.inputNames),1);
Interpolation = strings(numel(plan.inputNames),1);
for p = 1:numel(plan.inputNames)
    name = plan.inputNames{p};
    ReferenceValue(p) = patient.(name);
    FailureValue(p) = plan.values(targetSample,p);
    if FailureValue(p)>ReferenceValue(p)
        Direction(p) = "increase";
    elseif FailureValue(p)<ReferenceValue(p)
        Direction(p) = "decrease";
    else
        Direction(p) = "unchanged";
    end
    if strcmp(getSpec(spec,name).measure,'log')
        Interpolation(p) = "geometric/log-linear";
    else
        Interpolation(p) = "linear";
    end
end
EndpointTable = table(Parameter,ReferenceValue,FailureValue,Direction, ...
    Interpolation);

nRows = numel(pathPosition)*numel(plan.inputNames);
PathPosition = nan(nRows,1);
PathParameter = strings(nRows,1);
PathValue = nan(nRows,1);
PathFinalTumorCells = nan(nRows,1);
row = 0;
for q = 1:numel(pathPosition)
    patientUse = interpolatePatientToDraw( ...
        patient,plan,targetSample,pathPosition(q),spec);
    for p = 1:numel(plan.inputNames)
        row = row+1;
        PathPosition(row) = pathPosition(q);
        PathParameter(row) = string(plan.inputNames{p});
        PathValue(row) = patientUse.(plan.inputNames{p});
        PathFinalTumorCells(row) = finalTumor(q);
    end
end
PathTable = table(PathPosition,PathParameter,PathValue, ...
    PathFinalTumorCells);

Check = ["All path simulations completed"; ...
    "Surface finite and nonnegative"; ...
    "All surface day-360 values reconcile"; ...
    "Reference endpoint reconciles"; ...
    "Failure endpoint reconciles"; ...
    "Reference endpoint is controlled"; ...
    "Selected endpoint is detectable failure"];
Passed = [all(Completed);all(isfinite(tumor(:)) & tumor(:)>=0); ...
    surfaceEndpointDifference<=surfaceTolerance; ...
    referenceDifference<=referenceTolerance; ...
    failureDifference<=failureTolerance; ...
    finalTumor(1)<settings.detectionLimit; ...
    finalTumor(end)>=settings.detectionLimit];
Detail = [string(sum(Completed))+"/"+string(numel(Completed)); ...
    string(min(tumor(:)))+" to "+string(max(tumor(:)))+" cells"; ...
    string(surfaceEndpointDifference)+" cells maximum difference"; ...
    string(referenceDifference)+" cells difference"; ...
    string(failureDifference)+" cells difference"; ...
    string(finalTumor(1))+" cells";string(finalTumor(end))+" cells"];
AuditTable = table(Check,Passed,Detail);
end

function patientOut = interpolatePatientToDraw( ...
        patient,plan,targetSample,lambda,spec)
patientOut = patient;
if lambda==0
    return
end
for q = 1:numel(plan.inputNames)
    name = plan.inputNames{q};
    referenceValue = patient.(name);
    failureValue = plan.values(targetSample,q);
    if lambda==1
        value = failureValue;
    elseif strcmp(getSpec(spec,name).measure,'log') && ...
            referenceValue>0 && failureValue>0
        value = exp((1-lambda)*log(referenceValue) + ...
            lambda*log(failureValue));
    else
        value = (1-lambda)*referenceValue+lambda*failureValue;
    end
    patientOut.(name) = value;
end
end

function item = getSpec(spec,name)
index = find(strcmp({spec.name},name),1,'first');
if isempty(index)
    error('Unknown parameter %s.',name)
end
item = spec(index);
end

% =======================================================================
% FIGURE 8-SELECTED PATIENT-CASE PRESENTATION
% ========================================================================

function [fig,Audit] = makePatientCaseFigure( ...
        selectedMechanism,Tpatient,Surface,patient,viewAngles,settings)
pathPosition = double(Surface.PathPosition(:).');
timeGrid = double(Surface.TimeGrid(:));
Z = log10(double(Surface.Tumor).'+1);
[X,Y] = meshgrid(pathPosition,timeGrid);

zCapacity = log10(double(patient.k)+1);
zDetection = log10(settings.detectionLimit+1);
zOneCell = log10(settings.oneCellLimit+1);

fig = figure('Color','w','Units','pixels','Position',[40 40 1700 1100], ...
    'Renderer','opengl','Visible','off');
annotation(fig,'textbox',[0.080 0.920 0.03 0.05], ...
    'String','A','EdgeColor','none','FontName',settings.figureFont, ...
    'FontSize',20,'FontWeight','bold');
annotation(fig,'textbox',[0.080 0.625 0.03 0.05], ...
    'String','B','EdgeColor','none','FontName',settings.figureFont, ...
    'FontSize',20,'FontWeight','bold');

panelAPosition = [0.310 0.770 0.380 0.160];
axA = axes(fig,'Units','normalized','Position',panelAPosition);
hold(axA,'on'); box(axA,'on'); grid(axA,'on')
set(axA,'FontName',settings.figureFont,'FontSize',10, ...
    'FontWeight','bold','LineWidth',1.2,'YDir','reverse');

barBlue = [0.36 0.62 0.84];
barEdge = [0.20 0.35 0.50];
red = [0.85 0.05 0.05];
hPatientBars = barh(axA,Tpatient.MechanismNumber, ...
    Tpatient.FailurePercent,0.62, ...
    'FaceColor',barBlue,'EdgeColor',barEdge,'LineWidth',1.2);
selectedRow = Tpatient.MechanismNumber==selectedMechanism;
hSelectedOutline = barh(axA,Tpatient.MechanismNumber(selectedRow), ...
    Tpatient.FailurePercent(selectedRow),0.62,'FaceColor','none', ...
    'EdgeColor',red,'LineWidth',2.6);
for q = 1:height(Tpatient)
    if Tpatient.FailurePercent(q)<=85
        xText = Tpatient.FailurePercent(q)+1.8;
        alignment = 'left';
    else
        xText = Tpatient.FailurePercent(q)-1.8;
        alignment = 'right';
    end
    text(axA,xText,Tpatient.MechanismNumber(q),sprintf('%d/%d', ...
        Tpatient.FailureDraws(q),Tpatient.ValidDraws(q)), ...
        'FontName',settings.figureFont,'FontSize',11, ...
        'FontWeight','bold','HorizontalAlignment',alignment, ...
        'VerticalAlignment','middle','Clipping','on');
end
set(axA,'YTick',1:5,'YTickLabel',{'M1','M2','M3','M4','M5'});
xlabel(axA,'Failure-inducing perturbations (%)', ...
    'FontSize',12,'FontWeight','bold');
xlim(axA,[0 100]);

panelBPosition = [0.125 0.100 0.700 0.525];
axB = axes(fig,'Units','normalized','Position',panelBPosition);
hold(axB,'on'); box(axB,'on'); grid(axB,'on')
set(axB,'FontName',settings.figureFont,'FontSize',13, ...
    'FontWeight','bold','LineWidth',1.2);
colormap(axB,turbo(256));
surf(axB,X,Y,Z,'EdgeColor','none','FaceAlpha',0.97);
view(axB,viewAngles(1),viewAngles(2));
xlabel(axB,'Path position, s','FontSize',15,'FontWeight','bold');
ylabel(axB,'Time (days)','FontSize',15,'FontWeight','bold');
zlabel(axB,'log_{10}(T_s + T_u + 1)', ...
    'FontSize',15,'FontWeight','bold');
set(axB,'XTick',[0 0.25 0.50 0.75 1], ...
    'XTickLabel',{'0','0.25','0.50','0.75','1'}, ...
    'YTick',[0 120 240 360]);
xlim(axB,[0 1]); ylim(axB,[0 settings.tEnd]);

colorBarPosition = [0.855 0.135 0.017 0.405];
colorBar = colorbar(axB,'Position',colorBarPosition);
colorBar.Units = 'normalized';
colorBar.Label.String = '';
colorBar.FontName = settings.figureFont;
colorBar.FontSize = 12;

planeX = [0 1;0 1];
planeY = [0 0;settings.tEnd settings.tEnd];
surf(axB,planeX,planeY,zOneCell*ones(2), ...
    'FaceColor',[0.45 0.45 0.45],'FaceAlpha',0.16,'EdgeColor','none');
surf(axB,planeX,planeY,zDetection*ones(2), ...
    'FaceColor',[0.96 0.66 0.18],'FaceAlpha',0.22,'EdgeColor','none');
surf(axB,planeX,planeY,zCapacity*ones(2), ...
    'FaceColor',[0.78 0.22 0.22],'FaceAlpha',0.16,'EdgeColor','none');

orange = [0.82 0.44 0.00];
text(axB,0.90,180,zDetection+0.02,'Detection threshold', ...
    'Color',orange,'BackgroundColor','w','EdgeColor',orange, ...
    'Margin',3,'FontSize',10,'FontWeight','bold', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Clipping','off');
text(axB,0.84,300,zOneCell+0.03,'1-cell', ...
    'Color',[0.15 0.15 0.15],'BackgroundColor','w', ...
    'EdgeColor',[0.30 0.30 0.30],'Margin',3,'FontSize',10, ...
    'FontWeight','bold','HorizontalAlignment','center', ...
    'VerticalAlignment','middle');
text(axB,0.60,255,zCapacity+0.03,'Carrying capacity', ...
    'Color',[0.70 0.12 0.12],'BackgroundColor','w', ...
    'EdgeColor',[0.70 0.12 0.12],'Margin',3,'FontSize',10, ...
    'FontWeight','bold','HorizontalAlignment','center', ...
    'VerticalAlignment','middle');

zReference = Z(:,1);
zFailure = Z(:,end);
hReference = plot3(axB,zeros(size(timeGrid)),timeGrid,zReference, ...
    'k-','LineWidth',3.6);
hFailureTrajectory = customShortDashed3( ...
    axB,ones(size(timeGrid)),timeGrid,zFailure,red,3.0);
hFailureLegend = plot3(axB,nan,nan,nan,'--','Color',red, ...
    'LineWidth',3.2);
legendPosition = [0.625 0.630 0.330 0.060];
legendHandle = legend(axB,[hReference hFailureLegend], ...
    {'Reference (s=0)','Failure-inducing perturbation (s=1)'}, ...
    'Box','off','FontName',settings.figureFont,'FontSize',11, ...
    'Interpreter','none','AutoUpdate','off');
legendHandle.Units = 'normalized';
legendHandle.Position = legendPosition;

finiteZ = Z(isfinite(Z));
zUpper = max([max(finiteZ(:))+0.20,zCapacity+0.18,zDetection+0.18]);
zlim(axB,[0 zUpper]);
axA.Position = panelAPosition;
axB.Position = panelBPosition;
colorBar.Position = colorBarPosition;
legendHandle.Position = legendPosition;

zTicks = axB.ZTick;
zTickLabels = string(axB.ZTickLabel);
zeroZTick = find(abs(zTicks)<1e-12,1,'first');
if ~isempty(zeroZTick) && numel(zTickLabels)==numel(zTicks)
    zTickLabels(zeroZTick) = "";
    axB.ZTickLabel = cellstr(zTickLabels);
end

enableEditableFigureControls(fig,[axA;axB]);
set(fig,'Visible','on');
drawnow
axA.Position = panelAPosition;
axB.Position = panelBPosition;
colorBar.Position = colorBarPosition;
legendHandle.Position = legendPosition;

legendStrings = string(legendHandle.String(:));
legendInside = legendHandle.Position(1)>=0 && ...
    legendHandle.Position(2)>=0 && ...
    sum(legendHandle.Position([1 3]))<=1 && ...
    sum(legendHandle.Position([2 4]))<=1;
Check = ["Five mechanism bars are present"; ...
    "Selected mechanism is outlined"; ...
    "Reference and failure trajectories are present"; ...
    "Complete two-entry legend is retained inside canvas"];
Passed = [height(Tpatient)==5 && isgraphics(hPatientBars); ...
    sum(selectedRow)==1 && isgraphics(hSelectedOutline); ...
    all(isgraphics([hReference hFailureTrajectory hFailureLegend])); ...
    numel(legendStrings)==2 && ...
    legendStrings(1)=="Reference (s=0)" && ...
    legendStrings(2)=="Failure-inducing perturbation (s=1)" && ...
    legendInside];
Detail = ["5 mechanisms";"M"+string(selectedMechanism); ...
    "solid black reference; short-dashed red failure"; ...
    strjoin(legendStrings," | ")+"; inside canvas="+string(legendInside)];
Audit = table(Check,Passed,Detail);
if ~all(Passed)
    error('Figure 10 content audit failed: %s.', ...
        char(strjoin(Check(~Passed),'; ')))
end
end

function lineHandle = customShortDashed3(ax,x,y,z,lineColor,lineWidth)
x = x(:); y = y(:); z = z(:);
points = [x y z];
arc = [0;cumsum(sqrt(sum(diff(points,1,1).^2,2)))];
if arc(end)<=eps
    lineHandle = plot3(ax,x,y,z,'-','Color',lineColor, ...
        'LineWidth',lineWidth);
    return
end
arc = arc/arc(end);
arcDense = linspace(0,1,1000).';
xDense = interp1(arc,x,arcDense,'linear');
yDense = interp1(arc,y,arcDense,'linear');
zDense = interp1(arc,z,arcDense,'linear');
onMask = mod(arcDense,0.110)<=0.070;
xDense(~onMask) = nan;
yDense(~onMask) = nan;
zDense(~onMask) = nan;
lineHandle = plot3(ax,xDense,yDense,zDense,'-','Color',lineColor, ...
    'LineWidth',lineWidth);
end

function enableEditableFigureControls(fig,axesHandles)
set(fig,'WindowStyle','normal','MenuBar','figure','ToolBar','figure', ...
    'DockControls','on','HandleVisibility','on','NumberTitle','off');
for q = 1:numel(axesHandles)
    if ~isgraphics(axesHandles(q),'axes')
        continue
    end
    try
        enableDefaultInteractivity(axesHandles(q));
    catch
    end
    try
        axesHandles(q).Toolbar.Visible = 'on';
    catch
    end
end
end

function validateEditableFigureFile(figurePath,minimumAxes)
try
    probe = openfig(figurePath,'new','invisible');
catch ME
    error('Saved MATLAB figure cannot be reopened: %s',ME.message)
end
cleanup = onCleanup(@()close(probe)); %#ok<NASGU>
axesObjects = findall(probe,'Type','axes');
if numel(axesObjects)<minimumAxes
    error('Saved figure has %d editable axes; expected at least %d.', ...
        numel(axesObjects),minimumAxes)
end
if strcmpi(probe.ToolBar,'none') || strcmpi(probe.MenuBar,'none')
    error('Saved figure did not retain standard editing controls.')
end
fprintf('Validated editable MATLAB figure: %s (%d axes).\n', ...
    figurePath,numel(axesObjects))
end

% ========================================================================
% SELECTED-PATH ODE MODEL
% ========================================================================

function result = simulateScheduleWithAbsorption( ...
        y0,pars,settings,odeOpt,collectTrajectory)
starts = settings.treatmentTimes(settings.treatmentTimes<settings.tEnd);
ends = min(starts+settings.doseDuration,settings.tEnd);
breakPoints = sort(unique([0 settings.tEnd starts ends]));

result = struct('ok',false,'finalState',nan(7,1), ...
    'TsAbsorbed',false,'totalAbsorbed',false, ...
    'TsAbsorptionTime',nan,'totalAbsorptionTime',nan, ...
    'message',"not completed",'time',[],'state',[]);
yCurrent = y0(:);
TsIsAbsorbed = false;
totalIsAbsorbed = false;
TsAbsorptionTime = nan;
totalAbsorptionTime = nan;
trajectoryTime = [];
trajectoryState = [];

if ~hasAdmissibleInitialTumorState(yCurrent,settings.oneCellLimit)
    result.message = "invalid sub-one initial tumor state";
    return
end

for b = 1:numel(breakPoints)-1
    t0 = breakPoints(b);
    t1 = breakPoints(b+1);
    if t1<=t0
        continue
    end
    midpoint = 0.5*(t0+t1);
    mCurrent = activeInstillationRate( ...
        midpoint,starts,settings.doseDuration,pars.m);
    piece = integratePieceWithBoundaries( ...
        t0,t1,yCurrent,pars,mCurrent,settings.oneCellLimit,odeOpt, ...
        TsIsAbsorbed,totalIsAbsorbed,TsAbsorptionTime, ...
        totalAbsorptionTime,collectTrajectory);
    if ~piece.ok
        result.message = piece.message;
        return
    end
    yCurrent = piece.finalState;
    TsIsAbsorbed = piece.TsAbsorbed;
    totalIsAbsorbed = piece.totalAbsorbed;
    TsAbsorptionTime = piece.TsAbsorptionTime;
    totalAbsorptionTime = piece.totalAbsorptionTime;
    if collectTrajectory
        [trajectoryTime,trajectoryState] = appendTrajectorySegment( ...
            trajectoryTime,trajectoryState,piece.time,piece.state);
    end
end

if any(~isfinite(yCurrent)) || any(yCurrent<0)
    result.message = "nonfinite or negative final state";
    return
end
yCurrent = reconstructAbsorbingState( ...
    yCurrent,TsIsAbsorbed,totalIsAbsorbed);
result.ok = true;
result.finalState = yCurrent;
result.TsAbsorbed = TsIsAbsorbed;
result.totalAbsorbed = totalIsAbsorbed;
result.TsAbsorptionTime = TsAbsorptionTime;
result.totalAbsorptionTime = totalAbsorptionTime;
result.message = "ode15s completed day 360";
if collectTrajectory
    result.time = trajectoryTime;
    result.state = trajectoryState;
end
end

function tf = hasAdmissibleInitialTumorState(y,oneCellLimit)
y = y(:);
tf = numel(y)>=3 && all(isfinite(y(2:3))) && ...
    y(2)>=oneCellLimit && y(3)>=oneCellLimit && ...
    y(2)+y(3)>=oneCellLimit;
end

function piece = integratePieceWithBoundaries( ...
        t0,t1,y0,pars,mCurrent,oneCellLimit,odeOpt, ...
        TsIsAbsorbed,totalIsAbsorbed,TsAbsorptionTime, ...
        totalAbsorptionTime,collectTrajectory)
piece = struct('ok',false,'finalState',nan(7,1), ...
    'TsAbsorbed',TsIsAbsorbed,'totalAbsorbed',totalIsAbsorbed, ...
    'TsAbsorptionTime',TsAbsorptionTime, ...
    'totalAbsorptionTime',totalAbsorptionTime, ...
    'message',"integration failed",'time',[],'state',[]);
tCurrent = t0;
yCurrent = reconstructAbsorbingState( ...
    y0,TsIsAbsorbed,totalIsAbsorbed);
pieceTime = [];
pieceState = [];

for solverCall = 1:3
    if tCurrent>=t1-1e-12*max(1,abs(t1))
        piece.ok = true;
        piece.finalState = yCurrent;
        piece.message = "ode15s";
        piece.time = pieceTime;
        piece.state = pieceState;
        return
    end
    monitorTs = ~TsIsAbsorbed && ~totalIsAbsorbed;
    monitorTotal = ~totalIsAbsorbed;
    if monitorTs || monitorTotal
        eventOpt = odeset(odeOpt,'Events',@(t,y)oneCellBoundaryEvents( ...
            t,y,monitorTs,monitorTotal,oneCellLimit));
        [ok,t,Y,te,ye,ie,message] = solveSegment( ...
            tCurrent,t1,yCurrent,pars,mCurrent,eventOpt, ...
            TsIsAbsorbed,totalIsAbsorbed,true);
    else
        [ok,t,Y,te,ye,ie,message] = solveSegment( ...
            tCurrent,t1,yCurrent,pars,mCurrent,odeOpt, ...
            TsIsAbsorbed,totalIsAbsorbed,false);
    end
    if ~ok
        piece.message = message;
        return
    end
    if isempty(ie)
        if collectTrajectory
            [pieceTime,pieceState] = appendTrajectorySegment( ...
                pieceTime,pieceState,t,Y);
        end
        piece.ok = true;
        piece.finalState = reconstructAbsorbingState( ...
            Y(end,:).',TsIsAbsorbed,totalIsAbsorbed);
        piece.TsAbsorbed = TsIsAbsorbed;
        piece.totalAbsorbed = totalIsAbsorbed;
        piece.TsAbsorptionTime = TsAbsorptionTime;
        piece.totalAbsorptionTime = totalAbsorptionTime;
        piece.message = "ode15s";
        piece.time = pieceTime;
        piece.state = pieceState;
        return
    end

    eventTime = te(1);
    eventState = ye(1,:).';
    tolerance = 1e-10*max(1,abs(eventTime));
    eventsAtStop = ie(abs(te-eventTime)<=tolerance);
    if any(eventsAtStop==2)
        if ~TsIsAbsorbed
            TsAbsorptionTime = eventTime;
        end
        TsIsAbsorbed = true;
        totalIsAbsorbed = true;
        totalAbsorptionTime = eventTime;
    elseif any(eventsAtStop==1)
        TsIsAbsorbed = true;
        TsAbsorptionTime = eventTime;
        eventState = reconstructAbsorbingState( ...
            eventState,TsIsAbsorbed,totalIsAbsorbed);
        if eventState(2)+eventState(3)<oneCellLimit
            totalIsAbsorbed = true;
            totalAbsorptionTime = eventTime;
        end
    else
        piece.message = "unrecognized one-cell boundary event";
        return
    end
    yCurrent = reconstructAbsorbingState( ...
        eventState,TsIsAbsorbed,totalIsAbsorbed);
    if collectTrajectory
        [pieceTime,pieceState] = appendTrajectorySegment( ...
            pieceTime,pieceState,t,Y);
    end
    tCurrent = eventTime;
    piece.TsAbsorbed = TsIsAbsorbed;
    piece.totalAbsorbed = totalIsAbsorbed;
    piece.TsAbsorptionTime = TsAbsorptionTime;
    piece.totalAbsorptionTime = totalAbsorptionTime;
end
piece.message = "more than two boundary changes in one treatment piece";
end

function [timeOut,stateOut] = appendTrajectorySegment( ...
        timeIn,stateIn,timeNew,stateNew)
timeNew = timeNew(:);
if isempty(timeNew)
    timeOut = timeIn;
    stateOut = stateIn;
    return
end
if isempty(timeIn)
    timeOut = timeNew;
    stateOut = stateNew;
    return
end
sameBoundary = abs(timeNew(1)-timeIn(end)) <= ...
    1e-12*max(1,abs(timeIn(end)));
if sameBoundary
    timeNew(1) = [];
    stateNew(1,:) = [];
end
timeOut = [timeIn;timeNew]; %#ok<AGROW>
stateOut = [stateIn;stateNew]; %#ok<AGROW>
end

function [ok,t,Y,te,ye,ie,failureMessage] = solveSegment( ...
        t0,t1,y0,pars,mCurrent,options, ...
        TsIsAbsorbed,totalIsAbsorbed,withEvents)
ok = false;
t = []; Y = []; te = []; ye = []; ie = [];
failureMessage = "ode15s returned an incomplete or nonfinite segment";
try
    if withEvents
        [t,Y,te,ye,ie] = ode15s( ...
            @(t,y)mechanismODE(t,y,pars,mCurrent, ...
            TsIsAbsorbed,totalIsAbsorbed),[t0 t1],y0,options);
    else
        [t,Y] = ode15s( ...
            @(t,y)mechanismODE(t,y,pars,mCurrent, ...
            TsIsAbsorbed,totalIsAbsorbed),[t0 t1],y0,options);
    end
    reachedBoundaryOrEnd = ~isempty(ie) || ...
        (~isempty(t) && t(end)>=t1-1e-10*max(1,abs(t1)));
    if ~isempty(t) && ~isempty(Y) && reachedBoundaryOrEnd && ...
            all(isfinite(Y(:))) && all(Y(:)>=0)
        ok = true;
    end
catch ME
    failureMessage = "ode15s error: "+string(ME.message);
end
end

function [value,isterminal,direction] = oneCellBoundaryEvents( ...
        ~,y,monitorTs,monitorTotal,oneCellLimit)
value = ones(2,1);
if monitorTs
    value(1) = y(2)-oneCellLimit;
end
if monitorTotal
    value(2) = y(2)+y(3)-oneCellLimit;
end
isterminal = [1;1];
direction = [-1;-1];
end

function y = reconstructAbsorbingState(y,TsIsAbsorbed,totalIsAbsorbed)
y = y(:);
if totalIsAbsorbed
    y(2:3) = 0;
elseif TsIsAbsorbed
    y(2) = 0;
end
end

function mCurrent = activeInstillationRate(t,starts,tau,mDose)
if any(t>=starts & t<starts+tau)
    mCurrent = mDose;
else
    mCurrent = 0;
end
end

function dydt = mechanismODE(~,y,par,mCurrent, ...
        TsIsAbsorbed,totalIsAbsorbed)
M = y(1);
if totalIsAbsorbed
    Ts = 0;
    Tu = 0;
elseif TsIsAbsorbed
    Ts = 0;
    Tu = y(3);
else
    Ts = y(2);
    Tu = y(3);
end
Di = y(4);
Dm = y(5);
E = y(6);
R = y(7);

T = Ts+Tu;
F = M/(M+par.a);
capacityFactor = 1-T/par.k;
immuneAttenuation = exp(-(T/par.k)*(R/par.b)*(1-F));
signal1 = par.beta1*(Tu+par.theta1*Ts);
signal2 = par.beta2*(Tu+par.theta1*Ts);
maturation = par.p3*signal1/(signal1+par.h) + ...
    par.p4*F*signal2/(signal2+par.h);

dydt = zeros(7,1);
dydt(1) = -par.mu1*M+mCurrent;
dydt(2) = (par.r*Ts*(1-F)-par.alpha*Ts)*capacityFactor ...
    -Ts*((1-par.theta2)*par.p5*E*immuneAttenuation ...
    +par.p1*F+par.mus);
dydt(3) = (par.r*Tu*(1-F)+par.alpha*Ts)*capacityFactor ...
    -Tu*(par.p5*E*immuneAttenuation+par.p2*F+par.muu);
dydt(4) = par.d0-par.mu2*Di-Di*maturation;
dydt(5) = Di*maturation-par.mu2*Dm;
dydt(6) = par.gamma*Dm-par.p6*R*E-par.mu3*E;
dydt(7) = par.eta*Di*((par.d0/par.mu2-Di)/(par.d0/par.mu2)) ...
    *(1-F)-par.mu4*R-par.p7*R*F;
if TsIsAbsorbed
    dydt(2) = 0;
end
if totalIsAbsorbed
    dydt(2:3) = 0;
end
end

% ========================================================================
% FIGURE 8 EXPORT
% ========================================================================

function Audit = exportMemorySafeFigure( ...
        fig,outputFolder,stem,dpi,maximumPixels,pdfContent)
pngPath = fullfile(outputFolder,[stem '.png']);
pdfPath = fullfile(outputFolder,[stem '.pdf']);
drawnow

% Limit PNG dimensions using the paper size and requested resolution.
canvas = double(fig.Position(3:4));
aspectRatio = canvas(1)/max(canvas(2),1);
if ~isfinite(aspectRatio) || aspectRatio<=0
    error('Cannot export %s because its canvas aspect ratio is invalid.',stem)
end
maximumInches = maximumPixels/dpi;
if aspectRatio>=1
    paperWidth = maximumInches;
    paperHeight = maximumInches/aspectRatio;
else
    paperHeight = maximumInches;
    paperWidth = maximumInches*aspectRatio;
end

oldPaperUnits = fig.PaperUnits;
oldPaperMode = fig.PaperPositionMode;
oldPaperPosition = fig.PaperPosition;
oldPaperSize = fig.PaperSize;
cleanup = onCleanup(@()restorePaperProperties( ...
    fig,oldPaperUnits,oldPaperMode,oldPaperPosition,oldPaperSize)); %#ok<NASGU>
set(fig,'Color','w','InvertHardcopy','off','PaperUnits','inches', ...
    'PaperPositionMode','manual', ...
    'PaperPosition',[0 0 paperWidth paperHeight], ...
    'PaperSize',[paperWidth paperHeight]);

deleteIfFileExists(pngPath)
try
    print(fig,pngPath,'-dpng',sprintf('-r%d',dpi),'-opengl');
catch ME
    error(['Memory-safe PNG export failed for %s at a bounded maximum ' ...
        'dimension of %d pixels: %s'],stem,maximumPixels,ME.message)
end
[valid,detail,width,height,blackFraction,whiteCorners] = ...
    validatePngAcrossCanvas(pngPath,maximumPixels);
if ~valid
    deleteIfFileExists(pngPath)
    error('Saved PNG validation failed for %s: %s.',stem,char(detail))
end
fprintf('Validated PNG: %s (%s).\n',pngPath,char(detail))

deleteIfFileExists(pdfPath)
try
    if strcmp(pdfContent,'vector')
        exportgraphics(fig,pdfPath,'ContentType','vector', ...
            'BackgroundColor','white');
    else
       % Export the PDF as a 600-dpi image.

        print(fig,pdfPath,'-dpdf','-r600','-opengl','-bestfit');
    end
catch ME
    warning('PDF export failed for %s; .fig and validated PNG remain: %s', ...
        stem,ME.message)
end

Figure = string(stem);
PngPath = string(pngPath);
WidthPixels = width;
HeightPixels = height;
RequestedDPI = dpi;
MaximumDimensionPixels = maximumPixels;
BlackPixelPercent = 100*blackFraction;
WhiteCorners = whiteCorners;
Passed = valid;
Detail = detail;
Audit = table(Figure,PngPath,WidthPixels,HeightPixels,RequestedDPI, ...
    MaximumDimensionPixels,BlackPixelPercent,WhiteCorners,Passed,Detail);
end

function restorePaperProperties( ...
        fig,units,mode,position,paperSize)
if ~isgraphics(fig)
    return
end
try
    set(fig,'PaperUnits',units,'PaperPositionMode',mode, ...
        'PaperPosition',position,'PaperSize',paperSize);
catch
end
end

function [valid,detail,width,height,blackFraction,cornersWhite] = ...
        validatePngAcrossCanvas(path,maximumPixels)
valid = false;
detail = "file missing";
width = 0;
height = 0;
blackFraction = nan;
cornersWhite = false;
fileInfo = dir(path);
if numel(fileInfo)~=1 || fileInfo.bytes<=1024
    return
end
try
    info = imfinfo(path);
catch ME
    detail = "imfinfo failed: "+string(ME.message);
    return
end
width = double(info(1).Width);
height = double(info(1).Height);
if max(width,height)>maximumPixels+8 || min(width,height)<1000
    detail = "unexpected dimensions "+width+" x "+height;
    return
end

sampleLimit = 700;
rowStep = max(1,ceil(height/sampleLimit));
columnStep = max(1,ceil(width/sampleLimit));
try
   % Read and subsample the PNG before conversion to double.
    [imageData,colorMap,alphaData] = imread(path);
    rowIndex = 1:rowStep:height;
    columnIndex = 1:columnStep:width;
    if ndims(imageData)==2
        imageData = imageData(rowIndex,columnIndex);
    else
        imageData = imageData(rowIndex,columnIndex,:);
    end
    if ~isempty(alphaData)
        if ndims(alphaData)==2
            alphaData = alphaData(rowIndex,columnIndex);
        else
            alphaData = alphaData(rowIndex,columnIndex,:);
        end
    end
catch ME
    detail = "bounded imread failed: "+string(ME.message);
    return
end
if isempty(colorMap)
    originalClass = class(imageData);
    if ndims(imageData)==2
        rgb = repmat(imageData,1,1,3);
    else
        rgb = imageData(:,:,1:min(3,size(imageData,3)));
        if size(rgb,3)==1
            rgb = repmat(rgb,1,1,3);
        end
    end
    rgb = double(rgb);
    if isinteger(imageData)
        rgb = rgb/double(intmax(originalClass));
    elseif max(rgb(:))>1
        rgb = rgb/255;
    end
else
    rgb = ind2rgb(imageData,colorMap);
end
if ~isempty(alphaData)
    alphaClass = class(alphaData);
    alpha = double(alphaData);
    if isinteger(alphaData)
        alpha = alpha/double(intmax(alphaClass));
    elseif max(alpha(:))>1
        alpha = alpha/255;
    end
    if ndims(alpha)>2
        alpha = alpha(:,:,1);
    end
    rgb = rgb.*alpha+(1-alpha);
end

black = all(rgb<0.03,3);
white = all(rgb>0.97,3);
nonwhite = any(rgb<0.97,3);
blackFraction = mean(black(:));
whiteFraction = mean(white(:));
nonwhiteFraction = mean(nonwhite(:));
cornerSize = min([3,size(rgb,1),size(rgb,2)]);
cornerMeans = [mean(rgb(1:cornerSize,1:cornerSize,:),'all'), ...
    mean(rgb(1:cornerSize,end-cornerSize+1:end,:),'all'), ...
    mean(rgb(end-cornerSize+1:end,1:cornerSize,:),'all'), ...
    mean(rgb(end-cornerSize+1:end,end-cornerSize+1:end,:),'all')];
cornersWhite = all(cornerMeans>0.80);

[contentRows,contentColumns] = find(nonwhite);
if isempty(contentRows)
    contentWidthFraction = 0;
    contentHeightFraction = 0;
else
    contentWidthFraction = ...
        (max(contentColumns)-min(contentColumns)+1)/size(nonwhite,2);
    contentHeightFraction = ...
        (max(contentRows)-min(contentRows)+1)/size(nonwhite,1);
end
valid = blackFraction<0.05 && whiteFraction>0.25 && ...
    nonwhiteFraction>0.001 && cornersWhite && ...
    contentWidthFraction>0.55 && contentHeightFraction>0.50;
detail = string(sprintf([ ...
    '%d x %d pixels; %d bytes; sampled %d x %d; black %.3f%%; ' ...
    'white %.3f%%; content span %.1f%% x %.1f%%; white corners %d'], ...
    width,height,fileInfo.bytes,size(rgb,2),size(rgb,1), ...
    100*blackFraction,100*whiteFraction,100*contentWidthFraction, ...
    100*contentHeightFraction,cornersWhite));
end

function deleteIfFileExists(path)
if exist(path,'file')==2
    delete(path)
end
end
