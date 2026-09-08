%This is an analysis code and part of two codes for Figures 3-5.
% Virtual-patient analyses for the seven-state NMIBC/MMC
% model. All primary analyses use the continuous ODE system. No one-cell
% absorbing event, compartment floor, day-0 clamp, rounding or postprocessing
% reset is applied. Sub-one compartment values remain continuous modeled
% cell-equivalent burdens.
%
%
% The 0.5-domain grids are exact subsets of the extended grids. Identical LHS
% backgrounds are reused over every r-alpha coordinate, and overlapping
% coordinates are simulated once. The alpha=1/28 day^-1 wall is an orientation
% guide separating slower cancer-associated transition values from the
% normal-tissue-referenced 4-28-day band.
%
% STATE ORDER AND GOVERNING SYSTEM
%   y = [M; Ts; Tu; Di; Dm; E; R]
%   T = Ts+Tu; F = M/(M+a)
%   Q = exp[-(T/k)(R/b)(1-F)]
%   q1 = beta1(Tu+theta1 Ts); q2 = beta2(Tu+theta1 Ts)
%   A = p3 q1/(q1+h) + p4 F q2/(q2+h)
%
%   dM/dt  = -mu1 M + mCurrent
%   dTs/dt = (r Ts(1-F)-alpha Ts)(1-T/k)
%             -Ts[(1-theta2)p5 E Q+p1 F+mus]
%   dTu/dt = (r Tu(1-F)+alpha Ts)(1-T/k)
%             -Tu[p5 E Q+p2 F+muu]
%   dDi/dt = d0-mu2 Di-Di A
%   dDm/dt = Di A-mu2 Dm
%   dE/dt  = gamma Dm-p6 R E-mu3 E
%   dR/dt  = eta Di[(d0/mu2-Di)/(d0/mu2)](1-F)-mu4 R-p7 R F
%
% MODEL ASSUMPTIONS
%   - MMC-induced growth arrest multiplies r*Ts and r*Tu by (1-F).
%     The alpha*Ts differentiation flux is not multiplied by (1-F).
%   - The one-way differentiation-like Ts-to-Tu transition remains inside the
%     shared logistic factor; no Tu-to-Ts transition is introduced.
%   - dE/dt contains no MMC-dependent multiplier.
%   - theta2 is the protected Ts fraction; Ts killing uses 1-theta2.
%   - p1, mus, p2 and muu are four independently sampled inputs.
%   - T0 and p8 use independent LHS coordinates. Ts0=p8*T0 and
%     Tu0=(1-p8)*T0 are retained even when either compartment is below one.
%   - T0 is sampled conditionally below k.
%   - Positive ranges spanning at least 100-fold use logarithmic LHS mapping;
%     all other sampled ranges use linear mapping.
%
% INPUT BOUNDS AND SAMPLING
%   mu1 [5.34e-5,92.736] log; m [5,71856.28743] log
%   r [0.00625,0.5] linear in the non-factorial virtual-patient cohorts
%   k [9e7,2.89e11] log; a [7.8,75] linear
%   theta1 [0.158412888,0.35016835] linear
%   theta2 [0.05,1] linear; p5 [1.4e-6,6e-6] linear
%   p1 [0.76,1.07], mus [0.044,0.052], p2 [1.657,2.367],
%   muu [0.043,0.059], independently sampled and linear
%   mu2 [0.023,0.98] linear; beta1 [7e-4,0.195] log
%   beta2 [0.0178,0.406] linear; h [1,2.89e11] log
%   p3 [0.013,1.347] log; p4 [0,34.04181843] linear
%   gamma [0.0047,9.12] log; mu3 [0.0105,0.247] linear
%   eta [0.013,0.212] linear; mu4 [0.0204,0.888] linear
%   T0 [2,min(2.88737e11,k)] conditional log
%   p8 [4e-6,0.798] independent log; Di0 [108,10000] linear
%   fixed d0=1.032e5, p6=1.44e-5, p7=3.110210655,
%   b=395840.674352314.
%
% Twelve two-hour MMC sessions begin on days
% 0,7,14,21,28,35,70,98,126,154,182,210.
% States are evaluated at days 0,30,...,360.
% Monthly detectability is recorded separately at each post-baseline visit.
% Cumulative attainment is calculated by the separate figures script.

clearvars
clc
close all

overallTimer = tic;

%% ========================================================================
% COMPUTATIONAL SETTINGS (not biological parameters)
% ========================================================================

% Modes change only cohort/grid resolution. Biology and endpoints never change.
settings.runMode = 'paper'; % run test first; then select 'fast' or 'paper'
settings.randomSeed = 22;
settings.nParallelWorkers = 2;
settings.showFigures = false;
settings.saveFigures = false;
settings.confirmPaperRun = true;

switch lower(settings.runMode)
    case 'test'
        settings.nSizePerGroup = 4;
        settings.nP8PerGroup = 3;
        settings.nPersistence = 8;
        settings.nBackgrounds = 2;
        settings.nRBase = 3;
        settings.nRExtra = 2;
        settings.nAlphaBase = 5;
        settings.nAlphaExtra = 2;
        settings.requestParallel = false;
        settings.exportResolution = 150;
    case 'fast'
        settings.nSizePerGroup = 50;
        settings.nP8PerGroup = 30;
        settings.nPersistence = 100;
        settings.nBackgrounds = 10;
        settings.nRBase = 7;
        settings.nRExtra = 3;
        settings.nAlphaBase = 10;
        settings.nAlphaExtra = 3;
        settings.requestParallel = true;
        settings.exportResolution = 300;
    case 'paper'
        settings.nSizePerGroup = 1000;
        settings.nP8PerGroup = 1000;
        settings.nPersistence = 1000;
        settings.nBackgrounds = 50;
        settings.nRBase = 17;
        settings.nRExtra = 5;
        settings.nAlphaBase = 20;
        settings.nAlphaExtra = 5;
        settings.requestParallel = true;
        settings.exportResolution = 600;
    otherwise
        error('settings.runMode must be ''test'', ''fast'' or ''paper''.')
end

settings.tEnd = 360;
settings.visitTimes = 0:30:settings.tEnd;
settings.monthDays = 30;
settings.nMonths = 12;
settings.tau = 2/24;
settings.doseDays = [0 7 14 21 28 35 70 98 126 154 182 210];
settings.detectionDiameterMM = 1;
settings.sizeThresholdMM = 30;
settings.detectionLimit = diameterToCells(settings.detectionDiameterMM);
settings.sizeThreshold = diameterToCells(settings.sizeThresholdMM);

settings.odeAbsTol = 1e-9;
settings.odeRelTol = 1e-6;
settings.odeMaxStep = 1;

settings.normalAlphaLower = 1/28;
settings.normalAlphaUpper = 1/4;
settings.theoreticalAlphaMinimum = 1/280;
settings.rAlphaPrimaryUpper = 0.5;
settings.rAlphaExtendedUpper = 1.512;

settings.outputDir = ...
    'VP_MMC_Continuous_IndependentP8_Dual_rAlpha_Outputs';
settings.figureFiles = { ...
    'FigVP01_Continuous_Monthly_Below_Detection.pdf', ...
    'FigVP02_Continuous_Total_and_Ts_Burden.pdf', ...
    'FigVP03_rAlpha_Through_0p5.pdf', ...
    'FigVP03_rAlpha_Through_1p512.pdf'};

if ~exist(settings.outputDir,'dir')
    mkdir(settings.outputDir)
end

rng(settings.randomSeed,'twister')
settings.useParallel = startParallelIfPossible( ...
    settings.requestParallel,settings.nParallelWorkers);

%========================================================================
% MODEL INPUT SPECIFICATION
% ========================================================================

bounds = makeBounds();
[rGrid,alphaGrid,alphaRegion,primaryRMask,primaryAlphaMask] = ...
    makeRAlphaGrid(settings,bounds);
settings.plannedRAlphaSimulationCount = settings.nBackgrounds* ...
    numel(rGrid)*numel(alphaGrid);
settings.plannedSimulationCount = 2*settings.nSizePerGroup + ...
    3*settings.nP8PerGroup + settings.nPersistence + ...
    settings.plannedRAlphaSimulationCount;
if strcmpi(settings.runMode,'paper') && ~settings.confirmPaperRun
    error(['Paper mode requires %d full 360-day ODE simulations. ' ...
        'Set settings.confirmPaperRun=true only after the test run passes ' ...
        'and this workload is accepted.'],settings.plannedSimulationCount)
end
validateSpecification(settings,bounds,rGrid,alphaGrid,alphaRegion, ...
    primaryRMask,primaryAlphaMask)
printSpecification(settings,bounds,rGrid,alphaGrid,alphaRegion, ...
    primaryRMask,primaryAlphaMask)

odeOptions = odeset( ...
    'AbsTol',settings.odeAbsTol, ...
    'RelTol',settings.odeRelTol, ...
    'MaxStep',settings.odeMaxStep, ...
    'NonNegative',1:7);

%% ANALYSIS 1: INDEPENDENT, EQUAL-SIZED BASELINE T0 GROUPS
fprintf('\nANALYSIS 1/4: independent equal-sized initial-size cohorts ...\n')
[sizeGroups,sizeCoordinates] = makeIndependentSizeCohorts( ...
    settings.nSizePerGroup,bounds,settings);
sizeTrajectories = runSizeCohorts(sizeGroups,settings,odeOptions);
sizeAnalysis = analyzeInitialSize(sizeGroups,sizeCoordinates, ...
    sizeTrajectories,settings);

%% ANALYSIS 2: FIXED BASELINE p8 TERTILES IN DISTINCT VIRTUAL PATIENTS
fprintf('\nANALYSIS 2/4: equal-sized baseline-p8 tertiles ...\n')
rng(settings.randomSeed+2,'twister')
[p8Patients,p8Coordinates] = makeCancerCohort( ...
    3*settings.nP8PerGroup,bounds,settings);
[p8Groups,p8GroupIndices,p8Cutpoints] = splitBaselineP8Tertiles( ...
    p8Patients,settings.nP8PerGroup);
p8Trajectories = runPatientCohort(p8Patients,settings,odeOptions, ...
    'baseline-p8 cohort');
p8Analysis = analyzeInitialP8(p8Patients,p8Coordinates,p8Groups, ...
    p8GroupIndices,p8Cutpoints,p8Trajectories,settings);

%% ANALYSIS 3: TEMPORAL Ts PERSISTENCE IN THE CANCER-ALPHA RANGE
fprintf('\nANALYSIS 3/4: temporal Ts-persistence cohort ...\n')
rng(settings.randomSeed+3,'twister')
[persistencePatients,persistenceCoordinates] = makeCancerCohort( ...
    settings.nPersistence,bounds,settings);
persistenceTrajectories = runPatientCohort(persistencePatients, ...
    settings,odeOptions,'Ts-persistence cohort');
persistenceAnalysis = analyzePersistence(persistencePatients, ...
    persistenceCoordinates,persistenceTrajectories,settings);

%% ANALYSIS 4: CONTROLLED THREE-DIMENSIONAL r-ALPHA EXPERIMENT
fprintf('\nANALYSIS 4/4: controlled 3-D r-alpha experiment ...\n')
rng(settings.randomSeed+4,'twister')
[backgrounds,lhsCoordinates] = makeBackgroundCohort( ...
    settings.nBackgrounds,bounds);
validateBackgroundCohort(backgrounds,bounds)
rAlphaExtended = runRAlphaGrid(backgrounds,lhsCoordinates,rGrid,alphaGrid, ...
    alphaRegion,settings,odeOptions);
rAlphaPrimary = subsetRAlphaResult(rAlphaExtended,primaryRMask, ...
    primaryAlphaMask,settings);
validateAgreedResult(rAlphaPrimary)
validateAgreedResult(rAlphaExtended)
validateCompleteAnalysisContract(sizeAnalysis,p8Analysis, ...
    persistenceAnalysis,rAlphaPrimary,rAlphaExtended,settings)

if settings.showFigures || settings.saveFigures
    figureHandles = makeSelectedFigures(sizeAnalysis,p8Analysis, ...
        persistenceAnalysis,rAlphaPrimary,rAlphaExtended,settings);
end

writetable(sizeAnalysis.resultsTable,fullfile(settings.outputDir, ...
    'Initial_size_monthly_below_detection_proportions.csv'))
writetable(sizeAnalysis.day360BurdenTable,fullfile(settings.outputDir, ...
    'Initial_size_day360_tumor_burden.csv'))
writetable(p8Analysis.resultsTable,fullfile(settings.outputDir, ...
    'Initial_Ts_proportion_monthly_below_detection_and_burden.csv'))
writetable(persistenceAnalysis.resultsTable,fullfile(settings.outputDir, ...
    'Temporal_Ts_and_tumor_persistence.csv'))
writetable(rAlphaPrimary.responseTable,fullfile(settings.outputDir, ...
    'rAlpha_Through_0p5_day360_absolute_cell_counts.csv'))
writetable(rAlphaPrimary.convergenceTable,fullfile(settings.outputDir, ...
    'rAlpha_Through_0p5_resolution_diagnostics.csv'))
writetable(rAlphaExtended.responseTable,fullfile(settings.outputDir, ...
    'rAlpha_Through_1p512_day360_absolute_cell_counts.csv'))
writetable(rAlphaExtended.convergenceTable,fullfile(settings.outputDir, ...
    'rAlpha_Through_1p512_resolution_diagnostics.csv'))
parameterTable = makeParameterSpecificationTable(bounds,settings);
writetable(parameterTable,fullfile(settings.outputDir, ...
    'Parameter_bounds_sampling_and_fixed_values.csv'))
if settings.saveFigures
    exportSelectedFigures(figureHandles,settings)
end

softwareVersion = version;
matFile = fullfile(settings.outputDir, ...
    'VP_MMC_continuous_independent_p8_dual_rAlpha_results.mat');
save(matFile,'settings','bounds','sizeAnalysis','p8Analysis', ...
    'persistenceAnalysis', ...
    'rAlphaPrimary','rAlphaExtended','rGrid','alphaGrid','alphaRegion', ...
    'primaryRMask','primaryAlphaMask','softwareVersion','-v7.3')

fprintf('\nAll four agreed analysis designs completed successfully.\n')
fprintf('Output directory: %s\n',settings.outputDir)
fprintf('Summary tables and simulation results saved.\n')
fprintf('Total elapsed time: %s\n',formatElapsedTime(toc(overallTimer)))

%% ========================================================================
% LOCAL FUNCTIONS: MODEL INPUTS AND VALIDATION
% ========================================================================

function bounds = makeBounds()

    % MMC kinetics and tumor dynamics
    bounds.mu1    = makeBound(5.34e-5,      92.736,       'log');
    bounds.m      = makeBound(5,            71856.28743,  'log');
    bounds.r      = makeBound(6.25e-3,      0.5,          'linear');
    bounds.k      = makeBound(9.0e7,        2.89e11,      'log');
    bounds.a      = makeBound(7.8,          75,           'linear');

    % Tumor/immune coupling. theta2 is the protected Ts fraction.
    bounds.theta1 = makeBound(0.158412888,  0.35016835,   'linear');
    bounds.theta2 = makeBound(0.05,         1,            'linear');
    bounds.p5     = makeBound(1.4e-6,       6.0e-6,       'linear');

    % Four independent literature-derived tumor loss inputs.
    bounds.p1     = makeBound(0.76,         1.07,         'linear');
    bounds.mus    = makeBound(0.044,        0.052,        'linear');
    bounds.p2     = makeBound(1.657,        2.367,        'linear');
    bounds.muu    = makeBound(0.043,        0.059,        'linear');

    % Dendritic-cell, effector-cell and Treg dynamics
    bounds.mu2    = makeBound(0.023,        0.98,         'linear');
    bounds.beta1  = makeBound(7.0e-4,       0.195,        'log');
    bounds.beta2  = makeBound(0.0178,       0.406,        'linear');
    bounds.h      = makeBound(1,            2.89e11,      'log');
    bounds.p3     = makeBound(0.013,        1.347,        'log');
    bounds.p4     = makeBound(0,            34.04181843,  'linear');
    bounds.gamma  = makeBound(0.0047,       9.12,         'log');
    bounds.mu3    = makeBound(0.0105,       0.247,        'linear');
    bounds.eta    = makeBound(0.013,        0.212,        'linear');
    bounds.mu4    = makeBound(0.0204,       0.888,        'linear');

    % Fixed quantities
    bounds.d0.fixed = 1.032e5;
    bounds.p6.fixed = 1.44e-5;
    bounds.p7.fixed = 3.110210655;
    bounds.b.fixed  = 395840.674352314;

    % Initial-condition inputs. T0 remains conditional on k; p8 is independent.
    bounds.T0  = makeBound(2,       2.88737e11, 'conditional-log');
    bounds.p8  = makeBound(4e-6,    0.798,      'log');
    bounds.Di0 = makeBound(108,     10000,      'linear');
end

function b = makeBound(lowerValue,upperValue,samplingMode)
    b.lo = lowerValue;
    b.hi = upperValue;
    b.mode = samplingMode;
end

function names = backgroundFactorNames()
    % r and alpha are intentionally absent: the response grid assigns them.
    names = {'mu1','m','k','a','theta1','theta2','p5','p1','mus', ...
        'p2','muu','mu2','beta1','beta2','h','p3','p4','gamma', ...
        'mu3','eta','mu4','T0','p8','Di0'};
end

function [rGrid,alphaGrid,alphaRegion,primaryRMask,primaryAlphaMask] = ...
        makeRAlphaGrid(settings,bounds)

    % Construct one nested grid. The complete grid reaches 1.512 day^-1;
    % exact subsets reach 0.5 day^-1. The two figures therefore share every
    % overlapping simulation and differ only in displayed controlled domain.
    baseR = linspace(bounds.r.lo,settings.rAlphaPrimaryUpper, ...
        settings.nRBase);
    extraRWithBoundary = logspace(log10(settings.rAlphaPrimaryUpper), ...
        log10(settings.rAlphaExtendedUpper),settings.nRExtra+1);
    rGrid = unique([baseR extraRWithBoundary(2:end)]);

    baseAlpha = logspace(log10(settings.theoreticalAlphaMinimum), ...
        log10(settings.rAlphaPrimaryUpper),settings.nAlphaBase);
    baseAlpha = unique([baseAlpha settings.normalAlphaLower, ...
        settings.normalAlphaUpper,settings.rAlphaPrimaryUpper]);
    extraAlphaWithBoundary = logspace(log10(settings.rAlphaPrimaryUpper), ...
        log10(settings.rAlphaExtendedUpper),settings.nAlphaExtra+1);
    alphaGrid = unique([baseAlpha extraAlphaWithBoundary(2:end)]);

    primaryRMask = rGrid<=settings.rAlphaPrimaryUpper+1e-14;
    primaryAlphaMask = alphaGrid<=settings.rAlphaPrimaryUpper+1e-14;
    alphaRegion = strings(size(alphaGrid));
    alphaRegion(alphaGrid<settings.normalAlphaLower) = ...
        "cancer-associated slowing (>28-day reference)";
    alphaRegion(alphaGrid>=settings.normalAlphaLower & ...
        alphaGrid<=settings.normalAlphaUpper) = ...
        "normal-tissue-referenced characteristic transition (4-28 days)";
    alphaRegion(alphaGrid>settings.normalAlphaUpper) = ...
        "theoretical faster-than-4-day transition";
end

function validateSpecification(settings,bounds,rGrid,alphaGrid,alphaRegion, ...
        primaryRMask,primaryAlphaMask)

    validateModeSettings(settings)

% Validate parameter and initial-condition bounds.
    names = {'mu1','m','r','k','a','theta1','theta2','p5','p1','mus', ...
        'p2','muu','mu2','beta1','beta2','h','p3','p4','gamma', ...
        'mu3','eta','mu4','T0','p8','Di0'};

    expectedLower = [5.34e-5,5,6.25e-3,9.0e7,7.8, ...
        0.158412888,0.05,1.4e-6,0.76,0.044,1.657,0.043, ...
        0.023,7.0e-4,0.0178,1,0.013,0,0.0047,0.0105, ...
        0.013,0.0204,2,4e-6,108];

    expectedUpper = [92.736,71856.28743,0.5,2.89e11,75, ...
        0.35016835,1,6.0e-6,1.07,0.052,2.367,0.059, ...
        0.98,0.195,0.406,2.89e11,1.347,34.04181843,9.12, ...
        0.247,0.212,0.888,2.88737e11,0.798,10000];

    actualLower = zeros(size(expectedLower));
    actualUpper = zeros(size(expectedUpper));
    for j = 1:numel(names)
        actualLower(j) = bounds.(names{j}).lo;
        actualUpper(j) = bounds.(names{j}).hi;
    end
    if ~isequal(actualLower,expectedLower) || ...
            ~isequal(actualUpper,expectedUpper)
        error('One or more locked ranged-input values changed.')
    end

    fixedValues = [bounds.d0.fixed,bounds.p6.fixed, ...
        bounds.p7.fixed,bounds.b.fixed];
    expectedFixed = [1.032e5,1.44e-5,3.110210655,395840.674352314];
    if ~isequal(fixedValues,expectedFixed)
        error('One or more locked fixed values changed.')
    end

% Check sampling scales. Only T0 is conditional.
    conditionalNames = {'T0'};
    sampledNames = backgroundFactorNames();
    for j = 1:numel(sampledNames)
        name = sampledNames{j};
        if ismember(name,conditionalNames)
            continue
        end
        expectedMode = declaredScale(bounds.(name).lo,bounds.(name).hi);
        if ~strcmp(bounds.(name).mode,expectedMode)
            error('Sampling mode mismatch for %s: expected %s.', ...
                name,expectedMode)
        end
    end
    if ~strcmp(bounds.r.mode,'linear') || bounds.r.hi/bounds.r.lo~=80
        error('r must be linear on the agreed 80-fold range.')
    end

    % Check the supported normal alpha band and transparent simulated span.
    expectedAlphaMinimum = 1/280;
    if settings.normalAlphaLower~=1/28 || settings.normalAlphaUpper~=1/4
        error(['Normal-tissue-referenced alpha must correspond exactly to ' ...
            '4-28-day characteristic timescales.'])
    end
    if settings.theoreticalAlphaMinimum~=expectedAlphaMinimum || ...
            abs(alphaGrid(1)-expectedAlphaMinimum)>1e-14
        error('The theoretical alpha display minimum must be 1/280 day^-1.')
    end
    if any(alphaGrid<=0) || any(diff(alphaGrid)<=0)
        error('The alpha grid must be strictly positive and increasing.')
    end
    if sum(abs(alphaGrid-1/28)<1e-14)~=1
        error('alpha=1/28 must occur exactly once.')
    end
    if any(alphaGrid(alphaRegion== ...
            "cancer-associated slowing (>28-day reference)")>=1/28) || ...
       any(alphaGrid(alphaRegion== ...
            "normal-tissue-referenced characteristic transition (4-28 days)")<1/28) || ...
       any(alphaGrid(alphaRegion== ...
            "normal-tissue-referenced characteristic transition (4-28 days)")>1/4)
        error('The pathological and normal alpha regions overlap.')
    end
    if abs(rGrid(find(primaryRMask,1,'last'))-0.5)>1e-14 || ...
            abs(alphaGrid(find(primaryAlphaMask,1,'last'))-0.5)>1e-14
        error('The primary r-alpha domain must end at 0.5 day^-1 on both axes.')
    end
    if abs(rGrid(end)-1.512)>1e-14 || abs(alphaGrid(end)-1.512)>1e-14
        error('The extended r-alpha domain must end at 1.512 day^-1 on both axes.')
    end

    % Solver settings and treatment schedule.
    if settings.randomSeed~=22
        error('The locked random seed is 22.')
    end
    if settings.tEnd~=360 || settings.tau~=2/24 || ...
            ~isequal(settings.doseDays, ...
            [0 7 14 21 28 35 70 98 126 154 182 210])
        error('The endpoint, pulse duration or 12-session schedule changed.')
    end
    if ~isequal(settings.visitTimes,0:30:360)
        error('Monthly evaluation times must be 0,30,...,360 days.')
    end
    if settings.odeAbsTol~=1e-9 || settings.odeRelTol~=1e-6 || ...
            settings.odeMaxStep~=1
        error('The agreed ode15s controls changed.')
    end
    if numel(backgroundFactorNames())~=24
        error('The background LHS must contain 24 coordinates excluding r and alpha.')
    end

    % Independent-input guard for the four tumor loss parameters.
    if isequal([bounds.p1.lo bounds.p1.hi],[bounds.p2.lo bounds.p2.hi]) || ...
       isequal([bounds.mus.lo bounds.mus.hi],[bounds.muu.lo bounds.muu.hi])
        error('p1/p2 or mus/muu were incorrectly aliased.')
    end

    validateEquationRules(bounds)

    if rGrid(1)~=bounds.r.lo
        error('The controlled r grid does not begin at the agreed minimum.')
    end
end

function validateModeSettings(settings)

    if ~ismember(lower(settings.runMode),{'test','fast','paper'})
        error('Unknown run mode: %s.',settings.runMode)
    end
    counts = [settings.nSizePerGroup,settings.nP8PerGroup, ...
        settings.nPersistence,settings.nBackgrounds,settings.nRBase, ...
        settings.nRExtra,settings.nAlphaBase,settings.nAlphaExtra, ...
        settings.nParallelWorkers,settings.exportResolution];
    if any(~isfinite(counts)) || any(counts<1) || any(counts~=round(counts))
        error('All computational counts and resolutions must be positive integers.')
    end

    switch lower(settings.runMode)
        case 'test'
            expectedCounts = [4,3,8,2,3,2,5,2,2,150];
        case 'fast'
            expectedCounts = [50,30,100,10,7,3,10,3,2,300];
        case 'paper'
            expectedCounts = [1000,1000,1000,50,17,5,20,5,2,600];
    end
    if ~isequal(counts,expectedCounts)
        error('The computational counts do not match the agreed %s-mode design.', ...
            lower(settings.runMode))
    end
    if settings.plannedRAlphaSimulationCount<1 || ...
            settings.plannedSimulationCount~=2*settings.nSizePerGroup+ ...
            3*settings.nP8PerGroup+settings.nPersistence+ ...
            settings.plannedRAlphaSimulationCount
        error('The declared simulation count is inconsistent with the grids.')
    end
end

function validateEquationRules(bounds)

% Check modelRHS against the declared equations at a fixed test state. 

    p = midpointParameterSet(bounds);
    y = [10;2e6;3e6;500;250;100;50];

    mAudit = p.m;
    M = y(1); Ts = y(2); Tu = y(3);
    Di = y(4); Dm = y(5); E = y(6); R = y(7);
    T = Ts+Tu;
    F = M/(M+p.a);
    Q = exp(-((T/p.k)*(R/p.b)*(1-F)));
    q1 = p.beta1*(Tu+p.theta1*Ts);
    q2 = p.beta2*(Tu+p.theta1*Ts);
    A = p.p3*q1/(q1+p.h)+p.p4*F*q2/(q2+p.h);
    expected = [ ...
        -p.mu1*M+mAudit; ...
        (p.r*Ts*(1-F)-p.alpha*Ts)*(1-T/p.k) ...
            - Ts*((1-p.theta2)*p.p5*E*Q+p.p1*F+p.mus); ...
        (p.r*Tu*(1-F)+p.alpha*Ts)*(1-T/p.k) ...
            - Tu*(p.p5*E*Q+p.p2*F+p.muu); ...
        p.d0-p.mu2*Di-Di*A; ...
        Di*A-p.mu2*Dm; ...
        p.gamma*Dm-p.p6*R*E-p.mu3*E; ...
        p.eta*Di*((p.d0/p.mu2-Di)/(p.d0/p.mu2))*(1-F) ...
            - p.mu4*R-p.p7*R*F];
    actual = modelRHS(0,y,p,mAudit);
    equationTolerance = 1e-12*max(1,max(abs(expected)));
    if max(abs(actual-expected))>equationTolerance
        error('The implemented seven-state system differs from the declared system.')
    end

    % Check MMC independence of dE/dt and theta2-mediated protection.

    dAtM = modelRHS(0,y,p,0);
    yNoM = y;
    yNoM(1) = 0;
    dWithoutM = modelRHS(0,yNoM,p,0);
    if abs(dAtM(6)-dWithoutM(6))>1e-12
        error('dE/dt must not contain an MMC-dependent multiplier.')
    end

    pLow = p;
    pLow.theta2 = bounds.theta2.lo;
    pHigh = p;
    pHigh.theta2 = bounds.theta2.hi;
    dLow = modelRHS(0,y,pLow,0);
    dHigh = modelRHS(0,y,pHigh,0);
    if dHigh(2)<=dLow(2)
        error('theta2 protection is reversed; Ts killing must use 1-theta2.')
    end

end

function validateAgreedResult(result)

    expectedColumns = {'r_day_inverse','alpha_day_inverse', ...
        'DifferentiationLikeTransitionTime_days','AlphaRegion', ...
        'N_common_backgrounds', ...
        'MeanTotalTumorCellsDay360','MedianTotalTumorCellsDay360', ...
        'MeanTsCellsDay360','MedianTsCellsDay360', ...
        'MeanTuCellsDay360','MedianTuCellsDay360', ...
        'ReferenceMeanTotalTumorCellsAt28DayRate', ...
        'AbsoluteTotalTumorReductionFrom28DayRate', ...
        'RelativeTotalTumorReductionFrom28DayRate', ...
        'NormalReferenceReductionDomain'};
    if ~isequal(result.responseTable.Properties.VariableNames,expectedColumns)
        error('The output table contains an unagreed or missing endpoint.')
    end

    componentDifference = abs(result.meanTotalTumorCells- ...
        (result.meanTsCells+result.meanTuCells));
    allowedDifference = 1e-10*max(1,max(result.meanTotalTumorCells(:)));
    if max(componentDifference(:))>allowedDifference
        error('Total tumor burden must equal Ts+Tu at every reported point.')
    end

    referenceIndex = find(abs(result.alphaGrid-1/28)<1e-14);
    if numel(referenceIndex)~=1 || result.normalReferenceIndex~=referenceIndex
        error('The normal-reference boundary must occur exactly at alpha=1/28.')
    end
    expectedDomain = repmat(result.alphaGrid(:)>1/28 & ...
        result.alphaGrid(:)<=1/4, ...
        1,numel(result.rGrid)) & ...
        result.absoluteTotalTumorReduction>0;
    if ~isequal(result.normalReferenceReductionDomain,expectedDomain)
        error('The normal-reference reduction domain is inconsistent.')
    end
    if any(result.normalReferenceReductionDomain(referenceIndex,:)) || ...
            any(any(result.normalReferenceReductionDomain( ...
            1:referenceIndex-1,:)))
        error('The reduction domain must lie strictly above alpha=1/28.')
    end
    referenceDifference = abs(result.absoluteTotalTumorReduction( ...
        referenceIndex,:));
    referenceTolerance = 1e-12*max(1,max(result.referenceMeanTotalTumorCells));
    if any(referenceDifference>referenceTolerance)
        error('The reduction must be zero at the 28-day reference boundary.')
    end

    forbiddenResultFields = {'meanPatientTsFraction','meanTumorBurdenAUC', ...
        'detectionProbability','cumulativeIncidence','incidenceRateRatio'};
    if any(isfield(result,forbiddenResultFields))
        error('The result structure contains an explicitly excluded analysis.')
    end
end

function validateCompleteAnalysisContract(sizeAnalysis,p8Analysis, ...
        persistenceAnalysis,rAlphaPrimary,rAlphaExtended,settings)
    if height(sizeAnalysis.resultsTable)~=settings.nMonths || ...
            any(sizeAnalysis.resultsTable.N_Small~=settings.nSizePerGroup) || ...
            any(sizeAnalysis.resultsTable.N_Large~=settings.nSizePerGroup) || ...
            size(sizeAnalysis.monthlyDetectableStatus{1},2)~=settings.nMonths || ...
            size(sizeAnalysis.monthlyDetectableStatus{2},2)~=settings.nMonths
        error('The initial-size analysis is not the agreed balanced experiment.')
    end
    if any(p8Analysis.groupSizes~=settings.nP8PerGroup) || ...
            numel(p8Analysis.groupLabels)~=3 || ...
            numel(p8Analysis.monthlyDetectableStatus)~=3 || ...
            any(cellfun(@(status)size(status,2), ...
            p8Analysis.monthlyDetectableStatus)~=settings.nMonths)
        error('The baseline-p8 analysis must contain three equal fixed groups.')
    end
    if size(persistenceAnalysis.totalTumorCells,1)~=settings.nPersistence || ...
            size(persistenceAnalysis.totalTumorCells,2)~=13
        error('The temporal persistence cohort has the wrong size or time grid.')
    end
    rAlphaResults = {rAlphaPrimary,rAlphaExtended};
    for j = 1:2
        result = rAlphaResults{j};
        if result.endpointDay~=360 || any(result.alphaGrid<=0) || ...
                ~ismatrix(result.meanTotalTumorCells) || ...
                ~ismatrix(result.meanTsCells)
            error('An r-alpha result is not the agreed 3-D day-360 endpoint.')
        end
    end
    if abs(rAlphaPrimary.rGrid(end)-0.5)>1e-14 || ...
            abs(rAlphaPrimary.alphaGrid(end)-0.5)>1e-14
        error('The primary r-alpha figure must end at 0.5 on both axes.')
    end
    if abs(rAlphaExtended.rGrid(end)-1.512)>1e-14 || ...
            abs(rAlphaExtended.alphaGrid(end)-1.512)>1e-14
        error('The extended r-alpha figure must end at 1.512 on both axes.')
    end
    rMask = ismember(rAlphaExtended.rGrid,rAlphaPrimary.rGrid);
    alphaMask = ismember(rAlphaExtended.alphaGrid,rAlphaPrimary.alphaGrid);
    if ~isequal(rAlphaPrimary.meanTotalTumorCells, ...
            rAlphaExtended.meanTotalTumorCells(alphaMask,rMask)) || ...
       ~isequal(rAlphaPrimary.meanTsCells, ...
            rAlphaExtended.meanTsCells(alphaMask,rMask))
        error(['The 0.5-domain surfaces must be exact subsets of the ' ...
            '1.512-domain simulations.'])
    end
    forbiddenFields = {'untreated','founder','alphaZero', ...
        'subtraction','hazard','AUC','selectedTrajectories', ...
        'eventMonth','eventStatus','cumulativeIncidence','incidenceRate', ...
        'incidenceRateRatio','finalIRR','upperVsLowerIRR'};
    structures = {sizeAnalysis,p8Analysis,persistenceAnalysis, ...
        rAlphaPrimary,rAlphaExtended};
    for s = 1:numel(structures)
        if any(isfield(structures{s},forbiddenFields))
            error('An explicitly excluded analysis was added to the results.')
        end
    end
end

function p = midpointParameterSet(bounds)

    rangedNames = {'mu1','m','r','k','a','theta1','theta2','p5', ...
        'p1','mus','p2','muu','mu2','beta1','beta2','h','p3','p4', ...
        'gamma','mu3','eta','mu4'};
    for j = 1:numel(rangedNames)
        name = rangedNames{j};
        p.(name) = (bounds.(name).lo+bounds.(name).hi)/2;
    end
    p.alpha = 1/28;
    p.d0 = bounds.d0.fixed;
    p.p6 = bounds.p6.fixed;
    p.p7 = bounds.p7.fixed;
    p.b = bounds.b.fixed;
end

function printSpecification(settings,bounds,rGrid,alphaGrid,alphaRegion, ...
        primaryRMask,primaryAlphaMask)

    nPathological = sum(alphaRegion== ...
        "cancer-associated slowing (>28-day reference)");
    nNormal = sum(alphaRegion== ...
        "normal-tissue-referenced characteristic transition (4-28 days)");

    fprintf('\n============================================================\n')
    fprintf('AGREED VIRTUAL-PATIENT PUBLICATION ANALYSES\n')
    fprintf('============================================================\n')
    fprintf('Run mode: %s\n',upper(settings.runMode))
    fprintf('Initial-size experiment: %d independent patients per group\n', ...
        settings.nSizePerGroup)
    fprintf(['Initial-size endpoint: proportion currently below the 1-mm ' ...
        'threshold at each scheduled month, days 30-360 (day 0 excluded)\n'])
    fprintf(['Initial-Ts-proportion experiment: three fixed day-0 thirds, ' ...
        'n=%d each\n'], ...
        settings.nP8PerGroup)
    fprintf(['Every patient contributes one detectability status at every ' ...
        'post-baseline visit; baseline groups remain fixed\n'])
    fprintf('Temporal continuous-burden cancer-range cohort: N=%d\n', ...
        settings.nPersistence)
    fprintf('Temporal responses: continuous absolute Ts+Tu and Ts burdens\n')
    fprintf(['r-alpha experiment: %d independent 24-input backgrounds; ' ...
        'each is reused over the full grid\n'],settings.nBackgrounds)
    fprintf(['At each grid coordinate, mean and median cell counts are ' ...
        'calculated across those %d backgrounds\n'],settings.nBackgrounds)
    fprintf('Planned full ODE simulations: %d\n',settings.plannedSimulationCount)
    fprintf(['Primary r-alpha figure: both axes %.6g to %.6g day^-1 ' ...
        '(%d r x %d alpha values)\n'],rGrid(1),0.5, ...
        sum(primaryRMask),sum(primaryAlphaMask))
    fprintf(['Extended r-alpha figure: both axes %.6g to %.6g day^-1 ' ...
        '(%d r x %d alpha values)\n'],rGrid(1),1.512, ...
        numel(rGrid),numel(alphaGrid))
    fprintf(['  cancer-associated >28-day-reference region: %d values, ' ...
        'strictly below %.6g\n'], ...
        nPathological,settings.normalAlphaLower)
    fprintf(['  normal-tissue-referenced characteristic 4-28-day region: ' ...
        '%d values on ' ...
        '[%.6g, %.6g]\n'], ...
        nNormal,settings.normalAlphaLower,settings.normalAlphaUpper)
    fprintf('p1 [%.3g,%.3g], mus [%.3g,%.3g], p2 [%.3g,%.3g], muu [%.3g,%.3g]\n', ...
        bounds.p1.lo,bounds.p1.hi,bounds.mus.lo,bounds.mus.hi, ...
        bounds.p2.lo,bounds.p2.hi,bounds.muu.lo,bounds.muu.hi)
    fprintf('theta2: protected Ts fraction on [%.2f,%.2f]; killing uses 1-theta2\n', ...
        bounds.theta2.lo,bounds.theta2.hi)
    fprintf('MMC: 12 sessions, each %.3g day (2 hours)\n',settings.tau)
    fprintf(['Continuous ODE: no one-cell event, clamp, rounding or ' ...
        'postprocessing reset\n'])
    fprintf(['T0 and p8: independent LHS coordinates; sub-one Ts0 or Tu0 ' ...
        'is retained continuously\n'])
    fprintf('Figure export resolution: %d dpi\n',settings.exportResolution)
    fprintf(['Selected figures: monthly detectability by initial size and ' ...
        'initial Ts proportion; temporal continuous burdens; and two ' ...
        'separate two-panel mean 3-D r-alpha responses\n'])
    fprintf(['Rejected: every 2-D r-alpha display, alpha=0, treatment or ' ...
        'alpha=r curve/classification, alpha-minus-alpha0 subtraction, ' ...
        'founder, untreated, first-event, ' ...
        'cumulative-incidence, person-time, IRR, hazard and AUC analyses\n'])
    fprintf('Parallel requested: %d | parallel used: %d\n', ...
        settings.requestParallel,settings.useParallel)
    fprintf('Specification and equation-level audits passed.\n')
    fprintf('============================================================\n\n')
end

function parameterTable = makeParameterSpecificationTable(bounds,settings)

    names = {'mu1','m','r','k','a','alpha','theta1','theta2','p5', ...
        'p1','mus','p2','muu','mu2','beta1','beta2','h','p3','p4', ...
        'gamma','mu3','eta','mu4','T0','p8','Di0','d0','p6','p7','b'}.';
    units = {'day^-1','micromolar day^-1','day^-1','cells','micromolar', ...
        'day^-1','fraction','fraction','cells^-1 day^-1','day^-1', ...
        'day^-1','day^-1','day^-1','day^-1','dimensionless', ...
        'dimensionless','cells','day^-1','day^-1','day^-1','day^-1', ...
        'day^-1','day^-1','cells','fraction','cells','cells day^-1', ...
        'cells^-1 day^-1','day^-1','cells'}.';
    lower = nan(numel(names),1);
    upper = nan(numel(names),1);
    sampling = strings(numel(names),1);
    status = repmat("sampled",numel(names),1);
    for j = 1:numel(names)
        name = names{j};
        if strcmp(name,'alpha')
            lower(j) = settings.theoreticalAlphaMinimum;
            upper(j) = settings.normalAlphaLower;
            sampling(j) = "log (upper boundary excluded)";
            status(j) = "sampled cancer-associated cohort input";
        elseif isfield(bounds.(name),'fixed')
            lower(j) = bounds.(name).fixed;
            upper(j) = bounds.(name).fixed;
            sampling(j) = "fixed";
            status(j) = "fixed";
        else
            lower(j) = bounds.(name).lo;
            upper(j) = bounds.(name).hi;
            sampling(j) = string(bounds.(name).mode);
        end
    end
    parameterTable = table(string(names),lower,upper,sampling,string(units), ...
        status,'VariableNames',{'Input','Lower','Upper','Sampling','Units','Status'});
    controlledRows = table( ...
        ["r_alpha_primary";"alpha_primary"; ...
         "r_alpha_extended";"alpha_extended"], ...
        repmat(settings.theoreticalAlphaMinimum,4,1), ...
        [settings.rAlphaPrimaryUpper;settings.rAlphaPrimaryUpper; ...
         settings.rAlphaExtendedUpper;settings.rAlphaExtendedUpper], ...
        ["controlled nested grid";"controlled nested logarithmic grid"; ...
         "controlled nested grid";"controlled nested logarithmic grid"], ...
        repmat("day^-1",4,1), ...
        ["Figure 3 domain";"Figure 3 domain"; ...
         "extended Figure 3 domain";"extended Figure 3 domain"], ...
        'VariableNames',parameterTable.Properties.VariableNames);
    controlledRows.Lower([1 3]) = bounds.r.lo;
    parameterTable = [parameterTable;controlledRows];
end

function modeName = declaredScale(lowerValue,upperValue)
    if lowerValue>0 && upperValue/lowerValue>=100
        modeName = 'log';
    else
        modeName = 'linear';
    end
end

%% ========================================================================
% LOCAL FUNCTIONS: LHS BACKGROUND COHORT
% ========================================================================

function [backgrounds,U] = makeBackgroundCohort(nBackgrounds,bounds)

    factorCount = numel(backgroundFactorNames());
    U = simpleLHS(nBackgrounds,factorCount);
    backgrounds = cell(nBackgrounds,1);
    for i = 1:nBackgrounds
        backgrounds{i} = sampleBackground(U(i,:),bounds);
    end
end

function p = sampleBackground(u,bounds)

    index = 1;
    p.mu1    = sampleFromBound(u(index),bounds.mu1);    index=index+1;
    p.m      = sampleFromBound(u(index),bounds.m);      index=index+1;
    p.k      = sampleFromBound(u(index),bounds.k);      index=index+1;
    p.a      = sampleFromBound(u(index),bounds.a);      index=index+1;
    p.theta1 = sampleFromBound(u(index),bounds.theta1); index=index+1;
    p.theta2 = sampleFromBound(u(index),bounds.theta2); index=index+1;
    p.p5     = sampleFromBound(u(index),bounds.p5);     index=index+1;
    p.p1     = sampleFromBound(u(index),bounds.p1);     index=index+1;
    p.mus    = sampleFromBound(u(index),bounds.mus);    index=index+1;
    p.p2     = sampleFromBound(u(index),bounds.p2);     index=index+1;
    p.muu    = sampleFromBound(u(index),bounds.muu);    index=index+1;
    p.mu2    = sampleFromBound(u(index),bounds.mu2);    index=index+1;
    p.beta1  = sampleFromBound(u(index),bounds.beta1);  index=index+1;
    p.beta2  = sampleFromBound(u(index),bounds.beta2);  index=index+1;
    p.h      = sampleFromBound(u(index),bounds.h);      index=index+1;
    p.p3     = sampleFromBound(u(index),bounds.p3);     index=index+1;
    p.p4     = sampleFromBound(u(index),bounds.p4);     index=index+1;
    p.gamma  = sampleFromBound(u(index),bounds.gamma);  index=index+1;
    p.mu3    = sampleFromBound(u(index),bounds.mu3);    index=index+1;
    p.eta    = sampleFromBound(u(index),bounds.eta);    index=index+1;
    p.mu4    = sampleFromBound(u(index),bounds.mu4);    index=index+1;

    uT0 = u(index); index=index+1;
    uP8 = u(index); index=index+1;
    uDi0 = u(index);

    p.d0 = bounds.d0.fixed;
    p.p6 = bounds.p6.fixed;
    p.p7 = bounds.p7.fixed;
    p.b = bounds.b.fixed;

    % r and alpha are assigned later by the common factorial grid.
    p.r = NaN;
    p.alpha = NaN;

    % T0 uses log sampling on [2,min(2.88737e11,k)]. Because LHS points
    % lie strictly inside (0,1), the sampled T0 remains strictly below k.
    T0Upper = min(bounds.T0.hi,p.k);
    p.T0 = sampleLinearOrLog(uT0,bounds.T0.lo,T0Upper,'log');

    % p8 is an independent LHS input. Continuous sub-one initial compartment
    % burdens are retained; there is no conditional p8 interval or resampling.
    p.p8 = sampleFromBound(uP8,bounds.p8);
    p.Di0 = sampleFromBound(uDi0,bounds.Di0);

    p.Ts0 = p.p8*p.T0;
    p.Tu0 = (1-p.p8)*p.T0;
    p.y0 = [0;p.Ts0;p.Tu0;p.Di0;1;1;1];
end

function validateBackgroundCohort(backgrounds,bounds)

    T0 = cellfun(@(p)p.T0,backgrounds);
    k = cellfun(@(p)p.k,backgrounds);
    p8 = cellfun(@(p)p.p8,backgrounds);
    Ts0 = cellfun(@(p)p.Ts0,backgrounds);
    Tu0 = cellfun(@(p)p.Tu0,backgrounds);
    Dm0 = cellfun(@(p)p.y0(5),backgrounds);
    E0 = cellfun(@(p)p.y0(6),backgrounds);
    R0 = cellfun(@(p)p.y0(7),backgrounds);

    tolerance = 1e-10;
    if any(T0<bounds.T0.lo) || any(T0>=k) || any(T0>bounds.T0.hi)
        error('Conditional T0 sampling failed.')
    end
    if any(p8<bounds.p8.lo-tolerance) || ...
            any(p8>bounds.p8.hi+tolerance)
        error('Independent p8 sampling left the agreed global range.')
    end
    if any(Ts0<0) || any(Tu0<0)
        error('Continuous initial tumor-compartment burdens must be nonnegative.')
    end
    if any(abs(Ts0+Tu0-T0)>tolerance.*max(1,T0))
        error('Ts(0)+Tu(0) does not equal T0.')
    end
    if any(Dm0~=1) || any(E0~=1) || any(R0~=1)
        error('Dm(0), E(0) and R(0) must equal one.')
    end
end

function U = simpleLHS(n,d)

    U = zeros(n,d);
    for j = 1:d
        permutation = randperm(n).';
        U(:,j) = (permutation-1+rand(n,1))/n;
    end
end

function value = sampleFromBound(u,bound)
    value = sampleLinearOrLog(u,bound.lo,bound.hi,bound.mode);
end

function value = sampleLinearOrLog(u,lowerValue,upperValue,modeName)

    u = min(max(u,1e-12),1-1e-12);
    if upperValue<lowerValue
        error('Invalid sampling interval [%g,%g].',lowerValue,upperValue)
    elseif upperValue==lowerValue
        value = lowerValue;
    elseif strcmpi(modeName,'log')
        if lowerValue<=0
            error('Log sampling requires a strictly positive lower endpoint.')
        end
        value = 10^(log10(lowerValue)+u* ...
            (log10(upperValue)-log10(lowerValue)));
    elseif strcmpi(modeName,'linear')
        value = lowerValue+u*(upperValue-lowerValue);
    else
        error('Unknown sampling mode: %s.',modeName)
    end
end

function cells = diameterToCells(diameterMM)
    % The same declared geometric conversion is used for the 1-mm
    % detectability threshold and the 3-cm initial-size boundary: circular
    % footprint, 0.03-mm (three-cell) depth and 10^6 cells/mm^3.
    cells = pi*(diameterMM/2)^2*0.03*1e6;
end

%% ========================================================================
% LOCAL FUNCTIONS: INDEPENDENT CANCER-RANGE COHORTS
% ========================================================================

function [patients,U] = makeCancerCohort(n,bounds,settings)

    % Twenty-four background coordinates retain conditional T0 but use an
    % independent p8 coordinate. Two additional LHS coordinates assign r
    % and cancer-associated alpha. Normal-reference alpha values are used
    % only in the controlled r-alpha grid, never in these patient cohorts.
    [patients,backgroundU] = makeBackgroundCohort(n,bounds);
    kineticU = simpleLHS(n,2);
    alphaMinimum = settings.theoreticalAlphaMinimum;
    for i = 1:n
        p = patients{i};
        p.r = sampleFromBound(kineticU(i,1),bounds.r);
        p.alpha = sampleLinearOrLog(kineticU(i,2),alphaMinimum, ...
            settings.normalAlphaLower,'log');
        patients{i} = p;
    end
    U = [backgroundU kineticU];
    validateCancerCohort(patients,bounds,settings)
end

function [groups,coordinates] = makeIndependentSizeCohorts(n,bounds,settings)

    % The two conditions contain distinct LHS rows. A patient is sampled once
    % with one baseline T0; no background is evaluated under both conditions.
    rng(settings.randomSeed,'twister')
    [groups{1},coordinates{1}] = makeSizeCohort( ...
        n,bounds,settings,'small');
    rng(settings.randomSeed+1,'twister')
    [groups{2},coordinates{2}] = makeSizeCohort( ...
        n,bounds,settings,'large');
    validateIndependentSizeCohorts(groups,coordinates,bounds,settings)
end

function [patients,U] = makeSizeCohort(n,bounds,settings,condition)

    [patients,unusedU] = makeCancerCohort(n,bounds,settings);
    controlledU = simpleLHS(n,2);
    for i = 1:n
        p = patients{i};
        uT0 = controlledU(i,1);
        uP8 = controlledU(i,2);

        if strcmpi(condition,'small')
            lowerT0 = bounds.T0.lo;
            upperT0 = min(settings.sizeThreshold,p.k);
        elseif strcmpi(condition,'large')
            lowerT0 = settings.sizeThreshold;
            upperT0 = min(bounds.T0.hi,p.k);
        else
            error('Unknown initial-size condition: %s.',condition)
        end
        if upperT0<=lowerT0
            error('Patient-specific k does not permit the requested T0 condition.')
        end
        T0 = sampleLinearOrLog(uT0,lowerT0,upperT0,'log');
        p8 = sampleFromBound(uP8,bounds.p8);
        patients{i} = setInitialComposition(p,T0,p8);
    end
    % Remove the unused original T0/p8 coordinates. The final two columns are
    % the only coordinates used to construct this cohort's T0 and p8.
    U = [unusedU(:,1:21) unusedU(:,24:26) controlledU];
end

function p = setInitialComposition(p,T0,p8)
    p.T0 = T0;
    p.p8 = p8;
    p.Ts0 = p8*T0;
    p.Tu0 = (1-p8)*T0;
    p.y0 = [0;p.Ts0;p.Tu0;p.Di0;1;1;1];
end

function validateCancerCohort(patients,bounds,settings)
    validateBackgroundCohort(patients,bounds)
    r = cellfun(@(p)p.r,patients);
    alpha = cellfun(@(p)p.alpha,patients);
    alphaMinimum = settings.theoreticalAlphaMinimum;
    if any(r<bounds.r.lo) || any(r>bounds.r.hi)
        error('Cancer-cohort r sampling left the agreed range.')
    end
    if any(alpha<alphaMinimum) || any(alpha>=settings.normalAlphaLower) || ...
            any(alpha<=0)
        error('Cancer-cohort alpha must remain strictly below 1/28 day^-1.')
    end
end

function validateIndependentSizeCohorts(groups,coordinates,bounds,settings)
    if numel(groups)~=2 || numel(groups{1})~=numel(groups{2})
        error('The size experiment must have two exactly equal groups.')
    end
    n = numel(groups{1});
    for g = 1:2
        validateCancerCohort(groups{g},bounds,settings)
        for i = 1:n
            validateInitialParameterSet(groups{g}{i},bounds)
        end
    end
    if any(cellfun(@(p)p.T0,groups{1})>=settings.sizeThreshold) || ...
            any(cellfun(@(p)p.T0,groups{2})<settings.sizeThreshold)
            error('A T0 value lies in the wrong 3-cm condition.')
    end
    combinedCoordinates = [coordinates{1};coordinates{2}];
    if size(unique(combinedCoordinates,'rows'),1)~=2*n
        error('The initial-size groups contain a reused LHS patient row.')
    end
    fprintf(['Initial-size experiment validated: %d distinct virtual ' ...
        'patients in each fixed baseline group.\n'],n)
end

function validateInitialParameterSet(p,bounds)
    tolerance = 1e-9;
    if p.T0<bounds.T0.lo || p.T0>=p.k || p.T0>bounds.T0.hi
        error('A controlled initial T0 value is inadmissible.')
    end
    if p.Ts0<0 || p.Tu0<0
        error('Continuous initial tumor-compartment burdens must be nonnegative.')
    end
    if p.p8<bounds.p8.lo-tolerance || p.p8>bounds.p8.hi+tolerance
        error('A controlled p8 value left the agreed global range.')
    end
    if abs(p.Ts0+p.Tu0-p.T0)>tolerance*max(1,p.T0)
        error('Controlled Ts(0)+Tu(0) does not equal T0.')
    end
end

function [groups,groupIndices,cutpoints] = splitBaselineP8Tertiles( ...
        patients,nPerGroup)

    n = numel(patients);
    if n~=3*nPerGroup
        error('The p8 cohort must contain exactly three equal groups.')
    end
    p8 = cellfun(@(p)p.p8,patients);
    [sortedP8,order] = sort(p8,'ascend');
    groupIndices = {order(1:nPerGroup), ...
        order(nPerGroup+1:2*nPerGroup),order(2*nPerGroup+1:end)};
    groups = cell(1,3);
    for g = 1:3
        groups{g} = patients(groupIndices{g});
    end
    cutpoints = [sortedP8(nPerGroup),sortedP8(2*nPerGroup)];
    allIndices = sort(vertcat(groupIndices{:}));
    if ~isequal(allIndices,(1:n).') || ...
            any(cellfun(@numel,groupIndices)~=nPerGroup)
        error('The baseline-p8 tertiles are not disjoint and exactly equal.')
    end
    fprintf(['Baseline-p8 tertiles validated: %d distinct patients per ' ...
        'fixed baseline group.\n'],nPerGroup)
end

function cohorts = runSizeCohorts(groups,settings,odeOptions)
    cohorts = cell(1,2);
    cohorts{1} = runPatientCohort(groups{1},settings,odeOptions, ...
        'initial size <3 cm');
    cohorts{2} = runPatientCohort(groups{2},settings,odeOptions, ...
        'initial size >=3 cm');
end

function cohort = runPatientCohort(patients,settings,odeOptions,labelText)
    n = numel(patients);
    nVisits = numel(settings.visitTimes);
    states = nan(n,nVisits,7);
    solverOK = false(n,1);
    solverMessage = strings(n,1);
    reportEvery = max(1,ceil(n/20));
    progressTimer = tic;

    fprintf('Running %s: %d complete 360-day simulations ...\n',labelText,n)
    if settings.useParallel
        queue = parallel.pool.DataQueue;
        completed = 0;
        afterEach(queue,@updateProgress);
        parfor i = 1:n
            [yVisits,ok,message] = simulateAtVisits( ...
                patients{i},settings,odeOptions);
            states(i,:,:) = reshape(yVisits.',1,nVisits,7);
            solverOK(i) = ok;
            solverMessage(i) = message;
            send(queue,1)
        end
    else
        for i = 1:n
            [yVisits,ok,message] = simulateAtVisits( ...
                patients{i},settings,odeOptions);
            states(i,:,:) = reshape(yVisits.',1,nVisits,7);
            solverOK(i) = ok;
            solverMessage(i) = message;
            if i==1 || mod(i,reportEvery)==0 || i==n
                reportProgressLine(labelText,i,n,toc(progressTimer))
            end
        end
    end
    requireCompletePatientCohort(solverOK,solverMessage,labelText)
    cohort.states = states;
    cohort.solverOK = solverOK;
    cohort.solverMessage = solverMessage;

    function updateProgress(~)
        completed = completed+1;
        if completed==1 || mod(completed,reportEvery)==0 || completed==n
            reportProgressLine(labelText,completed,n,toc(progressTimer))
        end
    end
end

function requireCompletePatientCohort(ok,messages,labelText)
    if all(ok)
        return
    end
    failed = find(~ok);
    error('%s is incomplete: %d/%d failed; first failure: %s', ...
        labelText,numel(failed),numel(ok),messages(failed(1)))
end

%% ========================================================================
% LOCAL FUNCTIONS: INITIAL-SIZE, p8 AND Ts-PERSISTENCE ANALYSES
% ========================================================================

function analysis = analyzeInitialSize(groups,U,cohorts,settings)
    labels = {'Initial tumor size < 3 cm','Initial tumor size >= 3 cm'};
    totals = cell(1,2);
    monthlyDetectableStatus = cell(1,2);
    detectableCount = zeros(settings.nMonths,2);
    proportionDetectable = zeros(settings.nMonths,2);
    for g = 1:2
        totals{g} = cohorts{g}.states(:,:,2)+cohorts{g}.states(:,:,3);
        
        % Exclude baseline and evaluate detectability separately at each monthly visit.

        monthlyDetectableStatus{g} = ...
            totals{g}(:,2:end)>=settings.detectionLimit;
        detectableCount(:,g) = sum(monthlyDetectableStatus{g},1).';
        proportionDetectable(:,g) = detectableCount(:,g)/numel(groups{g});
        groupSummary(g) = summarizeBurden(totals{g}); %#ok<AGROW>
    end

    Day = (1:settings.nMonths).'*settings.monthDays;
    N_Small = repmat(numel(groups{1}),settings.nMonths,1);
    N_Large = repmat(numel(groups{2}),settings.nMonths,1);
    Small_DetectableCount = detectableCount(:,1);
    Small_BelowDetectionCount = N_Small-Small_DetectableCount;
    Small_ProportionDetectable = proportionDetectable(:,1);
    Small_ProportionBelowDetection = 1-Small_ProportionDetectable;
    Large_DetectableCount = detectableCount(:,2);
    Large_BelowDetectionCount = N_Large-Large_DetectableCount;
    Large_ProportionDetectable = proportionDetectable(:,2);
    Large_ProportionBelowDetection = 1-Large_ProportionDetectable;

    resultsTable = table(Day,N_Small,Small_DetectableCount, ...
        Small_BelowDetectionCount,Small_ProportionDetectable, ...
        Small_ProportionBelowDetection,N_Large, ...
        Large_DetectableCount,Large_BelowDetectionCount, ...
        Large_ProportionDetectable,Large_ProportionBelowDetection);

    analysis.labels = labels;
    analysis.coordinates = U;
    names = backgroundFactorNames();
    analysis.coordinateNames = [names(1:21),names(24),{'r','alpha', ...
        'T0_interval_coordinate','p8_coordinate'}];
    analysis.totalTumorCells = totals;
    analysis.groupSummary = groupSummary;
    analysis.monthlyDetectableStatus = monthlyDetectableStatus;
    analysis.proportionDetectable = proportionDetectable;
    analysis.resultsTable = resultsTable;
    analysis.initialT0 = {cellfun(@(p)p.T0,groups{1}), ...
        cellfun(@(p)p.T0,groups{2})};
    Group = string(labels(:));
    N = [numel(groups{1});numel(groups{2})];
    MeanTotalTumorCellsDay360 = [groupSummary(1).mean(end); ...
        groupSummary(2).mean(end)];
    Q25TotalTumorCellsDay360 = [groupSummary(1).q25(end); ...
        groupSummary(2).q25(end)];
    MedianTotalTumorCellsDay360 = [groupSummary(1).median(end); ...
        groupSummary(2).median(end)];
    Q75TotalTumorCellsDay360 = [groupSummary(1).q75(end); ...
        groupSummary(2).q75(end)];
    analysis.day360BurdenTable = table(Group,N,MeanTotalTumorCellsDay360, ...
        Q25TotalTumorCellsDay360,MedianTotalTumorCellsDay360, ...
        Q75TotalTumorCellsDay360);
    validateInitialSizeAnalysis(analysis,groups,settings)
end

function analysis = analyzeInitialP8(patients,U,groups,groupIndices, ...
        cutpoints,cohort,settings)

    labels = {'Lower third of initial Ts proportion', ...
        'Middle third of initial Ts proportion', ...
        'Upper third of initial Ts proportion'};
    allTotal = cohort.states(:,:,2)+cohort.states(:,:,3);
    nGroups = 3;
    nMonths = settings.nMonths;
    groupSummary = struct([]);
    monthlyDetectableStatus = cell(1,nGroups);
    proportionDetectable = nan(nGroups,nMonths);
    rows = nGroups*nMonths;
    Day = zeros(rows,1);
    InitialTsProportionGroup = strings(rows,1);
    N = zeros(rows,1);
    P8Minimum = zeros(rows,1);
    P8Maximum = zeros(rows,1);
    DetectableCount = zeros(rows,1);
    BelowDetectionCount = zeros(rows,1);
    ProportionDetectable = zeros(rows,1);
    ProportionBelowDetection = zeros(rows,1);
    MeanTotalTumorCells = zeros(rows,1);
    Q25TotalTumorCells = zeros(rows,1);
    MedianTotalTumorCells = zeros(rows,1);
    Q75TotalTumorCells = zeros(rows,1);

    row = 0;
    for g = 1:nGroups
        indices = groupIndices{g};
        groupTotal = allTotal(indices,:);
        groupP8 = cellfun(@(p)p.p8,groups{g});
        summary = summarizeBurden(groupTotal);
        summary.label = labels{g};
        summary.n = numel(indices);
        summary.p8Minimum = min(groupP8);
        summary.p8Maximum = max(groupP8);
        if g==1
            groupSummary = repmat(summary,nGroups,1);
        else
            groupSummary(g) = summary;
        end
        monthlyDetectableStatus{g} = ...
            groupTotal(:,2:end)>=settings.detectionLimit;
        proportionDetectable(g,:) = ...
            mean(monthlyDetectableStatus{g},1);
        for month = 1:nMonths
            row = row+1;
            Day(row) = month*settings.monthDays;
            InitialTsProportionGroup(row) = labels{g};
            N(row) = numel(indices);
            P8Minimum(row) = groupSummary(g).p8Minimum;
            P8Maximum(row) = groupSummary(g).p8Maximum;
            DetectableCount(row) = ...
                sum(monthlyDetectableStatus{g}(:,month));
            BelowDetectionCount(row) = N(row)-DetectableCount(row);
            ProportionDetectable(row) = proportionDetectable(g,month);
            ProportionBelowDetection(row) = 1-ProportionDetectable(row);
            MeanTotalTumorCells(row) = groupSummary(g).mean(month+1);
            Q25TotalTumorCells(row) = groupSummary(g).q25(month+1);
            MedianTotalTumorCells(row) = groupSummary(g).median(month+1);
            Q75TotalTumorCells(row) = groupSummary(g).q75(month+1);
        end
    end

    analysis.labels = labels;
    analysis.groupLabels = labels;
    analysis.groupSizes = cellfun(@numel,groupIndices);
    analysis.groupIndices = groupIndices;
    analysis.cutpoints = cutpoints;
    analysis.coordinates = U;
    analysis.coordinateNames = [backgroundFactorNames(),{'r','alpha'}];
    analysis.totalTumorCells = allTotal;
    analysis.groupSummary = groupSummary;
    analysis.monthlyDetectableStatus = monthlyDetectableStatus;
    analysis.proportionDetectable = proportionDetectable;
    analysis.resultsTable = table(Day,InitialTsProportionGroup,N, ...
        P8Minimum,P8Maximum,DetectableCount,BelowDetectionCount, ...
        ProportionDetectable,ProportionBelowDetection, ...
        MeanTotalTumorCells,Q25TotalTumorCells, ...
        MedianTotalTumorCells,Q75TotalTumorCells);
    validateInitialP8Analysis(analysis,patients,settings)
end

function summary = summarizeBurden(total)
    nVisits = size(total,2);
    summary.mean = mean(total,1);
    summary.q25 = zeros(1,nVisits);
    summary.median = zeros(1,nVisits);
    summary.q75 = zeros(1,nVisits);
    for v = 1:nVisits
        q = empiricalQuantile(total(:,v),[0.25 0.5 0.75]);
        summary.q25(v) = q(1);
        summary.median(v) = q(2);
        summary.q75(v) = q(3);
    end
end

function validateInitialSizeAnalysis(analysis,groups,settings)
    if numel(groups{1})~=numel(groups{2}) || ...
            any(analysis.initialT0{1}>=settings.sizeThreshold) || ...
            any(analysis.initialT0{2}<settings.sizeThreshold)
        error('The reported 3-cm experiment is not exactly balanced.')
    end
    if height(analysis.resultsTable)~=12 || ...
            ~isequal(analysis.resultsTable.Day,(30:30:360).')
        error('The size analysis must report all 12 scheduled monthly visits.')
    end
    for g = 1:2
        expectedStatus = ...
            analysis.totalTumorCells{g}(:,2:end)>=settings.detectionLimit;
        if ~isequal(analysis.monthlyDetectableStatus{g},expectedStatus)
            error('Initial-size monthly detectability was calculated incorrectly.')
        end
    end
    if any(analysis.resultsTable.Small_DetectableCount+ ...
            analysis.resultsTable.Small_BelowDetectionCount~= ...
            analysis.resultsTable.N_Small) || ...
       any(analysis.resultsTable.Large_DetectableCount+ ...
            analysis.resultsTable.Large_BelowDetectionCount~= ...
            analysis.resultsTable.N_Large)
        error('Every initial-size patient must be counted once at every month.')
    end
    expectedSmallProportion = ...
        analysis.resultsTable.Small_DetectableCount./ ...
        analysis.resultsTable.N_Small;
    expectedLargeProportion = ...
        analysis.resultsTable.Large_DetectableCount./ ...
        analysis.resultsTable.N_Large;
    if any(abs(analysis.resultsTable.Small_ProportionDetectable- ...
            expectedSmallProportion)>1e-14) || ...
       any(abs(analysis.resultsTable.Large_ProportionDetectable- ...
            expectedLargeProportion)>1e-14)
        error('Initial-size monthly detectable proportions are inconsistent.')
    end
    if any(abs(analysis.resultsTable.Small_ProportionBelowDetection- ...
            (1-expectedSmallProportion))>1e-14) || ...
       any(abs(analysis.resultsTable.Large_ProportionBelowDetection- ...
            (1-expectedLargeProportion))>1e-14)
        error('Initial-size monthly below-detection proportions are inconsistent.')
    end
end

function validateInitialP8Analysis(analysis,patients,settings)
    n = numel(patients);
    allIndices = sort(vertcat(analysis.groupIndices{:}));
    if ~isequal(allIndices,(1:n).') || ...
            any(analysis.groupSizes~=settings.nP8PerGroup)
        error('Baseline-p8 groups are not distinct and exactly equal.')
    end
    p8 = cellfun(@(p)p.p8,patients);
    groupMax = cellfun(@(idx)max(p8(idx)),analysis.groupIndices);
    groupMin = cellfun(@(idx)min(p8(idx)),analysis.groupIndices);
    if groupMax(1)>groupMin(2) || groupMax(2)>groupMin(3)
        error('Baseline-p8 tertiles overlap or were not sorted at baseline.')
    end
    if height(analysis.resultsTable)~=3*settings.nMonths
        error('The p8 analysis must report all 12 visits for three groups.')
    end
    for g = 1:3
        expectedStatus = analysis.totalTumorCells( ...
            analysis.groupIndices{g},2:end)>=settings.detectionLimit;
        if ~isequal(analysis.monthlyDetectableStatus{g},expectedStatus)
            error('Initial-Ts-proportion monthly detectability is incorrect.')
        end
    end
    if any(analysis.resultsTable.DetectableCount+ ...
            analysis.resultsTable.BelowDetectionCount~= ...
            analysis.resultsTable.N)
        error(['Every initial-Ts-proportion patient must be counted once ' ...
            'at every month.'])
    end
    expectedProportion = analysis.resultsTable.DetectableCount./ ...
        analysis.resultsTable.N;
    if any(abs(analysis.resultsTable.ProportionDetectable- ...
            expectedProportion)>1e-14)
        error(['Initial-Ts-proportion monthly detectable proportions are ' ...
            'inconsistent.'])
    end
    if any(abs(analysis.resultsTable.ProportionBelowDetection- ...
            (1-expectedProportion))>1e-14)
        error(['Initial-Ts-proportion monthly below-detection proportions ' ...
            'are inconsistent.'])
    end
end

function analysis = analyzePersistence(~,U,cohort,settings)
    total = cohort.states(:,:,2)+cohort.states(:,:,3);
    Ts = cohort.states(:,:,2);
    finalDetectable = total(:,end)>=settings.detectionLimit;
    groupMasks = {finalDetectable,~finalDetectable};
    labels = {'Detectable tumor at day 360', ...
        'Below 1-mm detection at day 360'};
    nGroups = numel(groupMasks);
    nVisits = numel(settings.visitTimes);

    rowCount = nGroups*nVisits;
    Day = zeros(rowCount,1);
    Day360Status = strings(rowCount,1);
    N = zeros(rowCount,1);
    MedianTotalTumorCells = nan(rowCount,1);
    Q25TotalTumorCells = nan(rowCount,1);
    Q75TotalTumorCells = nan(rowCount,1);
    MedianTsCells = nan(rowCount,1);
    Q25TsCells = nan(rowCount,1);
    Q75TsCells = nan(rowCount,1);

    groupSummary = repmat(struct(),nGroups,1);
    row = 0;
    for g = 1:nGroups
        mask = groupMasks{g};
        groupTotal = total(mask,:);
        groupTs = Ts(mask,:);
        groupSummary(g).label = labels{g};
        groupSummary(g).n = sum(mask);
        groupSummary(g).medianTotal = nan(1,nVisits);
        groupSummary(g).q25Total = nan(1,nVisits);
        groupSummary(g).q75Total = nan(1,nVisits);
        groupSummary(g).medianTs = nan(1,nVisits);
        groupSummary(g).q25Ts = nan(1,nVisits);
        groupSummary(g).q75Ts = nan(1,nVisits);
        for v = 1:nVisits
            row = row+1;
            Day(row) = settings.visitTimes(v);
            Day360Status(row) = labels{g};
            N(row) = sum(mask);
            if any(mask)
                totalValues = groupTotal(:,v);
                tsValues = groupTs(:,v);
                totalQuartiles = empiricalQuantile(totalValues,[0.25 0.5 0.75]);
                tsQuartiles = empiricalQuantile(tsValues,[0.25 0.5 0.75]);
                Q25TotalTumorCells(row) = totalQuartiles(1);
                MedianTotalTumorCells(row) = totalQuartiles(2);
                Q75TotalTumorCells(row) = totalQuartiles(3);
                Q25TsCells(row) = tsQuartiles(1);
                MedianTsCells(row) = tsQuartiles(2);
                Q75TsCells(row) = tsQuartiles(3);
                groupSummary(g).q25Total(v) = totalQuartiles(1);
                groupSummary(g).medianTotal(v) = totalQuartiles(2);
                groupSummary(g).q75Total(v) = totalQuartiles(3);
                groupSummary(g).q25Ts(v) = tsQuartiles(1);
                groupSummary(g).medianTs(v) = tsQuartiles(2);
                groupSummary(g).q75Ts(v) = tsQuartiles(3);
            end
        end
    end

    resultsTable = table(Day,Day360Status,N,MedianTotalTumorCells, ...
        Q25TotalTumorCells,Q75TotalTumorCells,MedianTsCells,Q25TsCells, ...
        Q75TsCells);

    analysis.visitTimes = settings.visitTimes;
    analysis.coordinates = U;
    analysis.coordinateNames = [backgroundFactorNames(),{'r','alpha'}];
    analysis.labels = labels;
    analysis.totalTumorCells = total;
    analysis.TsCells = Ts;
    analysis.finalDetectable = finalDetectable;
    analysis.groupSummary = groupSummary;
    analysis.resultsTable = resultsTable;
    validatePersistenceAnalysis(analysis,cohort,settings)
end

function q = empiricalQuantile(values,probabilities)
    values = sort(values(:));
    n = numel(values);
    q = nan(size(probabilities));
    if n==0
        return
    elseif n==1
        q = repmat(values,size(probabilities));
        return
    end
    positions = 1+(n-1)*probabilities;
    lower = floor(positions);
    upper = ceil(positions);
    weight = positions-lower;
    for j = 1:numel(probabilities)
        q(j) = values(lower(j))*(1-weight(j))+values(upper(j))*weight(j);
    end
end

function validatePersistenceAnalysis(analysis,cohort,settings)
    expectedTotal = cohort.states(:,:,2)+cohort.states(:,:,3);
    tolerance = 1e-10*max(1,max(expectedTotal(:)));
    if max(abs(analysis.totalTumorCells(:)-expectedTotal(:)))>tolerance
        error('Temporal total tumor burden must equal Ts+Tu.')
    end
    if size(analysis.totalTumorCells,2)~=13 || ...
            ~isequal(analysis.visitTimes,0:30:360)
        error('Temporal persistence must cover baseline plus 12 monthly visits.')
    end
    if settings.detectionLimit<=1
        error('The 1-mm detection threshold must remain above one cell-equivalent.')
    end
end

%% ========================================================================
% LOCAL FUNCTIONS: r-ALPHA GRID EXPERIMENT
% ========================================================================

function result = runRAlphaGrid(backgrounds,U,rGrid,alphaGrid, ...
        alphaRegion,settings,odeOptions)

    nBackgrounds = numel(backgrounds);
    nAlpha = numel(alphaGrid);
    nR = numel(rGrid);
    nTasks = nBackgrounds*nAlpha*nR;

    totalDay360Vector = nan(nTasks,1);
    tsDay360Vector = nan(nTasks,1);
    solverOK = false(nTasks,1);
    solverMessage = strings(nTasks,1);
    reportEvery = max(1,ceil(nTasks/100));
    progressTimer = tic;

    fprintf(['Running %d simulations on the common-background r-alpha ' ...
        'grid for the prespecified day-360 endpoint ...\n'],nTasks)

    if settings.useParallel
        queue = parallel.pool.DataQueue;
        completed = 0;
        afterEach(queue,@updateProgress);
        parfor task = 1:nTasks
            [backgroundIndex,alphaIndex,rIndex] = ind2sub( ...
                [nBackgrounds,nAlpha,nR],task);
            p = backgrounds{backgroundIndex};
            p.r = rGrid(rIndex);
            p.alpha = alphaGrid(alphaIndex);
            [yVisits,ok,message] = simulateAtVisits( ...
                p,settings,odeOptions);
            if ok
                tsDay360Vector(task) = yVisits(2,end);
                totalDay360Vector(task) = yVisits(2,end)+yVisits(3,end);
            end
            solverOK(task) = ok;
            solverMessage(task) = message;
            send(queue,1)
        end
    else
        for task = 1:nTasks
            [backgroundIndex,alphaIndex,rIndex] = ind2sub( ...
                [nBackgrounds,nAlpha,nR],task);
            p = backgrounds{backgroundIndex};
            p.r = rGrid(rIndex);
            p.alpha = alphaGrid(alphaIndex);
            [yVisits,ok,message] = simulateAtVisits( ...
                p,settings,odeOptions);
            if ok
                tsDay360Vector(task) = yVisits(2,end);
                totalDay360Vector(task) = yVisits(2,end)+yVisits(3,end);
            end
            solverOK(task) = ok;
            solverMessage(task) = message;
            if task==1 || mod(task,reportEvery)==0 || task==nTasks
                reportProgressLine('r-alpha grid',task,nTasks,toc(progressTimer))
            end
        end
    end

    requireCompleteCohort(solverOK,solverMessage)

    % Summary arrays are alpha x r. Every point uses the identical rows of
    % background inputs, so only r and alpha differ across the surface.
    totalDay360Cube = reshape(totalDay360Vector, ...
        [nBackgrounds,nAlpha,nR]);
    tsDay360Cube = reshape(tsDay360Vector, ...
        [nBackgrounds,nAlpha,nR]);
    meanDay360 = squeeze(mean(totalDay360Cube,1));
    meanTsDay360 = squeeze(mean(tsDay360Cube,1));
    meanTuDay360 = squeeze(mean(totalDay360Cube-tsDay360Cube,1));
    medianDay360 = squeeze(median(totalDay360Cube,1));
    medianTsDay360 = squeeze(median(tsDay360Cube,1));
    medianTuDay360 = squeeze(median(totalDay360Cube-tsDay360Cube,1));

    % Normal-reference response domain. The alpha=1/28 day^-1 grid row is
    % the slowest normal-tissue-referenced characteristic transition rate.
    % Because the identical parameter backgrounds are reused over the entire
    % grid, subtraction at fixed r is paired by construction. Positive values
    % mean that increasing alpha above the 28-day reference reduces mean
    % day-360 total-tumor burden. This is a descriptive model-response domain,
    % not a clinical threshold or a statistical-significance classification.
    normalReferenceIndex = find(abs(alphaGrid-settings.normalAlphaLower)<1e-14);
    if numel(normalReferenceIndex)~=1
        error('The alpha grid must contain exactly one 28-day reference row.')
    end
    referenceMeanTotalTumorCells = meanDay360(normalReferenceIndex,:);
    referenceMeanMatrix = repmat(referenceMeanTotalTumorCells,nAlpha,1);
    absoluteTotalTumorReduction = referenceMeanMatrix-meanDay360;
    relativeTotalTumorReduction = nan(nAlpha,nR);
    normalComparisonRows = alphaGrid(:)>=settings.normalAlphaLower & ...
        alphaGrid(:)<=settings.normalAlphaUpper;
    validReferenceColumns = referenceMeanTotalTumorCells>0;
    validComparison = repmat(normalComparisonRows,1,nR) & ...
        repmat(validReferenceColumns,nAlpha,1);
    relativeTotalTumorReduction(validComparison) = ...
        absoluteTotalTumorReduction(validComparison)./ ...
        referenceMeanMatrix(validComparison);
    normalReferenceReductionDomain = ...
        repmat(alphaGrid(:)>settings.normalAlphaLower & ...
        alphaGrid(:)<=settings.normalAlphaUpper,1,nR) & ...
        absoluteTotalTumorReduction>0;

    [alphaMatrix,rMatrix] = ndgrid(alphaGrid,rGrid);
    regionMatrix = repmat(alphaRegion(:),1,nR);
    transitionTimeDays = 1./alphaMatrix;

    responseTable = table( ...
        rMatrix(:),alphaMatrix(:),transitionTimeDays(:), ...
        regionMatrix(:),repmat(nBackgrounds,numel(rMatrix),1), ...
        meanDay360(:),medianDay360(:),meanTsDay360(:), ...
        medianTsDay360(:),meanTuDay360(:),medianTuDay360(:), ...
        referenceMeanMatrix(:),absoluteTotalTumorReduction(:), ...
        relativeTotalTumorReduction(:), ...
        normalReferenceReductionDomain(:), ...
        'VariableNames',{'r_day_inverse','alpha_day_inverse', ...
        'DifferentiationLikeTransitionTime_days','AlphaRegion', ...
        'N_common_backgrounds', ...
        'MeanTotalTumorCellsDay360','MedianTotalTumorCellsDay360', ...
        'MeanTsCellsDay360','MedianTsCellsDay360', ...
        'MeanTuCellsDay360','MedianTuCellsDay360', ...
        'ReferenceMeanTotalTumorCellsAt28DayRate', ...
        'AbsoluteTotalTumorReductionFrom28DayRate', ...
        'RelativeTotalTumorReductionFrom28DayRate', ...
        'NormalReferenceReductionDomain'});

    result.rGrid = rGrid;
    result.alphaGrid = alphaGrid;
    result.alphaRegion = alphaRegion;
    result.endpointDay = 360;
    result.nBackgrounds = nBackgrounds;
    result.meanTotalTumorCells = meanDay360;
    result.medianTotalTumorCells = medianDay360;
    result.meanTsCells = meanTsDay360;
    result.medianTsCells = medianTsDay360;
    result.meanTuCells = meanTuDay360;
    result.medianTuCells = medianTuDay360;
    result.normalReferenceIndex = normalReferenceIndex;
    result.referenceMeanTotalTumorCells = referenceMeanTotalTumorCells;
    result.absoluteTotalTumorReduction = absoluteTotalTumorReduction;
    result.relativeTotalTumorReduction = relativeTotalTumorReduction;
    result.normalReferenceReductionDomain = ...
        normalReferenceReductionDomain;
    result.totalDay360Cube = totalDay360Cube;
    result.tsDay360Cube = tsDay360Cube;
    result.responseTable = responseTable;
    result.backgroundFactorNames = backgroundFactorNames();
    result.lhsCoordinates = U;
    result.convergenceTable = makeRAlphaResolutionDiagnostics( ...
        result,rGrid,alphaGrid);

    fprintf(['r-alpha simulations completed. Mean day-360 total burden ' ...
        'across the grid: %.6e to %.6e cells.\n'], ...
        min(meanDay360(:)),max(meanDay360(:)))
    nEligibleNormalPoints = sum(alphaGrid>settings.normalAlphaLower & ...
        alphaGrid<=settings.normalAlphaUpper)*nR;
    fprintf(['Normal-reference reduction domain: %d/%d faster-normal ' ...
        'grid points have lower mean total burden than alpha=1/28 at ' ...
        'the same r.\n'],sum(normalReferenceReductionDomain(:)), ...
        nEligibleNormalPoints)

    function updateProgress(~)
        completed = completed+1;
        if completed==1 || mod(completed,reportEvery)==0 || ...
                completed==nTasks
            reportProgressLine('r-alpha grid',completed,nTasks, ...
                toc(progressTimer))
        end
    end
end

function subset = subsetRAlphaResult(fullResult,rMask,alphaMask,settings)
    % Extract the exact nested 0.5-domain result without rerunning any ODE.
    % The masks select coordinates from the already-completed extended grid.
    if numel(rMask)~=numel(fullResult.rGrid) || ...
            numel(alphaMask)~=numel(fullResult.alphaGrid)
        error('The nested r-alpha subset masks have incompatible sizes.')
    end
    if ~any(rMask) || ~any(alphaMask)
        error('The nested r-alpha subset cannot be empty.')
    end

    subset = fullResult;
    subset.rGrid = fullResult.rGrid(rMask);
    subset.alphaGrid = fullResult.alphaGrid(alphaMask);
    subset.alphaRegion = fullResult.alphaRegion(alphaMask);
    subset.meanTotalTumorCells = ...
        fullResult.meanTotalTumorCells(alphaMask,rMask);
    subset.medianTotalTumorCells = ...
        fullResult.medianTotalTumorCells(alphaMask,rMask);
    subset.meanTsCells = fullResult.meanTsCells(alphaMask,rMask);
    subset.medianTsCells = fullResult.medianTsCells(alphaMask,rMask);
    subset.meanTuCells = fullResult.meanTuCells(alphaMask,rMask);
    subset.medianTuCells = fullResult.medianTuCells(alphaMask,rMask);
    subset.referenceMeanTotalTumorCells = ...
        fullResult.referenceMeanTotalTumorCells(rMask);
    subset.absoluteTotalTumorReduction = ...
        fullResult.absoluteTotalTumorReduction(alphaMask,rMask);
    subset.relativeTotalTumorReduction = ...
        fullResult.relativeTotalTumorReduction(alphaMask,rMask);
    subset.normalReferenceReductionDomain = ...
        fullResult.normalReferenceReductionDomain(alphaMask,rMask);
    subset.totalDay360Cube = fullResult.totalDay360Cube(:,alphaMask,rMask);
    subset.tsDay360Cube = fullResult.tsDay360Cube(:,alphaMask,rMask);
    subset.normalReferenceIndex = find(abs( ...
        subset.alphaGrid-settings.normalAlphaLower)<1e-14);

    selectionMatrix = repmat(alphaMask(:),1,numel(fullResult.rGrid)) & ...
        repmat(rMask(:).',numel(fullResult.alphaGrid),1);
    subset.responseTable = fullResult.responseTable(selectionMatrix(:),:);
    subset.convergenceTable = makeRAlphaResolutionDiagnostics( ...
        subset,subset.rGrid,subset.alphaGrid);
end

function diagnostics = makeRAlphaResolutionDiagnostics(result,rGrid,alphaGrid)

    nBackgrounds = result.nBackgrounds;
    halfCount = max(1,floor(nBackgrounds/2));
    totalHalf = result.totalDay360Cube(1:halfCount,:,:);
    tsHalf = result.tsDay360Cube(1:halfCount,:,:);
    halfSurfaces = {squeeze(mean(totalHalf,1)), ...
        squeeze(median(totalHalf,1)),squeeze(mean(tsHalf,1)), ...
        squeeze(median(tsHalf,1))};
    fullSurfaces = {result.meanTotalTumorCells, ...
        result.medianTotalTumorCells,result.meanTsCells, ...
        result.medianTsCells};
    names = ["Mean total tumor cells";"Median total tumor cells"; ...
        "Mean Ts cells";"Median Ts cells"];

    coarseR = unique([1:2:numel(rGrid),numel(rGrid)]);
    coarseAlpha = unique([1:2:numel(alphaGrid),numel(alphaGrid)]);
    [queryR,queryLogAlpha] = meshgrid(rGrid,log10(alphaGrid));
    backgroundNRMSE = zeros(4,1);
    gridNRMSE = zeros(4,1);
    for j = 1:4
        fullSurface = fullSurfaces{j};
        backgroundNRMSE(j) = normalizedRMSE(halfSurfaces{j},fullSurface);
        coarseSurface = fullSurface(coarseAlpha,coarseR);
        interpolated = interp2(rGrid(coarseR),log10(alphaGrid(coarseAlpha)), ...
            coarseSurface,queryR,queryLogAlpha,'linear');
        gridNRMSE(j) = normalizedRMSE(interpolated,fullSurface);
    end

    Diagnostic = [repmat("First half of backgrounds versus full set",4,1); ...
        repmat("Coarse-grid interpolation versus full grid",4,1)];
    ResponseSummary = [names;names];
    NBackgroundsInReducedCalculation = [repmat(halfCount,4,1); ...
        repmat(nBackgrounds,4,1)];
    NormalizedRMSE = [backgroundNRMSE;gridNRMSE];
    diagnostics = table(Diagnostic,ResponseSummary, ...
        NBackgroundsInReducedCalculation,NormalizedRMSE);
end

function value = normalizedRMSE(candidate,reference)
    difference = candidate(:)-reference(:);
    scale = max(reference(:))-min(reference(:));
    if scale<=0
        scale = max(1,max(abs(reference(:))));
    end
    value = sqrt(mean(difference.^2))/scale;
end

function requireCompleteCohort(ok,messages)

    if all(ok)
        return
    end
    failed = find(~ok);
    firstFailure = failed(1);
    error(['The r-alpha grid is incomplete: %d/%d simulations failed. ' ...
        'First failure at linear task %d: %s'], ...
        numel(failed),numel(ok),firstFailure,messages(firstFailure))
end

%% ========================================================================
% LOCAL FUNCTIONS: PIECEWISE CONTINUOUS MMC INTEGRATION
% ========================================================================

function [yVisits,ok,message] = simulateAtVisits(p,settings,odeOptions)

    breakPoints = unique([0,settings.tEnd,settings.doseDays, ...
        settings.doseDays+settings.tau,settings.visitTimes]);
    breakPoints = sort(breakPoints( ...
        breakPoints>=0 & breakPoints<=settings.tEnd));

    nVisits = numel(settings.visitTimes);
    yVisits = nan(7,nVisits);
    yCurrent = p.y0(:);
    yVisits(:,settings.visitTimes==0) = yCurrent;
    ok = false;
    message = "continuous integration did not finish";

    for segment = 2:numel(breakPoints)
        t0 = breakPoints(segment-1);
        t1 = breakPoints(segment);
        if t1<=t0; continue; end

        mCurrent = activeInstillationRate( ...
            t0,settings.doseDays,settings.tau,p.m);
        try
            [tPiece,yPiece] = ode15s( ...
                @(time,state)modelRHS(time,state,p,mCurrent), ...
                [t0 t1],yCurrent,odeOptions);
        catch exception
            yVisits(:) = NaN;
            message = "ode15s error: "+string(exception.message);
            return
        end

        reachedEnd = ~isempty(tPiece) && ...
            tPiece(end)>=t1-1e-10*max(1,abs(t1));
        if isempty(yPiece) || ~reachedEnd || any(~isfinite(yPiece(:))) || ...
                any(yPiece(:)<-1e-8)
            yVisits(:) = NaN;
            message = "ode15s returned an incomplete or invalid solution";
            return
        end

        % Only solver-scale negative roundoff is corrected. This is not a
        % biological floor, one-cell rule or reporting reset.
        yCurrent = yPiece(end,:).';
        yCurrent(yCurrent<0 & yCurrent>-1e-8) = 0;

        visitIndex = find(abs(settings.visitTimes-t1)<= ...
            1e-10*max(1,abs(t1)),1,'first');
        if ~isempty(visitIndex); yVisits(:,visitIndex) = yCurrent; end
    end

    if any(~isfinite(yVisits(:))) || any(yVisits(:)<0)
        yVisits(:) = NaN;
        message = "nonfinite, negative or missing monthly state";
        return
    end
    ok = true;
    message = "continuous ode15s";
end

function mCurrent = activeInstillationRate(t,pulseStarts,tau,mDose)
    if tau>0 && any(t>=pulseStarts & t<pulseStarts+tau)
        mCurrent = mDose;
    else
        mCurrent = 0;
    end
end

function dydt = modelRHS(~,y,p,mCurrent)

    M = y(1);
    Ts = y(2);
    Tu = y(3);
    Di = y(4);
    Dm = y(5);
    E = y(6);
    R = y(7);

    T = Ts+Tu;
    F = M/(M+p.a);
    immuneAttenuation = exp(-((T/p.k)*(R/p.b)*(1-F)));

    signal1 = p.beta1*(Tu+p.theta1*Ts);
    signal2 = p.beta2*(Tu+p.theta1*Ts);
    maturation = p.p3*signal1/(signal1+p.h) ...
        + p.p4*F*signal2/(signal2+p.h);

    dydt = zeros(7,1);
    dydt(1) = -p.mu1*M+mCurrent;
    dydt(2) = (p.r*Ts*(1-F)-p.alpha*Ts)*(1-T/p.k) ...
        - Ts*((1-p.theta2)*p.p5*E*immuneAttenuation+p.p1*F+p.mus);
    dydt(3) = (p.r*Tu*(1-F)+p.alpha*Ts)*(1-T/p.k) ...
        - Tu*(p.p5*E*immuneAttenuation+p.p2*F+p.muu);
    dydt(4) = p.d0-p.mu2*Di-Di*maturation;
    dydt(5) = Di*maturation-p.mu2*Dm;
    dydt(6) = p.gamma*Dm-p.p6*R*E-p.mu3*E;
    dydt(7) = p.eta*Di*((p.d0/p.mu2-Di)/(p.d0/p.mu2))*(1-F) ...
        - p.mu4*R-p.p7*R*F;

end

%% ========================================================================
% LOCAL FUNCTIONS: FOUR SELECTED PUBLICATION FIGURES
% ========================================================================

function figures = makeSelectedFigures(sizeAnalysis,p8Analysis, ...
        persistenceAnalysis,rAlphaPrimary,rAlphaExtended,settings)
    figures = gobjects(4,1);
    figures(1) = plotMonthlyDetectabilityFigure( ...
        sizeAnalysis,p8Analysis,settings);
    figures(2) = plotPersistenceFigure(persistenceAnalysis,settings);
    figures(3) = plotRAlphaCombinedFigure(rAlphaPrimary,settings);
    figures(4) = plotRAlphaCombinedFigure(rAlphaExtended,settings);
end

function fig = plotMonthlyDetectabilityFigure( ...
        sizeAnalysis,p8Analysis,settings)
    days = (1:settings.nMonths)*settings.monthDays;
    sizeColors = [0.0000 0.4470 0.7410;0.8500 0.3250 0.0980];
    p8Colors = [0.00 0.52 0.32;0.49 0.18 0.56;0.93 0.55 0.02];

    belowSize = 1-sizeAnalysis.proportionDetectable;
    belowP8 = 1-p8Analysis.proportionDetectable;

    fig = publicationFigure('Monthly tumor burden below detection', ...
        [0.4 0.4 14.8 7.2],settings);
    layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','loose');
    title(layout,'Tumor burden below detection at monthly assessments', ...
        'FontName','Arial','FontWeight','bold','FontSize',17)

    axA = nexttile(layout); hold(axA,'on')
    for g = 1:2
        plot(axA,days,belowSize(:,g),'-o', ...
            'LineWidth',2.7,'MarkerSize',5.5,'Color',sizeColors(g,:), ...
            'MarkerFaceColor',sizeColors(g,:));
    end
    title(axA,'Initial tumor size')
    xlabel(axA,'Follow-up time [days]')
    ylabel(axA,'Proportion below the 1-mm detection threshold')
    legend(axA,{sprintf('Initial size < 3 cm (n=%d)', ...
        sizeAnalysis.resultsTable.N_Small(1)), ...
        sprintf('Initial size >= 3 cm (n=%d)', ...
        sizeAnalysis.resultsTable.N_Large(1))}, ...
        'Location','southeast','Box','off','Interpreter','none', ...
        'FontSize',12.5,'FontWeight','bold')

    axB = nexttile(layout); hold(axB,'on')
    p8Legend = cell(3,1);
    displayThirds = {'Lower third','Middle third','Upper third'};
    for g = 1:3
        plot(axB,days,belowP8(g,:),'-o', ...
            'LineWidth',2.5,'MarkerSize',5.2,'Color',p8Colors(g,:), ...
            'MarkerFaceColor',p8Colors(g,:));
        p8Legend{g} = sprintf('%s of initial T_s proportion (n=%d)', ...
            displayThirds{g},p8Analysis.groupSizes(g));
    end
    title(axB,'Initial tumor stem-cell proportion','Interpreter','tex')
    xlabel(axB,'Follow-up time [days]')
    legend(axB,p8Legend,'Location','southeast','Box','off', ...
        'Interpreter','tex','FontSize',12.5,'FontWeight','bold')

    yMaximum = 1;
    applyMonthlyDetectabilityAxes(axA,settings,yMaximum)
    applyMonthlyDetectabilityAxes(axB,settings,yMaximum)
    linkaxes([axA axB],'xy')
    drawnow
    addOutsidePanelLetter(fig,axA,'A')
    addOutsidePanelLetter(fig,axB,'B')
end

function applyMonthlyDetectabilityAxes(ax,settings,yMaximum)
    xlim(ax,[settings.monthDays settings.tEnd])
    xticks(ax,[30 60 120 180 240 300 360])
    ylim(ax,[0 yMaximum])
    yticks(ax,0:0.1:yMaximum)
    stylePublicationAxes(ax,'two-dimensional')
end

function fig = plotPersistenceFigure(analysis,settings)
    fig = publicationFigure('Continuous temporal tumor and Ts burdens', ...
        [0.25 0.25 15.4 7.3],settings);
    layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','loose');
    title(layout,'Continuous tumor burden according to day-360 detectability', ...
        'FontName','Arial','FontWeight','bold','FontSize',17, ...
        'Interpreter','tex')

    colors = [0.8500 0.3250 0.0980;0.0000 0.4470 0.7410];
    labels = cell(1,2);
    for g = 1:2
        labels{g} = sprintf('%s (n=%d)', ...
            conciseOutcomeLabel(analysis.groupSummary(g).label), ...
            analysis.groupSummary(g).n);
    end

    axA = nexttile(layout); hold(axA,'on')
    axB = nexttile(layout); hold(axB,'on')
    handles = gobjects(2,1);
    for g = 1:2
        s = analysis.groupSummary(g);
        handles(g) = plotMedianIQR(axA,analysis.visitTimes, ...
            s.medianTotal,s.q25Total,s.q75Total,colors(g,:));
        plotMedianIQR(axB,analysis.visitTimes, ...
            s.medianTs,s.q25Ts,s.q75Ts,colors(g,:));
    end
    yline(axA,log10(settings.detectionLimit+1),'k:', ...
        'LineWidth',1.5,'Label','1-mm detection threshold', ...
        'LabelHorizontalAlignment','left','HandleVisibility','off');
    xlabel(axA,'Follow-up time [days]')
    ylabel(axA,'log_{10}(T_s+T_u+1) [cells]','Interpreter','tex')
    legend(axA,handles,labels,'Location','northwest','Box','off', ...
        'Interpreter','none','FontSize',12.5,'FontWeight','bold')
    applyTemporalAxes(axA,settings)

    xlabel(axB,'Follow-up time [days]')
    ylabel(axB,'log_{10}(T_s+1) [cells]','Interpreter','tex')
    applyTemporalAxes(axB,settings)

    drawnow
    addOutsidePanelLetter(fig,axA,'A')
    addOutsidePanelLetter(fig,axB,'B')
end

function label = conciseOutcomeLabel(fullLabel)
    if contains(string(fullLabel),"Detectable tumor")
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
        'LineWidth',2.4,'MarkerSize',5,'MarkerFaceColor',color);
end

function applyTemporalAxes(ax,settings)
    xlim(ax,[0 settings.tEnd])
    xticks(ax,0:60:settings.tEnd)
    stylePublicationAxes(ax,'two-dimensional')
end

function fig = plotRAlphaCombinedFigure(result,settings)
    meanTotal = result.meanTotalTumorCells;
    meanTs = result.meanTsCells;
    maximumTotal = max(meanTotal(:));
    maximumTs = max(meanTs(:));
    if maximumTotal<=0; maximumTotal=1; end
    if maximumTs<=0; maximumTs=1; end

    domainUpper = result.rGrid(end);
    fig = publicationFigure(sprintf( ...
        'Three-dimensional r-alpha response through %.4g',domainUpper), ...
        [0.1 0.1 19.6 9.5],settings);
    layout = tiledlayout(fig,1,2,'TileSpacing','loose','Padding','loose');
    title(layout,sprintf( ...
        'Day-360 tumor burden (r and \\alpha through %.4g day^{-1})', ...
        domainUpper), ...
        'FontName','Arial','FontWeight','bold','FontSize',17, ...
        'Interpreter','tex')

    axA = nexttile(layout); hold(axA,'on')
    plotRAlphaSurfacePanel(axA,result,meanTotal,settings, ...
        'Mean total tumor-cell count at day 360 [cells]', ...
        'Mean total-tumor count',maximumTotal,true);
    cbTotal = colorbar(axA,'Location','eastoutside');
    cbTotal.Label.String = 'Mean total tumor-cell count [cells]';
    stylePublicationColorbar(cbTotal)

    axB = nexttile(layout); hold(axB,'on')
    plotRAlphaSurfacePanel(axB,result,meanTs,settings, ...
        'Mean T_s-cell count at day 360 [cells]', ...
        'Mean T_s-compartment count',maximumTs,false);
    cbTs = colorbar(axB,'Location','eastoutside');
    cbTs.Label.String = 'Mean T_s-cell count [cells]';
    stylePublicationColorbar(cbTs)

    view(axA,[42 27]); view(axB,[42 27]);
    linkedAxes = linkprop([axA axB],{'View','XLim','YLim'});
    setappdata(fig,'LinkedRAlphaAxes',linkedAxes)
    drawnow
    addOutsidePanelLetter(fig,axA,'A')
    addOutsidePanelLetter(fig,axB,'B')
end

function plotRAlphaSurfacePanel(ax,result,zCells,settings,zLabelText, ...
        titleText,maximumZ,showReductionDomain)
    [rMatrix,alphaMatrix] = meshgrid(result.rGrid,result.alphaGrid);
    surfaceHandle = surf(ax,rMatrix,alphaMatrix,zCells,zCells, ...
        'EdgeColor','none','FaceColor','interp');
    surfaceHandle.FaceAlpha = 0.98;

    % The visible wall marks alpha=1/28 day^-1, separating the theoretical
    % slower-than-normal region from the normal-tissue-referenced domain.
    % It is an orientation reference, not an outcome threshold.
    wallX = [result.rGrid(1) result.rGrid(end); ...
        result.rGrid(1) result.rGrid(end)];
    wallY = settings.normalAlphaLower*ones(2,2);
    wallZ = [0 0;1.04*maximumZ 1.04*maximumZ];
    surf(ax,wallX,wallY,wallZ,'FaceColor',[0.38 0.38 0.38], ...
        'FaceAlpha',0.16,'EdgeColor',[0.25 0.25 0.25], ...
        'EdgeAlpha',0.32,'HandleVisibility','off');

    if showReductionDomain
        addNormalReferenceReductionDomain(ax,result,settings)
    end

    xlabel(ax,'Tumor growth rate r [day^{-1}]','Interpreter','tex')
    ylabel(ax,{'Differentiation-like T_s-to-T_u transition', ...
        'rate \alpha [day^{-1}] (log scale)'},'Interpreter','tex')
    zlabel(ax,zLabelText,'Interpreter','tex')
    title(ax,titleText,'Interpreter','tex')
    ax.XScale = 'linear';
    ax.YScale = 'log';
    ax.XDir = 'normal';
    ax.YDir = 'normal';
    applyAlphaTicks(ax,result,settings)
    applyRTicks(ax,result)
    xlim(ax,[result.rGrid(1),result.rGrid(end)])
    ylim(ax,[result.alphaGrid(1),result.alphaGrid(end)])
    zlim(ax,[0,1.06*maximumZ])
    caxis(ax,[0,maximumZ])
    colormap(ax,vividSequentialMap(256))
    ax.Projection = 'orthographic';
    pbaspect(ax,[1.2 1 0.85])
    stylePublicationAxes(ax,'three-dimensional')
    ax.FontSize = 11.5;
end

function addNormalReferenceReductionDomain(ax,result,settings)
    % The green floor contains only coordinates within the
    % normal-tissue-referenced band, 1/28<alpha<=1/4, at which mean total
    % burden is lower than the paired alpha=1/28 mean at the same r.
    [rMatrix,alphaMatrix] = meshgrid(result.rGrid,result.alphaGrid);
    domainFloor = zeros(size(result.normalReferenceReductionDomain));
    domainFloor(~result.normalReferenceReductionDomain) = NaN;
    surf(ax,rMatrix,alphaMatrix,domainFloor, ...
        'FaceColor',[0.10 0.58 0.24],'FaceAlpha',0.30, ...
        'EdgeColor','none','HandleVisibility','off');
    if settings.normalAlphaLower~=1/28
        error('The displayed reduction domain requires alpha=1/28 reference.')
    end
end

function map = vividSequentialMap(numberColors)
    if exist('turbo','file')==2 || exist('turbo','builtin')==5
        map = turbo(numberColors);
    else
        warning(['This MATLAB release does not provide turbo. ' ...
            'Using parula without altering any data or analysis.'])
        map = parula(numberColors);
    end
end

function fig = publicationFigure(name,positionValue,settings)
    visibility = 'on';
    if ~settings.showFigures; visibility = 'off'; end
    fig = figure('Color','w','Name',name,'Visible',visibility, ...
        'Units','inches','Position',positionValue,'Renderer','painters');
    if ~strcmpi(settings.runMode,'paper')
        annotation(fig,'textbox',[0.012 0.965 0.34 0.024], ...
            'String',sprintf('%s MODE -- NOT PUBLICATION DATA', ...
            upper(settings.runMode)), ...
            'HorizontalAlignment','left','VerticalAlignment','middle', ...
            'EdgeColor','none','Color',[0.65 0 0], ...
            'FontName','Arial','FontSize',9.5,'FontWeight','bold');
    end
end

function stylePublicationAxes(ax,dimensionMode)
    grid(ax,'on')
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.08;
    ax.FontName = 'Arial';
    ax.FontSize = 13.5;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.05;
    ax.TickDir = 'out';
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontWeight = 'bold';
    ax.ZLabel.FontWeight = 'bold';
    ax.Title.FontWeight = 'bold';
    ax.XLabel.FontSize = 14.5;
    ax.YLabel.FontSize = 14.5;
    ax.ZLabel.FontSize = 13.5;
    ax.Title.FontSize = 15;
    if strcmpi(dimensionMode,'two-dimensional')
        box(ax,'off')
        ax.XAxisLocation = 'bottom';
        ax.YAxisLocation = 'left';
    else
        box(ax,'on')
    end
end

function stylePublicationColorbar(colorbarHandle)
    colorbarHandle.FontName = 'Arial';
    colorbarHandle.FontSize = 11.5;
    colorbarHandle.FontWeight = 'bold';
    colorbarHandle.Label.FontWeight = 'bold';
end

function applyAlphaTicks(ax,result,settings)
    ticks = unique([settings.theoreticalAlphaMinimum, ...
        settings.normalAlphaLower,settings.normalAlphaUpper, ...
        settings.rAlphaPrimaryUpper,result.alphaGrid(end)]);
    ticks = ticks(ticks>=result.alphaGrid(1) & ticks<=result.alphaGrid(end));
    yticks(ax,ticks)
    yticklabels(ax,alphaTickLabels(ticks))
end

function applyRTicks(ax,result)
    upper = result.rGrid(end);
    if upper<=0.5+1e-14
        ticks = [result.rGrid(1) 0.1 0.2 0.3 0.4 0.5];
    else
        ticks = [result.rGrid(1) 0.25 0.5 0.75 1 1.25 upper];
    end
    ticks = unique(ticks(ticks>=result.rGrid(1) & ticks<=upper));
    xticks(ax,ticks)
    xticklabels(ax,compose('%.4g',ticks))
end

function labels = alphaTickLabels(ticks)
    labels = strings(size(ticks));
    for j = 1:numel(ticks)
        if any(abs(ticks(j)-[1/280 1/28 1/4])<1e-12)
            reciprocal = round(1/ticks(j));
            labels(j) = sprintf('1/%d',reciprocal);
        else
            labels(j) = sprintf('%.4g',ticks(j));
        end
    end
end

function addOutsidePanelLetter(fig,ax,panelLetter)
    axesPosition = ax.Position;
    letterWidth = 0.032;
    letterHeight = 0.040;
    letterX = max(0.002,axesPosition(1)-0.060);
    letterY = min(0.965,axesPosition(2)+axesPosition(4)+0.020);
    annotation(fig,'textbox',[letterX letterY letterWidth letterHeight], ...
        'String',panelLetter,'EdgeColor','none','FitBoxToText','on', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Arial','FontSize',20,'FontWeight','bold', ...
        'Color',[0.08 0.08 0.08]);
end

function exportSelectedFigures(figures,settings)
    if numel(figures)~=numel(settings.figureFiles)
        error('Selected figure count does not match the declared filenames.')
    end
    for j = 1:numel(figures)
        outputFile = fullfile(settings.outputDir,settings.figureFiles{j});
        if exist(outputFile,'file')==2; delete(outputFile); end
        exportgraphics(figures(j),outputFile,'ContentType','vector', ...
            'Resolution',settings.exportResolution);
        fprintf('Saved selected figure %d/%d: %s\n',j,numel(figures), ...
            settings.figureFiles{j})
    end
end

%% ========================================================================
% LOCAL FUNCTIONS: PARALLEL EXECUTION AND PROGRESS REPORTING
% ========================================================================

function useParallel = startParallelIfPossible(requestParallel,nWorkers)

    useParallel = false;
    if ~requestParallel
        return
    end
    if exist('parpool','file')~=2 || exist('gcp','file')~=2
        warning('Parallel Computing Toolbox not found. Running serially.')
        return
    end

    try
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local',nWorkers);
        end
        useParallel = true;
    catch exception
        warning('Could not start the parallel pool; running serially: %s', ...
            exception.message)
    end
end

function reportProgressLine(label,completed,total,elapsedSeconds)

    fractionComplete = completed/total;
    if completed>0
        estimatedTotal = elapsedSeconds/fractionComplete;
        remainingSeconds = max(0,estimatedTotal-elapsedSeconds);
    else
        remainingSeconds = NaN;
    end
    fprintf('[%s] %6.2f%% (%d/%d) | elapsed %s | remaining %s\n', ...
        label,100*fractionComplete,completed,total, ...
        formatElapsedTime(elapsedSeconds),formatElapsedTime(remainingSeconds))
end

function textValue = formatElapsedTime(secondsValue)

    if ~isfinite(secondsValue)
        textValue = 'unknown';
        return
    end
    secondsValue = max(0,round(secondsValue));
    hours = floor(secondsValue/3600);
    minutes = floor(mod(secondsValue,3600)/60);
    seconds = mod(secondsValue,60);
    if hours>0
        textValue = sprintf('%dh %02dm %02ds',hours,minutes,seconds);
    elseif minutes>0
        textValue = sprintf('%dm %02ds',minutes,seconds);
    else
        textValue = sprintf('%ds',seconds);
    end
end
