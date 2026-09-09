% Variance-based global sensitivity analysis (GSA) for the 
% seven-state MMC/NMIBC model:
%       y = [M; Ts; Tu; Di; Dm; E; R].
%
% Scientific question
% -------------------
% Which uncertain kinetic parameters and initial conditions explain
% variation in modeled total tumor burden at day 360 under the 12-session
% MMC regimen?
%
% Primary output
% --------------
%       Y = log10(Ts(360) + Tu(360) + 1).
%
% Governing system 
% ---------------------------------------
% State order: y=[M;Ts;Tu;Di;Dm;E;R], T=Ts+Tu, F=M/(M+a),
% Q=exp[-(T/k)(R/b)(1-F)], q1=beta1(Tu+theta1*Ts),
% q2=beta2(Tu+theta1*Ts),
% A=p3*q1/(q1+h)+p4*F*q2/(q2+h).
%
% dM/dt  = -mu1*M + mCurrent
% dTs/dt = (r*Ts*(1-F)-alpha*Ts)(1-T/k)
%           -Ts[(1-theta2)*p5*E*Q+p1*F+mus]
% dTu/dt = (r*Tu*(1-F)+alpha*Ts)(1-T/k)
%           -Tu[p5*E*Q+p2*F+muu]
% dDi/dt = d0-mu2*Di-Di*A
% dDm/dt = Di*A-mu2*Dm
% dE/dt  = gamma*Dm-p6*R*E-mu3*E
% dR/dt  = eta*Di[(d0/mu2-Di)/(d0/mu2)](1-F)-mu4*R-p7*R*F
%
% theta1 occurs only in q1 and q2; theta2 is the protected Ts fraction;
% p8 is absent from the ODE; the Ts-to-Tu transition is one-way; no MMC
% factor multiplies dE/dt.
%
% Sobol inputs
% ------------
% The analysis always includes 26 ranged inputs: 23 kinetic parameters and
% T0, p8 and Di0. Fixed quantities are tau, d0, p6, p7, b, M0, Dm0, E0 and
% R0. T0 is mapped conditionally below the sampled carrying capacity k.
% p8 is sampled independently over its supported positive range; it is not
% conditioned on T0. Therefore, the continuous initial Ts or Tu cell-
% equivalent count may be below one. Such values are retained without
% rejection, resampling, rounding or absorption.
%
% The underlying Sobol coordinates are mutually independent. The physical
% T0 range depends on k. Therefore, indices for k and T0 quantify their
% corresponding independent generator coordinates under the declared
% conditional T0 mapping; the k coordinate also changes the admissible T0
% range. 
%
% Sampling rules and bounds
% -------------------------
% Exact numeric bounds are declared once in defineBounds(). Positive ordinary
% ranges spanning at least 100-fold are mapped log-uniformly; other ordinary
% ranges are uniform. Prespecified exceptions are: r is uniform in the
% primary GSA and log-uniform only in its robustness analysis; alpha is
% positive and log-uniform on [1/280,1/28] day^-1. T0 is conditionally
% log-uniform on [2,min(2.88737e11,k)] cells. p8 is independently
% log-uniform on [4e-6,0.798]. Di0 is uniform on [108,10000] cells.
% Broader controlled r-alpha grids belong to the separate response-surface
% experiment and are not GSA input distributions.
%
% Sampling-scale robustness
% -------------------------
% The primary analysis uses the specified distribution for every input. Because
% the r range spans 80-fold, r is uniform on its stated range in the primary
% analysis. A second complete Sobol analysis reuses the same unit designs but
% maps r log-uniformly, providing a stability check.
%
% Numerical and biological rules
% ------------------------------
% - ode15s; AbsTol=1e-9, RelTol=1e-6, MaxStep=1 day; NonNegative=1:7.
% - Continuous cell-equivalent ODE interpretation throughout the GSA.
% - No one-cell absorbing boundary and no positive-value clamp, reset,
%   rounding, rejection or resampling.
% - NonNegative is numerical protection only and is not an extinction rule.
% - p8 sets Ts0=p8*T0 and Tu0=(1-p8)*T0 and is absent from the ODE right-hand
%   side.
% - theta2 is the protected Ts fraction; effector killing uses 1-theta2.
% - p1, mus, p2 and muu are four independent literature-derived inputs.
% - No MMC multiplier appears in in dE/dt.
%
% Outputs
% -------
% - Raw MAT results and CSV tables for both sampling scenarios.
% - Sum of first-order indices, with a bootstrap 95% interval.
% - A diagnostic output-distribution/combined-index figure.
% - Figure 2, the section figure, in which S1 and ST are each independently
%   ranked and shown with paired-row bootstrap intervals.
% - A supplementary r-scale robustness figure and rank-comparison table.

clearvars
clc
close all
tic

% ========================================================================
% SETTINGS
% =======================================================================

settings = makeSettings();
rng(settings.randomSeed,'twister')

switch lower(settings.runMode)
    case "fast"
        settings.baseSampleSize = settings.fastBaseSampleSize;
        settings.nBootstrap = settings.fastBootstrap;
        settings.requestParallel = false;
        settings.exportResolution = 300;
    case "paper"
        settings.baseSampleSize = settings.paperBaseSampleSize;
        settings.nBootstrap = settings.paperBootstrap;
        settings.requestParallel = true;
        settings.exportResolution = 1200;
    otherwise
        error('Unknown runMode "%s". Use "fast" or "paper".',settings.runMode)
end

if settings.baseSampleSize < 2 || ...
        settings.baseSampleSize ~= floor(settings.baseSampleSize)
    error('The Sobol base sample size must be an integer of at least 2.')
end
if settings.nBootstrap < 1 || settings.nBootstrap ~= floor(settings.nBootstrap)
    error('The number of bootstrap replicates must be a positive integer.')
end

runStamp = datestr(now,'yyyymmdd_HHMMSS');
settings.outputDir = sprintf('GSA_Continuous_Day360_26Inputs_%s_%s', ...
    char(settings.runMode),runStamp);
if settings.saveOutputs && ~exist(settings.outputDir,'dir')
    mkdir(settings.outputDir)
end

bounds = defineBounds();
factors = defineFactors(bounds);
validateSpecification(bounds,factors,settings)

D = numel(factors);
N = settings.baseSampleSize;
useParallel = startParallelIfPossible( ...
    settings.requestParallel,settings.nParallelWorkers);

fprintf('\n============================================================\n')
fprintf(' MMC/NMIBC SOBOL GSA\n')
fprintf('Output: Y = log10(Ts(360) + Tu(360) + 1)\n')
fprintf('Run mode: %s\n',settings.runMode)
fprintf('Base sample size N: %d\n',N)
fprintf('Sobol factors D: %d (23 kinetic parameters + T0 + p8 + Di0)\n',D)
fprintf('Evaluations per sampling scenario: N*(D+2) = %d\n',N*(D+2))
fprintf('Two scenarios (primary and r-log robustness): %d evaluations\n', ...
    2*N*(D+2))
fprintf('Bootstrap replicates: %d\n',settings.nBootstrap)
fprintf('Solver: ode15s | AbsTol %.1e | RelTol %.1e | MaxStep %.1f day\n', ...
    settings.odeAbsTol,settings.odeRelTol,settings.odeMaxStep)
fprintf('MMC sessions (days): ')
fprintf('%g ',settings.doseDays)
fprintf('\n')
fprintf('Tumor-state interpretation: continuous cell equivalents; no one-cell absorption\n')
fprintf('Initial composition: p8 sampled independently of T0; sub-one Ts0/Tu0 retained\n')
fprintf('Parallel requested: %d | Parallel used: %d\n', ...
    settings.requestParallel,useParallel)
fprintf('============================================================\n\n')

samplingTable = makeSamplingTable(bounds,factors);
if settings.saveOutputs
    writetable(samplingTable,fullfile(settings.outputDir, ...
        'GSA_input_bounds_and_sampling_distributions.csv'))
end

% ========================================================================
% COMMON SALTELLI/JANSEN UNIT DESIGN
% ========================================================================

fprintf('Generating the common scrambled Sobol unit designs A and B ...\n')
[Aunit,Bunit,designMethod] = makeBaseUnitMatrices(N,D,settings);

% Use the same bootstrap rows in both sampling-scale scenarios so the
% comparison is not affected by different bootstrap draws.
bootstrapRows = randi(N,N,settings.nBootstrap);

% ========================================================================
% PRIMARY ANALYSIS:  SAMPLING DISTRIBUTIONS
% ========================================================================

primary = runSobolScenario(Aunit,Bunit,bootstrapRows,factors,bounds, ...
    settings,useParallel,'primary','linear');

% ========================================================================
% ROBUSTNESS ANALYSIS: CHANGE ONLY r TO LOG-UNIFORM
% ========================================================================

rLog = runSobolScenario(Aunit,Bunit,bootstrapRows,factors,bounds, ...
    settings,useParallel,'r_log_robustness','log');

% =======================================================================
% NUMERIC OUTPUTS AND ROBUSTNESS SUMMARY
% ========================================================================

primaryTable = makeResultsTable(primary.sobol,factors,'primary');
rLogTable = makeResultsTable(rLog.sobol,factors,'r_log_robustness');
summaryTable = makeScenarioSummaryTable(primary,rLog,N,D,designMethod);
robustnessTable = makeRobustnessTable(primary.sobol,rLog.sobol,factors);
primaryRanges = makePracticalRangeTable(primary.Aphysical, ...
    primary.Bphysical,factors,'primary');
rLogRanges = makePracticalRangeTable(rLog.Aphysical, ...
    rLog.Bphysical,factors,'r_log_robustness');

fprintf('\nPrimary analysis: top inputs by total-order index\n')
disp(primaryTable(1:min(15,height(primaryTable)),:))
fprintf('Primary sum(S1) = %.6f (bootstrap 95%% CI %.6f to %.6f)\n', ...
    primary.sobol.sumS1,primary.sobol.sumS1ciLow, ...
    primary.sobol.sumS1ciHigh)

fprintf('\nr-log robustness analysis: top inputs by total-order index\n')
disp(rLogTable(1:min(15,height(rLogTable)),:))
fprintf('r-log sum(S1) = %.6f (bootstrap 95%% CI %.6f to %.6f)\n', ...
    rLog.sobol.sumS1,rLog.sobol.sumS1ciLow,rLog.sobol.sumS1ciHigh)

fprintf('\nRank stability after changing only r to log-uniform:\n')
fprintf('  Spearman correlation of all ST ranks = %.4f\n', ...
    robustnessTable.Properties.UserData.spearmanST)
fprintf('  Overlap among the top %d ST inputs = %d\n', ...
    settings.topN,robustnessTable.Properties.UserData.topNOverlap)
rRow = robustnessTable.Factor=="r";
fprintf('  r rank: primary %d; r-log robustness %d\n', ...
    robustnessTable.Primary_ST_Rank(rRow), ...
    robustnessTable.Rlog_ST_Rank(rRow))

if settings.saveOutputs
    writetable(primaryTable,fullfile(settings.outputDir, ...
        'Sobol_indices_primary.csv'))
    writetable(rLogTable,fullfile(settings.outputDir, ...
        'Sobol_indices_r_log_robustness.csv'))
    writetable(summaryTable,fullfile(settings.outputDir, ...
        'Sobol_scenario_summary_including_sum_S1.csv'))
    writetable(robustnessTable,fullfile(settings.outputDir, ...
        'Sobol_r_sampling_scale_rank_stability.csv'))
    writetable(primaryRanges,fullfile(settings.outputDir, ...
        'Sobol_practical_ranges_primary.csv'))
    writetable(rLogRanges,fullfile(settings.outputDir, ...
        'Sobol_practical_ranges_r_log_robustness.csv'))

    save(fullfile(settings.outputDir,'Sobol_GSA_Day360_complete_results.mat'), ...
        'settings','bounds','factors','samplingTable','designMethod', ...
        'Aunit','Bunit','primary','rLog','primaryTable','rLogTable', ...
        'summaryTable','robustnessTable','primaryRanges','rLogRanges','-v7.3')
    writeDesignNotes(settings,designMethod,primary,rLog,robustnessTable)
end

% ========================================================================
% FIGURES
% ========================================================================

if settings.makeFigures
    makeDiagnosticSobolFigure(primary,factors,settings)
    makeSectionFigure2(primary.sobol,factors,settings)
    makeRScaleRobustnessFigure(primary.sobol,rLog.sobol,factors,settings, ...
        robustnessTable.Properties.UserData.spearmanST)
end

fprintf('\nCompleted successfully. Outputs saved in:\n%s\n',settings.outputDir)
fprintf('Total elapsed time: %s\n',formatElapsedTime(toc))

% =========================================================================
% LOCAL FUNCTIONS: SETTINGS, BOUNDS AND FACTORS
% ========================================================================

function settings = makeSettings()
    settings.runMode = "paper";       % "fast" tests code; "paper" produces final estimates
    settings.fastBaseSampleSize = 64;
    settings.paperBaseSampleSize = 1024;
    settings.fastBootstrap = 50;
    settings.paperBootstrap = 300;
    settings.randomSeed = 22;

    settings.requestParallel = true;
    settings.nParallelWorkers = 2;
    settings.requireSobolSequence = true;

    settings.tStart = 0;
    settings.tEnd = 360;
    settings.tau = 2/24;
    settings.doseDays = [0 7 14 21 28 35 70 98 126 154 182 210];
    settings.detectionLimit = pi*(1/2)^2*(3*0.01)*1e6;

    settings.odeAbsTol = 1e-9;
    settings.odeRelTol = 1e-6;
    settings.odeMaxStep = 1.0;

    settings.topN = 15;
    settings.figureFont = 'Arial';
    settings.mainFigureWidth = 13.2;
    settings.mainFigureHeight = 5.8;
    settings.saveOutputs = true;
    settings.makeFigures = false;
end

function B = defineBounds()
    B.mu1    = bound(5.34e-5, 92.736,       'log');
    B.m      = bound(5,       71856.28743,  'log');
    B.r      = bound(6.25e-3, 0.5,          'linear');
    B.k      = bound(9.0e7,   2.89e11,      'log');
    B.a      = bound(7.8,     75,           'linear');
    % The GSA uses the the theoretical cancer-associated slowing range.
    % The controlled r-alpha response-surface domain is a separate analysis.
    B.alpha  = bound(1/280,   1/28,         'log');
    B.theta1 = bound(0.158412888,0.35016835,'linear');
    B.theta2 = bound(0.05,    1,            'linear');
    B.p5     = bound(1.4e-6,  6.0e-6,       'linear');
    B.p1     = bound(0.76,    1.07,         'linear');
    B.mus    = bound(0.044,   0.052,        'linear');
    B.p2     = bound(1.657,   2.367,        'linear');
    B.muu    = bound(0.043,   0.059,        'linear');
    B.mu2    = bound(0.023,   0.98,         'linear');
    B.beta1  = bound(7.0e-4,  0.195,        'log');
    B.beta2  = bound(0.0178,  0.406,        'linear');
    B.h      = bound(1,       2.89e11,      'log');
    B.p3     = bound(0.013,   1.347,        'log');
    B.p4     = bound(0,       34.04181843,  'linear');
    B.gamma  = bound(0.0047,  9.12,         'log');
    B.mu3    = bound(0.0105,  0.247,        'linear');
    B.eta    = bound(0.013,   0.212,        'linear');
    B.mu4    = bound(0.0204,  0.888,        'linear');

    B.T0     = bound(2,       2.88737e11,   'conditional-log');
    B.p8     = bound(4e-6,    0.798,        'log');
    B.Di0    = bound(108,     10000,        'linear');

    B.tau.fixed = 2/24;
    B.d0.fixed = 1.032e5;
    B.p6.fixed = 1.44e-5;
    B.p7.fixed = 3.110210655;
    B.b.fixed = 395840.674352314;
    B.M0.fixed = 0;
    B.Dm0.fixed = 1;
    B.E0.fixed = 1;
    B.R0.fixed = 1;
end

function b = bound(lo,hi,modeName)
    b.lo = lo;
    b.hi = hi;
    b.mode = modeName;
end

function factors = defineFactors(B)
    names = {'mu1','m','r','k','a','alpha','theta1','theta2', ...
        'p5','p1','mus','p2','muu','mu2','beta1','beta2','h', ...
        'p3','p4','gamma','mu3','eta','mu4','T0','p8','Di0'};
    labels = {'\mu_1','m','r','k','a','\alpha','\theta_1','\theta_2', ...
        'p_5','p_1','\mu_s','p_2','\mu_u','\mu_2','\beta_1','\beta_2','h', ...
        'p_3','p_4','\gamma','\mu_3','\eta','\mu_4','T_0','p_8','D_i(0)'};

    factors = repmat(struct('name','','label','','role','','lo',NaN, ...
        'hi',NaN,'primaryMode','','robustMode','','constraint',''),numel(names),1);

    for j = 1:numel(names)
        nm = names{j};
        factors(j).name = nm;
        factors(j).label = labels{j};
        factors(j).lo = B.(nm).lo;
        factors(j).hi = B.(nm).hi;
        factors(j).primaryMode = B.(nm).mode;
        factors(j).robustMode = B.(nm).mode;
        if ismember(nm,{'T0','p8','Di0'})
            factors(j).role = 'initial condition';
        else
            factors(j).role = 'kinetic parameter';
        end
    end

    rID = find(strcmp(names,'r'),1,'first');
    factors(rID).robustMode = 'log';
    factors(strcmp(names,'T0')).constraint = ...
        'mapped to [2,min(2.88737e11,k)] with an open coordinate; T0<k';
    factors(strcmp(names,'p8')).constraint = ...
        'independent log-uniform fraction on [4e-6,0.798]; sub-one Ts0 or Tu0 retained';
    factors(strcmp(names,'Di0')).constraint = ...
        'sampled directly; D0=Di0+Dm0 and Dm0=1';
end

function validateSpecification(B,factors,settings)
    expectedNames = {'mu1','m','r','k','a','alpha','theta1','theta2', ...
        'p5','p1','mus','p2','muu','mu2','beta1','beta2','h', ...
        'p3','p4','gamma','mu3','eta','mu4','T0','p8','Di0'};
    if numel(factors) ~= 26 || ~isequal({factors.name},expectedNames)
        error('The GSA factor list must contain the specified 26 inputs in the declared order.')
    end
    if settings.tEnd ~= 360
        error('The  GSA endpoint is day 360.')
    end
    if abs(settings.tau-B.tau.fixed) > eps
        error('The MMC dwell time is inconsistent between settings and bounds.')
    end
    if ~isequal(settings.doseDays,[0 7 14 21 28 35 70 98 126 154 182 210])
        error('The 12-session MMC schedule does not match the administration days.')
    end

    % T0 is conditional and alpha is an explicitly specified log-scale
    % exception. r is the prespecified linear-primary/log-robustness input.
    specialNames = {'T0','alpha','r'};
    ordinary = setdiff(1:numel(factors),find(ismember({factors.name},specialNames)),'stable');
    for j = ordinary
        lo = factors(j).lo;
        hi = factors(j).hi;
        modeName = factors(j).primaryMode;
        required = declaredScale(lo,hi);
        if ~strcmp(modeName,required)
            error(['Primary sampling mismatch for %s: [%g,%g] requires %s ' ...
                'under the 100-fold rule, not %s.'], ...
                factors(j).name,lo,hi,required,modeName)
        end
    end
    if ~strcmp(B.r.mode,'linear') || ~strcmp(factors(strcmp({factors.name},'r')).robustMode,'log')
        error('r must be uniform in the primary analysis and log-uniform only in the robustness analysis.')
    end
    if abs(B.alpha.lo-1/280)>1e-14 || abs(B.alpha.hi-1/28)>1e-14 || ...
            ~strcmp(B.alpha.mode,'log')
        error('alpha must use the positive log range [1/280,1/28] day^-1.')
    end
    if B.theta2.lo~=0.05 || B.theta2.hi~=1 || ~strcmp(B.theta2.mode,'linear')
        error('theta2 must be the protected Ts fraction sampled linearly on [0.05,1].')
    end
    if ~strcmp(B.p8.mode,'log')
        error('p8 must be sampled independently and log-uniformly on its supported range.')
    end
    validateEquationRules(B)
    fprintf(['Specification audit passed: 26 factors, current bounds, current ' ...
        'sampling scales, continuous equations, current MMC schedule and ' ...
        'day-360 endpoint.\n'])
end

function modeName = declaredScale(lo,hi)
    if lo > 0 && hi/lo >= 100
        modeName = 'log';
    else
        modeName = 'linear';
    end
end

function validateEquationRules(B)
    % Independent equation-level tests guard against 
    % formulation errors. These tests do not change any simulated value.
    p.mu1 = sqrt(B.mu1.lo*B.mu1.hi);
    p.r = 0.2;
    p.k = 1e10;
    p.a = mean([B.a.lo B.a.hi]);
    p.alpha = sqrt(B.alpha.lo*B.alpha.hi);
    p.theta1 = mean([B.theta1.lo B.theta1.hi]);
    p.theta2 = 0.5;
    p.p5 = mean([B.p5.lo B.p5.hi]);
    p.p1 = mean([B.p1.lo B.p1.hi]);
    p.mus = mean([B.mus.lo B.mus.hi]);
    p.p2 = mean([B.p2.lo B.p2.hi]);
    p.muu = mean([B.muu.lo B.muu.hi]);
    p.mu2 = mean([B.mu2.lo B.mu2.hi]);
    p.beta1 = sqrt(B.beta1.lo*B.beta1.hi);
    p.beta2 = mean([B.beta2.lo B.beta2.hi]);
    p.h = 1e7;
    p.p3 = sqrt(B.p3.lo*B.p3.hi);
    p.p4 = mean([B.p4.lo B.p4.hi]);
    p.gamma = sqrt(B.gamma.lo*B.gamma.hi);
    p.mu3 = mean([B.mu3.lo B.mu3.hi]);
    p.eta = mean([B.eta.lo B.eta.hi]);
    p.mu4 = mean([B.mu4.lo B.mu4.hi]);
    p.d0 = B.d0.fixed;
    p.p6 = B.p6.fixed;
    p.p7 = B.p7.fixed;
    p.b = B.b.fixed;
    p.p8 = 0.1; % audit only: p8 must not enter the ODE right-hand side

    y = [20;2e6;4e6;2e4;3e3;5e3;2e3];

% Test states for checking that dE/dt has no direct dependence on M.
yNoDrug = y; yNoDrug(1) = 0;
yDrug = y;   yDrug(1) = 1e6;


    % theta2 is protection: increasing theta2 must reduce Ts effector loss,
    % and theta2 must not alter another state equation.
    pLow = p;  pLow.theta2 = B.theta2.lo;
    pHigh = p; pHigh.theta2 = B.theta2.hi;
    fLow = odeSystem(0,y,pLow,0);
    fHigh = odeSystem(0,y,pHigh,0);
    if ~(fHigh(2)>fLow(2))
        error('theta2 protection is reversed; Ts killing must use 1-theta2.')
    end
    assertNearlyEqual(fHigh([1 3:7]),fLow([1 3:7]), ...
        'theta2 appeared outside the Ts effector-killing term.')

    % theta1 belongs only to the DC maturation signal.
    pThetaA = p; pThetaA.theta1 = B.theta1.lo;
    pThetaB = p; pThetaB.theta1 = B.theta1.hi;
    fThetaA = odeSystem(0,y,pThetaA,0);
    fThetaB = odeSystem(0,y,pThetaB,0);
    assertNearlyEqual(fThetaA([1:3 6:7]),fThetaB([1:3 6:7]), ...
        'theta1 appeared outside the DC maturation equations.')
    if all(abs(fThetaA(4:5)-fThetaB(4:5))<1e-12)
        error('theta1 did not affect the DC maturation equations.')
    end

    % p8 is an initialization input only.
    pP8A = p; pP8A.p8 = B.p8.lo;
    pP8B = p; pP8B.p8 = B.p8.hi;
    assertNearlyEqual(odeSystem(0,y,pP8A,0),odeSystem(0,y,pP8B,0), ...
        'p8 appeared in the ODE right-hand side.')

    % No Tu-to-Ts transition is present.
    yNoTs = y; yNoTs(2) = 0;
    fNoTs = odeSystem(0,yNoTs,p,0);
    assertNearlyEqual(fNoTs(2),0,'A Tu-to-Ts source was introduced.')

    % dE/dt contains no MMC-dependent multiplier.
    fENoDrug = odeSystem(0,yNoDrug,p,0);
    fEDrug = odeSystem(0,yDrug,p,0);
    assertNearlyEqual(fENoDrug(6),fEDrug(6), ...
        'An MMC-dependent multiplier was introduced into dE/dt.')

    % At Di=d0/mu2 the  normalized Treg activation source vanishes.
    yHome = y;
    yHome(4) = p.d0/p.mu2;
    fHome = odeSystem(0,yHome,p,0);
    Fhome = yHome(1)/(yHome(1)+p.a);
    expectedR = -p.mu4*yHome(7)-p.p7*yHome(7)*Fhome;
    assertNearlyEqual(fHome(7),expectedR, ...
        'The normalized Treg activation source is not implemented.')

    if isequal([B.p1.lo B.p1.hi],[B.p2.lo B.p2.hi]) || ...
            isequal([B.mus.lo B.mus.hi],[B.muu.lo B.muu.hi])
        error('p1/p2 or mus/muu were incorrectly aliased.')
    end
end

function assertNearlyEqual(actual,expected,errorMessage)
    scale = max([1;abs(actual(:));abs(expected(:))]);
    if any(abs(actual(:)-expected(:))>1e-10*scale)
        error('%s',errorMessage)
    end
end

% ========================================================================
% LOCAL FUNCTIONS: DESIGN TRANSFORMATION AND MODEL EVALUATION
% ========================================================================

function [Aunit,Bunit,method] = makeBaseUnitMatrices(N,D,settings)
    if exist('sobolset','file') ~= 2
        if settings.requireSobolSequence
            error(['sobolset is unavailable. The paper run requires a scrambled ' ...
                'Sobol low-discrepancy design.'])
        end
        Aunit = rand(N,D);
        Bunit = rand(N,D);
        method = 'independent pseudorandom Monte Carlo';
        return
    end

    
    
    p = sobolset(2*D,'Skip',1000,'Leap',0);
    p = scramble(p,'MatousekAffineOwen');
    U = net(p,N);
    Aunit = makeOpenUnit(U(:,1:D));
    Bunit = makeOpenUnit(U(:,D+1:2*D));
    method = 'scrambled Sobol low-discrepancy sequence (MatousekAffineOwen)';
    fprintf('Design method: %s\n',method)
end

function U = makeOpenUnit(U)
    U = min(max(U,1e-12),1-1e-12);
end

function result = runSobolScenario(Aunit,Bunit,bootstrapRows,factors,B, ...
        settings,useParallel,scenarioName,rMode)
    N = size(Aunit,1);
    D = size(Aunit,2);

    fprintf('\n------------------------------------------------------------\n')
    fprintf('SCENARIO: %s\n',scenarioName)
    fprintf('r distribution: %s-uniform on [%g,%g]\n', ...
        rMode,B.r.lo,B.r.hi)
    fprintf('------------------------------------------------------------\n')

    Aphysical = transformUnitDesign(Aunit,factors,B,rMode);
    Bphysical = transformUnitDesign(Bunit,factors,B,rMode);

    fprintf('Evaluating A ...\n')
    [YA,okA,msgA] = evaluatePhysicalDesign(Aphysical,factors,B,settings, ...
        useParallel,[scenarioName ' / A']);
    requireCompleteDesign(okA,msgA,[scenarioName ' / A'])

    fprintf('Evaluating B ...\n')
    [YB,okB,msgB] = evaluatePhysicalDesign(Bphysical,factors,B,settings, ...
        useParallel,[scenarioName ' / B']);
    requireCompleteDesign(okB,msgB,[scenarioName ' / B'])

    YAB = nan(N,D);
    for j = 1:D
        fprintf('Hybrid %2d/%2d: A_B^(%s)\n',j,D,factors(j).name)
        UAB = Aunit;
        UAB(:,j) = Bunit(:,j);
        XAB = transformUnitDesign(UAB,factors,B,rMode);
        [YAB(:,j),okAB,msgAB] = evaluatePhysicalDesign(XAB,factors,B, ...
            settings,useParallel,[scenarioName ' / ' factors(j).name]);
        requireCompleteDesign(okAB,msgAB, ...
            [scenarioName ' / hybrid ' factors(j).name])
    end

    fprintf('Computing Saltelli first-order and Jansen total-order indices ...\n')
    sobol = computeSobolIndices(YA,YB,YAB,factors);
    fprintf('Bootstrapping %d paired-row replicates ...\n',size(bootstrapRows,2))
    sobol = addBootstrapIntervals(sobol,YA,YB,YAB,bootstrapRows,factors);

    result = struct();
    result.scenarioName = scenarioName;
    result.rMode = rMode;
    result.Aphysical = Aphysical;
    result.Bphysical = Bphysical;
    result.YA = YA;
    result.YB = YB;
    result.YAB = YAB;
    result.sobol = sobol;
end

function X = transformUnitDesign(U,factors,B,rMode)
    [N,D] = size(U);
    if D ~= numel(factors)
        error('Unit-design width does not match the factor table.')
    end
    X = nan(N,D);
    names = {factors.name};

    for i = 1:N
        u = U(i,:);

        % Map all ordinary factors first. r is the only scale changed by the
        % robustness scenario. T0 is mapped below after k. p8 is an ordinary
        % independent log-uniform input and is therefore mapped here.
        for j = 1:D
            nm = names{j};
            if strcmp(nm,'T0')
                continue
            end
            modeName = factors(j).primaryMode;
            if strcmp(nm,'r')
                modeName = rMode;
            end
            X(i,j) = sampleByMode(u(j),factors(j).lo,factors(j).hi,modeName);
        end

        kID = find(strcmp(names,'k'),1,'first');
        t0ID = find(strcmp(names,'T0'),1,'first');

        t0Upper = min(B.T0.hi,X(i,kID));
        X(i,t0ID) = sampleByMode(u(t0ID),B.T0.lo,t0Upper,'log');
        if X(i,t0ID) >= X(i,kID)
            error('Conditional mapping failure: T0 must be strictly below k.')
        end
    end
end

function x = sampleByMode(u,lo,hi,modeName)
    if ~(isfinite(u) && u > 0 && u < 1)
        error('Every Sobol coordinate must lie strictly inside (0,1).')
    end
    if hi < lo || ~isfinite(lo) || ~isfinite(hi)
        error('Invalid sampling interval [%g,%g].',lo,hi)
    end
    if hi == lo
        x = lo;
    elseif strcmpi(modeName,'log')
        if lo <= 0
            error('Log-uniform sampling requires a positive lower bound.')
        end
        x = 10^(log10(lo)+u*(log10(hi)-log10(lo)));
    elseif strcmpi(modeName,'linear')
        x = lo+u*(hi-lo);
    else
        error('Unknown sampling mode: %s',modeName)
    end
end

function [Y,ok,messages] = evaluatePhysicalDesign(X,factors,B,settings, ...
        useParallel,designLabel)
    N = size(X,1);
    Y = nan(N,1);
    ok = false(N,1);
    messages = strings(N,1);
    odeOpt = odeset('AbsTol',settings.odeAbsTol, ...
        'RelTol',settings.odeRelTol,'MaxStep',settings.odeMaxStep, ...
        'NonNegative',1:7);

    reportEvery = max(1,ceil(N/20));
    progressStart = tic;

    if useParallel
        q = parallel.pool.DataQueue;
        completed = 0;
        afterEach(q,@updateProgress);
        parfor i = 1:N
            p = physicalRowToPatient(X(i,:),factors,B);
            [Y(i),ok(i),messages(i)] = simulateModelOutput(p,settings,odeOpt);
            send(q,1)
        end
    else
        for i = 1:N
            p = physicalRowToPatient(X(i,:),factors,B);
            [Y(i),ok(i),messages(i)] = simulateModelOutput(p,settings,odeOpt);
            if i==1 || mod(i,reportEvery)==0 || i==N
                reportProgressLine(designLabel,i,N,toc(progressStart))
            end
        end
    end

    function updateProgress(~)
        completed = completed+1;
        if completed==1 || mod(completed,reportEvery)==0 || completed==N
            reportProgressLine(designLabel,completed,N,toc(progressStart))
        end
    end
end

function p = physicalRowToPatient(x,factors,B)
    p = struct();
    for j = 1:numel(factors)
        p.(factors(j).name) = x(j);
    end

    p.tau = B.tau.fixed;
    p.d0 = B.d0.fixed;
    p.p6 = B.p6.fixed;
    p.p7 = B.p7.fixed;
    p.b = B.b.fixed;
    p.M0 = B.M0.fixed;
    p.Dm0 = B.Dm0.fixed;
    p.E0 = B.E0.fixed;
    p.R0 = B.R0.fixed;

    p.Ts0 = p.p8*p.T0;
    p.Tu0 = (1-p.p8)*p.T0;
    if p.T0 >= p.k || p.Ts0 <= 0 || p.Tu0 <= 0
        error('Mapped physical row violates positivity or T0<k.')
    end
    p.y0 = [p.M0;p.Ts0;p.Tu0;p.Di0;p.Dm0;p.E0;p.R0];
end

function requireCompleteDesign(ok,messages,label)
    if all(ok)
        return
    end
    failed = find(~ok);
    nShow = min(10,numel(failed));
    fprintf(2,'\nFailed simulations in %s: %d\n',label,numel(failed))
    for q = 1:nShow
        i = failed(q);
        fprintf(2,'  row %d: %s\n',i,messages(i))
    end
    error(['The GSA requires a complete paired design. Publication output ' ...
        'was stopped rather than deleting failed Sobol rows.'])
end

function reportProgressLine(label,done,total,elapsedSeconds)
    fraction = done/total;
    if fraction > 0
        etaSeconds = elapsedSeconds*(1-fraction)/fraction;
    else
        etaSeconds = NaN;
    end
    fprintf('  %-28s %5d/%5d (%5.1f%%) | elapsed %s | ETA %s\n', ...
        label,done,total,100*fraction,formatElapsedTime(elapsedSeconds), ...
        formatElapsedTime(etaSeconds))
end

% ========================================================================
% LOCAL FUNCTIONS: CONTINUOUS ODE INTEGRATION
% ========================================================================

function [Y,solverOK,message] = simulateModelOutput(p,settings,odeOpt)
    Y = NaN;
    solverOK = false;
    message = "integration not completed";

    breaks = unique([settings.tStart settings.tEnd settings.doseDays ...
        settings.doseDays+p.tau]);
    breaks = breaks(breaks>=settings.tStart & breaks<=settings.tEnd);
    breaks = sort(breaks(:));

    yCurrent = p.y0(:);

    for bID = 2:numel(breaks)
        t0 = breaks(bID-1);
        t1 = breaks(bID);
        if t1 <= t0
            continue
        end
        mCurrent = activeInstillationRate(t0,settings.doseDays,p.tau,p.m);
        try
            [tSegment,ySegment] = ode15s( ...
                @(t,y) odeSystem(t,y,p,mCurrent),[t0 t1],yCurrent,odeOpt);
        catch ME
            message = "ode15s error: "+string(ME.message);
            return
        end
        if isempty(tSegment) || isempty(ySegment) || tSegment(end)<t1 || ...
                any(~isfinite(ySegment(:)))
            message = "ode15s returned an incomplete or nonfinite segment";
            return
        end
        yCurrent = ySegment(end,:).';
    end

    if any(~isfinite(yCurrent)) || any(yCurrent<0)
        message = "nonfinite or negative final state";
        return
    end

    Tfinal = yCurrent(2)+yCurrent(3);
    Y = log10(Tfinal+1);
    solverOK = isfinite(Y);
    if solverOK
        message = "continuous ode15s";
    else
        message = "nonfinite transformed output";
    end
end

function mCurrent = activeInstillationRate(t,pulseStarts,tau,mDose)
    if tau>0 && any(t>=pulseStarts & t<pulseStarts+tau)
        mCurrent = mDose;
    else
        mCurrent = 0;
    end
end

function dydt = odeSystem(~,y,p,mCurrent)
    M = y(1);
    Ts = y(2);
    Tu = y(3);
    Di = y(4);
    Dm = y(5);
    E = y(6);
    R = y(7);

    T = Ts+Tu;
    F = M/(M+p.a);
    immuneFactor = exp(-((T/p.k)*(R/p.b)*(1-F)));
    signal1 = p.beta1*(Tu+p.theta1*Ts);
    signal2 = p.beta2*(Tu+p.theta1*Ts);
    maturation = p.p3*signal1/(signal1+p.h) ...
        + p.p4*F*signal2/(signal2+p.h);

    dydt = zeros(7,1);
    dydt(1) = -p.mu1*M+mCurrent;
    dydt(2) = (p.r*Ts*(1-F)-p.alpha*Ts)*(1-T/p.k) ...
        - Ts*((1-p.theta2)*p.p5*E*immuneFactor+p.p1*F+p.mus);
    dydt(3) = (p.r*Tu*(1-F)+p.alpha*Ts)*(1-T/p.k) ...
        - Tu*(p.p5*E*immuneFactor+p.p2*F+p.muu);
    dydt(4) = p.d0-p.mu2*Di-Di*maturation;
    dydt(5) = Di*maturation-p.mu2*Dm;
    dydt(6) = p.gamma*Dm-p.p6*R*E-p.mu3*E;

    dydt(7) = p.eta*Di*((p.d0/p.mu2-Di)/(p.d0/p.mu2))*(1-F) ...
        -p.mu4*R-p.p7*R*F;
end

%========================================================================
% LOCAL FUNCTIONS: SOBOL ESTIMATORS AND TABLES
% ========================================================================

function sobol = computeSobolIndices(YA,YB,YAB,factors)
    Yall = [YA(:);YB(:)];
    VY = var(Yall,1);
    if ~(isfinite(VY) && VY>0)
        error('The output variance is zero or nonfinite.')
    end

    D = size(YAB,2);
    S1 = nan(D,1);
    ST = nan(D,1);
    for j = 1:D
        % A_B^(j) is A with coordinate j replaced by coordinate j from B.
        % Saltelli first-order estimator:
        S1(j) = mean(YB.*(YAB(:,j)-YA))/VY;
        % Jansen total-order estimator:
        ST(j) = mean((YA-YAB(:,j)).^2)/(2*VY);
    end

    [~,rankByST] = sort(ST,'descend');
    [~,rankByS1] = sort(S1,'descend');
    sobol.S1 = S1;
    sobol.ST = ST;
    sobol.outputVariance = VY;
    sobol.sumS1 = sum(S1);
    sobol.sumST = sum(ST);
    sobol.rankByST = rankByST;
    sobol.rankByS1 = rankByS1;
    sobol.factorNames = string({factors.name}).';
end

function sobol = addBootstrapIntervals(sobol,YA,YB,YAB,bootstrapRows,factors)
    nBoot = size(bootstrapRows,2);
    D = size(YAB,2);
    S1boot = nan(nBoot,D);
    STboot = nan(nBoot,D);
    for b = 1:nBoot
        idx = bootstrapRows(:,b);
        tmp = computeSobolIndices(YA(idx),YB(idx),YAB(idx,:),factors);
        S1boot(b,:) = tmp.S1.';
        STboot(b,:) = tmp.ST.';
    end

    S1ci = localPercentile(S1boot,[2.5 97.5]);
    STci = localPercentile(STboot,[2.5 97.5]);
    sumS1ci = localPercentile(sum(S1boot,2),[2.5 97.5]);
    sumSTci = localPercentile(sum(STboot,2),[2.5 97.5]);

    sobol.S1ciLow = S1ci(1,:).';
    sobol.S1ciHigh = S1ci(2,:).';
    sobol.STciLow = STci(1,:).';
    sobol.STciHigh = STci(2,:).';
    sobol.sumS1ciLow = sumS1ci(1);
    sobol.sumS1ciHigh = sumS1ci(2);
    sobol.sumSTciLow = sumSTci(1);
    sobol.sumSTciHigh = sumSTci(2);
    sobol.S1boot = S1boot;
    sobol.STboot = STboot;
end

function Q = localPercentile(X,pct)
    if isvector(X)
        X = X(:);
    end
    D = size(X,2);
    Q = nan(numel(pct),D);
    for j = 1:D
        x = sort(X(isfinite(X(:,j)),j));
        n = numel(x);
        if n==0
            continue
        end
        for pID = 1:numel(pct)
            pos = 1+(n-1)*pct(pID)/100;
            lo = floor(pos);
            hi = ceil(pos);
            if lo==hi
                Q(pID,j) = x(lo);
            else
                Q(pID,j) = x(lo)+(pos-lo)*(x(hi)-x(lo));
            end
        end
    end
end

function T = makeSamplingTable(B,factors)
    nRanged = numel(factors);
    fixedNames = {'tau','d0','p6','p7','b','M0','Dm0','E0','R0'};
    n = nRanged+numel(fixedNames);
    Factor = strings(n,1);
    Role = strings(n,1);
    LowerBound = nan(n,1);
    UpperBound = nan(n,1);
    Primary_Sobol_GSA = strings(n,1);
    Rlog_Robustness_GSA = strings(n,1);
    Constraint = strings(n,1);

    for j = 1:nRanged
        f = factors(j);
        Factor(j) = string(f.name);
        Role(j) = string(f.role);
        LowerBound(j) = f.lo;
        UpperBound(j) = f.hi;
        Primary_Sobol_GSA(j) = distributionLabel(f.primaryMode);
        Rlog_Robustness_GSA(j) = distributionLabel(f.robustMode);
        Constraint(j) = string(f.constraint);
    end
    for q = 1:numel(fixedNames)
        j = nRanged+q;
        nm = fixedNames{q};
        Factor(j) = string(nm);
        if ismember(nm,{'M0','Dm0','E0','R0'})
            Role(j) = "initial condition";
        else
            Role(j) = "fixed model quantity";
        end
        LowerBound(j) = B.(nm).fixed;
        UpperBound(j) = B.(nm).fixed;
        Primary_Sobol_GSA(j) = "fixed";
        Rlog_Robustness_GSA(j) = "fixed";
    end
    Constraint(Factor=="tau") = "2-hour MMC dwell";
    Constraint(ismember(Factor,["d0","p6","p7","b"])) = ...
        "fixed because no ranged value is used";
    Constraint(Factor=="M0") = "no MMC before the first instillation";
    Constraint(ismember(Factor,["Dm0","E0","R0"])) = "fixed at 1";

    T = table(Factor,Role,LowerBound,UpperBound, ...
        Primary_Sobol_GSA,Rlog_Robustness_GSA,Constraint);
end

function txt = distributionLabel(modeName)
    switch lower(modeName)
        case 'linear'
            txt = "uniform via Sobol coordinate";
        case 'log'
            txt = "log-uniform via Sobol coordinate";
        case 'conditional-log'
            txt = "conditional log-uniform via Sobol coordinate";
        otherwise
            txt = string(modeName);
    end
end

function T = makeResultsTable(sobol,factors,scenarioName)
    Factor = string({factors.name}).';
    Label = string({factors.label}).';
    Role = string({factors.role}).';
    LowerBound = [factors.lo].';
    UpperBound = [factors.hi].';
    Sampling = strings(numel(factors),1);
    for j = 1:numel(factors)
        if strcmp(scenarioName,'primary')
            modeName = factors(j).primaryMode;
        else
            modeName = factors(j).robustMode;
        end
        Sampling(j) = distributionLabel(modeName);
    end
    S1 = sobol.S1;
    S1_CI_Low = sobol.S1ciLow;
    S1_CI_High = sobol.S1ciHigh;
    ST = sobol.ST;
    ST_CI_Low = sobol.STciLow;
    ST_CI_High = sobol.STciHigh;
    ST_Rank = nan(numel(factors),1);
    S1_Rank = nan(numel(factors),1);
    ST_Rank(sobol.rankByST) = 1:numel(factors);
    S1_Rank(sobol.rankByS1) = 1:numel(factors);
    T = table(Factor,Label,Role,LowerBound,UpperBound,Sampling, ...
        S1,S1_CI_Low,S1_CI_High,S1_Rank,ST,ST_CI_Low,ST_CI_High,ST_Rank);
    T = sortrows(T,'ST_Rank','ascend');
end

function T = makeScenarioSummaryTable(primary,rLog,N,D,designMethod)
    Scenario = ["Primary: uniform r";"Robustness: r log-uniform"];
    R_Sampling = ["uniform";"log-uniform"];
    BaseSampleSize = repmat(N,2,1);
    NumberOfFactors = repmat(D,2,1);
    ModelEvaluations = repmat(N*(D+2),2,1);
    Design = repmat(string(designMethod),2,1);
    OutputVariance = [primary.sobol.outputVariance;rLog.sobol.outputVariance];
    Sum_S1 = [primary.sobol.sumS1;rLog.sobol.sumS1];
    Sum_S1_CI_Low = [primary.sobol.sumS1ciLow;rLog.sobol.sumS1ciLow];
    Sum_S1_CI_High = [primary.sobol.sumS1ciHigh;rLog.sobol.sumS1ciHigh];
    Sum_ST = [primary.sobol.sumST;rLog.sobol.sumST];
    Sum_ST_CI_Low = [primary.sobol.sumSTciLow;rLog.sobol.sumSTciLow];
    Sum_ST_CI_High = [primary.sobol.sumSTciHigh;rLog.sobol.sumSTciHigh];
    T = table(Scenario,R_Sampling,BaseSampleSize,NumberOfFactors, ...
        ModelEvaluations,Design,OutputVariance,Sum_S1,Sum_S1_CI_Low, ...
        Sum_S1_CI_High,Sum_ST,Sum_ST_CI_Low,Sum_ST_CI_High);
end

function T = makeRobustnessTable(primary,rLog,factors)
    D = numel(factors);
    Factor = string({factors.name}).';
    Label = string({factors.label}).';
    Primary_S1 = primary.S1;
    Rlog_S1 = rLog.S1;
    Primary_ST = primary.ST;
    Rlog_ST = rLog.ST;
    Primary_S1_Rank = nan(D,1);
    Rlog_S1_Rank = nan(D,1);
    Primary_ST_Rank = nan(D,1);
    Rlog_ST_Rank = nan(D,1);
    Primary_S1_Rank(primary.rankByS1) = 1:D;
    Rlog_S1_Rank(rLog.rankByS1) = 1:D;
    Primary_ST_Rank(primary.rankByST) = 1:D;
    Rlog_ST_Rank(rLog.rankByST) = 1:D;
    ST_Rank_Change = Rlog_ST_Rank-Primary_ST_Rank;
    T = table(Factor,Label,Primary_S1,Rlog_S1,Primary_S1_Rank, ...
        Rlog_S1_Rank,Primary_ST,Rlog_ST,Primary_ST_Rank, ...
        Rlog_ST_Rank,ST_Rank_Change);

    d = Primary_ST_Rank-Rlog_ST_Rank;
    spearmanST = 1-6*sum(d.^2)/(D*(D^2-1));
    topN = min(15,D);
    topNOverlap = numel(intersect(primary.rankByST(1:topN), ...
        rLog.rankByST(1:topN)));
    T.Properties.UserData = struct('spearmanST',spearmanST, ...
        'topNOverlap',topNOverlap);
end

function T = makePracticalRangeTable(Aphysical,Bphysical,factors,scenarioName)
    X = [Aphysical;Bphysical];
    Scenario = repmat(string(scenarioName),numel(factors),1);
    Factor = string({factors.name}).';
    PracticalMinimum = min(X,[],1).';
    PracticalMaximum = max(X,[],1).';
    T = table(Scenario,Factor,PracticalMinimum,PracticalMaximum);
end

% ========================================================================
% LOCAL FUNCTIONS: FIGURES
% ========================================================================

function makeDiagnosticSobolFigure(primary,factors,settings)
    Ydist = [primary.YA;primary.YB];
    fig = figure('Color','w','Name','Primary day-360 Sobol GSA');
    set(fig,'Units','inches','Position',[1 1 settings.mainFigureWidth ...
        settings.mainFigureHeight]);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    axA = nexttile(tl,1);
    histogram(axA,Ydist,32,'FaceColor',[0.35 0.35 0.35], ...
        'EdgeColor','none','FaceAlpha',0.88);
    hold(axA,'on')
    xRef = log10(settings.detectionLimit+1);
    yl = ylim(axA);
    plot(axA,[xRef xRef],yl,'--','Color',[0.1 0.1 0.1],'LineWidth',1.35)
    ylim(axA,yl)
    text(axA,xRef+0.08,yl(1)+0.82*range(yl),'1-mm detection threshold', ...
        'FontWeight','bold','FontSize',9.5,'BackgroundColor','w', ...
        'Margin',1,'Clipping','on')
    xlabel(axA,'log_{10}(T_s(360)+T_u(360)+1)','FontWeight','bold', ...
        'Interpreter','tex')
    ylabel(axA,'Number of simulations','FontWeight','bold')
    title(axA,'Day-360 GSA output','FontWeight','bold')
    formatAxes(axA,settings)
    addPanelLetter(axA,'A')

    axB = nexttile(tl,2);
    plotCombinedSobolPanel(axB,primary.sobol,factors,settings)
    title(axB,'Ranked sensitivity indices','FontWeight','bold')
    addPanelLetter(axB,'B')

    exportFigure(fig,settings,'FigS_GSA_Output_and_Combined_Indices')
end

function plotCombinedSobolPanel(ax,sobol,factors,settings)
    topN = min(settings.topN,numel(factors));
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
        drawHorizontalCI(ax,y(i)+0.18,STlo(i),SThi(i))
        drawHorizontalCI(ax,y(i)-0.18,S1lo(i),S1hi(i))
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
        'Location','northeast','Box','off','Interpreter','tex')
    formatAxes(ax,settings)
end

function makeSectionFigure2(sobol,factors,settings)
    fig = figure('Color','w','Name','Figure 2: global sensitivity analysis');
    set(fig,'Units','inches','Position',[1 1 13.4 6.2]);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    topN = min(settings.topN,numel(factors));
    xCommon = 1.12*max([sobol.S1ciHigh(:);sobol.STciHigh(:); ...
        sobol.S1(:);sobol.ST(:);0.05]);
    xCommon = min(1,max(0.05,xCommon));

    % Each panel is sorted by the index it displays. This 
    % removes the nonmonotonic S1 ordering produced when both panels share
    % the ST ranking.
    axA = nexttile(tl,1);
    idxS1 = flip(sobol.rankByS1(1:topN));
    plotSingleIndexPanel(axA,sobol.S1,sobol.S1ciLow,sobol.S1ciHigh, ...
        idxS1,factors,'First-order Sobol index, S_1', ...
        'Independent effects',[0.82 0.22 0.13],xCommon,settings)
    addPanelLetter(axA,'A')

    axB = nexttile(tl,2);
    idxST = flip(sobol.rankByST(1:topN));
    plotSingleIndexPanel(axB,sobol.ST,sobol.STciLow,sobol.STciHigh, ...
        idxST,factors,'Total-order Sobol index, S_T', ...
        'Total effects including interactions',[0.22 0.43 0.70], ...
        xCommon,settings)
    addPanelLetter(axB,'B')

    exportFigure(fig,settings,'Fig02_GSA_Day360_Sobol_S1_ST')
end

function plotSingleIndexPanel(ax,values,ciLow,ciHigh,idx,factors, ...
        xLabel,panelTitle,faceColor,xCommon,settings)
    y = 1:numel(idx);
    vals = max(values(idx),0);
    lo = max(ciLow(idx),0);
    hi = max(ciHigh(idx),0);
    hold(ax,'on')
    barh(ax,y,vals,0.68,'FaceColor',faceColor,'EdgeColor','none', ...
        'FaceAlpha',0.9)
    for i = 1:numel(idx)
        drawHorizontalCI(ax,y(i),lo(i),hi(i))
    end
    yticks(ax,y)
    yticklabels(ax,{factors(idx).label})
    set(ax,'TickLabelInterpreter','tex')
    xlabel(ax,xLabel,'FontWeight','bold','Interpreter','tex')
    ylabel(ax,'Model input','FontWeight','bold')
    title(ax,panelTitle,'FontWeight','bold')
    xlim(ax,[0 xCommon])
    ylim(ax,[0.35 numel(idx)+0.65])
    formatAxes(ax,settings)
end

function makeRScaleRobustnessFigure(primary,rLog,factors,settings,rho)
    D = numel(factors);
    topN = min(settings.topN,D);
    [~,rankByMax] = sort(max([primary.ST rLog.ST],[],2),'descend');
    idx = flip(rankByMax(1:topN));
    y = 1:topN;
    values = [max(primary.ST(idx),0) max(rLog.ST(idx),0)];

    fig = figure('Color','w','Name','Robustness to r sampling scale');
    set(fig,'Units','inches','Position',[1 1 8.2 6.2]);
    ax = axes(fig);
    barh(ax,y,values,0.75,'grouped','EdgeColor','none')
    yticks(ax,y)
    yticklabels(ax,{factors(idx).label})
    set(ax,'TickLabelInterpreter','tex')
    xlabel(ax,'Total-order Sobol index, S_T','FontWeight','bold', ...
        'Interpreter','tex')
    ylabel(ax,'Model input','FontWeight','bold')
    title(ax,sprintf('Robustness to the sampling distribution of r (rank \\rho = %.3f)',rho), ...
        'FontWeight','bold','Interpreter','tex')
    legend(ax,{'Primary: uniform r','Robustness: log-uniform r'}, ...
        'Location','northeast','Box','off')
    formatAxes(ax,settings)
    exportFigure(fig,settings,'FigS_GSA_r_Sampling_Scale_Robustness')
end

function drawHorizontalCI(ax,y,lo,hi)
    if ~(isfinite(lo) && isfinite(hi))
        return
    end
    if hi<lo
        tmp = lo;
        lo = hi;
        hi = tmp;
    end
    plot(ax,[lo hi],[y y],'-','Color',[0.05 0.05 0.05],'LineWidth',1.1)
    plot(ax,[lo lo],[y-0.11 y+0.11],'-','Color',[0.05 0.05 0.05],'LineWidth',1.1)
    plot(ax,[hi hi],[y-0.11 y+0.11],'-','Color',[0.05 0.05 0.05],'LineWidth',1.1)
end

function formatAxes(ax,settings)
    set(ax,'FontName',settings.figureFont,'FontSize',12.5, ...
        'LineWidth',1.25,'TickDir','out','Box','off')
    grid(ax,'off')
end

function addPanelLetter(ax,letterText)
    text(ax,-0.16,1.07,letterText,'Units','normalized','FontSize',20, ...
        'FontWeight','bold','HorizontalAlignment','left', ...
        'VerticalAlignment','bottom','Clipping','off')
end

function exportFigure(fig,settings,fileStem)
    if ~settings.saveOutputs
        return
    end
    pngFile = fullfile(settings.outputDir,[fileStem '.png']);
    pdfFile = fullfile(settings.outputDir,[fileStem '.pdf']);
    figFile = fullfile(settings.outputDir,[fileStem '.fig']);
    if exist('exportgraphics','file')==2
        exportgraphics(fig,pngFile,'Resolution',settings.exportResolution)
        exportgraphics(fig,pdfFile,'ContentType','vector')
    else
        print(fig,pngFile,'-dpng',sprintf('-r%d',settings.exportResolution))
        print(fig,pdfFile,'-dpdf','-painters')
    end
    savefig(fig,figFile)
end

% ========================================================================
% LOCAL FUNCTIONS: OUTPUT NOTES AND UTILITIES
% ========================================================================

function writeDesignNotes(settings,designMethod,primary,rLog,robustnessTable)
    fileName = fullfile(settings.outputDir,'GSA_design_and_interpretation_notes.txt');
    fid = fopen(fileName,'w');
    if fid<0
        warning('Could not write design-notes file.')
        return
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'Output: Y = log10(Ts(360)+Tu(360)+1)\n');
    fprintf(fid,'Design: %s\n',designMethod);
    fprintf(fid,'Factors: 26 (23 kinetic parameters plus T0, p8 and Di0)\n');
    fprintf(fid,'Primary r distribution: uniform\n');
    fprintf(fid,'Robustness r distribution: log-uniform\n');
    fprintf(fid,'Primary sum(S1): %.10g [%.10g, %.10g]\n', ...
        primary.sobol.sumS1,primary.sobol.sumS1ciLow,primary.sobol.sumS1ciHigh);
    fprintf(fid,'r-log sum(S1): %.10g [%.10g, %.10g]\n', ...
        rLog.sobol.sumS1,rLog.sobol.sumS1ciLow,rLog.sobol.sumS1ciHigh);
    fprintf(fid,'Spearman correlation of ST ranks: %.10g\n', ...
        robustnessTable.Properties.UserData.spearmanST);
    fprintf(fid,['Interpretation: classical Sobol indices apply to the independent ' ...
        'unit coordinates. T0 is conditionally mapped given k. Thus k and ' ...
        'T0 indices refer to their coordinates under the declared T0 ' ...
        'generator. p8 is independently log-uniform on [4e-6,0.798].\n']);
    fprintf(fid,['Tumor states are continuous cell equivalents. No one-cell ' ...
        'absorbing boundary or positive-value clamp is used in this GSA.\n']);
end

function useParallel = startParallelIfPossible(requestParallel,nWorkers)
    useParallel = false;
    if ~requestParallel
        return
    end
    if exist('parpool','file')~=2 || exist('gcp','file')~=2
        warning('Parallel Computing Toolbox was not found; running serially.')
        return
    end
    try
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local',nWorkers);
        end
        useParallel = true;
    catch ME
        warning('Could not start the parallel pool; running serially: %s',ME.message)
    end
end

function txt = formatElapsedTime(secondsValue)
    if ~isfinite(secondsValue)
        txt = '--:--:--';
        return
    end
    secondsValue = max(0,round(secondsValue));
    h = floor(secondsValue/3600);
    m = floor(mod(secondsValue,3600)/60);
    s = mod(secondsValue,60);
    txt = sprintf('%02d:%02d:%02d',h,m,s);
end
