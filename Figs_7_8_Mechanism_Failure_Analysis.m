% Mechanism-of-failure analysis for the seven-state NMIBC-MMC model.
% Generates the cohort, runs reference and perturbation simulations,
% and saves numerical results and audit tables for the figures scripts.
%
% Usage
% In the beginning, set runMode below to 'test', 'fast', or 'paper', then Run.
%
% Scientific scope
%   * Five prespecified biological mechanisms:
%       M1 = p1,p2
%       M2 = alpha,theta1,theta2
%       M3 = beta1,beta2,p3,p4
%       M4 = gamma,p5
%       M5 = eta,p6,p7; p6 and p7 remain fixed because only one supported
%            value is available for each, so only eta is varied.
%   * Twenty full-range, joint LHS perturbations are generated per
%     mechanism in paper mode. They are absolute parameter values, not
%     local +/- percentage changes. One saved design per mechanism is
%     applied to every virtual-patient background.
%   * One-cell absorption is used in this mechanism analysis only:
%       - downward entry of Ts into the sub-one region absorbs Ts at zero;
%       - downward entry of Ts+Tu into the sub-one region absorbs both.
%     Exactly one modeled cell is not automatically extinguished and may
%     regrow. Tu has no separate compartment-level absorbing boundary.
%   * Counts quantify model-encoded susceptibility conditional on the
%     mechanism definitions, bounds, sampling measures, and saved LHS
%     designs. They are not normalized rankings of biological importance.
%
% Reproducibility
%   A new cohort is generated in this file; no cohort, checkpoint, result,
%   or figure from the GSA or virtual-patient sections is loaded. Dedicated
%   seeds are declared in makeSettings() and exported with the results.
%   Integration is divided at MMC-session starts and ends. Monthly solver
%   boundaries are not introduced because this section classifies the
%   day-360 endpoint; the 3D display is evaluated on its declared time grid.
%   The separate figures scripts read the saved numerical results.
%
% Required MATLAB functionality
%   MATLAB R2020b or newer. Parallel execution is optional and requires
%   Parallel Computing Toolbox.

clearvars
clc
close all

%% USER SETTINGS
runMode = 'paper';  % 'test', 'fast', or 'paper'
runMode = validatestring(lower(string(runMode)),["test","fast","paper"]);

runTimer = tic;

scriptPath = [mfilename('fullpath') '.m'];
scriptFolder = fileparts(scriptPath);
settings = makeSettings(runMode,scriptFolder);
[spec,fixed] = makeInputSpecification();
mechanisms = makeMechanisms(spec);

validateSpecification(settings,spec,fixed,mechanisms);
runEquationUnitTests(spec,fixed,settings);

if exist(settings.outputDir,'dir')~=7
    mkdir(settings.outputDir)
end

fprintf('\n===============================================================\n')
fprintf('MMC MECHANISM-OF-FAILURE ANALYSIS\n')

fprintf('===============================================================\n')
fprintf('Mode: %s\n',settings.runMode)
fprintf('Virtual patients: %d\n',settings.nPatients)
fprintf('Joint perturbations per mechanism: %d\n',settings.nPerturbations)
fprintf('Absorption: Ts < 1 and total tumor < 1, mechanism section only\n')
fprintf('Solver: ode15s | AbsTol %.1e | RelTol %.1e | MaxStep %.3g day\n', ...
    settings.AbsTol,settings.RelTol,settings.MaxStep)
fprintf('Output directory: %s\n\n',settings.outputDir)
writeReproducibilityManifest(settings,spec,fixed,mechanisms)

% Sampling designs are generated before any parallel loop so serial and
% parallel execution receive identical virtual patients and perturbations.
[patients,cohortU] = generateVirtualPatients( ...
    settings.nPatients,spec,fixed,settings.seedCohort);
perturbationPlans = buildPerturbationPlans( ...
    mechanisms,spec,settings.nPerturbations,settings.seedPerturbations);

validateVirtualPatients(patients,spec,fixed);
validatePerturbationPlans(perturbationPlans,mechanisms,spec,settings);

inputTable = inputSpecificationTable(spec,fixed);
mechanismTable = mechanismDefinitionTable(mechanisms);
cohortTable = virtualPatientTable(patients,cohortU,spec);
perturbationTable = perturbationPlanTable(perturbationPlans,mechanisms);

writetable(inputTable,fullfile(settings.outputDir, ...
    'Input_specification_and_sampling_rules.csv'))
writetable(mechanismTable,fullfile(settings.outputDir, ...
    'Mechanism_definitions.csv'))
writetable(cohortTable,fullfile(settings.outputDir, ...
    'Virtual_patient_inputs.csv'))
writetable(perturbationTable,fullfile(settings.outputDir, ...
    'Mechanism_perturbation_design.csv'))

odeOpt = odeset('AbsTol',settings.AbsTol,'RelTol',settings.RelTol, ...
    'MaxStep',settings.MaxStep,'NonNegative',1:7);

simulation = runPrimarySimulations( ...
    patients,perturbationPlans,settings,fixed,odeOpt);

solverAuditTable = buildSolverAuditTable(simulation);
writetable(solverAuditTable,fullfile(settings.outputDir, ...
    'Numerical_completion_audit.csv'))
if settings.requireCompleteSimulations && ...
        (~all(simulation.ReferenceOK) || ~all(simulation.PerturbedOK(:)))
    save(fullfile(settings.outputDir, ...
        'INCOMPLETE_SIMULATION_DIAGNOSTIC.mat'),'simulation','settings', ...
        'patients','perturbationPlans','-v7.3')
    error(['At least one ODE simulation did not reach day 360. ' ...
        'Publication summaries and figures were not generated; inspect ' ...
        'INCOMPLETE_SIMULATION_DIAGNOSTIC.mat.'])
end

[summaryTable,patientMechanismTable,longTable,overlapTable] = ...
    summarizeMechanismTransitions( ...
    simulation,patients,mechanisms,perturbationPlans,settings);

bootstrapTable = patientClusterBootstrap( ...
    simulation,mechanisms,settings);
summaryTable = attachBootstrapIntervals(summaryTable,bootstrapTable);

if settings.runCohortStability
    cohortStabilityTable = cohortSizeStability( ...
        simulation,mechanisms,settings);
else
    cohortStabilityTable = table();
end

writetable(summaryTable,fullfile(settings.outputDir, ...
    'Mechanism_transition_summary.csv'))
writetable(patientMechanismTable,fullfile(settings.outputDir, ...
    'Patient_by_mechanism_summary.csv'))
writetable(longTable,fullfile(settings.outputDir, ...
    'Patient_mechanism_draw_long_table.csv'))
writetable(overlapTable,fullfile(settings.outputDir, ...
    'Patient_mechanism_overlap.csv'))
writetable(bootstrapTable,fullfile(settings.outputDir, ...
    'Patient_cluster_bootstrap_intervals.csv'))
if ~isempty(cohortStabilityTable)
    writetable(cohortStabilityTable,fullfile(settings.outputDir, ...
        'Cohort_size_stability.csv'))
end
validateCompletedAnalysis( ...
    simulation,summaryTable,patientMechanismTable,longTable, ...
    overlapTable,mechanisms,settings);


% Identify patient-mechanism pairs with control-to-failure transitions.
% Select a representative case using the rule in selectRepresentativeCases.
[candidateCaseTable,selectedCaseTable] = selectRepresentativeCases( ...
    simulation,patients,perturbationPlans,mechanisms,spec, ...
    summaryTable,settings);
writetable(candidateCaseTable,fullfile(settings.outputDir, ...
    'Candidate_patient_mechanism_cases.csv'))
writetable(selectedCaseTable,fullfile(settings.outputDir, ...
    'Selected_candidate_patient_cases.csv'))

% Save the numerical results for the separate figures scripts.
Results = struct();
Results.Settings = settings;
Results.InputSpecification = inputTable;
Results.Mechanisms = mechanismTable;
Results.VirtualPatients = cohortTable;
Results.PerturbationDesign = perturbationTable;
Results.Simulation = simulation;
Results.MechanismSummary = summaryTable;
Results.PatientMechanismSummary = patientMechanismTable;
Results.LongTable = longTable;
Results.PatientOverlap = overlapTable;
Results.Bootstrap = bootstrapTable;
Results.CohortStability = cohortStabilityTable;
Results.SolverAudit = solverAuditTable;
Results.CandidateCases = candidateCaseTable;
Results.SelectedCases = selectedCaseTable;
Results.FigureStyleAudit = table();
save(fullfile(settings.outputDir,'MMC_mechanism_failure_results.mat'), ...
    'Results','patients','cohortU','perturbationPlans','spec','fixed', ...
    'mechanisms','settings','-v7.3')

sourcePath = [mfilename('fullpath') '.m'];
if exist(sourcePath,'file')==2
    copyfile(sourcePath,fullfile(settings.outputDir, ...
        'MMC_Mechanism_Failure_NPJ_EXECUTED_SOURCE.m'))
end

printFinalSummary(summaryTable,simulation,settings,toc(runTimer));
validateOutputManifest(settings)
if ispc
    winopen(settings.outputDir)
end

%% ========================================================================
% SETTINGS
% ========================================================================

function settings = makeSettings(runMode,scriptFolder)
settings = struct();
settings.runMode = char(runMode);

% Version identifier used to validate saved results and checkpoints.
settings.codeVersion = 'MMC_MF_FINAL_CORRECTED_2026-08-23_B';

settings.tEnd = 360;
settings.doseDuration = 2/24;
settings.treatmentTimes = [0 7 14 21 28 35 70 98 126 154 182 210];
settings.oneCellLimit = 1;
settings.detectionLimit = diameterToCells(1);

settings.AbsTol = 1e-9;
settings.RelTol = 1e-6;
settings.MaxStep = 1;

% Independent seeds prevent changes in one sampling stage from silently
% changing another stage.
settings.seedCohort = 11;
settings.seedPerturbations = 22;
settings.seedBootstrap = 33;
settings.seedCohortStability = 44;

settings.figureFont = 'Arial';
settings.figureDPI = 1200;
settings.makeFigures = false;
settings.makePatientFigures = false;
settings.keepFiguresOpen = true;
settings.nSelectedPatients = 1;
settings.figureCodeVersion = ...
    'MMC_MF_EMBEDDED_STYLE_EDITABLE_VALIDATED_EXPORT_R6_2026-08-23';
settings.requireCompleteSimulations = true;
settings.resumeFromCheckpoint = true;
settings.checkpointBlockSize = 10;

switch settings.runMode
    case 'test'
        settings.nPatients = 4;
        settings.nPerturbations = 3;
        settings.nBootstrap = 50;
        settings.surfaceLevels = 21;
        settings.surfaceTimePoints = 91;
        settings.requestParallel = false;
        settings.runCohortStability = false;
        settings.cohortStabilitySizes = [];
        settings.nCohortStabilityReplicates = 0;
    case 'fast'
        settings.nPatients = 50;
        settings.nPerturbations = 20;
        settings.nBootstrap = 500;
        settings.surfaceLevels = 31;
        settings.surfaceTimePoints = 121;
        settings.requestParallel = false;
        settings.runCohortStability = false;
        settings.cohortStabilitySizes = [];
        settings.nCohortStabilityReplicates = 0;
    case 'paper'
        settings.nPatients = 500;
        settings.nPerturbations = 20;
        settings.nBootstrap = 2000;
        settings.surfaceLevels = 41;
        settings.surfaceTimePoints = 181;
        settings.requestParallel = true;
        settings.runCohortStability = true;
        settings.cohortStabilitySizes = [100 250 500];
        settings.nCohortStabilityReplicates = 250;
end

settings.outputDir = fullfile(scriptFolder,sprintf( ...
    'MMC_mechanism_failure_final_%s_results',settings.runMode));
settings.useParallel = startParallelIfPossible(settings.requestParallel);
end

function useParallel = startParallelIfPossible(requestParallel)
useParallel = false;
if ~requestParallel || isempty(ver('parallel'))
    return
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
    useParallel = true;
catch ME
    fprintf(2,'Parallel pool unavailable; continuing serially: %s\n',ME.message)
end
end

%% ========================================================================
%  INPUT SPECIFICATION
% ========================================================================

function [spec,fixed] = makeInputSpecification()
% Twenty-three ranged model parameters and three ranged initial-condition
% inputs: T(0), p8, and Di(0). Sampling measures are specified below.
% Alpha is the explicit positive log-scale
% exception to the general at-least-100-fold rule.

% The supported upper bounds of k and h are 2.89e11 cells. The slightly
% smaller value 2.88737e11 is a separate cap used only for the conditional
% T(0) interval.
T0max = 2.88737e11;
stateScaleMax = 2.89e11;
spec = struct('name',{},'label',{},'lo',{},'hi',{},'measure',{}, ...
    'unit',{},'role',{});

add = @(name,label,lo,hi,measure,unit,role) struct( ...
    'name',name,'label',label,'lo',lo,'hi',hi,'measure',measure, ...
    'unit',unit,'role',role);

spec(end+1) = add('mu1','\mu_1',5.34e-5,92.736,'log', ...
    'day^{-1}','background');
spec(end+1) = add('m','m',5,71856.28743,'log', ...
    'micromolar day^{-1}','background');
spec(end+1) = add('r','r',6.25e-3,0.5,'linear', ...
    'day^{-1}','background');
spec(end+1) = add('k','k',9.0e7,stateScaleMax,'log','cells','background');
spec(end+1) = add('a','a',7.8,75,'linear', ...
    'micromolar','background');

spec(end+1) = add('alpha','\alpha',1/280,1/28,'log', ...
    'day^{-1}','M2');
spec(end+1) = add('theta1','\theta_1',0.158412888,0.35016835, ...
    'linear','dimensionless','M2');
spec(end+1) = add('theta2','\theta_2',0.05,1, ...
    'linear','dimensionless','M2');

spec(end+1) = add('p5','p_5',1.4e-6,6.0e-6,'linear', ...
    'cell^{-1} day^{-1}','M4');
spec(end+1) = add('p1','p_1',0.76,1.07,'linear','day^{-1}','M1');
spec(end+1) = add('mus','\mu_s',0.044,0.052,'linear', ...
    'day^{-1}','background');
spec(end+1) = add('p2','p_2',1.657,2.367,'linear','day^{-1}','M1');
spec(end+1) = add('muu','\mu_u',0.043,0.059,'linear', ...
    'day^{-1}','background');

spec(end+1) = add('mu2','\mu_2',0.023,0.98,'linear', ...
    'day^{-1}','background');
spec(end+1) = add('beta1','\beta_1',7.0e-4,0.195,'log', ...
    'dimensionless','M3');
spec(end+1) = add('beta2','\beta_2',0.0178,0.406,'linear', ...
    'dimensionless','M3');
spec(end+1) = add('h','h',1,stateScaleMax,'log','cells','background');
spec(end+1) = add('p3','p_3',0.013,1.347,'log', ...
    'day^{-1}','M3');
spec(end+1) = add('p4','p_4',0,34.04181843,'linear', ...
    'day^{-1}','M3');

spec(end+1) = add('gamma','\gamma',0.0047,9.12,'log', ...
    'day^{-1}','M4');
spec(end+1) = add('mu3','\mu_3',0.0105,0.247,'linear', ...
    'day^{-1}','background');
spec(end+1) = add('eta','\eta',0.013,0.212,'linear', ...
    'day^{-1}','M5');
spec(end+1) = add('mu4','\mu_4',0.0204,0.888,'linear', ...
    'day^{-1}','background');

spec(end+1) = add('T0','T(0)',2,T0max,'conditional-log', ...
    'cells','initial condition');
spec(end+1) = add('p8','p_8',4e-6,0.798,'conditional', ...
    'fraction','initial condition');
spec(end+1) = add('Di0','D_i(0)',108,10000,'linear', ...
    'cells','initial condition');

fixed = struct();
fixed.d0 = 1.032e5;
fixed.p6 = 1.44e-5;
fixed.p7 = 3.110210655;
fixed.b = 395840.674352314;
fixed.tau = 2/24;
fixed.M0 = 0;
fixed.Dm0 = 1;
fixed.E0 = 1;
fixed.R0 = 1;
fixed.Tmax = T0max;
end

function mechanisms = makeMechanisms(spec)
mechanisms = struct('id',{},'title',{},'assignedInputs',{}, ...
    'variableInputs',{},'fixedInputs',{},'interpretation',{});

mechanisms(1) = makeMechanism('M1','Direct MMC effects', ...
    {'p1','p2'},{'p1','p2'},{}, ...
    'Direct MMC cytotoxicity in the two tumor compartments.');
mechanisms(2) = makeMechanism('M2','Tumor heterogeneity', ...
    {'alpha','theta1','theta2'},{'alpha','theta1','theta2'},{}, ...
    'Differentiation-like transfer and phenotype-dependent immune evasion.');
mechanisms(3) = makeMechanism('M3','ICD/DC maturation', ...
    {'beta1','beta2','p3','p4'},{'beta1','beta2','p3','p4'},{}, ...
    'Baseline and MMC-associated ICD-dependent DC maturation.');
mechanisms(4) = makeMechanism('M4','Antitumor immune response', ...
    {'gamma','p5'},{'gamma','p5'},{}, ...
    'Mature-DC-driven effector activation and tumor-cell killing.');
mechanisms(5) = makeMechanism('M5','Treg-mediated suppression', ...
    {'eta','p6','p7'},{'eta'},{'p6','p7'}, ...
    ['Treg activation, effector suppression, and MMC-dependent Treg ' ...
     'regulation; p6 and p7 remain fixed because only one supported ' ...
     'value is available for each.']);

known = {spec.name};
for g = 1:numel(mechanisms)
    if ~all(ismember(mechanisms(g).variableInputs,known))
        error('Mechanism %s contains an unknown ranged input.',mechanisms(g).id)
    end
end
end

function out = makeMechanism(id,titleText,assigned,variable,fixedInputs,note)
out = struct('id',id,'title',titleText,'assignedInputs',{assigned}, ...
    'variableInputs',{variable},'fixedInputs',{fixedInputs}, ...
    'interpretation',note);
end

%% ========================================================================
% SPECIFICATION VALIDATION
% ========================================================================

function validateSpecification(settings,spec,fixed,mechanisms)
if verLessThan('matlab','9.9')
    error('MATLAB R2020b or newer is required for this analysis.')
end
if numel(spec)~=26
    error('Exactly 26 ranged inputs are required; found %d.',numel(spec))
end
if numel(unique({spec.name}))~=26
    error('Input names must be unique.')
end

requiredNames = {'mu1','m','r','k','a','alpha','theta1','theta2', ...
    'p5','p1','mus','p2','muu','mu2','beta1','beta2','h','p3', ...
    'p4','gamma','mu3','eta','mu4','T0','p8','Di0'};
if ~isequal({spec.name},requiredNames)
    error('The ranged-input order differs from the specified 26-input order.')
end

expectedLower = [5.34e-5,5,6.25e-3,9.0e7,7.8,1/280, ...
    0.158412888,0.05,1.4e-6,0.76,0.044,1.657,0.043,0.023, ...
    7.0e-4,0.0178,1,0.013,0,0.0047,0.0105,0.013,0.0204, ...
    2,4e-6,108];
expectedUpper = [92.736,71856.28743,0.5,2.89e11,75,1/28, ...
    0.35016835,1,6.0e-6,1.07,0.052,2.367,0.059,0.98, ...
    0.195,0.406,2.89e11,1.347,34.04181843,9.12,0.247, ...
    0.212,0.888,2.88737e11,0.798,10000];
expectedMeasure = {'log','log','linear','log','linear','log', ...
    'linear','linear','linear','linear','linear','linear','linear', ...
    'linear','log','linear','log','log','linear','log','linear', ...
    'linear','linear','conditional-log','conditional','linear'};
actualLower = [spec.lo];
actualUpper = [spec.hi];
if ~isequal(actualLower,expectedLower) || ...
        ~isequal(actualUpper,expectedUpper) || ...
        ~isequal({spec.measure},expectedMeasure)
    error(['At least one input bound or computational sampling ' ...
        'measure differs from the mechanism specification.'])
end

for j = 1:numel(spec)
    if ~isfinite(spec(j).lo) || ~isfinite(spec(j).hi) || ...
            spec(j).hi<spec(j).lo
        error('Invalid bounds for %s.',spec(j).name)
    end
    if contains(spec(j).measure,'log') && spec(j).lo<=0
        error('Log-mapped input %s must have a positive lower bound.', ...
            spec(j).name)
    end
end

for j = 1:numel(spec)
    name = spec(j).name;
    if ismember(name,{'T0','p8','alpha'})
        continue
    end
    requiredMeasure = declaredScale(spec(j).lo,spec(j).hi);
    if ~strcmp(spec(j).measure,requiredMeasure)
        error(['Sampling measure mismatch for %s: [%g,%g] requires %s ' ...
            'under the specified range-width rule, not %s.'],name, ...
            spec(j).lo,spec(j).hi,requiredMeasure,spec(j).measure)
    end
end

assertExactBound(spec,'r',6.25e-3,0.5,'linear')
assertExactBound(spec,'alpha',1/280,1/28,'log')
assertExactBound(spec,'theta2',0.05,1,'linear')
assertExactBound(spec,'p8',4e-6,0.798,'conditional')
assertExactBound(spec,'T0',2,2.88737e11,'conditional-log')
assertExactBound(spec,'p1',0.76,1.07,'linear')
assertExactBound(spec,'p2',1.657,2.367,'linear')
assertExactBound(spec,'mus',0.044,0.052,'linear')
assertExactBound(spec,'muu',0.043,0.059,'linear')

if fixed.d0~=1.032e5 || fixed.p6~=1.44e-5 || ...
        fixed.p7~=3.110210655 || fixed.b~=395840.674352314 || ...
        fixed.tau~=settings.doseDuration || ...
        fixed.M0~=0 || fixed.Dm0~=1 || fixed.E0~=1 || fixed.R0~=1 || ...
        fixed.Tmax~=2.88737e11
    error('One or more fixed biological quantities differ from specification.')
end
if settings.oneCellLimit~=1
    error('The absorbing boundary must be one modeled cell.')
end
if settings.detectionLimit~=diameterToCells(1)
    error('The detection limit must use the specified 1-mm conversion.')
end
if ~isequal(settings.treatmentTimes, ...
        [0 7 14 21 28 35 70 98 126 154 182 210])
    error('Treatment start times differ from the 12-session schedule.')
end
if settings.doseDuration~=2/24 || settings.tEnd~=360
    error('Treatment duration or simulation endpoint is incorrect.')
end
if settings.AbsTol~=1e-9 || settings.RelTol~=1e-6 || ...
        settings.MaxStep~=1
    error('The final ode15s tolerance or maximum-step setting changed.')
end
if settings.seedCohort~=11 || settings.seedPerturbations~=22 || ...
        settings.seedBootstrap~=33 || settings.seedCohortStability~=44
    error('A prespecified mechanism-analysis random seed changed.')
end
if activeInstillationRate(0,settings.treatmentTimes, ...
        settings.doseDuration,7)~=7 || ...
        activeInstillationRate(settings.doseDuration/2, ...
        settings.treatmentTimes,settings.doseDuration,7)~=7 || ...
        activeInstillationRate(settings.doseDuration, ...
        settings.treatmentTimes,settings.doseDuration,7)~=0 || ...
        activeInstillationRate(1,settings.treatmentTimes, ...
        settings.doseDuration,7)~=0
    error('Piecewise MMC input does not switch at the declared boundaries.')
end
if strcmp(settings.runMode,'paper') && ...
        (settings.nPatients~=500 || settings.nPerturbations~=20)
    error('Paper mode must use 500 patients and 20 perturbations/mechanism.')
end
if settings.nSelectedPatients~=1 || settings.surfaceLevels<2 || ...
        settings.surfaceTimePoints<2
    error('The primary analysis must select exactly one patient case.')
end

expectedID = {'M1','M2','M3','M4','M5'};
expectedTitle = {'Direct MMC effects','Tumor heterogeneity', ...
    'ICD/DC maturation','Antitumor immune response', ...
    'Treg-mediated suppression'};
expectedAssigned = {{'p1','p2'}, ...
    {'alpha','theta1','theta2'}, ...
    {'beta1','beta2','p3','p4'}, ...
    {'gamma','p5'}, {'eta','p6','p7'}};
expectedVariable = {{'p1','p2'}, ...
    {'alpha','theta1','theta2'}, ...
    {'beta1','beta2','p3','p4'}, ...
    {'gamma','p5'}, {'eta'}};
expectedFixed = {{},{},{},{},{'p6','p7'}};
for g = 1:5
    if ~strcmp(mechanisms(g).id,expectedID{g}) || ...
            ~strcmp(mechanisms(g).title,expectedTitle{g}) || ...
            ~isequal(mechanisms(g).assignedInputs,expectedAssigned{g}) || ...
            ~isequal(mechanisms(g).variableInputs,expectedVariable{g}) || ...
            ~isequal(mechanisms(g).fixedInputs,expectedFixed{g})
        error('The prespecified definition of M%d was altered.',g)
    end
end
end

function assertExactBound(spec,name,lo,hi,measure)
s = getSpec(spec,name);
if s.lo~=lo || s.hi~=hi || ~strcmp(s.measure,measure)
    error('Bound or sampling-measure mismatch for %s.',name)
end
end

function runEquationUnitTests(spec,fixed,settings)
% Test the model equations and one-cell boundary rules independently of the sampled cohort.
p = midpointParameterStruct(spec,fixed);
y = [10;1000;2000;5000;20;30;40];

d1 = mechanismODE(0,y,p,0,false,false);
yHighM = y;
yHighM(1) = 1e9;
d2 = mechanismODE(0,yHighM,p,0,false,false);

% Reconstruct all seven equations independently at the audit state.
M = y(1); Ts = y(2); Tu = y(3); Di = y(4);
Dm = y(5); E = y(6); R = y(7); T = Ts+Tu;
F = M/(M+p.a);
Q = exp(-(T/p.k)*(R/p.b)*(1-F));
q1 = p.beta1*(Tu+p.theta1*Ts);
q2 = p.beta2*(Tu+p.theta1*Ts);
maturation = p.p3*q1/(q1+p.h)+p.p4*F*q2/(q2+p.h);


expectedAll = zeros(7,1);
expectedAll(1) = -p.mu1*M;
expectedAll(2) = (p.r*Ts*(1-F)-p.alpha*Ts)*(1-T/p.k) ...
    -Ts*((1-p.theta2)*p.p5*E*Q+p.p1*F+p.mus);
expectedAll(3) = (p.r*Tu*(1-F)+p.alpha*Ts)*(1-T/p.k) ...
    -Tu*(p.p5*E*Q+p.p2*F+p.muu);
expectedAll(4) = p.d0-p.mu2*Di-Di*maturation;
expectedAll(5) = Di*maturation-p.mu2*Dm;
expectedAll(6) = p.gamma*Dm-p.p6*R*E-p.mu3*E;
expectedAll(7) = ...
    p.eta*Di*((p.d0/p.mu2-Di)/(p.d0/p.mu2))*(1-F) ...
    -p.mu4*R-p.p7*R*F;
assertVectorClose(d1,expectedAll, ...
    'The implemented seven-state system differs from the declared system.')

% E has no direct MMC multiplier: with the same Dm,E,R, its derivative is
% identical when only M changes.
if abs(d1(6)-d2(6))>1e-12*max(1,abs(d1(6)))
    error('Equation audit failed: MMC has entered the E equation.')
end


pGrowth = p;
pGrowth.alpha = 0;
pGrowth.p1 = 0;
pGrowth.p2 = 0;
pGrowth.p5 = 0;
pGrowth.mus = 0;
pGrowth.muu = 0;
dGrowthLowM = mechanismODE(0,y,pGrowth,0,false,false);
dGrowthHighM = mechanismODE(0,yHighM,pGrowth,0,false,false);
if abs(dGrowthLowM(2)-dGrowthHighM(2))<1e-12 || ...
        abs(dGrowthLowM(3)-dGrowthHighM(3))<1e-12
    error('Equation audit failed: tumor MMC growth-arrest factor is absent.')
end

% Differentiation is a balanced Ts-to-Tu transfer inside the common
% capacity factor. With all other tumor terms disabled, total tumor change
% must be zero.
pTransfer = pGrowth;
pTransfer.r = 0;
pTransfer.alpha = p.alpha;
dTransfer = mechanismODE(0,y,pTransfer,0,false,false);
if abs(dTransfer(2)+dTransfer(3))>1e-10*max(1,abs(dTransfer(2)))
    error('Equation audit failed: differentiation transfer is unbalanced.')
end

% theta2 is the protected fraction, so susceptible Ts killing is scaled by
% 1-theta2. Check the full expected Ts derivative explicitly.
expectedTs = (p.r*y(2)*(1-F)-p.alpha*y(2))*(1-T/p.k) ...
    -y(2)*((1-p.theta2)*p.p5*y(6)*Q+p.p1*F+p.mus);
if abs(d1(2)-expectedTs)>1e-10*max(1,abs(expectedTs))
    error('Equation audit failed: theta2 or suppression scaling is wrong.')
end

% Check the normalized Treg source and both Treg loss terms.
expectedR = p.eta*y(4)*((p.d0/p.mu2-y(4))/(p.d0/p.mu2)) ...
    *(1-F)-p.mu4*y(7)-p.p7*y(7)*F;
if abs(d1(7)-expectedR)>1e-10*max(1,abs(expectedR))
    error('Equation audit failed: Treg equation differs from Eq. 7.')
end

% theta1 occurs only in the two DC-maturation signals.
pTheta = p;
pTheta.theta1 = getSpec(spec,'theta1').hi;
dTheta = mechanismODE(0,y,pTheta,0,false,false);
assertVectorClose(dTheta([1:3 6:7]),d1([1:3 6:7]), ...
    'theta1 appeared outside the Di/Dm maturation equations.')
if all(abs(dTheta(4:5)-d1(4:5))<1e-12)
    error('theta1 did not affect the Di/Dm maturation equations.')
end

% p8 sets the initial tumor composition and does not appear explicitly in the ODEs.
pP8a = p; pP8a.p8 = getSpec(spec,'p8').lo;
pP8b = p; pP8b.p8 = getSpec(spec,'p8').hi;
assertVectorClose(mechanismODE(0,y,pP8a,0,false,false), ...
    mechanismODE(0,y,pP8b,0,false,false), ...
    'p8 appeared in the ODE right-hand side.')

% The model includes no Tu-to-Ts transition.
yNoTs = y;
yNoTs(2) = 0;
dNoTs = mechanismODE(0,yNoTs,p,0,false,false);
if dNoTs(2)~=0
    error('A Tu-to-Ts transition was introduced into the tumor equations.')
end

% At zero tumor burden, maturation must vanish.
yNoTumor = y;
yNoTumor(2:3) = 0;
d0 = mechanismODE(0,yNoTumor,p,0,false,false);
if abs(d0(5)+p.mu2*yNoTumor(5))>1e-10
    error('Equation audit failed: DC maturation persists without tumor.')
end

% An absorbed state must have zero tumor derivatives.
dTsAbs = mechanismODE(0,y,p,0,true,false);
dTotAbs = mechanismODE(0,y,p,0,true,true);
if dTsAbs(2)~=0 || any(dTotAbs(2:3)~=0)
    error('Equation audit failed: absorbed tumor state can regrow.')
end

[eventValue,isTerminal,direction] = oneCellBoundaryEvents( ...
    0,[0;1;4;1;1;1;1],true,true,1);
if eventValue(1)~=0 || eventValue(2)~=4 || ...
        any(isTerminal~=[1;1]) || any(direction~=[-1;-1])
    error('Boundary audit failed: Ts event is not a downward one-cell event.')
end
[eventValue,~,~] = oneCellBoundaryEvents( ...
    0,[0;0.4;0.6;1;1;1;1],false,true,1);
if eventValue(2)~=0
    error('Boundary audit failed: total-tumor one-cell surface is wrong.')
end
[eventValue,~,~] = oneCellBoundaryEvents( ...
    0,[0;2;0.5;1;1;1;1],true,true,1);
if any(eventValue==0)
    error('Boundary audit failed: a separate Tu threshold was introduced.')
end
yTsAbs = reconstructAbsorbingState([0;0.5;10;1;1;1;1],true,false);
yTotalAbs = reconstructAbsorbingState([0;0.5;0.4;1;1;1;1],true,true);
if yTsAbs(2)~=0 || yTsAbs(3)~=10 || any(yTotalAbs(2:3)~=0)
    error('Boundary audit failed: absorbed-state reconstruction is wrong.')
end

% Exactly one modeled cell is a valid initial state. This section-specific
% cohort admits neither a sub-one Ts nor a sub-one Tu initial compartment;
% absorption is applied only after initialization and only on a downward
% crossing. Tu still has no separate dynamic absorbing boundary.
if ~hasAdmissibleInitialTumorState([0;1;1;1;1;1;1], ...
        settings.oneCellLimit) || ...
        hasAdmissibleInitialTumorState([0;0.999;2;1;1;1;1], ...
        settings.oneCellLimit) || ...
        hasAdmissibleInitialTumorState([0;2;0.999;1;1;1;1], ...
        settings.oneCellLimit) || ...
        hasAdmissibleInitialTumorState([0;0.4;0.5;1;1;1;1], ...
        settings.oneCellLimit)
    error(['Boundary audit failed: exactly-one and sub-one initial-state ' ...
        'rules are not implemented as declared.'])
end

% A one-cell Ts state with a positive local derivative must be allowed to
% grow; it must not be reset merely because it is exactly on the threshold.
% The positive derivative and the downward-only event direction jointly
% test this rule without starting ode15s exactly on a root surface, whose
% initial-root reporting can differ between MATLAB releases.
yOne = [0;1;1;fixed.d0/p.mu2;1;1;1];
dOne = mechanismODE(0,yOne,p,0,false,false);
if dOne(2)<=0
    error('Boundary audit setup did not give positive one-cell Ts growth.')
end

% After Ts absorption, a single remaining Tu cell is likewise retained and
% can regrow. Tu has no separate compartment-level absorbing boundary; only
% the total-tumor boundary is monitored in this state.
yOneTu = [0;0;1;fixed.d0/p.mu2;1;1;1];
dOneTu = mechanismODE(0,yOneTu,p,0,true,false);
if dOneTu(3)<=0
    error('Boundary audit setup did not give positive one-cell Tu growth.')
end

fprintf('Equation-level unit tests passed.\n')
end

function assertVectorClose(actual,expected,errorMessage)
tolerance = 1e-11*max(1,max(abs(expected(:))));
if any(~isfinite(actual(:))) || any(~isfinite(expected(:))) || ...
        max(abs(actual(:)-expected(:)))>tolerance
    error('%s',errorMessage)
end
end

function p = midpointParameterStruct(spec,fixed)
p = fixed;
for j = 1:numel(spec)
    if ismember(spec(j).name,{'T0','p8','Di0'})
        continue
    end
    if contains(spec(j).measure,'log')
        p.(spec(j).name) = sqrt(spec(j).lo*spec(j).hi);
    else
        p.(spec(j).name) = 0.5*(spec(j).lo+spec(j).hi);
    end
end
end

%% ========================================================================
% LHS AND VIRTUAL-PATIENT GENERATION
% ========================================================================

function [patients,U] = generateVirtualPatients(N,spec,fixed,seed)
U = makeLHS(N,numel(spec),seed);
patients = cell(N,1);
for n = 1:N
    patients{n} = decodeVirtualPatient(U(n,:),spec,fixed);
end
end

function patient = decodeVirtualPatient(u,spec,fixed)
patient = fixed;

% Ordinary model inputs are mapped first so k is available for T(0).
for j = 1:numel(spec)
    name = spec(j).name;
    if ismember(name,{'T0','p8','Di0'})
        continue
    end
    patient.(name) = mapUnitCoordinate(u(j),spec(j).lo, ...
        spec(j).hi,spec(j).measure);
end

jT0 = find(strcmp({spec.name},'T0'));
jP8 = find(strcmp({spec.name},'p8'));
jDi0 = find(strcmp({spec.name},'Di0'));

sT0 = spec(jT0);
T0hi = min(sT0.hi,patient.k);
if T0hi<=sT0.lo
    error('No admissible T(0) interval for k=%g.',patient.k)
end
patient.T0 = mapUnitCoordinate(u(jT0),sT0.lo,T0hi,'log');
patient.T0ConditionalUpper = T0hi;

% Conditional p8 bounds ensure Ts(0) >= 1 and Tu(0) >= 1.
sP8 = spec(jP8);
p8lo = max(sP8.lo,1/patient.T0);
p8hi = min(sP8.hi,1-1/patient.T0);
if p8hi<p8lo
    error('No admissible p8 interval for T(0)=%g.',patient.T0)
elseif p8hi==p8lo
    p8Measure = 'linear';
    patient.p8 = p8lo;
else
    p8Measure = declaredScale(p8lo,p8hi);
    patient.p8 = mapUnitCoordinate(u(jP8),p8lo,p8hi,p8Measure);
end
patient.p8ConditionalLower = p8lo;
patient.p8ConditionalUpper = p8hi;
patient.p8ConditionalMeasure = p8Measure;

patient.Di0 = mapUnitCoordinate(u(jDi0),spec(jDi0).lo, ...
    spec(jDi0).hi,spec(jDi0).measure);
patient.Ts0 = patient.p8*patient.T0;
patient.Tu0 = (1-patient.p8)*patient.T0;
end

function plans = buildPerturbationPlans(mechanisms,spec,nSamples,baseSeed)
G = numel(mechanisms);
plans = repmat(struct('mechanismID','', 'inputNames',{{}}, ...
    'unitCoordinates',[],'values',[],'measures',{{}}),G,1);

for g = 1:G
    names = mechanisms(g).variableInputs;
    U = makeLHS(nSamples,numel(names),baseSeed+g-1);
    values = nan(size(U));
    measures = cell(1,numel(names));
    for q = 1:numel(names)
        s = getSpec(spec,names{q});
        values(:,q) = arrayfun(@(x)mapUnitCoordinate( ...
            x,s.lo,s.hi,s.measure),U(:,q));
        measures{q} = s.measure;
    end
    plans(g).mechanismID = mechanisms(g).id;
    plans(g).inputNames = names;
    plans(g).unitCoordinates = U;
    plans(g).values = values;
    plans(g).measures = measures;
end
end

function U = makeLHS(N,D,seed)
stream = RandStream('mt19937ar','Seed',seed);
U = nan(N,D);
for d = 1:D
    permutation = randperm(stream,N).';
    withinStratum = rand(stream,N,1);
    U(:,d) = (permutation-withinStratum)/N;
end
if any(U(:)<=0 | U(:)>=1)
    error('LHS coordinates must lie strictly inside the unit interval.')
end
end

function value = mapUnitCoordinate(u,lo,hi,measure)
u = min(max(u,eps),1-eps);
switch measure
    case 'linear'
        value = lo+u*(hi-lo);
    case 'log'
        value = exp(log(lo)+u*(log(hi)-log(lo)));
    otherwise
        error('Unsupported direct mapping measure: %s.',measure)
end
end

function measure = declaredScale(lo,hi)
if lo>0 && hi/lo>=100
    measure = 'log';
else
    measure = 'linear';
end
end

function s = getSpec(spec,name)
idx = find(strcmp({spec.name},name),1,'first');
if isempty(idx)
    error('Unknown input: %s.',name)
end
s = spec(idx);
end

function validateVirtualPatients(patients,spec,fixed)
N = numel(patients);
for n = 1:N
    p = patients{n};
    if p.T0<2 || p.T0>=min(fixed.Tmax,p.k)
        error('Virtual patient %d violates 2 <= T(0) < min(Tmax,k).',n)
    end
    if p.Ts0<1-1e-10 || p.Tu0<1-1e-10
        error('Virtual patient %d has a sub-one initial tumor compartment.',n)
    end
    if p.T0ConditionalUpper~=min(fixed.Tmax,p.k) || ...
            p.p8ConditionalLower~=max(getSpec(spec,'p8').lo,1/p.T0) || ...
            p.p8ConditionalUpper~=min(getSpec(spec,'p8').hi,1-1/p.T0) || ...
            ~ismember(p.p8ConditionalMeasure,{'linear','log'})
        error('Virtual patient %d has undocumented conditional mappings.',n)
    end
    if abs(p.Ts0+p.Tu0-p.T0)>1e-10*max(1,p.T0)
        error('Virtual patient %d has inconsistent initial tumor counts.',n)
    end
    if p.p8<getSpec(spec,'p8').lo || p.p8>getSpec(spec,'p8').hi
        error('Virtual patient %d violates the global p8 bounds.',n)
    end
    if p.Di0<getSpec(spec,'Di0').lo || p.Di0>getSpec(spec,'Di0').hi
        error('Virtual patient %d violates the Di(0) bounds.',n)
    end
end
fprintf('Virtual-patient validation passed for %d patients.\n',N)
end

function validatePerturbationPlans(plans,mechanisms,spec,settings)
if numel(plans)~=5
    error('Exactly five mechanism plans are required.')
end
for g = 1:5
    if ~strcmp(plans(g).mechanismID,mechanisms(g).id) || ...
            ~isequal(plans(g).inputNames,mechanisms(g).variableInputs)
        error('Perturbation-plan allocation mismatch for %s.',mechanisms(g).id)
    end
    if size(plans(g).values,1)~=settings.nPerturbations
        error('Incorrect number of perturbations for %s.',mechanisms(g).id)
    end
    for q = 1:numel(plans(g).inputNames)
        s = getSpec(spec,plans(g).inputNames{q});
        v = plans(g).values(:,q);
        if any(v<s.lo | v>s.hi)
            error('Perturbation values for %s exceed their bounds.',s.name)
        end
    end
end
fprintf('Mechanism perturbation-plan validation passed.\n')
end

%% ========================================================================
% PRIMARY SIMULATIONS
% ========================================================================

function simulation = runPrimarySimulations( ...
        patients,plans,settings,fixed,odeOpt)
N = numel(patients);
G = numel(plans);
S = settings.nPerturbations;

ReferenceFinalState = nan(N,7);
ReferenceOK = false(N,1);
ReferenceTsAbsorbed = false(N,1);
ReferenceTotalAbsorbed = false(N,1);
ReferenceTsAbsorptionTime = nan(N,1);
ReferenceTotalAbsorptionTime = nan(N,1);
ReferenceMessage = strings(N,1);

PerturbedFinalState = nan(N,G,S,7);
PerturbedOK = false(N,G,S);
PerturbedTsAbsorbed = false(N,G,S);
PerturbedTotalAbsorbed = false(N,G,S);
PerturbedTsAbsorptionTime = nan(N,G,S);
PerturbedTotalAbsorptionTime = nan(N,G,S);
PerturbedMessage = strings(N,G,S);
CompletedPatients = false(N,1);

checkpointFile = fullfile(settings.outputDir, ...
    'PRIMARY_SIMULATIONS_CHECKPOINT_MECHANISM_FINAL.mat');
if settings.resumeFromCheckpoint && exist(checkpointFile,'file')==2
    checkpoint = load(checkpointFile);
    requiredFields = {'CheckpointCodeVersion','CheckpointPatients', ...
        'CheckpointPlans','CheckpointSettings','CompletedPatients', ...
        'ReferenceFinalState','ReferenceOK','ReferenceTsAbsorbed', ...
        'ReferenceTotalAbsorbed','ReferenceTsAbsorptionTime', ...
        'ReferenceTotalAbsorptionTime','ReferenceMessage', ...
        'PerturbedFinalState','PerturbedOK','PerturbedTsAbsorbed', ...
        'PerturbedTotalAbsorbed','PerturbedTsAbsorptionTime', ...
        'PerturbedTotalAbsorptionTime','PerturbedMessage'};
    if ~all(isfield(checkpoint,requiredFields))
        error(['The mechanism-analysis checkpoint is incomplete and will ' ...
            'not be used. ' ...
            'Move it out of the output directory before restarting.'])
    end
    if ~strcmp(checkpoint.CheckpointCodeVersion,settings.codeVersion) || ...
            ~isequaln(checkpoint.CheckpointPatients,patients) || ...
            ~isequaln(checkpoint.CheckpointPlans,plans) || ...
            ~isequaln(checkpoint.CheckpointSettings, ...
            checkpointScientificSettings(settings))
        error(['The existing checkpoint is incompatible with the current ' ...
            'equations, bounds, sampling designs, or numerical rules. It ' ...
            'was not loaded; move it out of the output directory.'])
    end

    CompletedPatients = checkpoint.CompletedPatients;
    ReferenceFinalState = checkpoint.ReferenceFinalState;
    ReferenceOK = checkpoint.ReferenceOK;
    ReferenceTsAbsorbed = checkpoint.ReferenceTsAbsorbed;
    ReferenceTotalAbsorbed = checkpoint.ReferenceTotalAbsorbed;
    ReferenceTsAbsorptionTime = checkpoint.ReferenceTsAbsorptionTime;
    ReferenceTotalAbsorptionTime = ...
        checkpoint.ReferenceTotalAbsorptionTime;
    ReferenceMessage = checkpoint.ReferenceMessage;
    PerturbedFinalState = checkpoint.PerturbedFinalState;
    PerturbedOK = checkpoint.PerturbedOK;
    PerturbedTsAbsorbed = checkpoint.PerturbedTsAbsorbed;
    PerturbedTotalAbsorbed = checkpoint.PerturbedTotalAbsorbed;
    PerturbedTsAbsorptionTime = checkpoint.PerturbedTsAbsorptionTime;
    PerturbedTotalAbsorptionTime = ...
        checkpoint.PerturbedTotalAbsorptionTime;
    PerturbedMessage = checkpoint.PerturbedMessage;
    fprintf('Resuming compatible checkpoint: %d/%d patients complete.\n', ...
        sum(CompletedPatients),N)
end

fprintf('\nRunning %d reference and %d perturbed simulations...\n', ...
    N,N*G*S)
progressStart = tic;

remainingPatients = find(~CompletedPatients);
blockSize = max(1,settings.checkpointBlockSize);
for blockStart = 1:blockSize:numel(remainingPatients)
    block = remainingPatients(blockStart: ...
        min(blockStart+blockSize-1,numel(remainingPatients)));
    blockOutput = cell(numel(block),1);

    if settings.useParallel
        parfor q = 1:numel(block)
            n = block(q);
            blockOutput{q} = runOnePatient( ...
                patients{n},plans,settings,fixed,odeOpt);
        end
    else
        for q = 1:numel(block)
            n = block(q);
            blockOutput{q} = runOnePatient( ...
                patients{n},plans,settings,fixed,odeOpt);
        end
    end

    for q = 1:numel(block)
        n = block(q);
        out = blockOutput{q};
        ReferenceFinalState(n,:) = out.referenceFinalState;
        ReferenceOK(n) = out.referenceOK;
        ReferenceTsAbsorbed(n) = out.referenceTsAbsorbed;
        ReferenceTotalAbsorbed(n) = out.referenceTotalAbsorbed;
        ReferenceTsAbsorptionTime(n) = out.referenceTsAbsorptionTime;
        ReferenceTotalAbsorptionTime(n) = out.referenceTotalAbsorptionTime;
        ReferenceMessage(n) = out.referenceMessage;

        PerturbedFinalState(n,:,:,:) = out.perturbedFinalState;
        PerturbedOK(n,:,:) = out.perturbedOK;
        PerturbedTsAbsorbed(n,:,:) = out.perturbedTsAbsorbed;
        PerturbedTotalAbsorbed(n,:,:) = out.perturbedTotalAbsorbed;
        PerturbedTsAbsorptionTime(n,:,:) = out.perturbedTsAbsorptionTime;
        PerturbedTotalAbsorptionTime(n,:,:) = ...
            out.perturbedTotalAbsorptionTime;
        PerturbedMessage(n,:,:) = out.perturbedMessage;
        CompletedPatients(n) = true;
    end

    CheckpointCodeVersion = settings.codeVersion; %#ok<NASGU>
    CheckpointPatients = patients; %#ok<NASGU>
    CheckpointPlans = plans; %#ok<NASGU>
    CheckpointSettings = checkpointScientificSettings(settings); %#ok<NASGU>
    save(checkpointFile,'CheckpointCodeVersion','CheckpointPatients', ...
        'CheckpointPlans','CheckpointSettings','CompletedPatients', ...
        'ReferenceFinalState','ReferenceOK','ReferenceTsAbsorbed', ...
        'ReferenceTotalAbsorbed','ReferenceTsAbsorptionTime', ...
        'ReferenceTotalAbsorptionTime','ReferenceMessage', ...
        'PerturbedFinalState','PerturbedOK','PerturbedTsAbsorbed', ...
        'PerturbedTotalAbsorbed','PerturbedTsAbsorptionTime', ...
        'PerturbedTotalAbsorptionTime','PerturbedMessage','-v7.3')

    printProgress('patients',sum(CompletedPatients),N,toc(progressStart))
end

ReferenceFinalTumor = ReferenceFinalState(:,2)+ReferenceFinalState(:,3);
PerturbedFinalTumor = PerturbedFinalState(:,:,:,2) + ...
    PerturbedFinalState(:,:,:,3);
ReferenceControlled = ReferenceOK & ...
    ReferenceFinalTumor<settings.detectionLimit;

simulation = struct();
simulation.ReferenceFinalState = ReferenceFinalState;
simulation.ReferenceFinalTumor = ReferenceFinalTumor;
simulation.ReferenceOK = ReferenceOK;
simulation.ReferenceControlled = ReferenceControlled;
simulation.ReferenceTsAbsorbed = ReferenceTsAbsorbed;
simulation.ReferenceTotalAbsorbed = ReferenceTotalAbsorbed;
simulation.ReferenceTsAbsorptionTime = ReferenceTsAbsorptionTime;
simulation.ReferenceTotalAbsorptionTime = ReferenceTotalAbsorptionTime;
simulation.ReferenceMessage = ReferenceMessage;
simulation.PerturbedFinalState = PerturbedFinalState;
simulation.PerturbedFinalTumor = PerturbedFinalTumor;
simulation.PerturbedOK = PerturbedOK;
simulation.PerturbedTsAbsorbed = PerturbedTsAbsorbed;
simulation.PerturbedTotalAbsorbed = PerturbedTotalAbsorbed;
simulation.PerturbedTsAbsorptionTime = PerturbedTsAbsorptionTime;
simulation.PerturbedTotalAbsorptionTime = PerturbedTotalAbsorptionTime;
simulation.PerturbedMessage = PerturbedMessage;

fprintf('Primary simulations completed in %.1f min.\n',toc(progressStart)/60)
end

function core = checkpointScientificSettings(settings)
fields = {'codeVersion','runMode','nPatients','nPerturbations', ...
    'tEnd','doseDuration','treatmentTimes','oneCellLimit', ...
    'detectionLimit','AbsTol','RelTol','MaxStep','seedCohort', ...
    'seedPerturbations'};
core = struct();
for j = 1:numel(fields)
    core.(fields{j}) = settings.(fields{j});
end
end

function out = runOnePatient(patient,plans,settings,fixed,odeOpt)
G = numel(plans);
S = settings.nPerturbations;

out.referenceFinalState = nan(1,7);
out.referenceOK = false;
out.referenceTsAbsorbed = false;
out.referenceTotalAbsorbed = false;
out.referenceTsAbsorptionTime = nan;
out.referenceTotalAbsorptionTime = nan;
out.referenceMessage = "not run";
out.perturbedFinalState = nan(1,G,S,7);
out.perturbedOK = false(1,G,S);
out.perturbedTsAbsorbed = false(1,G,S);
out.perturbedTotalAbsorbed = false(1,G,S);
out.perturbedTsAbsorptionTime = nan(1,G,S);
out.perturbedTotalAbsorptionTime = nan(1,G,S);
out.perturbedMessage = strings(1,G,S);

ref = simulateEndpoint(patient,settings,fixed,odeOpt);
out.referenceFinalState = ref.finalState.';
out.referenceOK = ref.ok;
out.referenceTsAbsorbed = ref.TsAbsorbed;
out.referenceTotalAbsorbed = ref.totalAbsorbed;
out.referenceTsAbsorptionTime = ref.TsAbsorptionTime;
out.referenceTotalAbsorptionTime = ref.totalAbsorptionTime;
out.referenceMessage = ref.message;

for g = 1:G
    for s = 1:S
        perturbedPatient = applyPerturbation(patient,plans(g),s);
        sim = simulateEndpoint(perturbedPatient,settings,fixed,odeOpt);
        out.perturbedFinalState(1,g,s,:) = ...
            reshape(sim.finalState,1,1,1,7);
        out.perturbedOK(1,g,s) = sim.ok;
        out.perturbedTsAbsorbed(1,g,s) = sim.TsAbsorbed;
        out.perturbedTotalAbsorbed(1,g,s) = sim.totalAbsorbed;
        out.perturbedTsAbsorptionTime(1,g,s) = sim.TsAbsorptionTime;
        out.perturbedTotalAbsorptionTime(1,g,s) = ...
            sim.totalAbsorptionTime;
        out.perturbedMessage(1,g,s) = sim.message;
    end
end
end

function patientOut = applyPerturbation(patientIn,plan,sampleIndex)
patientOut = patientIn;
for q = 1:numel(plan.inputNames)
    patientOut.(plan.inputNames{q}) = plan.values(sampleIndex,q);
end
end

function result = simulateEndpoint(patient,settings,fixed,odeOpt)
y0 = [fixed.M0;patient.Ts0;patient.Tu0;patient.Di0; ...
    fixed.Dm0;fixed.E0;fixed.R0];

result = simulateScheduleWithAbsorption( ...
    y0,patient,settings,odeOpt,false);
end

function result = simulateScheduleWithAbsorption( ...
        y0,pars,settings,odeOpt,collectTrajectory)
if nargin<5
    collectTrajectory = false;
end
starts = settings.treatmentTimes(settings.treatmentTimes<settings.tEnd);
ends = min(starts+settings.doseDuration,settings.tEnd);
breakPoints = unique([0 settings.tEnd starts ends]);
breakPoints = sort(breakPoints);

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

% The cohort generator prevents sub-one initial tumor compartments. These
% guards stop silently inconsistent externally supplied initial states.
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

% At most two biological boundary changes are possible: Ts absorption and
% total-tumor absorption. The third call completes the remaining interval.
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
    sameTimeTolerance = 1e-10*max(1,abs(eventTime));
    eventsAtStop = ie(abs(te-eventTime)<=sameTimeTolerance);

    % Total-tumor extinction has priority when both events are simultaneous.
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
        % Retain the pre-absorption event state (one cell) for accurate
        % interpolation up to the crossing. The trajectory consumer applies
        % the recorded absorbing state at and after the event time.
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
t = [];
Y = [];
te = [];
ye = [];
ie = [];
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

%% ========================================================================
% SEVEN-STATE MODEL FROM THE MANUSCRIPT FORMULATION
% ========================================================================

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

function printProgress(labelText,completed,total,elapsedSeconds)
fraction = completed/max(total,1);
if completed>0
    remainingSeconds = elapsedSeconds*(1-fraction)/fraction;
else
    remainingSeconds = nan;
end
fprintf('%s: %d/%d (%.1f%%), elapsed %.1f min, estimated remaining %.1f min\n', ...
    labelText,completed,total,100*fraction,elapsedSeconds/60, ...
    remainingSeconds/60)
end

%% ========================================================================
% TRANSITION SUMMARIES AND LONG-FORM OUTPUTS
% ========================================================================

function [summaryTable,patientTable,longTable,overlapTable] = ...
        summarizeMechanismTransitions( ...
        simulation,patients,mechanisms,plans,settings)
N = numel(patients);
G = numel(mechanisms);
S = settings.nPerturbations;
eligible = simulation.ReferenceControlled(:);
nEligible = sum(eligible);

MechanismID = string({mechanisms.id}).';
Mechanism = string({mechanisms.title}).';
NumberVaried = zeros(G,1);
EligibleReferencePatients = repmat(nEligible,G,1);
ValidPerturbations = zeros(G,1);
TransitionEvents = zeros(G,1);
TransitionPercent = nan(G,1);
AffectedPatients = zeros(G,1);
AffectedPatientPercent = nan(G,1);
ShareOfAllTransitionEventsPercent = nan(G,1);
ValidNoAbsorptionDraws = zeros(G,1);
ValidTsOnlyAbsorptionDraws = zeros(G,1);
ValidTotalExtinctionDraws = zeros(G,1);

affectedMatrix = false(N,G);

for g = 1:G
    NumberVaried(g) = numel(plans(g).inputNames);
    ok = reshape(simulation.PerturbedOK(:,g,:),N,S);
    tumor = reshape(simulation.PerturbedFinalTumor(:,g,:),N,S);
    valid = eligible & ok;
    transitions = valid & tumor>=settings.detectionLimit;
    affected = any(transitions,2);
    affectedMatrix(:,g) = affected;

    ValidPerturbations(g) = sum(valid(:));
    TransitionEvents(g) = sum(transitions(:));
    if ValidPerturbations(g)>0
        TransitionPercent(g) = ...
            100*TransitionEvents(g)/ValidPerturbations(g);
    end
    AffectedPatients(g) = sum(affected);
    if nEligible>0
        AffectedPatientPercent(g) = 100*AffectedPatients(g)/nEligible;
    end

    tsAbs = reshape(simulation.PerturbedTsAbsorbed(:,g,:),N,S);
    totAbs = reshape(simulation.PerturbedTotalAbsorbed(:,g,:),N,S);
    ValidNoAbsorptionDraws(g) = ...
        sum(valid(:) & ~tsAbs(:) & ~totAbs(:));
    ValidTsOnlyAbsorptionDraws(g) = ...
        sum(valid(:) & tsAbs(:) & ~totAbs(:));
    ValidTotalExtinctionDraws(g) = sum(valid(:) & totAbs(:));
    if ValidNoAbsorptionDraws(g)+ValidTsOnlyAbsorptionDraws(g)+ ...
            ValidTotalExtinctionDraws(g)~=ValidPerturbations(g)
        error('Absorption categories do not reconcile for %s.', ...
            mechanisms(g).id)
    end
end

totalTransitions = sum(TransitionEvents);
if totalTransitions>0
    ShareOfAllTransitionEventsPercent = ...
        100*TransitionEvents/totalTransitions;
end

summaryTable = table(MechanismID,Mechanism,NumberVaried, ...
    EligibleReferencePatients,ValidPerturbations,TransitionEvents, ...
    TransitionPercent,AffectedPatients,AffectedPatientPercent, ...
    ShareOfAllTransitionEventsPercent,ValidNoAbsorptionDraws, ...
    ValidTsOnlyAbsorptionDraws,ValidTotalExtinctionDraws);

% Patient-by-mechanism table. Every eligible patient remains visible even
% when no transition occurs.
nPatientRows = N*G;
VirtualPatient = repelem((1:N).',G);
MechanismIDP = repmat(MechanismID,N,1);
ReferenceValid = repelem(simulation.ReferenceOK(:),G);
ReferenceControlled = repelem(eligible,G);
ReferenceTumorCells = repelem(simulation.ReferenceFinalTumor(:),G);
ValidDraws = zeros(nPatientRows,1);
TransitionDraws = zeros(nPatientRows,1);
TransitionPercentWithinPatient = nan(nPatientRows,1);
AnyTransition = false(nPatientRows,1);
NoAbsorptionDraws = zeros(nPatientRows,1);
TsOnlyAbsorptionDraws = zeros(nPatientRows,1);
TotalExtinctionDraws = zeros(nPatientRows,1);

row = 0;
for n = 1:N
    for g = 1:G
        row = row+1;
        ok = reshape(simulation.PerturbedOK(n,g,:),S,1);
        tumor = reshape(simulation.PerturbedFinalTumor(n,g,:),S,1);
        valid = eligible(n) & ok;
        transition = valid & tumor>=settings.detectionLimit;
        ValidDraws(row) = sum(valid);
        TransitionDraws(row) = sum(transition);
        if ValidDraws(row)>0
            TransitionPercentWithinPatient(row) = ...
                100*TransitionDraws(row)/ValidDraws(row);
        end
        AnyTransition(row) = any(transition);
        tsAbs = reshape(simulation.PerturbedTsAbsorbed(n,g,:),S,1);
        totalAbs = reshape(simulation.PerturbedTotalAbsorbed(n,g,:),S,1);
        NoAbsorptionDraws(row) = sum(valid & ~tsAbs & ~totalAbs);
        TsOnlyAbsorptionDraws(row) = sum(valid & tsAbs & ~totalAbs);
        TotalExtinctionDraws(row) = sum(valid & totalAbs);
    end
end

patientTable = table(VirtualPatient,MechanismIDP,ReferenceValid, ...
    ReferenceControlled,ReferenceTumorCells,ValidDraws,TransitionDraws, ...
    TransitionPercentWithinPatient,AnyTransition,NoAbsorptionDraws, ...
    TsOnlyAbsorptionDraws,TotalExtinctionDraws, ...
    'VariableNames',{'VirtualPatient','MechanismID','ReferenceValid', ...
    'ReferenceControlled','ReferenceTumorCells','ValidDraws', ...
    'TransitionDraws','TransitionPercentWithinPatient','AnyTransition', ...
    'NoAbsorptionDraws','TsOnlyAbsorptionDraws', ...
    'TotalExtinctionDraws'});

% Long-form table for descriptive summaries, with one row per patient,
% mechanism and perturbation draw.
nRows = N*G*S;
VirtualPatientL = zeros(nRows,1);
MechanismIDL = strings(nRows,1);
PerturbationDraw = zeros(nRows,1);
ReferenceValidL = false(nRows,1);
ReferenceControlledL = false(nRows,1);
ReferenceTumorCellsL = nan(nRows,1);
PerturbedValid = false(nRows,1);
PerturbedTsCells = nan(nRows,1);
PerturbedTuCells = nan(nRows,1);
PerturbedTumorCells = nan(nRows,1);
ControlToFailure = false(nRows,1);
TsAbsorbed = false(nRows,1);
TotalTumorAbsorbed = false(nRows,1);
TsAbsorptionTime = nan(nRows,1);
TotalAbsorptionTime = nan(nRows,1);
SolverMessage = strings(nRows,1);

row = 0;
for n = 1:N
    for g = 1:G
        for s = 1:S
            row = row+1;
            VirtualPatientL(row) = n;
            MechanismIDL(row) = mechanisms(g).id;
            PerturbationDraw(row) = s;
            ReferenceValidL(row) = simulation.ReferenceOK(n);
            ReferenceControlledL(row) = eligible(n);
            ReferenceTumorCellsL(row) = simulation.ReferenceFinalTumor(n);
            PerturbedValid(row) = simulation.PerturbedOK(n,g,s);
            PerturbedTsCells(row) = ...
                simulation.PerturbedFinalState(n,g,s,2);
            PerturbedTuCells(row) = ...
                simulation.PerturbedFinalState(n,g,s,3);
            PerturbedTumorCells(row) = ...
                simulation.PerturbedFinalTumor(n,g,s);
            ControlToFailure(row) = eligible(n) && ...
                PerturbedValid(row) && ...
                PerturbedTumorCells(row)>=settings.detectionLimit;
            TsAbsorbed(row) = ...
                simulation.PerturbedTsAbsorbed(n,g,s);
            TotalTumorAbsorbed(row) = ...
                simulation.PerturbedTotalAbsorbed(n,g,s);
            TsAbsorptionTime(row) = ...
                simulation.PerturbedTsAbsorptionTime(n,g,s);
            TotalAbsorptionTime(row) = ...
                simulation.PerturbedTotalAbsorptionTime(n,g,s);
            SolverMessage(row) = simulation.PerturbedMessage(n,g,s);
        end
    end
end

longTable = table(VirtualPatientL,MechanismIDL,PerturbationDraw, ...
    ReferenceValidL,ReferenceControlledL,ReferenceTumorCellsL, ...
    PerturbedValid,PerturbedTsCells,PerturbedTuCells, ...
    PerturbedTumorCells,ControlToFailure,TsAbsorbed, ...
    TotalTumorAbsorbed,TsAbsorptionTime,TotalAbsorptionTime, ...
    SolverMessage, ...
    'VariableNames',{'VirtualPatient','MechanismID','PerturbationDraw', ...
    'ReferenceValid','ReferenceControlled','ReferenceTumorCells', ...
    'PerturbedValid','PerturbedTsCells','PerturbedTuCells', ...
    'PerturbedTumorCells','ControlToFailure','TsAbsorbed', ...
    'TotalTumorAbsorbed','TsAbsorptionTime','TotalAbsorptionTime', ...
    'SolverMessage'});

% Wide patient overlap table in prespecified M1-M5 order.
overlapTable = table((1:N).',simulation.ReferenceOK(:),eligible, ...
    simulation.ReferenceFinalTumor(:),affectedMatrix(:,1), ...
    affectedMatrix(:,2),affectedMatrix(:,3),affectedMatrix(:,4), ...
    affectedMatrix(:,5),sum(affectedMatrix,2), ...
    'VariableNames',{'VirtualPatient','ReferenceValid', ...
    'ReferenceControlled','ReferenceTumorCells','AffectedByM1', ...
    'AffectedByM2','AffectedByM3','AffectedByM4','AffectedByM5', ...
    'NumberOfMechanismsWithTransition'});
end

%% ========================================================================
% PATIENT-CLUSTERED NUMERICAL RESAMPLING
% ========================================================================

function bootstrapTable = patientClusterBootstrap( ...
        simulation,mechanisms,settings)
eligible = find(simulation.ReferenceControlled);
G = numel(mechanisms);
B = settings.nBootstrap;

MechanismID = string({mechanisms.id}).';
EventRateLower95 = nan(G,1);
EventRateUpper95 = nan(G,1);
AffectedPatientRateLower95 = nan(G,1);
AffectedPatientRateUpper95 = nan(G,1);
BootstrapReplicates = repmat(B,G,1);

if isempty(eligible) || B<1
    bootstrapTable = table(MechanismID,BootstrapReplicates, ...
        EventRateLower95,EventRateUpper95, ...
        AffectedPatientRateLower95,AffectedPatientRateUpper95);
    return
end

stream = RandStream('mt19937ar','Seed',settings.seedBootstrap);
eventRate = nan(B,G);
patientRate = nan(B,G);

fprintf('Running %d patient-cluster bootstrap replicates...\n',B)
for b = 1:B
    sampled = eligible(randi(stream,numel(eligible),numel(eligible),1));
    for g = 1:G
        ok = reshape(simulation.PerturbedOK(sampled,g,:), ...
            numel(sampled),settings.nPerturbations);
        tumor = reshape(simulation.PerturbedFinalTumor(sampled,g,:), ...
            numel(sampled),settings.nPerturbations);
        valid = ok & isfinite(tumor);
        transition = valid & tumor>=settings.detectionLimit;
        eventRate(b,g) = 100*sum(transition(:))/max(sum(valid(:)),1);
        patientRate(b,g) = 100*sum(any(transition,2))/numel(sampled);
    end
end

for g = 1:G
    EventRateLower95(g) = localPercentile(eventRate(:,g),2.5);
    EventRateUpper95(g) = localPercentile(eventRate(:,g),97.5);
    AffectedPatientRateLower95(g) = ...
        localPercentile(patientRate(:,g),2.5);
    AffectedPatientRateUpper95(g) = ...
        localPercentile(patientRate(:,g),97.5);
end

bootstrapTable = table(MechanismID,BootstrapReplicates, ...
    EventRateLower95,EventRateUpper95, ...
    AffectedPatientRateLower95,AffectedPatientRateUpper95);
end

function summaryTable = attachBootstrapIntervals(summaryTable,bootstrapTable)
summaryTable.EventRateLower95 = nan(height(summaryTable),1);
summaryTable.EventRateUpper95 = nan(height(summaryTable),1);
summaryTable.AffectedPatientRateLower95 = nan(height(summaryTable),1);
summaryTable.AffectedPatientRateUpper95 = nan(height(summaryTable),1);

for g = 1:height(summaryTable)
    idx = find(bootstrapTable.MechanismID==summaryTable.MechanismID(g), ...
        1,'first');
    if isempty(idx)
        continue
    end
    summaryTable.EventRateLower95(g) = ...
        bootstrapTable.EventRateLower95(idx);
    summaryTable.EventRateUpper95(g) = ...
        bootstrapTable.EventRateUpper95(idx);
    summaryTable.AffectedPatientRateLower95(g) = ...
        bootstrapTable.AffectedPatientRateLower95(idx);
    summaryTable.AffectedPatientRateUpper95(g) = ...
        bootstrapTable.AffectedPatientRateUpper95(idx);
end
end

function stabilityTable = cohortSizeStability( ...
        simulation,mechanisms,settings)
eligibleMask = simulation.ReferenceControlled(:);
NTotal = numel(eligibleMask);
G = numel(mechanisms);
sizes = settings.cohortStabilitySizes;
sizes = sizes(sizes<=NTotal);
R = settings.nCohortStabilityReplicates;

NVirtualPatients = zeros(0,1);
MechanismID = strings(0,1);
Replicates = zeros(0,1);
MedianEligiblePatients = zeros(0,1);
MedianTransitionPercent = zeros(0,1);
Lower95 = zeros(0,1);
Upper95 = zeros(0,1);

if ~any(eligibleMask) || isempty(sizes) || R<1
    stabilityTable = table(NVirtualPatients,MechanismID,Replicates, ...
        MedianEligiblePatients,MedianTransitionPercent,Lower95,Upper95);
    return
end

stream = RandStream('mt19937ar','Seed',settings.seedCohortStability);
for nTarget = sizes
    estimates = nan(R,G);
    eligibleCounts = zeros(R,1);
    for r = 1:R
        subset = randperm(stream,NTotal,nTarget);
        eligible = subset(eligibleMask(subset));
        eligibleCounts(r) = numel(eligible);
        if isempty(eligible)
            continue
        end
        for g = 1:G
            ok = reshape(simulation.PerturbedOK(eligible,g,:), ...
                numel(eligible),settings.nPerturbations);
            tumor = reshape(simulation.PerturbedFinalTumor(eligible,g,:), ...
                numel(eligible),settings.nPerturbations);
            valid = ok & isfinite(tumor);
            transition = valid & tumor>=settings.detectionLimit;
            estimates(r,g) = 100*sum(transition(:))/max(sum(valid(:)),1);
        end
    end
    for g = 1:G
        NVirtualPatients(end+1,1) = nTarget; %#ok<AGROW>
        MechanismID(end+1,1) = mechanisms(g).id; %#ok<AGROW>
        Replicates(end+1,1) = R; %#ok<AGROW>
        MedianEligiblePatients(end+1,1) = ...
            localPercentile(eligibleCounts,50); %#ok<AGROW>
        MedianTransitionPercent(end+1,1) = ...
            localPercentile(estimates(:,g),50); %#ok<AGROW>
        Lower95(end+1,1) = ...
            localPercentile(estimates(:,g),2.5); %#ok<AGROW>
        Upper95(end+1,1) = ...
            localPercentile(estimates(:,g),97.5); %#ok<AGROW>
    end
end

stabilityTable = table(NVirtualPatients,MechanismID,Replicates, ...
    MedianEligiblePatients,MedianTransitionPercent,Lower95,Upper95);
end

function q = localPercentile(x,p)
x = sort(x(isfinite(x)));
if isempty(x)
    q = nan;
    return
end
if numel(x)==1
    q = x;
    return
end
position = 1+(numel(x)-1)*(p/100);
lo = floor(position);
hi = ceil(position);
if lo==hi
    q = x(lo);
else
    q = x(lo)+(position-lo)*(x(hi)-x(lo));
end
end

%% ========================================================================
% AUDIT TABLES
% ========================================================================

function writeReproducibilityManifest(settings,spec,fixed,mechanisms)
pathValue = fullfile(settings.outputDir, ...
    'Mechanism_analysis_reproducibility_manifest.txt');
fid = fopen(pathValue,'w');
if fid<0
    error('Could not create the reproducibility manifest.')
end
fileCleanup = onCleanup(@()fclose(fid));
try
    fprintf(fid,'MMC mechanism-of-failure analysis\n');
    fprintf(fid,'Code version: %s\n',settings.codeVersion);
    fprintf(fid,'Run mode: %s\n',settings.runMode);
    fprintf(fid,'New section-specific virtual patients: %d\n', ...
        settings.nPatients);
    fprintf(fid,'LHS strata per cohort input: %d\n',settings.nPatients);
    fprintf(fid,'Joint perturbations per mechanism: %d\n', ...
        settings.nPerturbations);
    fprintf(fid,['Perturbations: absolute full-range joint LHS values; ' ...
        'not local percentage changes. Each mechanism-specific design is ' ...
        'applied unchanged to every virtual patient.\n']);
    fprintf(fid,'Cohort seed: %d\n',settings.seedCohort);
    fprintf(fid,'Perturbation base seed: %d (M1--M5 use %d--%d)\n', ...
        settings.seedPerturbations,settings.seedPerturbations, ...
        settings.seedPerturbations+4);
    fprintf(fid,'Patient-cluster bootstrap seed: %d\n', ...
        settings.seedBootstrap);
    fprintf(fid,['Bootstrap: virtual patients are resampled as clusters; ' ...
        'intervals are numerical cohort-resampling summaries, not ' ...
        'biological or clinical confidence intervals.\n']);
    fprintf(fid,'Cohort-stability seed: %d\n', ...
        settings.seedCohortStability);
    fprintf(fid,'Endpoint: day %g\n',settings.tEnd);
    fprintf(fid,'Detection limit: %.15g modeled cells\n', ...
        settings.detectionLimit);
    fprintf(fid,['One-cell rules: mechanism section only; exactly one cell ' ...
        'is retained unless the solution crosses downward into the sub-one ' ...
        'region and may regrow; Tu has no separate compartment-level ' ...
        'absorbing boundary.\n']);
    fprintf(fid,['Initial composition: T(0) is log-mapped on ' ...
        '[2,min(2.88737e11,k)); p8 is mapped on the patient-specific ' ...
        'intersection [4e-6,0.798] with Ts(0)>=1 and Tu(0)>=1, using ' ...
        'the prespecified range-width scale rule.\n']);
    fprintf(fid,['Tumor equations: 1-F(M) retained in both proliferation ' ...
        'terms; no MMC factor in dE/dt; Ts-to-Tu only.\n']);
    fprintf(fid,'ode15s: AbsTol %.1e, RelTol %.1e, MaxStep %g day\n', ...
        settings.AbsTol,settings.RelTol,settings.MaxStep);
    fprintf(fid,'Treatment starts:');
    fprintf(fid,' %g',settings.treatmentTimes);
    fprintf(fid,' days; duration %.15g day.\n',settings.doseDuration);
    fprintf(fid,['Integration boundaries: MMC-session starts and ends ' ...
        'only; no monthly solver boundaries are required for the day-360 ' ...
        'mechanism endpoint. The 3D case is evaluated on its declared ' ...
        'time grid.\n']);
    fprintf(fid,'Ranged inputs: %d; fixed M0=%g, Dm0=%g, E0=%g, R0=%g.\n', ...
        numel(spec),fixed.M0,fixed.Dm0,fixed.E0,fixed.R0);
    fprintf(fid,'Mechanism allocations:\n');
    for g = 1:numel(mechanisms)
        fprintf(fid,'  %s varied: %s',mechanisms(g).id, ...
            strjoin(mechanisms(g).variableInputs,', '));
        if ~isempty(mechanisms(g).fixedInputs)
            fprintf(fid,'; biologically assigned but fixed: %s', ...
                strjoin(mechanisms(g).fixedInputs,', '));
        end
        fprintf(fid,'\n');
    end
    clear fileCleanup
catch ME
    rethrow(ME)
end
end

function T = inputSpecificationTable(spec,fixed)
Input = strings(0,1);
Label = strings(0,1);
LowerBound = zeros(0,1);
UpperBound = zeros(0,1);
SamplingMeasure = strings(0,1);
Unit = strings(0,1);
Role = strings(0,1);
Status = strings(0,1);

for j = 1:numel(spec)
    Input(end+1,1) = spec(j).name; %#ok<AGROW>
    Label(end+1,1) = spec(j).label; %#ok<AGROW>
    LowerBound(end+1,1) = spec(j).lo; %#ok<AGROW>
    UpperBound(end+1,1) = spec(j).hi; %#ok<AGROW>
    SamplingMeasure(end+1,1) = spec(j).measure; %#ok<AGROW>
    Unit(end+1,1) = spec(j).unit; %#ok<AGROW>
    Role(end+1,1) = spec(j).role; %#ok<AGROW>
    Status(end+1,1) = "ranged"; %#ok<AGROW>
end

fixedNames = {'d0','p6','p7','b','tau','M0','Dm0','E0','R0'};
fixedLabels = {'d_0','p_6','p_7','b','\tau','M(0)','D_m(0)', ...
    'E(0)','R(0)'};
fixedUnits = {'cells day^{-1}','cell^{-1} day^{-1}', ...
    'day^{-1}','cells','day','micromolar','cells','cells','cells'};
fixedRoles = {'source','M5 fixed','M5 fixed','immune scaling', ...
    'MMC-session duration','initial condition','initial condition', ...
    'initial condition','initial condition'};

for j = 1:numel(fixedNames)
    Input(end+1,1) = fixedNames{j}; %#ok<AGROW>
    Label(end+1,1) = fixedLabels{j}; %#ok<AGROW>
    LowerBound(end+1,1) = fixed.(fixedNames{j}); %#ok<AGROW>
    UpperBound(end+1,1) = fixed.(fixedNames{j}); %#ok<AGROW>
    SamplingMeasure(end+1,1) = "fixed"; %#ok<AGROW>
    Unit(end+1,1) = fixedUnits{j}; %#ok<AGROW>
    Role(end+1,1) = fixedRoles{j}; %#ok<AGROW>
    Status(end+1,1) = "fixed"; %#ok<AGROW>
end

T = table(Input,Label,LowerBound,UpperBound,SamplingMeasure, ...
    Unit,Role,Status);
end

function T = mechanismDefinitionTable(mechanisms)
MechanismID = string({mechanisms.id}).';
Mechanism = string({mechanisms.title}).';
AssignedParameters = strings(numel(mechanisms),1);
VariedParameters = strings(numel(mechanisms),1);
FixedAssignedParameters = strings(numel(mechanisms),1);
NumberAssigned = zeros(numel(mechanisms),1);
NumberVaried = zeros(numel(mechanisms),1);
Interpretation = string({mechanisms.interpretation}).';

for g = 1:numel(mechanisms)
    AssignedParameters(g) = strjoin(string(mechanisms(g).assignedInputs),', ');
    VariedParameters(g) = strjoin(string(mechanisms(g).variableInputs),', ');
    FixedAssignedParameters(g) = ...
        strjoin(string(mechanisms(g).fixedInputs),', ');
    NumberAssigned(g) = numel(mechanisms(g).assignedInputs);
    NumberVaried(g) = numel(mechanisms(g).variableInputs);
end

T = table(MechanismID,Mechanism,AssignedParameters,VariedParameters, ...
    FixedAssignedParameters,NumberAssigned,NumberVaried,Interpretation);
end

function T = virtualPatientTable(patients,U,spec)
N = numel(patients);
VirtualPatient = (1:N).';
T = table(VirtualPatient);

for j = 1:numel(spec)
    name = spec(j).name;
    values = zeros(N,1);
    for n = 1:N
        values(n) = patients{n}.(name);
    end
    T.(matlab.lang.makeValidName(name)) = values;
    T.([matlab.lang.makeValidName(name) '_unit']) = U(:,j);
end

T.Ts0 = cellfun(@(p)p.Ts0,patients);
T.Tu0 = cellfun(@(p)p.Tu0,patients);
T.T0ConditionalUpper = cellfun(@(p)p.T0ConditionalUpper,patients);
T.p8ConditionalLower = cellfun(@(p)p.p8ConditionalLower,patients);
T.p8ConditionalUpper = cellfun(@(p)p.p8ConditionalUpper,patients);
T.p8ConditionalMeasure = string(cellfun( ...
    @(p)p.p8ConditionalMeasure,patients,'UniformOutput',false));
end

function T = perturbationPlanTable(plans,mechanisms)
MechanismID = strings(0,1);
PerturbationDraw = zeros(0,1);
Parameter = strings(0,1);
UnitCoordinate = zeros(0,1);
MappedValue = zeros(0,1);
SamplingMeasure = strings(0,1);
MechanismDimension = zeros(0,1);

for g = 1:numel(plans)
    for s = 1:size(plans(g).values,1)
        for q = 1:numel(plans(g).inputNames)
            MechanismID(end+1,1) = mechanisms(g).id; %#ok<AGROW>
            PerturbationDraw(end+1,1) = s; %#ok<AGROW>
            Parameter(end+1,1) = plans(g).inputNames{q}; %#ok<AGROW>
            UnitCoordinate(end+1,1) = ...
                plans(g).unitCoordinates(s,q); %#ok<AGROW>
            MappedValue(end+1,1) = plans(g).values(s,q); %#ok<AGROW>
            SamplingMeasure(end+1,1) = plans(g).measures{q}; %#ok<AGROW>
            MechanismDimension(end+1,1) = ...
                numel(plans(g).inputNames); %#ok<AGROW>
        end
    end
end

T = table(MechanismID,PerturbationDraw,Parameter,UnitCoordinate, ...
    MappedValue,SamplingMeasure,MechanismDimension);
end

function T = buildSolverAuditTable(simulation)
ReferenceSimulations = numel(simulation.ReferenceOK);
ReferenceCompleted = sum(simulation.ReferenceOK);
ReferenceInvalid = ReferenceSimulations-ReferenceCompleted;
PerturbedSimulations = numel(simulation.PerturbedOK);
PerturbedCompleted = sum(simulation.PerturbedOK(:));
PerturbedInvalid = PerturbedSimulations-PerturbedCompleted;
ReferenceNoAbsorption = sum(simulation.ReferenceOK & ...
    ~simulation.ReferenceTsAbsorbed & ~simulation.ReferenceTotalAbsorbed);
ReferenceTsOnlyAbsorbed = sum(simulation.ReferenceOK & ...
    simulation.ReferenceTsAbsorbed & ~simulation.ReferenceTotalAbsorbed);
ReferenceTotalAbsorbed = sum(simulation.ReferenceTotalAbsorbed);
PerturbedNoAbsorption = sum(simulation.PerturbedOK(:) & ...
    ~simulation.PerturbedTsAbsorbed(:) & ...
    ~simulation.PerturbedTotalAbsorbed(:));
PerturbedTsOnlyAbsorbed = sum(simulation.PerturbedOK(:) & ...
    simulation.PerturbedTsAbsorbed(:) & ...
    ~simulation.PerturbedTotalAbsorbed(:));
PerturbedTotalAbsorbed = sum(simulation.PerturbedTotalAbsorbed(:));

T = table(ReferenceSimulations,ReferenceCompleted,ReferenceInvalid, ...
    PerturbedSimulations,PerturbedCompleted,PerturbedInvalid, ...
    ReferenceNoAbsorption,ReferenceTsOnlyAbsorbed, ...
    ReferenceTotalAbsorbed,PerturbedNoAbsorption, ...
    PerturbedTsOnlyAbsorbed,PerturbedTotalAbsorbed);

if ReferenceNoAbsorption+ReferenceTsOnlyAbsorbed+ ...
        ReferenceTotalAbsorbed~=ReferenceCompleted || ...
        PerturbedNoAbsorption+PerturbedTsOnlyAbsorbed+ ...
        PerturbedTotalAbsorbed~=PerturbedCompleted
    error('Mutually exclusive absorption-audit categories do not reconcile.')
end
end

%% ========================================================================
% REPRESENTATIVE PATIENT CASE SELECTION
% ========================================================================

function [Candidate,Selected] = selectRepresentativeCases( ...
        simulation,patients,plans,mechanisms,spec,S,settings)
% Enumerate every reference-controlled patient-mechanism pair with at least
% one actual control-to-failure draw. The displayed patient is selected by a
% prespecified patient-level rule: maximize the patient's leading failure
% count out of the mechanism draws (20 in paper mode); if tied, minimize the RMS
% range-normalized input distance of an actual failing draw from that
% patient's reference vector. Remaining exact ties use patient number and
% then the prespecified M1-M5 order. Day-360 output distance is never used.

N = numel(patients);
G = numel(mechanisms);
Sdraw = settings.nPerturbations;

MechanismNumber = zeros(0,1);
MechanismID = strings(0,1);
Mechanism = strings(0,1);
VirtualPatient = zeros(0,1);
FailureSample = zeros(0,1);
FailureDraws = zeros(0,1);
ValidDraws = zeros(0,1);
FailurePercent = zeros(0,1);
NormalizedDistance = zeros(0,1);
ReferenceTumorCells = zeros(0,1);
FailureTumorCells = zeros(0,1);
AggregateMechanismEvents = zeros(0,1);

for g = 1:G
    ok = reshape(simulation.PerturbedOK(:,g,:),N,Sdraw);
    tumor = reshape(simulation.PerturbedFinalTumor(:,g,:),N,Sdraw);
    for n = find(simulation.ReferenceControlled(:)).'
        validSamples = find(ok(n,:) & isfinite(tumor(n,:)));
        failingSamples = find(ok(n,:) & isfinite(tumor(n,:)) & ...
            tumor(n,:)>=settings.detectionLimit);
        if isempty(failingSamples)
            continue
        end
        [closestSample,closestDistance] = closestFailureDraw( ...
            patients{n},plans(g),failingSamples,spec);

        MechanismNumber(end+1,1) = g; %#ok<AGROW>
        MechanismID(end+1,1) = string(mechanisms(g).id); %#ok<AGROW>
        Mechanism(end+1,1) = string(mechanisms(g).title); %#ok<AGROW>
        VirtualPatient(end+1,1) = n; %#ok<AGROW>
        FailureSample(end+1,1) = closestSample; %#ok<AGROW>
        FailureDraws(end+1,1) = numel(failingSamples); %#ok<AGROW>
        ValidDraws(end+1,1) = numel(validSamples); %#ok<AGROW>
        FailurePercent(end+1,1) = ...
            100*numel(failingSamples)/numel(validSamples); %#ok<AGROW>
        NormalizedDistance(end+1,1) = closestDistance; %#ok<AGROW>
        ReferenceTumorCells(end+1,1) = ...
            simulation.ReferenceFinalTumor(n); %#ok<AGROW>
        FailureTumorCells(end+1,1) = tumor(n,closestSample); %#ok<AGROW>
        AggregateMechanismEvents(end+1,1) = ...
            S.TransitionEvents(g); %#ok<AGROW>
    end
end

Candidate = table(MechanismNumber,MechanismID,Mechanism,VirtualPatient, ...
    FailureSample,FailureDraws,ValidDraws,FailurePercent, ...
    NormalizedDistance,ReferenceTumorCells,FailureTumorCells, ...
    AggregateMechanismEvents);

Candidate.PatientMechanismBreadth = zeros(height(Candidate),1);
Candidate.PatientTotalFailureDraws = zeros(height(Candidate),1);
Candidate.PatientLeadingFailureDraws = zeros(height(Candidate),1);
for n = unique(Candidate.VirtualPatient).'
    mask = Candidate.VirtualPatient==n;
    Candidate.PatientMechanismBreadth(mask) = ...
        numel(unique(Candidate.MechanismNumber(mask)));
    Candidate.PatientTotalFailureDraws(mask) = ...
        sum(Candidate.FailureDraws(mask));
    Candidate.PatientLeadingFailureDraws(mask) = ...
        max(Candidate.FailureDraws(mask));
end
Candidate.CandidateRow = (1:height(Candidate)).';

Selected = Candidate([],:);
Selected.Rank = zeros(0,1);
Selected.SelectionCriterion = strings(0,1);
if isempty(Candidate)
    fprintf(['No reference-controlled patient had an actual ' ...
        'failure-producing perturbation; no 3D case was fabricated.\n'])
    return
end

patientRepresentativeRows = zeros(0,1);
for n = reshape(unique(Candidate.VirtualPatient),1,[])
    patientRows = find(Candidate.VirtualPatient==n);
    leadingCount = max(Candidate.FailureDraws(patientRows));
    leadingRows = patientRows(Candidate.FailureDraws(patientRows)==leadingCount);
    scoreWithinPatient = [Candidate.NormalizedDistance(leadingRows), ...
        Candidate.MechanismNumber(leadingRows)];
    [~,withinOrder] = sortrows(scoreWithinPatient,[1 2]);
    patientRepresentativeRows(end+1,1) = ...
        leadingRows(withinOrder(1)); %#ok<AGROW>
end

selectionScore = [ ...
    -Candidate.PatientLeadingFailureDraws(patientRepresentativeRows), ...
    Candidate.NormalizedDistance(patientRepresentativeRows), ...
    Candidate.VirtualPatient(patientRepresentativeRows), ...
    Candidate.MechanismNumber(patientRepresentativeRows)];
[~,selectionOrder] = sortrows(selectionScore,1:size(selectionScore,2));
selectedRow = patientRepresentativeRows(selectionOrder(1));

Selected = Candidate(selectedRow,:);
Selected.Rank = 1;
Selected.SelectionCriterion = ...
    "largest patient-leading failure count; then smallest RMS " + ...
    "range-normalized failing-input distance; then patient number; " + ...
    "then M1--M5 order";
Selected = movevars(Selected,{'Rank','SelectionCriterion'},'Before',1);

fprintf('Qualifying patient--mechanism candidate pairs: %d.\n', ...
    height(Candidate))
fprintf('Primary 3D patient case selected: VP %d, %s, %d/%d failures.\n', ...
    Selected.VirtualPatient,Selected.MechanismID, ...
    Selected.FailureDraws,Selected.ValidDraws)
end

function [sampleIndex,distance] = closestFailureDraw( ...
        patient,plan,failingSamples,spec)
distance = inf;
sampleIndex = failingSamples(1);
for s = reshape(failingSamples,1,[])
    currentDistance = normalizedPerturbationDistance( ...
        patient,plan,s,spec);
    if currentDistance<distance
        distance = currentDistance;
        sampleIndex = s;
    end
end
end

function distance = normalizedPerturbationDistance( ...
        patient,plan,sampleIndex,spec)
component = zeros(1,numel(plan.inputNames));
for q = 1:numel(plan.inputNames)
    name = plan.inputNames{q};
    s = getSpec(spec,name);
    referenceValue = patient.(name);
    perturbedValue = plan.values(sampleIndex,q);
    if strcmp(s.measure,'log') && referenceValue>0 && perturbedValue>0
        component(q) = (log(perturbedValue)-log(referenceValue))/ ...
            max(log(s.hi)-log(s.lo),eps);
    elseif s.hi>s.lo
        component(q) = (perturbedValue-referenceValue)/(s.hi-s.lo);
    end
end
% RMS normalization makes this tie-breaker comparable across mechanisms
% containing different numbers of varied inputs. It does not alter any
% transition count or mechanism definition.
distance = sqrt(mean(component.^2));
end

function validateCompletedAnalysis( ...
        simulation,summaryTable,patientTable,longTable, ...
        overlapTable,mechanisms,settings)
N = numel(simulation.ReferenceOK);
G = numel(mechanisms);
S = settings.nPerturbations;

if height(summaryTable)~=G
    error('Summary table must contain exactly five mechanism rows.')
end
if ~isequal(summaryTable.MechanismID,string({mechanisms.id}).')
    error('Summary rows are not in prespecified M1--M5 order.')
end
if height(patientTable)~=N*G
    error('Patient-by-mechanism table has an incorrect number of rows.')
end
if height(longTable)~=N*G*S
    error('Long-form table has an incorrect number of rows.')
end
if height(overlapTable)~=N
    error('Patient-overlap table has an incorrect number of rows.')
end

if any(simulation.ReferenceOK & ...
        (~isfinite(simulation.ReferenceFinalTumor) | ...
        simulation.ReferenceFinalTumor<0))
    error('An accepted reference simulation has an invalid endpoint.')
end
acceptedPerturbed = simulation.PerturbedOK;
if any(acceptedPerturbed(:) & ...
        (~isfinite(simulation.PerturbedFinalTumor(:)) | ...
        simulation.PerturbedFinalTumor(:)<0))
    error('An accepted perturbed simulation has an invalid endpoint.')
end

% Absorbed states must remain exactly zero at day 360.
refTs = simulation.ReferenceFinalState(:,2);
refTu = simulation.ReferenceFinalState(:,3);
if any(simulation.ReferenceTsAbsorbed & refTs~=0)
    error('A reference Ts-absorbed state is nonzero at day 360.')
end
if any(simulation.ReferenceTotalAbsorbed & (refTs~=0 | refTu~=0))
    error('A reference total-absorbed state is nonzero at day 360.')
end

pertTs = simulation.PerturbedFinalState(:,:,:,2);
pertTu = simulation.PerturbedFinalState(:,:,:,3);
if any(simulation.PerturbedTsAbsorbed(:) & pertTs(:)~=0)
    error('A perturbed Ts-absorbed state is nonzero at day 360.')
end
if any(simulation.PerturbedTotalAbsorbed(:) & ...
        (pertTs(:)~=0 | pertTu(:)~=0))
    error('A perturbed total-absorbed state is nonzero at day 360.')
end

% Reconcile all published counts with the long-form table.
for g = 1:G
    mask = longTable.MechanismID==mechanisms(g).id;
    valid = mask & longTable.ReferenceControlled & ...
        longTable.PerturbedValid;
    transitions = valid & longTable.ControlToFailure;
    affected = unique(longTable.VirtualPatient(transitions));

    if sum(valid)~=summaryTable.ValidPerturbations(g)
        error('Valid-perturbation denominator mismatch for %s.', ...
            mechanisms(g).id)
    end
    if sum(transitions)~=summaryTable.TransitionEvents(g)
        error('Transition-event mismatch for %s.',mechanisms(g).id)
    end
    if numel(affected)~=summaryTable.AffectedPatients(g)
        error('Distinct-patient mismatch for %s.',mechanisms(g).id)
    end

    patientMask = patientTable.MechanismID==mechanisms(g).id;
    if sum(patientTable.ValidDraws(patientMask))~= ...
            summaryTable.ValidPerturbations(g) || ...
            sum(patientTable.TransitionDraws(patientMask))~= ...
            summaryTable.TransitionEvents(g) || ...
            sum(patientTable.AnyTransition(patientMask))~= ...
            summaryTable.AffectedPatients(g)
        error('Patient-level table does not reconcile for %s.', ...
            mechanisms(g).id)
    end

    overlapName = sprintf('AffectedByM%d',g);
    expectedOverlap = patientTable.AnyTransition(patientMask);
    if ~isequal(overlapTable.(overlapName),expectedOverlap)
        error('Patient-overlap flags do not reconcile for %s.', ...
            mechanisms(g).id)
    end
end

if any(summaryTable.ValidPerturbations~= ...
        sum(simulation.ReferenceControlled)*S)
    error('Valid denominators do not equal eligible patients times draws.')
end
if ~isequal(overlapTable.NumberOfMechanismsWithTransition, ...
        sum([overlapTable.AffectedByM1,overlapTable.AffectedByM2, ...
        overlapTable.AffectedByM3,overlapTable.AffectedByM4, ...
        overlapTable.AffectedByM5],2))
    error('Patient-overlap mechanism totals do not reconcile.')
end

if any(summaryTable.TransitionEvents>summaryTable.ValidPerturbations)
    error('A transition count exceeds its valid perturbation denominator.')
end
if any(summaryTable.AffectedPatients> ...
        summaryTable.EligibleReferencePatients)
    error('An affected-patient count exceeds the eligible cohort.')
end
if sum(summaryTable.TransitionEvents)>0 && ...
        abs(sum(summaryTable.ShareOfAllTransitionEventsPercent)-100)>1e-8
    error('Transition-event shares do not sum to 100%%.')
end

fprintf('Completed-analysis reconciliation passed.\n')
end

function printFinalSummary(S,simulation,settings,elapsedSeconds)
fprintf('\n===============================================================\n')
fprintf('ANALYSIS COMPLETED AND RECONCILED\n')
fprintf('===============================================================\n')
fprintf('Reference simulations completed: %d/%d\n', ...
    sum(simulation.ReferenceOK),numel(simulation.ReferenceOK))
fprintf('Reference-controlled virtual patients: %d/%d\n', ...
    sum(simulation.ReferenceControlled),numel(simulation.ReferenceControlled))
fprintf('Total control-to-failure events: %d\n',sum(S.TransitionEvents))
fprintf('Distinct patients with >=1 transition: %d\n', ...
    countDistinctAffectedPatients(simulation,settings))
fprintf('\nMechanism  Events/valid  Event rate  Affected VPs/eligible\n')
for g = 1:height(S)
    fprintf('%-4s       %4d/%-5d   %7.3f%%      %3d/%-3d\n', ...
        char(S.MechanismID(g)),S.TransitionEvents(g), ...
        S.ValidPerturbations(g),S.TransitionPercent(g), ...
        S.AffectedPatients(g),S.EligibleReferencePatients(g))
end
fprintf('\nElapsed time: %.1f min\n',elapsedSeconds/60)
fprintf('Outputs: %s\n',settings.outputDir)
fprintf(['Interpretation: conditional model-encoded susceptibility; ' ...
    'not a normalized ranking of biological importance.\n'])
end

function n = countDistinctAffectedPatients(simulation,settings)
eligible = simulation.ReferenceControlled(:);
N = numel(eligible);
G = size(simulation.PerturbedOK,2);
affected = false(N,1);
for g = 1:G
    ok = reshape(simulation.PerturbedOK(:,g,:), ...
        N,settings.nPerturbations);
    tumor = reshape(simulation.PerturbedFinalTumor(:,g,:), ...
        N,settings.nPerturbations);
    affected = affected | any(eligible & ok & ...
        tumor>=settings.detectionLimit,2);
end
n = sum(affected);
end

function validateOutputManifest(settings)
required = { ...
    'Mechanism_analysis_reproducibility_manifest.txt', ...
    'Input_specification_and_sampling_rules.csv', ...
    'Mechanism_definitions.csv', ...
    'Virtual_patient_inputs.csv', ...
    'Mechanism_perturbation_design.csv', ...
    'Numerical_completion_audit.csv', ...
    'Mechanism_transition_summary.csv', ...
    'Patient_by_mechanism_summary.csv', ...
    'Patient_mechanism_draw_long_table.csv', ...
    'Patient_mechanism_overlap.csv', ...
    'Patient_cluster_bootstrap_intervals.csv', ...
    'Candidate_patient_mechanism_cases.csv', ...
    'Selected_candidate_patient_cases.csv', ...
    'MMC_mechanism_failure_results.mat', ...
    'MMC_Mechanism_Failure_NPJ_EXECUTED_SOURCE.m'};

if settings.runCohortStability
    required = [required,{'Cohort_size_stability.csv'}]; %#ok<AGROW>
end

File = string(required(:));
Exists = false(numel(required),1);
Bytes = zeros(numel(required),1);
for q = 1:numel(required)
    pathValue = fullfile(settings.outputDir,required{q});
    info = dir(pathValue);
    Exists(q) = numel(info)==1 && info.bytes>0;
    if Exists(q)
        Bytes(q) = info.bytes;
    end
end
manifest = table(File,Exists,Bytes);
writetable(manifest,fullfile(settings.outputDir,'Output_file_manifest.csv'))
if ~all(Exists)
    error('Output-file audit failed; missing or empty: %s.', ...
        strjoin(cellstr(File(~Exists)),', '))
end
fprintf('Output-file audit passed: %d required files exist and are nonempty.\n', ...
    numel(required))
end

function cells = diameterToCells(diameterMM)
cells = pi*(diameterMM/2)^2*(3*0.01)*1e6;
end
