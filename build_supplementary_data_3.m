%% File locations
% Build one Supplementary Data 3 workbook from the saved CSV and MAT files.
% README text, formulas and formatting are specified below.
% No existing workbook or template is needed.
% Requires MATLAB R2020b or later and Microsoft Excel on Windows.

resultsFolder = 'C:\Users\marom\Downloads\MMC_mechanism_failure_final_paper_results';
figuresFolder = fullfile(resultsFolder, 'publication_figures');

outputFile = fullfile(resultsFolder, ...
    'Supplementary_Data_3_MMC_treatment_failure_analysis.xlsx');

if ~isfolder(resultsFolder)
    error('The results folder does not exist: %s', resultsFolder);
end
if isfile(outputFile)
    error('The output already exists. Change outputFile above or rename the existing output.');
end
%% Read the CSV files
% Each row lists a worksheet, its input filename, and its expected size.
csvFiles = {
    'Mechanism_Summary',      'Mechanism_transition_summary.csv',                5, 17;
    'Bootstrap_Intervals',    'Patient_cluster_bootstrap_intervals.csv',         5,  6;
    'Patient_by_Mechanism',   'Patient_by_mechanism_summary.csv',              2500, 12;
    'Patient_Overlap',        'Patient_mechanism_overlap.csv',                  500, 10;
    'Mechanism_Definitions',  'Mechanism_definitions.csv',                       5,  8;
    'Perturbation_Design',    'Mechanism_perturbation_design.csv',              240,  7;
    'Draw_Level_Results',     'Patient_mechanism_draw_long_table.csv',        50000, 16;
    'Virtual_Patient_Inputs', 'Virtual_patient_inputs.csv',                    500, 59;
    'Candidate_Cases',        'Candidate_patient_mechanism_cases.csv',           47, 16;
    'Selected_Case',          'Selected_candidate_patient_cases.csv',            1, 18;
    'Case_Screen',     'Fig10_Selected_Patient_VP0127_M4_patient_mechanism_screen.csv',          5, 6;
    'Case_Parameters', 'Fig10_Selected_Patient_VP0127_M4_reference_and_failure_parameters.csv', 2, 5;
    'Case_Path',       'Fig10_Selected_Patient_VP0127_M4_continuous_parameter_path.csv',       82, 4;
    'Input_Rules',            'Input_specification_and_sampling_rules.csv',     35,  8;
    'Numerical_Audit',        'Numerical_completion_audit.csv',                  1, 12;
    'Group_Size_Stability',   'Cohort_size_stability.csv',                      15,  7
};

sheetNames = [csvFiles(:,1); {'Case_Surface'}];
sheetData = cell(17, 1);
for k = 1:size(csvFiles, 1)
    filename = csvFiles{k, 2};
    inputFile = fullfile(resultsFolder, filename);
    if startsWith(filename, 'Fig10_')
        inputFile = fullfile(figuresFolder, filename);
    end
    assert(isfile(inputFile), 'Missing input file: %s', inputFile);

    % Read every row as text, including the header. The CSV layout is known,
    % so MATLAB does not need to guess the header or the first data row.
    expectedRows = csvFiles{k, 3};
    expectedColumns = csvFiles{k, 4};
    options = delimitedTextImportOptions('NumVariables', expectedColumns);
    options.Delimiter = ',';
    options.DataLines = [1 Inf];
    options.VariableNamesLine = 0;
    options.VariableUnitsLine = 0;
    options.VariableDescriptionsLine = 0;
    options.ConsecutiveDelimitersRule = 'split';
    options.ExtraColumnsRule = 'error';
    options.EmptyLineRule = 'skip';
    options.Encoding = 'UTF-8';
    options = setvartype(options, options.VariableNames, 'string');
    options = setvaropts(options, options.VariableNames, ...
        'WhitespaceRule', 'preserve', 'QuoteRule', 'remove', 'TreatAsMissing', {});
    T = readtable(inputFile, options);
    textValues = string(T{:,:});
    textValues(ismissing(textValues)) = "";

    assert(size(textValues,1)==expectedRows+1 && size(textValues,2)==expectedColumns, ...
        ['%s: expected %d data rows and %d columns; read %d data rows ' ...
         'and %d columns. The first row is the CSV header.'], ...
        inputFile, expectedRows, expectedColumns, ...
        size(textValues,1)-1, size(textValues,2));

    headers = cellstr(textValues(1,:));
    values = cellstr(textValues(2:end,:));
    % Convert finite numbers only; keep the literal entry NaN as text.
    numbers = str2double(textValues(2:end,:));
    numericFields = isfinite(numbers);
    values(numericFields) = num2cell(numbers(numericFields));
    sheetData{k} = [headers; values];
    fprintf('Read %s: %d data rows, %d columns.\n', filename, expectedRows, expectedColumns);
end

%% Read the saved case surface
% Either MAT encoding is accepted, in either of the two input folders.
surfaceFiles = {
    fullfile(resultsFolder, 'Fig10_Selected_Patient_VP0127_M4_surface_data_v7.mat');
    fullfile(figuresFolder, 'Fig10_Selected_Patient_VP0127_M4_surface_data_v7.mat');
    fullfile(resultsFolder, 'Fig10_Selected_Patient_VP0127_M4_surface_data.mat');
    fullfile(figuresFolder, 'Fig10_Selected_Patient_VP0127_M4_surface_data.mat')
};
surfaceIndex = find(isfile(string(surfaceFiles)), 1);
assert(~isempty(surfaceIndex), 'The saved VP127/M4 surface MAT-file is missing.');
saved = load(surfaceFiles{surfaceIndex}, 'Surface');
S = saved.Surface;

positions = S.PathPosition(:);
times = S.TimeGrid(:);
assert(numel(positions)==41 && isequal(times,(0:2:360).') && ...
    isequal(size(S.Tumor),[41,181]) && S.TargetSample==16, ...
    'The saved surface does not match the audited VP127/M4 case.');

nTimes = numel(times);
surfaceValues = [
    repelem(positions, nTimes), ...
    repmat(times, numel(positions), 1), ...
    reshape(S.Tumor.', [], 1), ...
    repelem(S.FinalTumor(:), nTimes), ...
    repelem(double(S.TsAbsorbed(:)), nTimes), ...
    repelem(double(S.TotalTumorAbsorbed(:)), nTimes), ...
    repelem(S.TsAbsorptionTime(:), nTimes), ...
    repelem(S.TotalTumorAbsorptionTime(:), nTimes), ...
    repmat(double(S.TargetSample), numel(positions)*nTimes, 1)
];
surfaceHeaders = {'PathPosition','TimeDays','TumorCells','PathFinalTumorCells', ...
    'TsAbsorbed','TotalTumorAbsorbed','TsAbsorptionTimeDays', ...
    'TotalTumorAbsorptionTimeDays','TargetPerturbationDraw'};
surfaceCells = num2cell(surfaceValues);
surfaceCells(isnan(surfaceValues)) = {''};
sheetData{17} = [surfaceHeaders; surfaceCells];

%% Check the completed run
summary = sheetData{1};
eligibleColumn = strcmp(summary(1,:), 'EligibleReferencePatients');
transitionColumn = strcmp(summary(1,:), 'TransitionEvents');
assert(all(cell2mat(summary(2:end,eligibleColumn))==322) && ...
    sum(cell2mat(summary(2:end,transitionColumn)))==253, ...
    'The input files do not match the run described in the audited README.');

selectedCase = sheetData{10};
caseFields = {'VirtualPatient','MechanismNumber','FailureSample','FailureDraws','ValidDraws'};
expectedCase = [127, 4, 16, 17, 20];
for k = 1:numel(caseFields)
    column = strcmp(selectedCase(1,:), caseFields{k});
    assert(isequal(selectedCase{2,column}, expectedCase(k)), ...
        'The selected case differs from the audited case.');
end

pathData = sheetData{13};
pathPositions = cell2mat(pathData(2:end,strcmp(pathData(1,:), 'PathPosition')));
pathFinal = cell2mat(pathData(2:end,strcmp(pathData(1,:), 'PathFinalTumorCells')));
for k = 1:numel(positions)
    rows = abs(pathPositions-positions(k)) < 1e-12;
    assert(any(rows) && all(abs(pathFinal(rows)-S.FinalTumor(k)) <= ...
        1e-11*max(1,abs(S.FinalTumor(k)))), 'The case path and surface do not agree.');
end

%% README
% Cell positions and wording follow the audited workbook.
readme = repmat({''}, 40, 8);
readme{1,1} = 'Supplementary Data 3 | Mechanism-of-failure analysis';
readme{3,1} = 'Paper-mode source data for the virtual-patient analysis of model-encoded routes to MMC treatment failure';
readme{5,1} = 'Analysis item';
readme{5,2} = 'Value';
readme{5,4} = 'Interpretation and scope';
readme{6,1} = 'Reference simulations';
readme{6,2} = 500;
readme{6,4} = 'This analysis tests whether changes confined to one of five predefined model mechanisms convert a reference-controlled virtual patient to detectable failure at day 360.';
readme{7,1} = 'Perturbed simulations';
readme{7,2} = 50000;
readme{7,4} = 'M4 denotes the antitumor immune response: gamma governs mature-DC-driven effector-cell activation, and p5 scales effector-mediated tumor-cell killing.';
readme{8,1} = 'Reference-controlled virtual patients';
readme{8,2} = 322;
readme{8,4} = 'Event counts report outcomes under the specified parameter ranges and sampling rules.';
readme{9,1} = 'Control-to-failure transitions';
readme{9,2} = 253;
readme{9,4} = 'Bootstrap limits summarize numerical patient-cluster resampling and are not biological, clinical, or patient-population confidence intervals.';
readme{10,1} = 'Distinct affected virtual patients';
readme{10,2} = 21;
readme{11,1} = 'Selected virtual patient';
readme{11,2} = 127;
readme{12,1} = 'Selected mechanism';
readme{12,2} = 'M4';
readme{15,1} = 'Worksheet';
readme{15,2} = 'Data class';
readme{15,3} = 'Rows';
readme{15,4} = 'Description';
readme{15,5} = 'Source output';
readme{16,1} = 'Mechanism_Summary';
readme{16,2} = 'processed';
readme{16,3} = 5;
readme{16,4} = 'Mechanism-level event and affected-patient summaries, including bootstrap limits.';
readme{16,5} = 'Mechanism_transition_summary.csv';
readme{17,1} = 'Bootstrap_Intervals';
readme{17,2} = 'processed';
readme{17,3} = 5;
readme{17,4} = 'Patient-cluster bootstrap intervals from 2,000 replicates.';
readme{17,5} = 'Patient_cluster_bootstrap_intervals.csv';
readme{18,1} = 'Patient_by_Mechanism';
readme{18,2} = 'processed';
readme{18,3} = 2500;
readme{18,4} = 'One row per virtual-patient--mechanism pair.';
readme{18,5} = 'Patient_by_mechanism_summary.csv';
readme{19,1} = 'Patient_Overlap';
readme{19,2} = 'processed';
readme{19,3} = 500;
readme{19,4} = 'Mechanism overlap and number of mechanisms producing transitions in each patient.';
readme{19,5} = 'Patient_mechanism_overlap.csv';
readme{20,1} = 'Mechanism_Definitions';
readme{20,2} = 'specification';
readme{20,3} = 5;
readme{20,4} = 'Biological mechanism definitions and assigned, varied, and fixed parameters.';
readme{20,5} = 'Mechanism_definitions.csv';
readme{21,1} = 'Perturbation_Design';
readme{21,2} = 'raw';
readme{21,3} = 240;
readme{21,4} = 'The 20 joint LHS vectors for M1--M5, including unit coordinates and mappings.';
readme{21,5} = 'Mechanism_perturbation_design.csv';
readme{22,1} = 'Draw_Level_Results';
readme{22,2} = 'raw';
readme{22,3} = 50000;
readme{22,4} = 'Complete 50,000-row mechanism-perturbation output table.';
readme{22,5} = 'Patient_mechanism_draw_long_table.csv';
readme{23,1} = 'Virtual_Patient_Inputs';
readme{23,2} = 'raw';
readme{23,3} = 500;
readme{23,4} = 'Input values and unit coordinates for the 500-patient mechanism-analysis group.';
readme{23,5} = 'Virtual_patient_inputs.csv';
readme{24,1} = 'Candidate_Cases';
readme{24,2} = 'processed';
readme{24,3} = 47;
readme{24,4} = 'All qualifying patient--mechanism cases considered for the individual example.';
readme{24,5} = 'Candidate_patient_mechanism_cases.csv';
readme{25,1} = 'Selected_Case';
readme{25,2} = 'processed';
readme{25,3} = 1;
readme{25,4} = 'Prespecified selection audit for virtual patient 127 and M4.';
readme{25,5} = 'Selected_candidate_patient_cases.csv';
readme{26,1} = 'Case_Screen';
readme{26,2} = 'figure source';
readme{26,3} = 5;
readme{26,4} = 'The five-mechanism screen for virtual patient 127.';
readme{26,5} = 'Fig10_Selected_Patient_VP0127_M4_patient_mechanism_screen.csv';
readme{27,1} = 'Case_Parameters';
readme{27,2} = 'figure source';
readme{27,3} = 2;
readme{27,4} = 'Reference and selected failure-inducing M4 parameter values.';
readme{27,5} = 'Fig10_Selected_Patient_VP0127_M4_reference_and_failure_parameters.csv';
readme{28,1} = 'Case_Path';
readme{28,2} = 'figure source';
readme{28,3} = 82;
readme{28,4} = 'Continuous M4 path coordinates and day-360 tumor burden for Figure 10B.';
readme{28,5} = 'Fig10_Selected_Patient_VP0127_M4_continuous_parameter_path.csv';
readme{29,1} = 'Input_Rules';
readme{29,2} = 'specification';
readme{29,3} = 35;
readme{29,4} = 'Mechanism-analysis input bounds, units, roles, and computational sampling measures.';
readme{29,5} = 'Input_specification_and_sampling_rules.csv';
readme{30,1} = 'Numerical_Audit';
readme{30,2} = 'audit';
readme{30,3} = 1;
readme{30,4} = 'Simulation-completion and one-cell-rule audit totals.';
readme{30,5} = 'Numerical_completion_audit.csv';
readme{31,1} = 'Group_Size_Stability';
readme{31,2} = 'audit';
readme{31,3} = 15;
readme{31,4} = 'Random without-replacement group-size resampling summaries produced by the paper-mode code.';
readme{31,5} = 'Cohort_size_stability.csv';
readme{32,1} = 'Case_Surface';
readme{32,2} = 'figure source';
readme{32,3} = 7421;
readme{32,4} = 'Time-resolved tumor burden for all 41 positions along the selected VP127/M4 parameter path and all 181 saved time points.';
readme{32,5} = 'Fig10_Selected_Patient_VP0127_M4_surface_data.mat';
readme{34,1} = 'Technical notes';
readme{35,1} = 'Units';
readme{35,4} = 'Headers containing Percent, Lower95, or Upper95 are stored in percentage-point units exactly as exported by MATLAB.';
readme{36,1} = 'Perturbation design';
readme{36,4} = 'The same 20 absolute joint vectors for a mechanism were applied to every virtual patient.';
readme{37,1} = 'Missing values';
readme{37,4} = 'Blank cells or NaN in absorption-time fields mean that the corresponding absorbing event was not observed; they do not represent time zero.';
readme{38,1} = 'Group-size stability';
readme{38,4} = 'Random samples without replacement were used. The N=500 entries contain the complete mechanism-analysis group, so their lower and upper limits equal the point estimate.';
readme{39,1} = 'Surface source';
readme{39,4} = 'Case_Surface was imported from a MATLAB -v7 conversion of the saved v7.3 surface file; only file encoding changed, not numerical values.';
readme{40,1} = 'Data preservation';
readme{40,4} = 'All other worksheets preserve the generated values exactly. No simulation outcomes were recalculated or modified during this revision. The one-cell absorbing rules apply only to this mechanism-of-failure analysis.';

readmeFormulas = {
    'B6', '=''Numerical_Audit''!A2';
    'B7', '=''Numerical_Audit''!D2';
    'B8', '=''Mechanism_Summary''!D2';
    'B9', '=SUM(''Mechanism_Summary''!F2:F6)';
    'B10', '=COUNTIF(''Patient_Overlap''!J2:J501,">0")';
    'B11', '=''Selected_Case''!F2';
    'B12', '=''Selected_Case''!D2';
};

%% Worksheet widths and number formats
columnWidths = {
    [38 22 16 92 58 15 15 15]; % README
    [18 30 15 22 21 19 18 19 18 18 22 22 22 18 18 18 18]; % Mechanism_Summary
    [18 22 18 18 18 18]; % Bootstrap_Intervals
    [15 18 17 22 18 13 18 18 16 20 22 22]; % Patient_by_Mechanism
    [15 17 22 18 15 15 15 15 15 22]; % Patient_Overlap
    [18 30 30 30 30 17 15 52]; % Mechanism_Definitions
    [18 15 18 18 18 22 21]; % Perturbation_Design
    [15 18 15 17 22 18 17 18 18 18 19 13 21 18 18 52]; % Draw_Level_Results
    [15 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 22 12 12 21 21 21 22]; % Virtual_Patient_Inputs
    [15 18 30 15 15 15 13 18 18 18 18 22 22 22 22 15]; % Candidate_Cases
    [12 52 15 18 30 15 15 15 13 18 18 18 18 22 22 22 22 15]; % Selected_Case
    [15 18 30 15 13 18]; % Case_Screen
    [18 18 18 22 22]; % Case_Parameters
    [15 18 18 18]; % Case_Path
    [18 18 18 18 22 22 22 22]; % Input_Rules
    [22 21 19 22 21 19 22 22 22 22 22 22]; % Numerical_Audit
    [15 18 13 22 18 18 18]; % Group_Size_Stability
    [16 12 22 22 15 24 24 24 24]; % Case_Surface
};

% Worksheet, data range, Excel number format.
numberFormats = {
    'Mechanism_Summary', 'C2:C6', '0';
    'Mechanism_Summary', 'D2:D6,F2:F6,H2:H6,J2:M6', '#,##0';
    'Mechanism_Summary', 'G2:G6,I2:I6,N2:Q6', '0.000';
    'Bootstrap_Intervals', 'B2:B6', '0';
    'Bootstrap_Intervals', 'C2:F6', '0.000';
    'Patient_by_Mechanism', 'A2:A2501,C2:D2501,I2:I2501', '0';
    'Patient_by_Mechanism', 'E2:E2501', '0.000000E+00';
    'Patient_by_Mechanism', 'F2:G2501,J2:L2501', '#,##0';
    'Patient_by_Mechanism', 'H2:H2501', '0.000';
    'Patient_Overlap', 'A2:C501,E2:I501', '0';
    'Patient_Overlap', 'D2:D501', '0.000000E+00';
    'Patient_Overlap', 'J2:J501', '#,##0';
    'Mechanism_Definitions', 'F2:G6', '0';
    'Perturbation_Design', 'B2:B241,G2:G241', '0';
    'Perturbation_Design', 'D2:D241', '0.000000';
    'Perturbation_Design', 'E2:E241', '0.000000E+00';
    'Draw_Level_Results', 'A2:A50001,C2:E50001,G2:G50001,K2:M50001', '0';
    'Draw_Level_Results', 'F2:F50001,H2:J50001', '0.000000E+00';
    'Draw_Level_Results', 'N2:O50001', '#,##0';
    'Virtual_Patient_Inputs', 'A2:A501', '0';
    'Virtual_Patient_Inputs', 'C2:C501,E2:E501,G2:G501,I2:I501,K2:K501,M2:M501,O2:O501,Q2:Q501,S2:S501,U2:U501,W2:W501,Y2:Y501,AA2:AA501,AC2:AC501,AE2:AE501,AG2:AG501,AI2:AI501,AK2:AK501,AM2:AM501,AO2:AO501,AQ2:AQ501,AS2:AS501,AU2:AU501,AW2:AW501,AY2:AY501,BA2:BA501', '0.000000';
    'Candidate_Cases', 'A2:A48,D2:D48,P2:P48', '0';
    'Candidate_Cases', 'F2:G48,L2:O48', '#,##0';
    'Candidate_Cases', 'H2:H48', '0.000';
    'Candidate_Cases', 'I2:I48', '0.000000';
    'Candidate_Cases', 'J2:K48', '0.000000E+00';
    'Selected_Case', 'A2,C2,F2,R2', '0';
    'Selected_Case', 'H2:I2,N2:Q2', '#,##0';
    'Selected_Case', 'J2', '0.000';
    'Selected_Case', 'K2', '0.000000';
    'Selected_Case', 'L2:M2', '0.000000E+00';
    'Case_Screen', 'A2:A6', '0';
    'Case_Screen', 'D2:E6', '#,##0';
    'Case_Screen', 'F2:F6', '0.000';
    'Case_Parameters', 'B2:C3', '0.000000E+00';
    'Case_Path', 'A2:A83', '0.000000';
    'Case_Path', 'C2:D83', '0.000000E+00';
    'Input_Rules', 'C2:D36', '0.000000E+00';
    'Numerical_Audit', 'A2:G2,J2', '#,##0';
    'Group_Size_Stability', 'A2:A16,D2:D16', '#,##0';
    'Group_Size_Stability', 'C2:C16', '0';
    'Group_Size_Stability', 'E2:G16', '0.000';
    'Case_Surface', 'A2:A7422', '0.000000';
    'Case_Surface', 'B2:B7422,E2:F7422,I2:I7422', '0';
    'Case_Surface', 'C2:D7422', '0.000000E+00';
    'Case_Surface', 'G2:H7422', '0.0';
};

%% README formatting
% Range, font, point size, bold, italic, RGB text color.
readmeFonts = {
    'A1:H2', 'Aptos Display', 18, true, false, [255 255 255];
    'A3:E3', 'Aptos', 11, false, true, [69 90 100];
    'F3:H3', 'Aptos', 11, false, true, [31 41 55];
    'A5:B5,D5:H5,A15:E15,A34:H34', 'Aptos', 10, true, false, [255 255 255];
    'A6:B12,F6:H12,B16:E31,F35:H39', 'Aptos', 10, false, false, [31 41 55];
    'D6:E12,B35:E39', 'Aptos', 10, false, false, [38 50 56];
    'A16:A31', 'Aptos', 10, true, false, [31 78 61];
    'A32,A40', 'Carlito', 11, true, false, [31 90 71];
    'B32:E32,B40:E40', 'Carlito', 11, false, false, [38 50 56];
    'A35:A39', 'Aptos', 10, true, false, [31 90 71];
};

% Range and RGB background color.
readmeFills = {
    'A1:E1', [31 90 71];
    'F1:H2,A2:E2,A15:E15', [31 78 61];
    'A3:E3', [226 240 217];
    'F3:H3', [226 240 217];
    'A5:B5,F5:H5,F34:H34', [112 173 71];
    'D5:E5,A34:E34', [112 173 71];
    'A6:A12,F6:H12,F35:H39', [243 248 239];
    'D6:E12,A32:E32,A35:E40', [243 248 239];
};

% Range, horizontal alignment, vertical alignment, wrap text.
readmeAlignment = {
    'A1:H2', 'left', 'center', false;
    'A3:H3', 'general', 'center', true;
    'A5:B12,F5:H5,A32:C32,F34:H34', 'general', 'bottom', false;
    'D5:E5,A34:E34', 'general', 'center', false;
    'D6:H12,A16:E31,A35:H39,A40:E40', 'general', 'top', true;
    'A15:E15', 'center', 'top', true;
    'D32:E32', 'general', 'bottom', true;
};

%% Write and format one Excel workbook
% Excel constants are named here to keep the formatting commands readable.
xlCenter = -4108;
xlLeft = -4131;
xlTop = -4160;
xlBottom = -4107;
xlGeneral = 1;
xlContinuous = 1;
xlThin = 2;
xlMedium = -4138;
xlEdgeLeft = 7;
xlEdgeTop = 8;
xlEdgeBottom = 9;
xlEdgeRight = 10;
xlInsideHorizontal = 12;
rgbWeights = [1; 256; 65536];

temporaryFile = [tempname(resultsFolder) '.xlsx'];
excel = [];
book = [];
try
    excel = actxserver('Excel.Application');
    excel.Visible = false;
    excel.DisplayAlerts = false;
    excel.ScreenUpdating = false;

    writecell(readme, temporaryFile, 'Sheet', 'README', ...
        'UseExcel', false, 'AutoFitWidth', false);
    for k = 1:numel(sheetNames)
        fprintf('Writing %s...\n', sheetNames{k});
        writecell(sheetData{k}, temporaryFile, 'Sheet', sheetNames{k}, ...
            'UseExcel', false, 'AutoFitWidth', false);
    end
    book = excel.Workbooks.Open(temporaryFile);
    assert(book.Worksheets.Count==18, 'The workbook must contain 18 worksheets.');

    % Normal font determines how Excel measures column widths.
    normalStyle = book.Styles.Item(1); % Built-in Normal style
    normalStyle.Font.Name = 'Carlito';
    normalStyle.Font.Size = 11;
    themeColors = [0 0 0; 255 255 255; 14 40 65; 232 232 232;
        21 96 130; 233 113 50; 25 107 36; 15 158 213;
        160 43 147; 78 167 46; 70 120 134; 150 96 125];
    for k = 1:12
        color = book.Theme.ThemeColorScheme.Colors(k);
        color.RGB = themeColors(k,:) * rgbWeights;
    end

    allNames = [{'README'}; sheetNames];
    lastColumns = cell(18,1);
    for k = 1:18
        ws = book.Worksheets.Item(allNames{k});
        ws.Activate;
        excel.ActiveWindow.DisplayGridlines = strcmp(allNames{k}, 'Draw_Level_Results');
        for c = 1:numel(columnWidths{k})
            n = c;
            label = '';
            while n > 0
                label = [char(65+mod(n-1,26)) label];
                n = floor((n-1)/26);
            end
            ws.Range([label ':' label]).ColumnWidth = columnWidths{k}(c);
        end
        lastColumns{k} = label;
    end

    % Data sheets: existing green headers, table stripes and number formats.
    for k = 1:17
        ws = book.Worksheets.Item(sheetNames{k});
        lastRow = size(sheetData{k},1);
        lastColumn = lastColumns{k+1};
        fullRange = ws.Range(sprintf('A1:%s%d', lastColumn, lastRow));
        body = ws.Range(sprintf('A2:%s%d', lastColumn, lastRow));
        header = ws.Range(['A1:' lastColumn '1']);
        fullRange.RowHeight = 15;

        if k~=7 && k~=17
            excelTable = ws.ListObjects.Add(1, fullRange, [], 1);
            excelTable.Name = sprintf('SD3Table%d', k);
            excelTable.TableStyle = 'TableStyleMedium4';
            excelTable.ShowTableStyleRowStripes = true;
        end
        body.Font.Name = 'Aptos';
        body.Font.Size = 10;
        body.Font.Color = [31 41 55] * rgbWeights;
        body.VerticalAlignment = xlCenter;
        if k==7, body.Font.Size = 9; end

        header.Font.Name = 'Aptos';
        header.Font.Size = 10;
        header.Font.Bold = true;
        header.Font.Color = [255 255 255] * rgbWeights;
        header.Interior.Color = [31 78 61] * rgbWeights;
        header.HorizontalAlignment = xlCenter;
        header.VerticalAlignment = xlCenter;
        header.WrapText = true;
        header.RowHeight = 32;
        for edge = [xlEdgeLeft xlEdgeTop xlEdgeBottom xlEdgeRight]
            border = header.Borders.Item(edge);
            border.LineStyle = xlContinuous;
            border.Weight = xlMedium;
            border.Color = [21 55 43] * rgbWeights;
        end

        if k==17
            body.Font.Name = 'Carlito';
            body.Font.Size = 11;
            body.Font.Color = [38 50 56] * rgbWeights;
            header.Font.Name = 'Carlito';
            header.Font.Size = 11;
            header.Interior.Color = [31 90 71] * rgbWeights;
            header.RowHeight = 34;
            header.Borders.LineStyle = xlContinuous;
            header.Borders.Weight = xlThin;
            header.Borders.Color = [40 216 91] * rgbWeights;
            for edge = [xlInsideHorizontal xlEdgeBottom]
                border = body.Borders.Item(edge);
                border.LineStyle = xlContinuous;
                border.Weight = xlThin;
                border.Color = [40 216 91] * rgbWeights;
            end
            banding = body.FormatConditions.Add(2, [], '=MOD(ROW(),2)=0');
            banding.Interior.Color = [198 239 206] * rgbWeights;
        end
    end
    for k = 1:size(numberFormats,1)
        ws = book.Worksheets.Item(numberFormats{k,1});
        ws.Range(numberFormats{k,2}).NumberFormat = numberFormats{k,3};
    end

    % README: original text, seven formulas, row heights and formatting.
    ws = book.Worksheets.Item('README');
    ws.Range('A1:H40').RowHeight = 15;
    ws.Range('A1').RowHeight = 34;
    ws.Range('A3').RowHeight = 52;
    ws.Range('A6:A9').RowHeight = 34;
    ws.Range('A16:A32').RowHeight = 30;
    ws.Range('A35:A40').RowHeight = 34;
    for k = 1:size(readmeFonts,1)
        area = ws.Range(readmeFonts{k,1});
        area.Font.Name = readmeFonts{k,2};
        area.Font.Size = readmeFonts{k,3};
        area.Font.Bold = readmeFonts{k,4};
        area.Font.Italic = readmeFonts{k,5};
        area.Font.Color = readmeFonts{k,6} * rgbWeights;
    end
    for k = 1:size(readmeFills,1)
        ws.Range(readmeFills{k,1}).Interior.Color = readmeFills{k,2} * rgbWeights;
    end
    for k = 1:size(readmeAlignment,1)
        area = ws.Range(readmeAlignment{k,1});
        switch readmeAlignment{k,2}
            case 'left', area.HorizontalAlignment = xlLeft;
            case 'center', area.HorizontalAlignment = xlCenter;
            otherwise, area.HorizontalAlignment = xlGeneral;
        end
        switch readmeAlignment{k,3}
            case 'top', area.VerticalAlignment = xlTop;
            case 'center', area.VerticalAlignment = xlCenter;
            otherwise, area.VerticalAlignment = xlBottom;
        end
        area.WrapText = readmeAlignment{k,4};
    end
    ws.Range('B6:B11,C16:C32').NumberFormat = '#,##0';

    % Border rectangles: range, edges, weight and RGB color.
    readmeBorders = {
        'A1:H2',   [7 8 9 10],       xlMedium, [21 55 43];
        'A5:B12',  [7 8 9 10 11 12], xlThin,   [198 224 180];
        'D6:E12',  [7 8 9 10],       xlThin,   [198 224 180];
        'F6:H12',  [8 9 10],         xlThin,   [198 224 180];
        'A15:E31', [7 8 9 10 11 12], xlThin,   [217 234 211];
        'A32:E32', [7 8 9 10 11 12], xlThin,   [198 224 180];
        'A35:E40', [7 8 9 10 11 12], xlThin,   [198 224 180];
        'F35:H39', [8 9 10],         xlThin,   [198 224 180]
    };
    for k = 1:size(readmeBorders,1)
        area = ws.Range(readmeBorders{k,1});
        for edge = readmeBorders{k,2}
            border = area.Borders.Item(edge);
            border.LineStyle = xlContinuous;
            border.Weight = readmeBorders{k,3};
            border.Color = readmeBorders{k,4} * rgbWeights;
        end
    end
    for k = 1:size(readmeFormulas,1)
        ws.Range(readmeFormulas{k,1}).Formula = readmeFormulas{k,2};
    end
    ws.Activate;
    excel.Calculate;
    book.Save;
    book.Close(false);
    book = [];
    excel.Quit;
    delete(excel);
    excel = [];
    assert(~isfile(outputFile), 'The output already exists and was not overwritten.');
    movefile(temporaryFile, outputFile);
catch exception
    if ~isempty(book)
        try, book.Close(false); catch, end
    end
    if ~isempty(excel)
        try, excel.Quit; delete(excel); catch, end
    end
    if isfile(temporaryFile), delete(temporaryFile); end
    rethrow(exception);
end

fprintf('Created: %s\n', outputFile);
