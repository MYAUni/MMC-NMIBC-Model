%% Build Supplementary Data 2 from the saved GSA results
% Run this script in MATLAB R2020b or later, with Microsoft Excel on Windows.
% It writes one workbook containing a README and six source-data worksheets.
% Running it again replaces that output workbook. Close the workbook first.
% The CSV files are only read; no simulations or bootstrap samples are rerun.

%% File locations
resultsFolder = ...
    'C:\Users\marom\Desktop\Final Codes\Paper codes\19.8\GSA_Continuous_Day360_26Inputs_paper_20260821_174247';
convergenceFolder = fullfile(resultsFolder, ...
    'GSA_Postprocessed_Publication_20260822_221826');
outputFile = fullfile(resultsFolder, ...
    'Supplementary_Data_2_GSA_Sobol_Results.xlsx');

% Worksheet name, CSV filename, and folder. Upload numbers are not used.
sources = {
    'Primary Sobol',       'Sobol_indices_primary.csv',                        resultsFolder;
    'r-log robustness',    'Sobol_indices_r_log_robustness.csv',               resultsFolder;
    'Scenario summary',    'Sobol_scenario_summary_including_sum_S1.csv',      resultsFolder;
    'Rank stability',      'Sobol_r_sampling_scale_rank_stability.csv',       resultsFolder;
    'Convergence summary', 'Sobol_nested_sample_convergence_summary.csv',     convergenceFolder;
    'Convergence indices', 'Sobol_nested_sample_convergence.csv',             convergenceFolder
};

%% Read the six CSV tables
% Specify the header and first data row explicitly to avoid import guessing.
textColumns = {'Factor','Label','Role','Sampling','Scenario','R_Sampling','Design'};
tables = cell(size(sources,1),1);

for s = 1:size(sources,1)
    inputFile = fullfile(sources{s,3}, sources{s,2});
    if ~isfile(inputFile)
        error('Missing source file: %s', inputFile);
    end
    [fid, message] = fopen(inputFile, 'rt', 'n', 'UTF-8');
    if fid < 0
        error('Cannot read %s: %s', inputFile, message);
    end
    headerLine = fgetl(fid);
    fclose(fid);
    headers = strsplit(erase(headerLine, char(65279)), ',');

    options = delimitedTextImportOptions('NumVariables', numel(headers));
    options.Delimiter = ',';
    options.DataLines = [2 Inf];
    options.VariableNamesLine = 0;
    options.VariableNames = headers;
    options.Encoding = 'UTF-8';
    options.ConsecutiveDelimitersRule = 'split';
    options.ExtraColumnsRule = 'error';
    options.EmptyLineRule = 'skip';
    options.ImportErrorRule = 'error';
    options = setvartype(options, headers, 'double');
    textNames = headers(ismember(headers, textColumns));
    options = setvartype(options, textNames, 'char');
    options = setvaropts(options, textNames, ...
        'WhitespaceRule', 'preserve', 'QuoteRule', 'remove');
    tables{s} = readtable(inputFile, options);

    values = tables{s}{:,~ismember(headers,textColumns)};
    if isempty(values) || any(~isfinite(values(:)))
        error('Missing or nonnumeric result values in %s.', inputFile);
    end
    fprintf('Read %s: %d rows.\n', sources{s,2}, height(tables{s}));
end

%% Check that the main and convergence tables describe the same results
primary = tables{1};
summary = tables{3};
rankStability = tables{4};
convergenceSummary = tables{5};
convergence = tables{6};
scenarioNames = {'Primary_uniform_r','Robustness_log_uniform_r'};
samplingNames = {'uniform','log-uniform'};
rankPrefixes = {'Primary_','Rlog_'};
indexFields = {'S1','ST','S1_Rank','ST_Rank'};
convergenceFields = {'S1_Raw','ST_Raw','S1_Rank','ST_Rank'};
tolerance = 1e-10;

for s = 1:2
    T = tables{s};
    scenarioRow = summary(strcmp(summary.R_Sampling, samplingNames{s}),:);
    assert(height(scenarioRow)==1, 'Missing or duplicate scenario summary.');
    assert(height(T)==scenarioRow.NumberOfFactors && ...
        numel(unique(T.Factor))==height(T), 'The factor lists are inconsistent.');

    C = convergence(strcmp(convergence.Scenario,scenarioNames{s}) & ...
        convergence.BaseSampleSize==scenarioRow.BaseSampleSize,:);
    [found, order] = ismember(T.Factor, C.Factor);
    assert(height(C)==height(T) && all(found) && ...
        numel(unique(C.Factor))==height(C), ...
        'The final convergence table has a different factor list.');
    difference = T{:,indexFields} - C{order,convergenceFields};
    assert(all(abs(difference(:))<tolerance), ...
        'The convergence results do not match the main Sobol results.');

    [found, order] = ismember(T.Factor, rankStability.Factor);
    assert(height(rankStability)==height(T) && all(found), ...
        'The rank-stability table has a different factor list.');
    rankFields = strcat(rankPrefixes{s}, indexFields);
    difference = T{:,indexFields} - rankStability{order,rankFields};
    assert(all(abs(difference(:))<tolerance), ...
        'The rank-stability values do not match the Sobol tables.');

    assert(abs(sum(T.S1)-scenarioRow.Sum_S1)<tolerance && ...
        abs(sum(T.ST)-scenarioRow.Sum_ST)<tolerance, ...
        'The scenario summary does not agree with the individual indices.');
end

for r = 1:height(convergenceSummary)
    S = convergenceSummary(r,:);
    C = convergence(strcmp(convergence.Scenario,S.Scenario{1}) & ...
        convergence.BaseSampleSize==S.BaseSampleSize,:);
    assert(height(C)==height(primary) && ...
        abs(sum(C.S1_Raw)-S.Sum_S1)<tolerance && ...
        abs(sum(C.ST_Raw)-S.Sum_ST)<tolerance, ...
        'The detailed and summary convergence tables disagree.');
end
fprintf('The source tables agree.\n');

%% README
readme = {
    'Supplementary Data 2 | Sobol global sensitivity analysis', '', '';
    '', '', '';
    'Item', 'Specification', '';
    'Purpose', 'Numerical source data for the main and supplementary GSA results.', '';
    'Model output', 'Y = log10[Ts(360) + Tu(360) + 1]. Time is in days.', '';
    'Design', sprintf('%s; %d inputs; N = %d; %d model evaluations per scenario.', ...
        summary.Design{1}, height(primary), summary.BaseSampleSize(1), ...
        summary.ModelEvaluations(1)), '';
    'Sampling interpretation', 'Uniform and log-uniform mappings specify computational sampling over supported bounds, not biological or patient-population distributions.', '';
    'Primary r sampling', 'Uniform over 0.00625-0.5 day^-1.', '';
    'Robustness analysis', 'The same design coordinates and bounds are used; only the r mapping changes to log-uniform.', '';
    'Sobol estimators', 'First-order S1: Saltelli estimator. Total-order ST: Jansen estimator.', '';
    'Bootstrap intervals', '95% percentile intervals from 300 paired-row bootstrap resamples (2.5th and 97.5th percentiles).', '';
    'Interval interpretation', 'Numerical resampling intervals conditional on the bounds, sampling measures and evaluated design. These are not biological, clinical or patient-level confidence intervals.', '';
    'Numerical values', 'Source estimates and interval endpoints, including negative values, are preserved. Excel number formats affect display only.', '';
    'State interpretation', 'Tumor states are continuous cell equivalents. No one-cell absorbing boundary or positive-value clamp is used in this GSA.', '';
    'Numerical solver', 'MATLAB ode15s; RelTol = 1e-6; AbsTol = 1e-9; MaxStep = 1 day; NonNegative = 1:7.', '';
    'MMC schedule', 'Twelve sessions, each with a 2-hour dwell. Integration is split at every session start and end; m is constant during each dwell and zero between sessions.', '';
    'Randomization', 'rng(22,''twister''); scrambled Sobol sequence with Skip = 1000 and Leap = 0.', '';
    'Convergence interpretation', 'Nested prefixes of one scrambled Sobol design at N = 128, 256, 512 and 1024. This assesses numerical stability, not uncertainty across independent scrambles.', '';
    'Rank-change convention', 'ST_Rank_Change = robustness rank minus primary rank. A negative value indicates movement toward a higher rank.', '';
    'Conditional input mapping', 'Sobol indices refer to independent unit coordinates. T0 is mapped conditionally on k; the k and T0 indices therefore depend on this declared generator. p8 is sampled independently.', '';
    'Labels and units', 'Factor identifies the model input; Label contains source LaTeX notation. Sobol indices and ranks are dimensionless. Input units and supporting evidence are documented in Supplementary Data 1.', '';
    '', '', '';
    'Workbook sheet', 'Source file', 'Content'
};
descriptions = {
    'Primary indices, bootstrap intervals and ranks';
    'Indices, intervals and ranks with log-uniform r';
    'Scenario settings, output variance and summed indices';
    'Primary-versus-robustness indices and ranks';
    'Nested-sample convergence summaries';
    'All individual indices and ranks at each nested sample size'
};
readme = [readme; sources(:,1:2), descriptions];

%% Write and format one workbook
% Build a temporary workbook, then replace the output after it is saved.
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
    for s = 1:size(sources,1)
        writetable(tables{s}, temporaryFile, 'Sheet', sources{s,1}, ...
            'UseExcel', false, 'AutoFitWidth', false);
    end
    book = excel.Workbooks.Open(temporaryFile);
    assert(book.Worksheets.Count==7, 'Expected seven worksheets.');

    % Plain, consistent formatting. Raw source column names are unchanged.
    headerColor = 31 + 78*256 + 61*65536;
    sheetNames = [{'README'}; sources(:,1)];
    for s = 1:numel(sheetNames)
        ws = book.Worksheets.Item(sheetNames{s});
        ws.Activate;
        excel.ActiveWindow.DisplayGridlines = false;
        excel.ActiveWindow.Zoom = 85;
        ws.UsedRange.Font.Name = 'Arial';
        ws.UsedRange.Font.Size = 11;
        ws.UsedRange.VerticalAlignment = -4160; % Excel: top
        ws.UsedRange.WrapText = true;

        if s==1
            ws.Range('A:A').ColumnWidth = 29;
            ws.Range('B:B').ColumnWidth = 76;
            ws.Range('C:C').ColumnWidth = 45;
            ws.Range('A1:C1').Merge;
            ws.Range('A1').Font.Size = 16;
            ws.Range('A1:C1').Font.Bold = true;
            ws.Range('A1:C1').Interior.Color = headerColor;
            ws.Range('A1:C1').Font.Color = 16777215; % White
            for row = [3 23]
                header = ws.Range(sprintf('A%d:C%d', row, row));
                header.Font.Bold = true;
                header.Interior.Color = headerColor;
                header.Font.Color = 16777215;
            end
            ws.Range('A4:A21').Font.Bold = true;
            ws.UsedRange.Rows.AutoFit;
            ws.Range('A1:C1').RowHeight = 32;
        else
            T = tables{s-1};
            lastRow = height(T)+1;
            lastColumn = char('A'+width(T)-1); % These tables have at most 14 columns.
            fullRange = ws.Range(sprintf('A1:%s%d', lastColumn, lastRow));
            excelTable = ws.ListObjects.Add(1, fullRange, [], 1);
            excelTable.Name = sprintf('SD2Table%d', s-1);
            excelTable.TableStyle = 'TableStyleLight1';
            excelTable.ShowTableStyleRowStripes = true;
            names = T.Properties.VariableNames;
            for c = 1:width(T)
                letter = char('A'+c-1);
                column = ws.Range([letter ':' letter]);
                body = ws.Range(sprintf('%s2:%s%d', letter, letter, lastRow));
                name = names{c};
                if ismember(name,textColumns)
                    column.ColumnWidth = 28;
                    if ismember(name,{'Factor','Label'})
                        column.ColumnWidth = 13;
                    elseif ismember(name,{'Sampling','Design'})
                        column.ColumnWidth = 38;
                    end
                else
                    column.ColumnWidth = 22;
                    body.WrapText = false;
                    body.NumberFormat = '0.000000';
                    if contains(name,'Rank') || ismember(name, ...
                            {'BaseSampleSize','NumberOfFactors','ModelEvaluations', ...
                             'Top5_ST_Overlap_with_Final_N'})
                        body.NumberFormat = '0';
                    end
                    if contains(name,'Spearman')
                        body.NumberFormat = '0.000000';
                    end
                    if ismember(name,{'LowerBound','UpperBound'})
                        body.NumberFormat = '0.000000E+00';
                    end
                end
            end
            header = ws.Range(['A1:' lastColumn '1']);
            header.Font.Bold = true;
            header.Font.Color = 16777215;
            header.Interior.Color = headerColor;
            header.WrapText = true;
            ws.UsedRange.Rows.AutoFit;
            header.RowHeight = max(66, header.RowHeight);
            ws.Range('B2').Select;
            excel.ActiveWindow.SplitRow = 1;
            excel.ActiveWindow.SplitColumn = 1;
            excel.ActiveWindow.FreezePanes = true;
        end
        ws.Range('A1').Select;
    end

    book.Worksheets.Item('README').Activate;
    book.Save;
    book.Close(false);
    book = [];
    excel.Quit;
    delete(excel);
    excel = [];

    [saved, message] = movefile(temporaryFile, outputFile, 'f');
    if ~saved
        error('Cannot save %s. Close it in Excel and run again. %s', ...
            outputFile, message);
    end
catch exception
    if ~isempty(book)
        try
            book.Close(false);
        catch
        end
    end
    if ~isempty(excel)
        try
            excel.Quit;
            delete(excel);
        catch
        end
    end
    if isfile(temporaryFile)
        delete(temporaryFile);
    end
    rethrow(exception);
end

fprintf('\nCreated one workbook:\n%s\n', outputFile);
