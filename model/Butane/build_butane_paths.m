function paths = build_butane_paths(datapath, iteration)
%BUILD_BUTANE_PATHS Construct full file paths for Butane artifacts.
%   PATHS = BUILD_BUTANE_PATHS(DATAPATH) returns a struct whose fields map
%   artifact names to their corresponding files under DATAPATH. The fields
%   include:
%       chart, chart_part, chart_plot, TranM, TranM_plot,
%       chart_relearn, well, FPT, diary, diary_FPT, Butane_analysis.
%
%   PATHS = BUILD_BUTANE_PATHS(DATAPATH, ITERATION) appends the ITERATION
%   identifier to the iteration-specific artifacts: chart, chart_part,
%   TranM, diary, and FPT. All other artifact paths remain unchanged so the
%   function can be used for both per-iteration and global files.
%
%   ITERATION can be numeric, char, or string. Empty ITERATION values are
%   treated the same as omitting the argument.

if nargin < 2 || isempty(iteration)
    iterationSuffix = '';
else
    if isnumeric(iteration)
        if ~isscalar(iteration)
            error('Iteration must be a scalar value when numeric.');
        end
        iterationSuffix = num2str(iteration);
    else
        iterationSuffix = char(iteration);
    end
end

iterationFields = struct( ...
    'chart',           true, ...
    'chart_part',      true, ...
    'chart_plot',      false, ...
    'TranM',           true, ...
    'TranM_plot',      false, ...
    'chart_relearn',   false, ...
    'well',            false, ...
    'FPT',             true, ...
    'diary',           true, ...
    'diary_FPT',       false, ...
    'Butane_analysis', false ...
);

paths = struct();
paths.chart           = makePath('chart', '.mat',         iterationFields.chart);
paths.chart_part      = makePath('chart_part', '.mat',    iterationFields.chart_part);
paths.chart_plot      = makePath('chart_plot', '.mat',    iterationFields.chart_plot);
paths.TranM           = makePath('TranM', '.mat',         iterationFields.TranM);
paths.TranM_plot      = makePath('TranM_plot', '.mat',    iterationFields.TranM_plot);
paths.chart_relearn   = makePath('chart_relearn', '.mat', iterationFields.chart_relearn);
paths.well            = makePath('well', '.mat',          iterationFields.well);
paths.FPT             = makePath('Butane_FPT', '.mat',    iterationFields.FPT);
paths.diary           = makePath('diary', '.txt',         iterationFields.diary);
paths.diary_FPT       = makePath('diary_FPT', '.txt',     iterationFields.diary_FPT);
paths.Butane_analysis = makePath('Butane_analysis', '.mat', iterationFields.Butane_analysis);

    function filepath = makePath(baseName, extension, useIteration)
        suffix = '';
        if useIteration
            suffix = iterationSuffix;
        end
        filepath = [datapath, baseName, suffix, extension];
    end
end
