function tbl = getHouseholdData(fNameData)

fprintf('Reading raw houshold composition data... ')

opts = detectImportOptions(fNameData);

% Specify variable names and types
ageGroupVarNames = {'V1_0_9_years', 'V2_10_19_years', 'V3_20_29_years', 'V4_30_39_years', 'V5_40_49_years', 'V6_50_59_years', 'V7_60_69_years', 'V8_70__years'};        % column names in the household compposition data 
nAgeGroupsHH = length(ageGroupVarNames);
opts = setvartype(opts, ageGroupVarNames, {'double'});

% Read table
tbl = readtable(fNameData, opts);

% Discard TOTALS row
tbl = tbl(1:end-1, :);

%  Make the first nAgeGroups columns of tbl into a matrix representing frequency of each group within a given household
tbl.HHfreq = table2array(tbl(:, 1:nAgeGroupsHH));       

% Remove the columns that went into HHfreq: 
tbl = removevars(tbl, ageGroupVarNames);

% Add total household size variable
tbl.totHHsize = sum(tbl.HHfreq, 2);

fprintf('done\n')