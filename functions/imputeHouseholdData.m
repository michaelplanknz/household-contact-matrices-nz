function x = imputeHouseholdData(tbl, targetPop, displayFlag)

% Run the impuation procedure on the household compsition data in table
% tbl, with target population vector targetPop
% The output x is the vector of household type counts (i.e. the new version of tbl.Count) after imputation

% Minimum allowed count for a household type
minFreq = 4;

% Maximum allowed number of iterations of the imputation process
maxIter = 1e6;

% Relative tolerance for discrepancy between househol population and target
% population 
fTol = 0.01;

% Count vector in the raw data
x = tbl.Count;
nHouseTypes = length(x);

% Number of households to add at each step of the imputation procedure
% (default 1)
nToAdd = 1;

% Initial population size by age according to the household composition data
popSizeCurrent = sum(tbl.HHfreq.*x, 1)';

ii = 1;
done = false;
% Main imputation iteration loop
while ~done && ii <= maxIter
    popDf = targetPop - popSizeCurrent;         % discrepancy vector between the current and target population size 
    ck = (tbl.HHfreq * popDf)./tbl.totHHsize;          % calculate ck for each household type according to how similar it is to the pop discrepancy vector (normalised by household size)

    weights = tbl.Count.*abs(ck);               % calculate sampling weights by multiplying ck by the count for each household type in the raw data tbl.Count

    % Sample a household type to add/remove
    k = randsample(nHouseTypes, nToAdd, true, weights); 
    x(k) = max(minFreq, x(k)+sign(weights(k)) );

    % Current population size by age according to the household composition data
    popSizeCurrent = sum(tbl.HHfreq.*x, 1)';

    % Check convergence criterion
    done = norm(popSizeCurrent - targetPop)/norm(targetPop) < fTol ;
    ii = ii+1;
end


if displayFlag
    % Calculate summary statisticvs on the numbert of individuals, households,
    % and distinct household types added/removed:
    Dx = x-tbl.Count;
    nAdded = sum(tbl.totHHsize(Dx > 0).*Dx(Dx > 0));
    nHousesAdded = sum(Dx(Dx > 0));
    nTypesAdded = sum(Dx > 0);
    nRemoved = -sum(tbl.totHHsize(Dx < 0).*Dx(Dx < 0));
    nHousesRemoved = -sum(Dx(Dx < 0));
    nTypesRemoved = sum(Dx < 0);

    fprintf('Imputation added %i individuals in %i households of %i distinct types and removed %i individuals in %i households of %i distinct types\n', nAdded, nHousesAdded, nTypesAdded, nRemoved, nHousesRemoved, nTypesRemoved)
end






