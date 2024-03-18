function  [t, S, E, I, R, synPop] = runABM(synPop, Alpha, Cn, par)

% Run the ABM on a synthetic population defined in synPop with a household contact daily infection rate Alpha, non-home infectoin per unit time matrix Cn, and parameters par

popSizeTot = height(synPop);
[nAgeGroups, ~] = size(Cn);
assert(nAgeGroups == max(synPop.ageGroup))      % check number of age groups in synPop is the same as the dimension of the non-household contact matrix

nHouses = max(synPop.houseID);

% Record indices of individuals in each group to enable more efficient lookup
indAge = cell(nAgeGroups, 1);
for iAge = 1:nAgeGroups
     indAge{iAge} = find(synPop.ageGroup == iAge);
end


% Set up table variables and compartment total variables to record epidemic
synPop.state = zeros(popSizeTot, 1);        % 0 for S, 1 for E, 2 for I, 3 for R
synPop.tExp = nan(popSizeTot, 1);
synPop.tInf = nan(popSizeTot, 1);
synPop.tRec = nan(popSizeTot, 1);
S = zeros(par.tMax+1, nAgeGroups);
E = zeros(par.tMax+1, nAgeGroups);
I = zeros(par.tMax+1, nAgeGroups);
R = zeros(par.tMax+1, nAgeGroups);



% set up vector of time steps
t = 0:1:par.tMax;


% Choose initial exposed population at random and record compartment variables at t=0
indExp = randsample(popSizeTot, round(par.Efrac*popSizeTot), false);
synPop.state(indExp) = 1;
S(1, :) = histcounts(synPop.ageGroup(synPop.state == 0), 1:nAgeGroups+1); 
E(1, :) = histcounts(synPop.ageGroup(synPop.state == 1), 1:nAgeGroups+1); 
I(1, :) = histcounts(synPop.ageGroup(synPop.state == 2), 1:nAgeGroups+1); 
R(1, :) = histcounts(synPop.ageGroup(synPop.state == 3), 1:nAgeGroups+1); 

iStep = 1;
while iStep <= par.tMax && sum(E(iStep, :)+I(iStep, :), 2) > 0

    % Simulate household infections, with probability Alpha for each susceptible contact
    % (If there are multiple active houses in a household, the household members will be listed multiple times in contIDs and will become infected if at least one of the infection 'attempts' is successful
    contIDs = getInfCont(synPop);
    contInfIDs = contIDs( synPop.state(contIDs) == 0 & rand(length(contIDs), 1) < 1-exp(-Alpha));

    % Simulate non-household infections
    % Number of attempted infections by age group, aggregated across all source cases:
    expAttemptedInfections = sum(Cn(synPop.ageGroup(synPop.state == 2) , :), 1);
    nAttemptedInfections = poissrnd(expAttemptedInfections);
    % Make a vector of the IDs of target infections in each age group by randomly sampling from all individuals in that age group
    indTarget = cell(nAgeGroups, 1);
    for iAge = 1:nAgeGroups
        indTarget{iAge} = randsample( indAge{iAge}, nAttemptedInfections(iAge) );
    end
    % concatenate into a single of vector of target IDs
    indTargetAll = cat(1, indTarget{1:nAgeGroups});
    % Any susceptible target IDs become infected
    outInfIDs = indTargetAll(synPop.state(indTargetAll) == 0);

    % Simulate which individuals will enter the infectius or the recovered state
    onsetFlag = synPop.state == 1 & rand(popSizeTot, 1) < par.Gamma;
    recoverFlag = synPop.state == 2 & rand(popSizeTot, 1) < par.Mu;

    % Update state variables and record event times 
    allInfIDs = union(contInfIDs, outInfIDs);
    synPop.state(allInfIDs) = 1;
    synPop.state(onsetFlag) = 2;
    synPop.state(recoverFlag) = 3;
    synPop.tExp(allInfIDs) = t(iStep+1);
    synPop.tInf(onsetFlag) =  t(iStep+1);
    synPop.tRec(recoverFlag) = t(iStep+1);

    % Calculate compatment totals for each age group at the current time step:
    S(iStep+1, :) = histcounts(synPop.ageGroup(synPop.state == 0), 1:nAgeGroups+1); 
    E(iStep+1, :) = histcounts(synPop.ageGroup(synPop.state == 1), 1:nAgeGroups+1); 
    I(iStep+1, :) = histcounts(synPop.ageGroup(synPop.state == 2), 1:nAgeGroups+1); 
    R(iStep+1, :) = histcounts(synPop.ageGroup(synPop.state == 3), 1:nAgeGroups+1); 

    iStep = iStep+1;

end

% If epidemic finished early, copy final states into remaining rows of the compartment matrices
nFill = par.tMax+1-iStep;
S(iStep+1:end, :) = repmat(S(iStep, :), nFill, 1);
E(iStep+1:end, :) = repmat(E(iStep, :), nFill, 1);
I(iStep+1:end, :) = repmat(I(iStep, :), nFill, 1);
R(iStep+1:end, :) = repmat(R(iStep, :), nFill, 1);



