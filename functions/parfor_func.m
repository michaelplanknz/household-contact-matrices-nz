function [results_Prem, results_HH10, results_ABM] = parfor_func(a_h, a_n, Ph5, C10, Pn5, Pn10, par5, par10, tblRaw, popSize10, nReps, Alpha )

%  Contents of parfor loop in a function to ensure correct parallelisation

[nAgeGroupsHH, ~] = size(C10);
maxHHsize = max(tblRaw.totHHsize);

% Run compartment-based model with Prem matrix
results_Prem = getODEResults(  a_h*Ph5, a_n*Pn5, par5 );
results_Prem.I = results_Prem.I(:, 1:2:15) + results_Prem.I(:, 2:2:16);     % aggregate results from 5-year bands into 10-year bands 
nInf = results_Prem.finalSize.*par5.popSize';
results_Prem.finalSize = (nInf(1:2:15)+nInf(2:2:16)) ./ par10.popSize';       

% Run compartment-based model with household composotion matrix on 10-year age groups
results_HH10 = getODEResults(  a_h*C10, a_n*Pn10, par10 );

% Run realisations of ABM
if nReps > 0
    fs_age = zeros(nReps, nAgeGroupsHH);
    fs_size = zeros(nReps, maxHHsize);
    for iRep = 1:nReps
        tblInd = tblRaw;
        tblInd.Count = imputeHouseholdData(tblRaw, popSize10, 0);
        synPopInd = makeSynPop(tblInd);
        resultsTemp  = getABMResults( synPopInd, a_h, a_n*Pn10, par10 );
        fs_age(iRep, :) = resultsTemp.finalSize;
        fs_size(iRep, :) = resultsTemp.finalSize_byHH;
    end

    % Save the detailed results structure of the ABM for a single realisation
    results_ABM = resultsTemp;
    
    % For final size results, append the summary quantiles across all realistions
    results_ABM.finalSize_med = median(fs_age, 1);
    results_ABM.finalSize_lo  = quantile(fs_age, Alpha/2, 1);
    results_ABM.finalSize_hi = quantile(fs_age, 1-Alpha/2, 1);
    results_ABM.finalSize_byHH_med  = median(fs_size, 1);
    results_ABM.finalSize_byHH_lo  = quantile(fs_size, Alpha/2, 1);
    results_ABM.finalSize_byHH_hi  = quantile(fs_size, 1-Alpha/2, 1);
else
    results_ABM = nan;
end