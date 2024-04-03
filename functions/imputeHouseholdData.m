function x = imputeHouseholdData(tbl, targetPop)

%fprintf('Running imputation... ')

minFreq = 4;


maxIter = 1e6;
fTol = 0.01;

x = tbl.Count;
nHouseTypes = length(x);
nToAdd = 1;


popSizeCurrent = sum(tbl.HHfreq.*x, 1)';
ii = 1;
done = false;
while ~done && ii <= maxIter
    popDf = targetPop - popSizeCurrent;
    ck = (tbl.HHfreq * popDf);                               % cosine angle between the pop discrepancy vector and the age composition vector for each household type
    

     weights = tbl.Count.*abs(ck);
     k = randsample(nHouseTypes, nToAdd, true, weights); 

    x(k) = max(minFreq, x(k)+sign(weights(k)) );

    

    popSizeCurrent = sum(tbl.HHfreq.*x, 1)';

    %fprintf('i = %i, total pop = %i \n', ii, sum(popSizeCurrent));

    done = norm(popSizeCurrent - targetPop)/norm(targetPop) < fTol ;
    ii = ii+1;

end

Dx = x-tbl.Count;
nAdded = sum(tbl.totHHsize(Dx > 0).*Dx(Dx > 0));
nHousesAdded = sum(Dx(Dx > 0));
nTypesAdded = sum(Dx > 0);
nRemoved = -sum(tbl.totHHsize(Dx < 0).*Dx(Dx < 0));
nHousesRemoved = -sum(Dx(Dx < 0));
nTypesRemoved = sum(Dx < 0);

%fprintf('complete with %i individuals in %i households of %i distinct types added and %i individuals in %i households of %i distinct types removed\n', nAdded, nHousesAdded, nTypesAdded, nRemoved, nHousesRemoved, nTypesRemoved)







