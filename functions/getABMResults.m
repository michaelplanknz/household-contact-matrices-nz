function results = getABMResults(synPop, ah, Mn, par)   

% Run agent-based model using synthetic pop
[t, ~, ~, I, R, synPopOut] = runABM(synPop, ah, Mn, par);

% Save required results in output structure
results.t = t;
results.I = I;

% Final epidemic size in each age group:
results.finalSize = R(end, :)./par.popSize';  

maxHHsize = max(synPop.houseSize);
HHsize = 1:maxHHsize;       % vector of household sizes
% Make a binary matrix whose (i,j) element is 1 iff individual i is in a household size of size j:
HHsize_indicator = (synPopOut.houseSize == HHsize);
% Make a binary matrix whose (i,j) element is 1 iff individual i is in a household size of size j and individual i got infected:
HHsize_inf_indicator = (synPopOut.state > 0).*(synPopOut.houseSize == HHsize);

% Final epidemic size by household size:
results.finalSize_byHH = sum(HHsize_inf_indicator)./sum(HHsize_indicator);


