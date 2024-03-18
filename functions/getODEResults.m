function results = getODEResults(Mh, Mn, par)        

% overal contact matrix
M = Mh+Mn;

% Get dominant eigenvalue and eigenvector
[v, l] = getDomEV(M'/par.Mu);

% calculate pseudo-eigenvalues (projection of NGM_h*v onto v where v is the dominant eigenvector of NGM_tot) of the household andn non-household contributions
results.R0 = l;
results.Rh = (Mh'*v)' * v / par.Mu / norm(v)^2;
results.Rn = (Mn'*v)' * v / par.Mu / norm(v)^2;

% Run ODE model with the specified overall contact matrix and parameter set
[t, ~, ~, I, R] = runODEModel(M, par);

% Save required results in output structure
results.t = t;
results.I = I;        

% Final epidemic size in each age group
results.finalSize = R(end, :)./par.popSize';
