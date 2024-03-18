function dxdt = mySEIRmodel(t, x, C, par)

% C is the average infections per day matrix 
% such that C_ij/Mu where 1/Mu is the average infectious period is the avg number of people in age group i infected by one
% case in age group j in a fully susceptible pop

xMat = reshape(x, par.nAgeGroups, 3);
S = xMat(:, 1);
E = xMat(:, 2);
I = xMat(:, 3);



dSdt = - (C') * I .* (S./par.popSize);
dEdt = (C') * I .* (S./par.popSize)  -  par.Gamma*E;
dIdt = par.Gamma*E - par.Mu*I;

dxdt = [dSdt; dEdt; dIdt];


