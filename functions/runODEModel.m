function [t, S, E, I, R] = runODEModel(C, par)

% Function to run the epidemic ODE model and return the vector of times and
% matrix of solution values
% C is the average infections per day matrix 
% such that C_ij/Mu where 1/Mu is the average infectious period is the avg number of people in age group i infected by one
% case in age group j in a fully susceptible pop


[par.nAgeGroups, ~] = size(C);

% output solution at daily time points up to tMax
tSpan = 0:1:par.tMax;

% initial condition
E0 = par.Efrac*par.popSize;
S0 = par.popSize - E0;
I0 = zeros(par.nAgeGroups, 1);
IC = [S0; E0; I0];

% Call ODE solver
[t, Y] = ode45(@(t, x)mySEIRmodel(t, x, C, par), tSpan, IC);

S = Y(:, 1:par.nAgeGroups);
E = Y(:, par.nAgeGroups+1:2*par.nAgeGroups);
I = Y(:, 2*par.nAgeGroups+1:end);
R = par.popSize' - (S+E+I);



