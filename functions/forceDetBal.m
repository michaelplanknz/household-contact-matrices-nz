function C_detBal = forceDetBal(C, popSize)

% Function to force a contact matrix C to satisfy the detailed balance equations with respect to a specified population size vector popSize
% Note popSize should be a *column* vector containing the population size in each group

C_detBal = 0.5 * (C + (popSize'./popSize).*C');  
