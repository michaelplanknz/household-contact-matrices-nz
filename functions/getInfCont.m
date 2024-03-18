function contIDs = getInfCont(synPop)

% Return indices of susceptible household contacts of currently infectious individuals
% Those in a household with n infectious people will have their index listed n times
% Includes the infectious individuals themselves

contIDs = cat(2, synPop.neighbourIDs{synPop.state == 2})';

