function [C, popSize] = makeHouseholdContMat(tbl)

% Make household contact matrix from raw household composition data in the table bl
% tbl.HHfreq shoule be a n x m array where n is the number of household types and m is the number of age groups
% The (i,j) element of tbl.HHfreq is the number of individuals of age group j in households on type i
% tbl.Counts should be a n x 1 vector containing the number of households of each type

[nRows, nAgeGroupsHH] = size(tbl.HHfreq);

% Find the total size of the population represented in the data in tbl
popSize = sum(tbl.HHfreq.*tbl.Count, 1)';       

fprintf('Making contact matrix from household composition data for %i individuals in %i houses... ', sum(popSize), sum(tbl.Count))

% Construct household contact matrix from household composition data
idM = eye(nAgeGroupsHH);
nCont = zeros(nAgeGroupsHH);
for iAge = 1:nAgeGroupsHH
   % This sums across households the number of contacts that people in age group i have with people in each other age group
   % The ith row in the sum is the product of count, the number of people in age group i in that household, and the number of contacts they have in each other age gorup (subtracting 1 from the ith column to avoid counting self as a contact)
   nCont(iAge, :) = sum( tbl.Count .* tbl.HHfreq(:, iAge) .* (tbl.HHfreq-repmat(idM(iAge, :), nRows, 1))  , 1);
end
% contacts per person is the total number of contacts i has with j (nCont_ij) divided by the total number of people in i (pop_i), so divide each row of matrix by the corresponding pop size
C = nCont./popSize;

fprintf('done\n')
