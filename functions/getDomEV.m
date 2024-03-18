function [v, a] = getDomEV(A)

% Get the dominant eigenvalue a and corresponding eigenvector v of a matrix A 
% Dominant means largest real part 

[vu, a] = eigs(A, 1, 'largestreal');
v = vu/sum(vu);





