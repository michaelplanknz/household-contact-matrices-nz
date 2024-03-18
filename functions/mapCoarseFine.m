function [Afine, Bcoarse, popSizeCoarse] = mapCoarseFine(Acoarse, Bfine, popSizeFine)

% function to map a representation of contact matrix on coarse age groups onto fine age groups using a corresponding fine-scale representation

popSizeCoarse = popSizeFine(1:2:end-1)+popSizeFine(2:2:end);

% start by constructing a coarse verison of B by combining element in 2x2 blocks
Bw = popSizeFine.*Bfine;                % multiply each row of B by the corresponding age dist
Bt = Bw(:, 1:2:end-1)+Bw(:, 2:2:end);   % combine elements in 2x2 blocks
Bt2 = Bt(1:2:end-1, :)+Bt(2:2:end, :);
Bcoarse = Bt2 ./ popSizeCoarse;

% Now simply scale each element of Bfine by the corresponding ratio of Acoarse to Bcoarse
scaleFactor = repelem(Acoarse./Bcoarse, 2, 2);
Afine = Bfine.*scaleFactor;




