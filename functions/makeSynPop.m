function [synPop, gr] = makeSynPop(tbl)

%fprintf('Making synthetic population... ')

HHsizeDist = groupsummary(tbl, "totHHsize", "sum", "Count");  


[nHouseTypes, nAgeGroups] = size(tbl.HHfreq);
popSize = sum(tbl.HHfreq.*tbl.Count)'; 
popSizeTot = sum(popSize);


synPop.personID = (1:popSizeTot)';
synPop.ageGroup = nan(popSizeTot, 1);
synPop.houseID = nan(popSizeTot, 1);
synPop.houseSize = nan(popSizeTot, 1);
synPop.neighbourIDs = cell(popSizeTot, 1);

nEdges = sum( HHsizeDist.totHHsize.*(HHsizeDist.totHHsize-1)/2 .* HHsizeDist.sum_Count);
EndNodes = nan(nEdges, 2);

edgeCount = 1;
houseCount = 1;
indCount = 1;
for iHouseType = 1:nHouseTypes
    nHouses = tbl.Count(iHouseType);        % number of houses of this type

    houseIDs = houseCount:houseCount+nHouses-1;       % vector of household IDs for household of the current type
    ageGroup1 = repelem( (1:nAgeGroups)', tbl.HHfreq(iHouseType, :) );  % vector of age groups of people in the current household type
    houseSize = length(ageGroup1);                                  % total household size of the current type

    nNew = houseSize*nHouses;                         % total number of people in households on the current type

    % columns for entry into synPop table
    synPop.ageGroup(indCount:indCount+nNew-1) = repmat(ageGroup1, nHouses, 1);
    synPop.houseID(indCount:indCount+nNew-1) = repelem(houseIDs, houseSize);
    synPop.houseSize(indCount:indCount+nNew-1) = houseSize;

    % Make a lookup table of the IDs of household members (neighbours) to improve efficiency of epidemic simulation algorithm
    neighbourIDmat1 = 0:houseSize-1;     % list of members of household of this type with person IDs 0:houseSize-1
    neighbourIDmat = repmat(neighbourIDmat1, nHouses*houseSize, 1) + indCount + repelem((0:houseSize:nHouses*houseSize-1)', houseSize, 1);      % matrix of neighbour IDs for all individuals being added at this stage
    synPop.neighbourIDs(indCount:indCount+nNew-1, 1) = num2cell(neighbourIDmat, 2);                                                              % store neighbour IDs in synPop as a cell array (due to different number of neighbours for each individual)

    if houseSize > 1    
        edges1 = nchoosek( 0:houseSize-1, 2);      % list of edges for a single household of this type with person IDs 0:houseSize-1
        nEdges1 = houseSize*(houseSize-1)/2;
        offset = repelem((indCount:houseSize:indCount+(nHouses-1)*houseSize)', nEdges1);
        edges = repmat(edges1, nHouses, 1) + offset;      % repeat edge matrix with an offset for the person IDs 
        nEdges = nEdges1*nHouses;
        EndNodes(edgeCount:edgeCount+nEdges-1, :) = edges;
        edgeCount = edgeCount+nEdges;
    end
    
    houseCount = houseCount+nHouses;
    indCount = indCount+nNew;
end

nHouses = houseCount - 1;

synPop = struct2table(synPop);

% Make a graph object wthat is a collectiopn of complete graphs, one for each household, that are disconnnected from one another
EdgeTable = table(EndNodes);
gr = graph(EdgeTable, synPop);

%fprintf('done: %i individuals in %i households\n', height(synPop), nHouses )


