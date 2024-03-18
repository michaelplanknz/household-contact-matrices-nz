clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of realisations of the synthetic pop and ABM to run:
nReps = 0;
Alpha = 0.05;           % Alpha level for confidence intervals

% Numerical parameters - range of values for ah and an
aRateHouse = (0.025:0.01:0.065)';        % daily infection rate for household contacts
aRateNonhouse = (0.025:0.01:0.065)';        % daily infection rate for non-household contacts

% Vectors representing the lower boundaries of 5-year and 10-year age groups
age5 = 0:5:75;
age10 = 0:10:70;

% Add directories containing data files and Matlab functions to the path
addpath('data')
addpath('functions')

% For reproducibility
rng(707349);

% Input data file names
fNameData = 'Table_1_rounded.csv';
fNamePremHome = 'nz_contacts_home.xlsx';
fNamePremSchool = 'nz_contacts_school.xlsx';
fNamePremWork = 'nz_contacts_work.xlsx';
fNamePremOther = 'nz_contacts_other.xlsx';
fNamePop = 'pop_size_by_age_2018census_usually_resident.csv';

% Output folder and file names
savFlag = 0;                % set to 1 to write outputs or 0 otherwise
savFolder = "results/";
fNameContMat = "household_matrix";
fNamePopSize = "popsize";
fNameSynPop = "synpop";




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read in Prem matrices for NZ 
P5_home = readmatrix(fNamePremHome);
P5_school = readmatrix(fNamePremSchool);
P5_work = readmatrix(fNamePremWork);
P5_other = readmatrix(fNamePremOther);
[nAgeGroupsMat, ~] = size(P5_home);

% Read in pop size data
popData = readtable(fNamePop);

% Make a vector of pop sizes by age (5 year bands & 10 year bands)
popSizeCensus5 = popData.Population;
popSizeCensus5 = [popSizeCensus5(1:nAgeGroupsMat-1); sum(popSizeCensus5(nAgeGroupsMat:end))];                            % Merge last few age groups into one if required
popSizeCensus10 = popSizeCensus5(1:2:end-1)+popSizeCensus5(2:2:end);

% Read household data
tblRaw = getHouseholdData(fNameData);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make contact matrix from the raw data
[Craw, popSizeRaw] = makeHouseholdContMat(tblRaw);

% Impute data for missing part of pop by adjusting frequency of household types
x = imputeHouseholdData(tblRaw, popSizeCensus10);
tbl = tblRaw;
tbl.Count = x;

% Make contact matrix from imputed data
[C10, popSize10] = makeHouseholdContMat(tbl);

% Map the represented population in 10 year age bands to 5 year age bands
scaleFactor = repelem(popSize10./popSizeCensus10, 2);
popSize5 = round(popSizeCensus5 .* scaleFactor);

% Version of Prem matrices satisfying detailed balance equations (for the represented HH pop)
Ph5 = forceDetBal(P5_home, popSize5);
Pn5 = forceDetBal(P5_school+P5_work+P5_other, popSize5);

% Map contact matrix derived from 10-year household composition data onto a finer (5 year age bands) matrix using the Prem matrix to weight contributions within finer age groups
% (function also outputs a coarse version of the Prem matrix, Ph10)
C5 = mapCoarseFine(C10, Ph5, popSize5);

% Get a coarse version (Pn10) of the Prem non-household contact matrix to use in the ABM
[~, Pn10, ~] = mapCoarseFine(C10, Pn5, popSize5);

% Calculate dominant eigenvectors of various home contact matrices
[evCraw, ~] = getDomEV(Craw');
[evC5, dC5] = getDomEV(C5');
[evC10, ~] = getDomEV(C10');
[evPh5, dPh5] = getDomEV(Ph5');


% Scale Prem matrices so the dominant eigvenvalue of Prem household matrix matches that of the household composition data matrix
Ph5 = Ph5* (dC5/dPh5);
Pn5 = Pn5 * (dC5/dPh5);
Pn10 = Pn10 * (dC5/dPh5);

% Make synthetic population for simulating the ABM
synPop = makeSynPop(tbl);





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run epidemic models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running epidemics\n')

% Get parameter values and append the pop size in 5-year or 10-year age bands to the parameter structure
par5 = getPar();
par10 = par5;
par5.popSize = popSize5;
par10.popSize = popSize10;

nChange1 = length(aRateHouse);
nChange2 = length(aRateNonhouse);
aRateHouseList = repmat(aRateHouse, nChange2, 1);
aRateNonhouseList = repelem(aRateNonhouse, nChange1, 1);

[~, nAgeGroupsHH] = size(tbl.HHfreq);
parfor iChangeAll = 1:nChange1*nChange2
    [results_Prem(iChangeAll), results_HH5(iChangeAll), results_HH10(iChangeAll), results_ABM(iChangeAll)] = parfor_func(aRateHouseList(iChangeAll), aRateNonhouseList(iChangeAll), Ph5, C5, C10, Pn5, Pn10, par5, par10, tblRaw, popSizeCensus10, nReps, Alpha );
    fprintf('   done %i/%i\n', iChangeAll, nChange1*nChange2)
end
% Reshape results sturcutre arrays into 2D arrays corresponding to the two varied parameters:
results_Prem = reshape(results_Prem, nChange2, nChange1)';
results_HH5  = reshape(results_HH5, nChange2, nChange1)';
results_HH10 = reshape(results_HH10, nChange2, nChange1)';
results_ABM  = reshape(results_ABM, nChange2, nChange1)';






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save outputs and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a table with population size after imputation to save alongside
% household contact matrices
popSizeSav = table(age5', popSize5, 'VariableNames', {'age', 'population'});

% Make a table with household size distribution data (raw and imputed)
HHsizeDistRaw = groupsummary(tblRaw, "totHHsize", "sum", "Count");        
HHsizeDistRaw.p = HHsizeDistRaw.sum_Count/sum(HHsizeDistRaw.sum_Count);
HHsizeDist = groupsummary(tbl, "totHHsize", "sum", "Count");        
HHsizeDist.p = HHsizeDist.sum_Count/sum(HHsizeDist.sum_Count);


% Save contact matrices and synthetic populatoin as .csv files and all
% results in .mat file
if savFlag
    writematrix(C10, savFolder+fNameContMat+"_10.csv");
    writematrix(C5, savFolder+fNameContMat+"_5.csv");
    writetable(popSizeSav, savFolder+fNamePopSize+".csv");
    writetable(synPop, savFolder+fNameSynPop+".csv");
    save(savFolder+"results.mat");
end



% Plot graphs
h = figure(1);
h.Position =  [560   296   856   652];
subplot(2, 2, 1);
plot(age10+5, popSizeRaw/1e6, 'o-', age10+5, popSize10/1e6, 'o-', age10+5, popSizeCensus10/1e6, 'o-')
xlabel('age (years)')
ylabel('population size (millions)')
legend('raw', 'imputed', 'census pop', 'Location', 'southwest')
title('(a)')
subplot(2, 2, 2);
plot(HHsizeDistRaw.totHHsize, HHsizeDistRaw.p, 'o-', HHsizeDist.totHHsize, HHsizeDist.p, 'o-')
xlabel('household size')
ylabel('proportion of households')
legend('raw','imputed')
title('(b)')
subplot(2, 2, 3)
histogram(log10(tblRaw.Count));
hold on
histogram(log10(tbl.Count));
h = gca;
h.XTickLabel = 10.^(1:5);
xlabel('frequency of household type')
ylabel('number of household types')
legend('raw', 'imputed')
title('(c)')
subplot(2, 2, 4);
plot(age10+5, sum(Craw, 2), 'o-', age10+5, sum(C10, 2), 'o-', age5+2.5, sum(C5, 2), 'o-', age5+2.5, sum(Ph5, 2), 'o-')
xlabel('age (years)')
ylabel('avg household contacts per person')
legend('raw(10)', 'imputed(10)', 'imputed(5)', 'Prem(5)', 'Location', 'southwest')
ylim([0 inf])
title('(d)')
if savFlag
    saveas(h, savFolder+'Fig1.png');
end



h = figure(2);
tl = tiledlayout(2, 2, 'Padding', 'Compact');
nexttile;
imagesc(age10, age10, Craw); 
clim([0 1.2])
title('(a) raw')
nexttile;
imagesc(age5, age5, C10);
clim([0 1.2])
title('(b) imputed')
nexttile;
imagesc(age10, age10, C5);
clim([0 1.2])
title('(c) imputed(5)')
nexttile;
imagesc(age5, age5, Ph5); 
clim([0 1.2])
title('(d) Prem')
cb = colorbar;
cb.Layout.Tile = 'east';
xlabel(tl, 'age of contact (years)')
ylabel(tl, 'age of individual (years)')
if savFlag
    saveas(h, savFolder+'Fig2.png');
end



colOrd = colororder;
colOrd = [colOrd; [0 0 0]];

h = figure(3);
h.Position = [    89          87        1040         763];
set(gcf, 'DefaultAxesColorOrder', colOrd);
set(gcf, 'DefaultAxesLineStyleOrder', {'-','--',':'});
subplot(2, 2, 1)
plot(results_Prem(1, 1).t, results_Prem(1, 1).I./par5.popSize')
ylim([0 0.02])
xlabel('t (days)')
ylabel('proportion infectious')
l = legend(string(0:5:75));
l.Position = [    0.3897    0.5974    0.0625    0.3506];
title('(a) Prem-ODE(5)')
subplot(2, 2, 2)
plot(results_HH5(1, 1).t, results_HH5(1, 1).I./par5.popSize')
ylim([0 0.02])
xlabel('t (days)')
ylabel('proportion infectious')
title('(b) Household-ODE(5)')
subplot(2, 2, 3)
set(gcf, 'DefaultAxesColorOrder', colOrd(1:2:end, :));
set(gcf, 'DefaultAxesLineStyleOrder', {'-','--',':'});
plot(results_HH10(1, 1).t, results_HH10(1, 1).I./par10.popSize')
ylim([0 0.02])
xlabel('t (days)')
ylabel('proportion infectious')
legend(string(0:10:70))
title('(b) Household-ODE(10)')
subplot(2, 2, 4)
plot(results_ABM(1, 1).t, results_ABM(1, 1).I'./par10.popSize )
ylim([0 0.02])
xlabel('t (days)')
ylabel('proportion infectious')
title('(d) Household-ABM(10)')
if savFlag
    saveas(h, savFolder+'Fig3.png');
end

R0_Prem = reshape([results_Prem.R0], nChange1, nChange2);
R0_HH5 = reshape([results_HH5.R0], nChange1, nChange2);
R0_HH10 = reshape([results_HH10.R0], nChange1, nChange2);

h = figure(4);
h.Position = [ 1          41        1680         933];
tl = tiledlayout(nChange2, nChange1, 'Padding', 'Compact');
for iChange2 = 1:nChange2
    for iChange1 = 1:nChange1
        nexttile
        plot(age5+2.5, results_Prem(iChange1, iChange2).finalSize, '.-', age5+2.5, results_HH5(iChange1, iChange2).finalSize, '.-', age10+5, results_HH10(iChange1, iChange2).finalSize, '.-', age10+5, results_ABM(iChange1, iChange2).finalSize, '.-')
        grid on
        if iChange1 == 1 && iChange2 == 1
            legend('Prem-ODE(5)', 'Household-ODE(5)', 'Household-ODE(10)', 'Household-ABM(10)')
        end
       ylim([0 1])
       xlabel(tl, 'age (years)');
       ylabel(tl, 'proportion infected')
       subtitle(sprintf('R_0 = %.2f ', R0_HH10(iChange1, iChange2) ));
       if iChange2 == 1
            title(sprintf('a_h = %.3f', aRateHouse(iChange1)));
       end
       if iChange1 == 1
            text(-20, 0.25, sprintf('a_n = %.3f', aRateNonhouse(iChange2)), 'FontWeight', 'bold', 'Rotation', 90)
       end
    end
end
if savFlag
    saveas(h, savFolder+'Fig4.png');
end


h = figure(5);
h.Position = [ 1          41        1680         933];
tl = tiledlayout(nChange2, nChange1, 'Padding', 'Compact');
for iChange2 = 1:nChange2
    for iChange1 = 1:nChange1
       nexttile
       med = results_ABM(iChange1, iChange2).finalSize_med;
       neg = results_ABM(iChange1, iChange2).finalSize_med - results_ABM(iChange1, iChange2).finalSize_lo;
       pos = results_ABM(iChange1, iChange2).finalSize_hi - results_ABM(iChange1, iChange2).finalSize_med;
       errorbar(age10+5, med, neg, pos, '.-')
       grid on
       ylim([0 1])
       xlabel(tl, 'age (years)');
       ylabel(tl, 'proportion infected')
       subtitle(sprintf('R_0 = %.2f ', R0_HH10(iChange1, iChange2) ));
       if iChange2 == 1
            title(sprintf('a_h = %.3f', aRateHouse(iChange1)));
       end
       if iChange1 == 1
            text(-20, 0.25, sprintf('a_n = %.3f', aRateNonhouse(iChange2)), 'FontWeight', 'bold', 'Rotation', 90)
       end
    end
end
if savFlag
    saveas(h, savFolder+'FigS1.png');
end


maxHHsize = max(tbl.totHHsize);
h = figure(6);
h.Position = [ 1          41        1680         933];
tl = tiledlayout(nChange2, nChange1, 'Padding', 'Compact');
for iChange2 = 1:nChange2
    for iChange1 = 1:nChange1
       nexttile
       med = results_ABM(iChange1, iChange2).finalSize_byHH_med;
       neg = results_ABM(iChange1, iChange2).finalSize_byHH_med - results_ABM(iChange1, iChange2).finalSize_byHH_lo;
       pos = results_ABM(iChange1, iChange2).finalSize_byHH_hi - results_ABM(iChange1, iChange2).finalSize_byHH_med;
       errorbar(1:maxHHsize, med, neg, pos, '.-')
       grid on
       ylim([0 1])
       xlim([1 maxHHsize])
       xlabel(tl, 'household size');
       ylabel(tl, 'proportion infected')
       subtitle(sprintf('R_0 = %.2f ', R0_HH10(iChange1, iChange2) ));
       if iChange2 == 1
            title(sprintf('a_h = %.3f', aRateHouse(iChange1)));
       end
       if iChange1 == 1
            text(-2.8, 0.25, sprintf('a_n = %.3f', aRateNonhouse(iChange2)), 'FontWeight', 'bold', 'Rotation', 90)
       end
    end
end
if savFlag
    saveas(h, savFolder+'FigS2.png');
end


