function x = imputeHouseholdData(tbl, targetPop)

%fprintf('Running imputation... ')

minFreq = 4;


maxIter = 1e6;
fTol = 0.01;

nM = tbl.HHfreq;
x = tbl.Count;
nHouseTypes = length(x);
nToAdd = 1;

%vn = vecnorm(nM, 2, 2);
popSizeCurrent = sum(nM.*x, 1)';
ii = 1;
done = false;
while ~done & ii <= maxIter
    popDf = targetPop - popSizeCurrent;
    %cosAngle = (nM * popDf) ./ (vn * norm(popDf) );       % cosine angle between the pop discrepancy vector and the age composition vector for each household type
    cosAngle = (nM * popDf);                               % cosine angle between the pop discrepancy vector and the age composition vector for each household type
    
    % sample a household type to add, weighted by (original frequency * (1+cosAngle))
%     weights = tbl.Count.*(1+cosAngle);
%     k = randsample(nHouseTypes, nToAdd, true, weights); 

     weights = tbl.Count.*abs(cosAngle);
     k = randsample(nHouseTypes, nToAdd, true, weights); 

    x(k) = max(minFreq, x(k)+sign(weights(k)) );

    

    popSizeCurrent = sum(nM.*x, 1)';

    %fprintf('i = %i, total pop = %i \n', ii, sum(popSizeCurrent));

    done = norm(popSizeCurrent - targetPop)/norm(targetPop) < fTol ;
    ii = ii+1;

end

Dx = x-tbl.Count;
nAdded = sum(tbl.totHHsize(Dx > 0).*Dx(Dx > 0));
nHousesAdded = sum(Dx(Dx > 0));
nTypesAdded = sum(Dx > 0);
nRemoved = -sum(tbl.totHHsize(Dx < 0).*Dx(Dx < 0));
nHousesRemoved = -sum(Dx(Dx < 0));
nTypesRemoved = sum(Dx < 0);

%fprintf('complete with %i individuals in %i households of %i distinct types added and %i individuals in %i households of %i distinct types removed\n', nAdded, nHousesAdded, nTypesAdded, nRemoved, nHousesRemoved, nTypesRemoved)







% Second option hard to get a balance between staying close to original x and apporixmating the target pop dist 


% Alpha = 5e-4;
% 
% nRows = height(tbl);
% 
% %opts=optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 1e6, 'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-11);
% opts = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 1e6);
% 
% myF1 = @(xx)( Alpha*norm(xx-tbl.Count) );
% myF2 = @(xx)( norm( sum(tbl.HHfreq.*xx, 1)'/sum(sum(tbl.HHfreq.*xx)) - targetPop/sum(targetPop)   ) );
% myF = @(xx)(myF1(xx)+myF2(xx));
% 
% x0 = tbl.Count.*(1+0.2*(rand(nRows, 1)-0.5));
% x = fmincon(myF, x0, [], [], [], [], minFreq*ones(size(x0)), [], [], opts);
% 
% fprintf('Initial    f = %.3e + %.3e = %.3e\n', myF1(tbl.Count), myF2(tbl.Count), myF(tbl.Count))
% fprintf('Randomized f = %.3e + %.3e = %.3e\n', myF1(x0), myF2(x0), myF(x0))
% fprintf('Final      f = %.3e + %.3e = %.3e\n', myF1(x), myF2(x), myF(x))
% fprintf('Difference from initial 2-norm = %.3e, inf-norm = %.3e\n', norm(x-tbl.Count), norm(x-tbl.Count, inf) );
% 
% popDists = [sum(tbl.HHfreq.*tbl.Count, 1)'/sum(sum(tbl.HHfreq.*tbl.Count)), sum(tbl.HHfreq.*x0, 1)'/sum(sum(tbl.HHfreq.*x0)), sum(tbl.HHfreq.*x, 1)'/sum(sum(tbl.HHfreq.*x)), targetPop/sum(targetPop)]
% 
% figure;
% subplot(1, 2, 1)
% plot(1:8, popDists, 'o-' )
% legend('original', 'randomized', 'final', 'target')
% subplot(1, 2, 2)
% histogram( log10(tbl.Count) )
% hold on
% histogram( log10(x) )
% 










% First option - move along a line that reduces pop discrepancyt faster
% Converges to correct pop dist but changes the distribution of household type frequencies (becomes more normal as opposed to right skewed)


% maxIter = 1e3;
% xTol = 1e-3;
% fTol = 10;
% 
% nM = tbl.HHfreq;
% x = tbl.Count;
% 
% opts = optimoptions('fmincon', 'display', 'off');
% 
% ii = 1;
% done = false;
% while ~done & ii <= maxIter
%     popSizeCurrent = sum(nM.*x, 1)';
%     popDf = targetPop - popSizeCurrent;
%     cosAngle = (nM * popDf) ./ (vecnorm(nM, 2, 2) * norm(popDf) );       % cosine angle between the pop discrepancy vector and the age composition vector for each household type
%     %dotProd = nM * popDf;       % cosine angle between the pop discrepancy vector and the age composition vector for each household type
%     
% 
%     objFn = @(dx)( norm(nM'*max(minFreq, x+dx.*cosAngle) - targetPop) );      % try to find the x such that changing f to f+c*dotProd gives the best match with the target pop
%     IC = norm(popDf)/norm(cosAngle);
%     dx = fmincon(objFn, IC, [], [], [], [], [], [], [], opts);
% 
%     xOld = x;
%     x = max(minFreq, x+dx.*cosAngle);
% 
%     done = norm(x-xOld)/norm(xOld) < xTol & norm(popDf) < fTol;
% 
%     ii = ii+1
% end
% 
% popDists = [sum(tbl.HHfreq.*tbl.Count, 1)'/sum(sum(tbl.HHfreq.*tbl.Count)), sum(tbl.HHfreq.*x, 1)'/sum(sum(tbl.HHfreq.*x)), targetPop/sum(targetPop)]
% 
% figure;
% subplot(1, 2, 1)
% plot(1:8, popDists, 'o-' )
% legend('original',  'final', 'target')
% subplot(1, 2, 2)
% histogram( log10(tbl.Count) )
% hold on
% histogram( log10(x) )
% 

