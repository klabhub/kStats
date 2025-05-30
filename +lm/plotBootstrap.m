function plotBootstrap(results)
% Plot the results of a bootstrap analysis.
% 
% The input struct is the output of the lm.bootstrap function.
% See Also lm.bootstrap
arguments
    results (1,1) struct
end

FE = results.m.fixedEffects;
nrFixedEffects =numel(FE);
clf;
layout=  tiledlayout("flow");

% Show a histogram for the estimates of each fixed effect, with 
% a thick line representing the original model estimate.
for f=1:nrFixedEffects
    ax(f) = nexttile;
    histogram(results.fe.all(f,:),'Normalization','probability');
    hold on
    plot(FE(f)*[1 1],ylim,'k','LineWidth',2);
    legStr = ["Simulated Sets" "Orignal FE"];
    if results.pv.mode =="TYPE-I"
        plot(results.fe.simulated(f)*[1 1],ylim,'r','LineWidth',2);
        legStr = [legStr "Null FE"];            %#ok<AGROW>
    end

    plot([0 0],ylim,'k','LineWidth',0.5)
    title (sprintf('%s: %.3G CI [%.3G %.3G]',results.m.CoefficientNames{f},results.fe.mean(f),results.fe.ci(f,1),results.fe.ci(f,2)) ,"Interpreter","none");
    xlabel 'Coefficient'
    ylabel 'Probability'
    xlim(max(abs(xlim))*[-1 1])
    legend(legStr)
end

% Show the log likelihood distrivution
ax(nrFixedEffects+1) = nexttile;ax(nrFixedEffects+1);
histogram(results.ll,'Normalization','probability')
hold on
plot(results.m.ModelCriterion.LogLikelihood*[1 1],ylim,'k','LineWidth',2);
title (sprintf('%s: %.3G CI [%.3G %.3G]','Log Likelihood:',mean(results.ll,"omitnan"),prctile(results.ll,2.5),prctile(results.ll,97.5)) ,"Interpreter","none");
xlabel 'Log Likehood'
ylabel 'Probability'

% Show a panel with the residuals and the estimated distributions
if ismember(results.pv.mode,["TYPE-I" "TYPE-II"])
    % Show the residuals and how they are fit by the kernel density
    ax(nrFixedEffects+2) =nexttile;
    maxResidual = prctile(abs(results.m.residuals),97.5);
    nrBins = min(20,round(numel(results.m.residuals)/10));
    x= linspace(-maxResidual,maxResidual,nrBins);
    hold on
    R= results.m.residuals;
    % Match histogram and density colors. Histogram handle does not
    % conatin the color that was used (only 'auto')
    axColors = ax(nrFixedEffects+2).ColorOrder;
    nrColors = size(axColors,1);
    nrUHeteroBins = numel(results.uHeteroBins);
    for b=1:nrUHeteroBins
        thisGroup = results.responseGroupingIx==results.uHeteroBins(b);
        h = histogram(R(thisGroup),x,'Normalization','Probability');
        thisColor = axColors(mod(b-1,nrColors)+1,:);
        h.FaceColor = thisColor;
        thisPdf =pdf(results.noiseDistribution{b},x);
        plot(x,thisPdf./sum(thisPdf),'LineWidth',2,'color',thisColor)
    end
    xlim(maxResidual*[-1 1])
    xlabel 'Residual'
    ylabel 'Probability'
    legend('Residuals','Kernel Density Estimate')
end

% Generate a title as a summary
switch (results.pv.mode)
    case 'TYPE-I'
        % Show type-I error probability per fixed effect
        str = strjoin(string(results.m.CoefficientNames)' + ": p=" + string(results.fe.pTypeI),' ');
    case 'TYPE-II'
        % Show type-II  error probability per fixed effect
        str = strjoin(string(results.m.CoefficientNames)' + ": p=" + string(results.fe.pTypeII),' ');
    case 'RESAMPLE'
        % Show stdev across samples as a percentage of the mean fixed effect
        str = strjoin(string(results.m.CoefficientNames)' + ": std =" + string(round(100*results.fe.std./abs(results.fe.mean))) + "%",' ');
end
title(layout,[results.pv.mode + " analysis" ; str],'Interpreter','none')
end