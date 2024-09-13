function [lmeWithRandomSlope,stats,reTable] = randomSlope(lme,pv)
% Add a slope for a factor as a random effect for a group, evaluate whether
% that random slope effect is significant, and evaluate whether that factor (at a given level)
% is significantly different for each element in the group.
% For example if the lme is 'output ~ input + (1|subject)' but I want
% to investigate whether adding a random slope effect input improves the
% model:
% randomSlope(lme,pv.factor= "input", pv.group = "subject")
% will then return the LMM for 'output ~ input + (1|subject) + (-1 + input|subject)'
% and return stats for the comparison of this model with the original
% (outpuf of linearmixedmodel/compare).
% To also generate a table of random effects and evaluate whether
% individual subjects had significant overall effects of input, specify
% the level of the factor at which you want to evaluate this. 
% 
% 
% Note
% The reTable always has the output from linearmixedmodel/randomEffects, which
% assesses whether a random effect is different from zero (pValue, tStat,
% DF). If a level is provided in the call to this function, this table also
% contains columns that determine whether the (total= fixed+random) slope effect 
%  at the given level was significantly different from zero for each level of the group
% (p,F,df1,df2) (Output of coefTes, using a 'REContrast' input).
%
% BK - June 2024
arguments (Input)
    lme (1,1)                       % The original LMM
    pv.factor (1,1) string          % The factor in the LMM for which to add a slope random effect
    pv.group (1,1) string           % Grouping variable for the random effect
    pv.level (1,1) = NaN            % Level of the factor at which to evaluate significance for each level of the grouping variable
    
    pv.plot (1,1) logical = false   % Show a histogram of effects
    pv.alpha (1,1) double = 0.05    % Alpha level (plot only)
end
arguments (Output)
    lmeWithRandomSlope (1,1)        % The model with the added slope
    stats table                     % Comparison of the extended model with the original (output of compare)    
    reTable table                   % Table with estimates of the random effects for the group                                                       
end

% Construct the new formula and compute the model
newFormula = sprintf([char(lme.Formula) ' + (-1 + %s|%s)'], pv.factor , pv.group);
lmeWithRandomSlope= fitlme(lme.Variables,newFormula,'Exclude',lme.ObservationInfo.Excluded,'DummyVarCoding',lm.dummyVarCoding(lme));
% Compare the new with the old model
if nargout>1 || pv.plot
    stats = dataset2table(compare(lme,lmeWithRandomSlope));
end
% Determine random effects
if nargout >2 || pv.plot
    [~,~,reTable] = randomEffects(lmeWithRandomSlope);
    reTable =convertvars(dataset2table(reTable),@iscellstr,"string");
end

% Evaluate significance of the factor at the specified level
if  (isstring(pv.level) && pv.level~="") || (isnumeric(pv.level) && ~isnan(pv.level))
    % Find the fixed effect
    if lme.VariableInfo{pv.factor,"IsCategorical"}
        level = pv.factor+ "_" + pv.level;
        isFe  = lmeWithRandomSlope.Coefficients.Name == level;
        fe = lmeWithRandomSlope.Coefficients.Estimate(isFe);
        feContrast = double(isFe);
    else
        isFe  = lmeWithRandomSlope.Coefficients.Name == pv.factor;
        fe = lmeWithRandomSlope.Coefficients.Estimate(isFe)*pv.level;
        feContrast = double(isFe)*pv.level;
    end
    levels = reTable.Level(reTable.Group==pv.group);
    nrLevels = numel(levels);
    tIndividual = table('size',[nrLevels 5],'VariableNames',{'Level','p','F','df1','df2'},'VariableTypes',{'string', 'double','double','double','double'});
    for lvl =1:nrLevels
        reContrast = double(reTable.Level == levels(lvl));
        tIndividual.Level(lvl) = levels(lvl);
        [tIndividual.p(lvl),tIndividual.F(lvl), tIndividual.df1(lvl),tIndividual.df2(lvl)] = coefTest(lmeWithRandomSlope,feContrast',0,'REContrast',reContrast');
    end
    reTable = reTable(reTable.Group==pv.group,:);  % Keep only the specifeid group
    reTable = innerjoin(reTable,tIndividual);
    individualSignificance = true;
    estimate= reTable.Estimate +fe;
else
    reTable = reTable(reTable.Group==pv.group,:);
    estimate= reTable.Estimate; % Just show the random effects below
    individualSignificance = false;
end

if pv.plot   
    [allC,edges]= histcounts(estimate); 
    centers = edges(1:end-1) + (edges(2:end)-edges(1:end-1))/2;
    bar(centers',allC','b');
    xlabel(pv.factor + " Random Effect")      
    str= sprintf('AIC %d (Chi^2 (%d) = %3.3f, p= %3.3g)',round(diff(stats.AIC)),stats.deltaDF(2),stats.LRStat(2),stats.pValue(2));
    if individualSignificance    
        isSignificant = (reTable.p<pv.alpha);
        sigC = histcounts(estimate(isSignificant),edges); 
        hold on
        h= bar(centers',sigC','r');
        h.FaceAlpha=0.8;        
        xlabel(pv.factor + " Fixed+Random Effect")
        legend('All',sprintf('p<%.3f',pv.alpha));
        str = [str sprintf('\n %d/%d (%.2f%%) significant',sum(isSignificant),nrLevels,100*mean(isSignificant))];
    end
    ylabel("#" + pv.group)
    title (str)
    hold off
end

