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
%
% To also generate a table of random effects and evaluate whether
% individual subjects had significant overall effects of input, specify
% the contrast that you want to evaluate.
% For instance
% randomSlope(lme,pv.factor= "input", pv.group = "subject",contrast =
% {"A","B"})
% Will determine the contrast between input ="A" and input ="B" that
% includes random + fixed effects, for each subject (And assess
% significance for each subject).
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
    pv.contrast (1,:) cell ={}      % Categorical level of factor at which to evaluate the contrast and  significance for each level of the grouping variable
    pv.exclude (:,1) logical =[] 

    pv.plot (1,1) logical = false   % Show a histogram of effects
    pv.alpha (1,1) double = 0.05    % Alpha level (plot only)
    pv.tail (1,1) string = "both"
end
arguments (Output)
    lmeWithRandomSlope (1,1)        % The model with the added slope
    stats table                     % Comparison of the extended model with the original (output of compare)
    reTable table                   % Table with estimates of the random effects for the group
end

% Construct the new formula and compute the model
newFormula = sprintf([char(lme.Formula) ' + (-1 + %s|%s)'], pv.factor , pv.group);
dvCoding = lm.dummyVarCoding(lme);

exclude = lme.ObservationInfo.Excluded;
if ~isempty(pv.exclude)
    exclude = exclude | pv.exclude;
end
lmeWithRandomSlope= fitlme(lme.Variables,newFormula,'Exclude',exclude,'DummyVarCoding',dvCoding);
% Compare the new with the old model
if nargout>1 || pv.plot
    stats = dataset2table(compare(lme,lmeWithRandomSlope));
end
% Determine random effects
if nargout >2 || pv.plot
    [~,~,reTable] = randomEffects(lmeWithRandomSlope);
    reTable =convertvars(dataset2table(reTable),@iscellstr,"string");
    isSlopeRe = reTable.Group==pv.group & (reTable.Name == pv.factor | startsWith(reTable.Name,pv.factor +"_"));
    grpLevels = reTable.Level(isSlopeRe);
    nrGrpLevels = numel(grpLevels);
end


if  ~isempty(pv.contrast)
    % Determine the fixed effect contrast.
    feContrast = lm.contrast(lmeWithRandomSlope,{pv.factor,pv.contrast{1}},{pv.factor,pv.contrast{2}});

    % Evaluate significance of the factor at the specified level
    thisVarInfo = lmeWithRandomSlope.VariableInfo(strcmp(pv.factor,lmeWithRandomSlope.VariableInfo.Row),:);

    % the reTable will have entries for all levels, except the
    % first(reference coded) or last (effects coded).
    if strcmpi(dvCoding,"REFERENCE")
        inTable = thisVarInfo.Range{1}(2:end);
    else
        inTable = thisVarInfo.Range{1}(1:end-1);
    end



   
    tIndividual = table('size',[nrGrpLevels 6],'VariableNames',{'Contrast','Level','p','F','df1','df2'},'VariableTypes',{'double','string', 'double','double','double','double'});
    for grpLvl =1:nrGrpLevels
        stayThisGrp = reTable.Level == grpLevels(grpLvl);
        % Build the RE contrast
        if ~ismember(pv.contrast{1},inTable)
            % Level 1 is not in the model - construct from others
            if dvCoding=="reference"
                lvl1 = zeros(height(reTable),1);
            else
                lvl1 = -stayThisGrp;
            end
        else
            % level 1 is in the model
            lvl1 = double(stayThisGrp & ismember(reTable.Name,pv.factor +"_" + pv.contrast{1}));
        end

        if ~ismember(pv.contrast{2},inTable)
            % Level 2 is not in the model - construct from others
            if dvCoding=="reference"
                lvl2 = zeros(height(reTable),1);
            else
                lvl2 = -stayThisGrp;
            end
        else
            % level 2 is in the model
            lvl2 = double(stayThisGrp & ismember(reTable.Name,pv.factor +"_" + pv.contrast{2}));
        end
        reContrast = (lvl1-lvl2);
        tIndividual.Level(grpLvl) = grpLevels(grpLvl);
        [tIndividual.p(grpLvl),tIndividual.F(grpLvl), tIndividual.df1(grpLvl),tIndividual.df2(grpLvl)] = coefTest(lmeWithRandomSlope,feContrast,0,'REContrast',reContrast');
        tIndividual.Contrast(grpLvl) = reContrast'*reTable.Estimate +feContrast*lmeWithRandomSlope.Coefficients.Estimate;
    end
    name = pv.factor + ":" + pv.contrast{1} + "-" + pv.contrast{2};
    tIndividual= addvars(tIndividual,repmat(name,height(tIndividual),1),'NewVariableNames','ContrastName');
    reTable = innerjoin(reTable,tIndividual);

    isSlopeRe = reTable.Group==pv.group & (reTable.Name == pv.factor | startsWith(reTable.Name,pv.factor +"_")); 
    reTable = reTable(isSlopeRe,:); % Return only the slope effects that were added.
   
    estimate = reTable.Contrast;  % Show the estimated contrasts in the plot below
    p = reTable.p;
    xLbl  = name + " Fixed +";
else    

    isSlopeRe = reTable.Group==pv.group & (reTable.Name == pv.factor | startsWith(reTable.Name,pv.factor +"_")); 
    reTable = reTable(isSlopeRe,:); % Return only the slope effects that were added.
    estimate = reTable.Estimate ;  % Show the random effects in the plot below
    p = reTable.pValue;
    xLbl  = pv.factor;
end

 switch lower(pv.tail)
     case "both"
         isSignificant = p <pv.alpha;
     case "left"
         isSignificant  =estimate<0 & p < pv.alpha;         
     case "right"
         isSignificant = estimate>0 & p < pv.alpha;
 end

if pv.plot
    [allC,edges]= histcounts(estimate);
    centers = edges(1:end-1) + (edges(2:end)-edges(1:end-1))/2;
    bar(centers',allC','b');
    xlabel(xLbl + " Random Effect",'Interpreter','none')
    str= sprintf('AIC %d (Chi^2 (%d) = %3.3f, p= %3.3g)',round(diff(stats.AIC)),stats.deltaDF(2),stats.LRStat(2),stats.pValue(2));
    sigC = histcounts(estimate(isSignificant),edges);
    hold on
    h= bar(centers',sigC','r');
    h.FaceAlpha=0.8;

    legend('All',sprintf('p<%.3f',pv.alpha));
    str = [str sprintf('\n %d/%d (%.2f%%) significant',sum(isSignificant),nrGrpLevels,100*mean(isSignificant))];
    ylabel("#" + pv.group)
    title (str,'Interpreter','none')
    hold off
end

