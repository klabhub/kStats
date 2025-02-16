function [effectsT,ciT] = plotPerSubject(m,pv)
% For a given linear model based on one or more subjects, refit the model
% per subject and show the results to get an idea of variability in the
% sample.
%
% INPUT
% glm -  a (generalized) linear mixed model
%
% OUTPUT
% effects  - table with fixed effects. First row is the group, the other
%               rows are the individual subjects
% ci - table with confidence intervals.
%
% BK - Feb 2020.
arguments
    m (1,1)  % Linear model
    pv.showHistogram (1,1) logical = false  % Histogram of effects across subjects
    pv.showLine (1,1)logical = true         % Line graph of effects for each subject
    pv.NUMPRECISION (1,1) double {mustBeInteger,mustBeNonnegative} =3; % num2str for effects
    pv.effects (1,:) string = ""            % Which effects to show. Defaults to all except intercept.Set to "*" to include the intercept.
end

dummyVarCoding = lm.dummyVarCoding(m);

%% Createa a table of fixed effects and confidence intervals
fe = m.fixedEffects;
feNames = m.CoefficientNames;
effectsT = table(fe,'RowNames',feNames,'VariableNames',{'Group'});
ci = [m.Coefficients.Lower m.Coefficients.Upper];
ciT  = table(ci,'RowNames',feNames,'VariableNames',{'Group'});


if pv.effects ==""
    pv.effects = feNames(2:end); % Exclude intercept
elseif pv.effects == "*"
    pv.effects = feNames; % Include intercept
end
%% Now fit each subject separately and add to the tables

T = m.Variables(:,m.VariableInfo.InModel | ismember(m.VariableNames,m.Formula.ResponseName));
subjects = unique(T{:,m.Formula.GroupingVariableNames{:}}); % only keep subjects who're relevant to this condition
formula = m.Formula.char;
for s=subjects'
    try
        % Extract the relevant subset of data for this subjects
        thisT = T(~m.ObservationInfo.Excluded & T.subject==s,:);
        % Refit for this subject
        if isa(m,'LinearMixedModel')
            thisGlm = fitlme(thisT,formula,'FitMethod',m.FitMethod,'DummyVarCoding',dummyVarCoding) ;
        else
            lastwarn('')
            thisGlm = fitglme(thisT,formula,'FitMethod',m.FitMethod,'Distribution',...
                m.Distribution,'Link',m.Link,'DummyVarCoding',dummyVarCoding) ;
            [msg,id] = lastwarn;
            if strcmpi(id,'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_PLUnableToConverge')
                lastwarn('')
                error(msg);
            end
        end
        
        fe = thisGlm.fixedEffects;       
        thisFeNames = thisGlm.CoefficientNames;
        % Some fe in the per-subject fit may be ordered differently; make
        % sure they are placed in the correct row of the effects table.
        [tf,order] =ismember(thisFeNames,feNames);
        assert(all(tf),'Missing FE in the per-subject fit for %s.',s);
        effectsT = [effectsT   table(fe(order),'VariableNames',string(s))]; %#ok<AGROW>               
        ci = [thisGlm.Coefficients.Lower thisGlm.Coefficients.Upper];
        ciT  = [ciT table(ci(order,:),'VariableNames',string(s))];%#ok<AGROW>        
    catch me
        fprintf('perSubject lmm for %s failed on %s (%s)\n',formula,s,me.message)
        continue
    end
end

%% Visualize the main group and individual results
nrSubjects = width(effectsT)-1;
% Restrict effects 
keep = ismember(effectsT.Properties.RowNames,pv.effects);
effectsT  =effectsT(keep,:);
ciT = ciT(keep,:);

nrEffects = height(effectsT);
nrRows = sum(pv.showHistogram+pv.showLine);
for e =1:nrEffects
    if pv.showHistogram
        subplot(nrRows,nrEffects,e);
        % Top row shows histogram of effects
        histogram(effectsT{e,2:end});
    end

    if pv.showLine
        subplot(nrRows,nrEffects,e+nrEffects*pv.showHistogram);
        % Bottom row shows line plots with CI per subject.
        line(reshape(ciT{e,2:end},[2 nrSubjects]),repmat(1:nrSubjects,[2 1]),'Color','k')
        hold on
        plot(effectsT{e,2:end},1:nrSubjects,'k.')
        line([0 0],[0 nrSubjects])
        % Show the group effect as a red line
        line(ciT{e,1},0.5+nrSubjects/2*[1 1],'Color','r','LineWidth',2);
        ylim([0 nrSubjects]);
        set(gca,'yTickLabels',{})
        if e==1
            ylabel 'Subject #'
            set(gca,'yTick',1:nrSubjects,'ytickLabel',subjects)
        end
    end

    if prod(sign(ciT{e,1})) >0 % CI both <0 or both >0; significant at alpha level.
        sigStr = '(*)';
    else
        sigStr = '';
    end
    t=  title([effectsT.Properties.RowNames{e} ': ' num2str(effectsT{e,1},pv.NUMPRECISION) ' ' sigStr]);
    t.Interpreter = 'None';
    xlabel 'Effect'
    hold off
end

