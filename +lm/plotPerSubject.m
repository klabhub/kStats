function [effectsTable,ciTable,pTable,models] = plotPerSubject(m,pv)
% For a linear model based on one or more groups (e.g., subjects), refit the model
% per level of the group (i.e., per subject) and show the results for a subset of model coefficients
% to get an idea of  variability across the different levels.
%
% INPUT
% glm -  a (generalized) linear mixed model
%
% OUTPUT
% effectsTable  - table with fixed effects. First column is the group, the other
%               columns are the individual subjects
% ciTable - table with confidence intervals.
% pTable - table with p-values
% models - Cell array with the models per subject.
%
% BK - Feb 2020.
arguments
    m (1,1)  % Linear model
    pv.coefficients (1,:) string = ""        % Show only these coefficients  (Defaults to all model coefficients)
    pv.zscore (1,:) string = ""              % Z-score these variables in the model for each level of the grouping variable. Defaults to none.
    pv.showHistogram (1,1) logical = false  % Histogram of effects across subjects
    pv.showLine (1,1)logical = true         % Line graph of effects for each subject
    pv.showGroup (1,1) logical = true       % Add the group level result as a red line
    pv.NUMPRECISION (1,1) double {mustBeInteger,mustBeNonnegative} =3; % num2str for coefficients
    pv.alpha (1,1) double = 0.05;
    pv.groupingVariable (1,:) string = ""  % The name of the variable that defines the grouping. Defaults to all grouping variables.    
end
assert(all(pv.coefficients=="") || all(pv.coefficients=="*") || all(ismember(pv.coefficients,m.Coefficients.Name)),"pv.coefficients should list the names of coefficients in the model");
assert(all(pv.zscore =="") || all(ismember(pv.zscore,m.CoefficientNames)),"pv.zscore should list the names of coefficients in the model");
%% Extract the data from the full model,
T = m.Variables(~m.ObservationInfo.Excluded,m.VariableInfo.InModel | ismember(m.VariableNames,m.Formula.ResponseName));
dummyVarCoding = lm.dummyVarCoding(m);
% Determine which variable represents the grouping (e.g. subject)
if pv.groupingVariable ==""
    % Use all
    pv.groupingVariable = string(m.Formula.GroupingVariableNames{:});
end
groupLevels = unique(T(:,pv.groupingVariable),"rows");
formula = m.Formula.char;

%% Createa the output tables and add the results for the full sample
feNames = m.Coefficients.Name;
fe = m.fixedEffects;
effectsTable = table(fe,'RowNames',feNames,'VariableNames',{'Group'});
ci = [m.Coefficients.Lower m.Coefficients.Upper];
ciTable  = table(ci,'RowNames',feNames,'VariableNames',{'Group'});
pTable  = table(m.Coefficients.pValue,'RowNames',feNames,'VariableNames',{'Group'});

%% Now fit each group level separately and add to the tables
nrGroupLevels = height(groupLevels);
if nargout>3
    models= cell(1,nrGroupLevels);
end
groupNames= convertvars(groupLevels,@(x)(~isstring(x)),'string');
if width(groupNames)>1
    groupNames= join(groupNames{:,:},"/"); % Join columns as string array
else
    groupNames = groupNames{:,1}; % Extract as string array
end
for grpCntr=1:nrGroupLevels    
    try
        thisGroupName = groupNames(grpCntr);
        % Extract the relevant subset of data for this subjects
        thisT = innerjoin(T,groupLevels(grpCntr,:));        
        if pv.zscore ~=""
            for factor = pv.zscore
                thisT.(factor) = zscore(thisT.(factor));
            end
        end
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
        % Make sure Fe are placed in the correct row of the effects table.
        [tf,order] =ismember(thisFeNames,feNames);
        assert(all(tf),'Missing FE in the per-level fit for %s.',thisGroupName);
        effectsTable = [effectsTable   table(fe(order),'VariableNames',thisGroupName)]; %#ok<AGROW>               
        ci = [thisGlm.Coefficients.Lower thisGlm.Coefficients.Upper];
        ciTable  = [ciTable table(ci(order,:),'VariableNames',thisGroupName)];%#ok<AGROW>
        pTable = [pTable table(thisGlm.Coefficients.pValue,'VariableNames',thisGroupName)]; %#ok<AGROW>        
    catch me
        fprintf('perSubject lmm for group %s failed (%s) - %s\n',thisGroupName,me.message,formula)
        continue
    end
end

%%
if pv.coefficients == ""
    % All except intercept
    pv.coefficients = string(feNames(2:end));
elseif pv.coefficients == "*"
    % Include intercept
    pv.coefficients = string(feNames);
end
keepCoeffs = ismember(effectsTable.Properties.RowNames,pv.coefficients);
effectsTable =effectsTable(keepCoeffs,:);
ciTable = ciTable(keepCoeffs,:);
pTable =pTable(keepCoeffs,:);
feNames= feNames(keepCoeffs);
%% Visualize the main group and individual results
if pv.showHistogram || pv.showLine
    clf;
    nrGroups = width(effectsTable)-1;
    nrEffects = height(effectsTable);
    nrRows = sum(pv.showHistogram+pv.showLine);
    nrCols = nrEffects;
    for e =1:nrEffects
        if pTable{e,1} <pv.alpha
            sigStr = '(*)';
        else
            sigStr = '';
        end
        titleStr = [feNames{e} ': ' num2str(effectsTable{e,1},pv.NUMPRECISION) ' ' sigStr] ;
        titleStr = char(titleStr,sprintf('Same sign: %d/%d ',sum(sign(effectsTable{e,1})==sign(effectsTable{e,2:end})),nrGroups));
        titleStr = char(titleStr,sprintf('Significant (alpha = %.3f): %d/%d',pv.alpha, sum(pTable{e,2:end}<pv.alpha),nrGroups));

        if pv.showHistogram
            h1= subplot(nrRows,nrCols,e);
            % Top row shows histogram of effects
            histogram(effectsTable{e,2:end});
            title(titleStr,'Interpreter','None');
            xlabel([feNames{e} ' Effect'],'Interpreter','none')
        end

        if pv.showLine
            h2= subplot(nrRows,nrCols,e+nrCols*pv.showHistogram);
            % Bottom row shows line plots with CI per group level (subject)
            line(reshape(ciTable{e,2:end},[2 nrGroups]),repmat(1:nrGroups,[2 1]),'Color','k')
            hold on
            plot(effectsTable{e,2:end},1:nrGroups,'k.')
            line([0 0],[0 nrGroups])
            if pv.showGroup
            % Show the group effect as a red line in the middle
            y = nrGroups/2 + mod(nrGroups+1,2)/2;
            plot(effectsTable{e,1},y,'r*')
            line(ciTable{e,1},y*[1 1],'Color','r','LineWidth',2);
            end
            ylim([0 nrGroups]);
            set(gca,'yTickLabels',{})
            if e==1
                ylabel(thisGroupName + "#")
                set(gca,'yTick',1:nrGroups,'ytickLabel',groupNames)
            end
            if ~pv.showHistogram
                title(titleStr,'Interpreter','None');
            end
            xlabel([feNames{e} ' Effect'],'Interpreter','none')
        end

        if pv.showHistogram && pv.showLine
            linkaxes([h1 h2],'x')
        end
    end
end

