function [prediction,low,high,conditions] = plotPredict(m,pv)
%plotPredict - Plot predictions for all categorical predictor combinations
%   [PREDICTION,LOW,HIGH,CONDITIONS] = plotPredict(M) generates a table of
%   all combinations of categorical fixed-effect predictors, predicts the
%   response for every row, and plots the predictions with confidence
%   intervals.
%
%   [...] = plotPredict(M,ContinuousValues=VALUES) supplies one row of
%   values for every continuous fixed-effect predictor. VALUES must be a
%   table with one scalar variable for each continuous predictor.
%
%   [...] = plotPredict(...,sigOnly=true) plots and returns only predictions
%   whose confidence intervals exclude zero.
%
%   [...] = plotPredict(...,sigDiffOnly=true) plots and returns only
%   predictions whose confidence intervals do not overlap the intercept
%   confidence interval.
%
%   [...] = plotPredict(...,plot2D=true) displays predictions for a model
%   with two categorical fixed predictors as a two-dimensional grid.
%
%   See also predict, bar, errorbar
arguments
    m (1,1)
    pv.continuousValues table = table()
    pv.sigOnly (1,1) logical = false
    pv.sigDiffOnly (1,1) logical = false
    pv.plot2D (1,1) logical = false
    pv.conditional (1,1) logical = false
    pv.alpha (1,1) double {mustBeGreaterThanOrEqual(pv.alpha,0), ...
        mustBeLessThanOrEqual(pv.alpha,1)} = 0.05
end

predictorNames = string(m.Formula.FELinearFormula.PredictorNames);
[conditions,labels] = predictionTable(m,predictorNames,pv.continuousValues);
[prediction,ci] = predict(m,conditions,'Conditional',pv.conditional, ...
    'Alpha',pv.alpha);
low = ci(:,1);
high = ci(:,2);
alphaData = ones(size(prediction));

if pv.sigOnly
    if pv.plot2D
        alphaData = min(alphaData,significanceAlpha(low,high,0,0));
    else
        keep = low > 0 | high < 0;
        prediction = prediction(keep);
        low = low(keep);
        high = high(keep);
        conditions = conditions(keep,:);
        labels = labels(keep);
    end
end

if pv.sigDiffOnly
    [interceptLow,interceptHigh] = interceptInterval(m);
    if pv.plot2D
        alphaData = min(alphaData,significanceAlpha(low,high, ...
            interceptLow,interceptHigh));
    else
        keep = high < interceptLow | low > interceptHigh;
        prediction = prediction(keep);
        low = low(keep);
        high = high(keep);
        conditions = conditions(keep,:);
        labels = labels(keep);
    end
end

if pv.plot2D
    plotPredictionGrid(m,prediction,low,high,conditions,alphaData);
else
    set(gcf,'defaultTextInterpreter','none')
    set(gcf,'defaultLegendInterpreter','none')
    set(gcf,'defaultAxesTickLabelInterpreter','none')
    x = 1:numel(prediction);
    bar(x,prediction,'FaceColor','w','EdgeColor','k');
    hold on
    errorbar(x,prediction,prediction-low,high-prediction,'.', ...
        'CapSize',0,'Color','k');
    ylabel(m.ResponseName)
    xlabel('Prediction conditions')
    set(gca,'XTick',x,'XTickLabel',labels,'XtickLabelRotation',45)
end
end

function plotPredictionGrid(m,prediction,~,~,conditions,alphaData)
predictorNames = string(m.Formula.FELinearFormula.PredictorNames);
variableNames = string(m.VariableInfo.Properties.RowNames);
variableIx = zeros(1,numel(predictorNames));
for i = 1:numel(predictorNames)
    variableIx(i) = find(variableNames == predictorNames(i),1);
end
if numel(predictorNames) ~= 2 || ...
        ~all(m.VariableInfo.IsCategorical(variableIx))
    error('lm:plotPredict:Plot2DRequiresTwoCategoricalPredictors', ...
        'plot2D requires exactly two categorical fixed predictors.');
end

levels1 = variableLevels(m,variableIx(1));
levels2 = variableLevels(m,variableIx(2));
grid = NaN(numel(levels1),numel(levels2));
alphaGrid = NaN(size(grid));
for row = 1:height(conditions)
    level1 = string(conditions.(predictorNames(1))(row));
    level2 = string(conditions.(predictorNames(2))(row));
    i = find(levels1 == level1,1);
    j = find(levels2 == level2,1);
    grid(i,j) = prediction(row);
    alphaGrid(i,j) = alphaData(row);
end

set(gcf,'defaultTextInterpreter','none')
set(gcf,'defaultColorbarTickLabelInterpreter','none')
imageHandle = imagesc(grid);
imageHandle.AlphaData = alphaGrid;
imageHandle.AlphaDataMapping = 'none';
axis tight
colorbar
xlabel(predictorNames(2))
ylabel(predictorNames(1))
set(gca,'XTick',1:numel(levels2),'XTickLabel',levels2, ...
    'YTick',1:numel(levels1),'YTickLabel',levels1)
title(m.ResponseName)
end

function alpha = significanceAlpha(low,high,nullLow,nullHigh)
% Map separation from a null interval to opacity.
distance = max(max(nullLow-high,low-nullHigh),0);
errorBarLength = max(high-low,eps);
alpha = 0.2 + 0.8*min(1,distance./errorBarLength);
end

function [conditions,labels] = predictionTable(m,predictorNames,continuous)
if isempty(predictorNames)
    conditions = table('Size',[1 0],'VariableTypes',{}, ...
        'VariableNames',{});
    labels = "(Intercept)";
    return
end

variableNames = string(m.VariableInfo.Properties.RowNames);
values = cell(1,numel(predictorNames));
isCategorical = false(1,numel(predictorNames));
categoriesForPredictor = cell(1,numel(predictorNames));
for predictorIx = 1:numel(predictorNames)
    variableIx = find(variableNames == predictorNames(predictorIx),1);
    if isempty(variableIx)
        error('lm:plotPredict:UnknownPredictor', ...
            'Predictor %s was not found in the fitted model.', ...
            predictorNames(predictorIx));
    end

    isCategorical(predictorIx) = m.VariableInfo.IsCategorical(variableIx);
    if isCategorical(predictorIx)
        levels = variableLevels(m,variableIx);
        categoriesForPredictor{predictorIx} = levels;
        values{predictorIx} = levels;
    else
        if height(continuous) ~= 1 || ...
                ~ismember(predictorNames(predictorIx),string(continuous.Properties.VariableNames))
            error('lm:plotPredict:ContinuousValueRequired', ...
                ['Continuous predictor %s requires one scalar value in ' ...
                'ContinuousValues.'],predictorNames(predictorIx));
        end
        value = continuous.(predictorNames(predictorIx));
        if ~isscalar(value)
            error('lm:plotPredict:ContinuousValueRequired', ...
                'Continuous predictor %s must have one scalar value.', ...
                predictorNames(predictorIx));
        end
        values{predictorIx} = value;
    end
end

grid = cell(1,numel(values));
gridInputs = values;
for predictorIx = 1:numel(values)
    if ~isCategorical(predictorIx)
        gridInputs{predictorIx} = values{predictorIx};
    end
end
[grid{:}] = ndgrid(gridInputs{:});
nRows = numel(grid{1});
conditions = table();
for predictorIx = 1:numel(predictorNames)
    name = predictorNames(predictorIx);
    if isCategorical(predictorIx)
        levels = categoriesForPredictor{predictorIx};
        observed = m.Variables.(name);
        conditions.(name) = matchingClass(string(grid{predictorIx}(:)), ...
            observed,levels);
    else
        conditions.(name) = repmat(grid{predictorIx}(1),nRows,1);
    end
end

% LinearMixedModel and GeneralizedLinearMixedModel validate grouping
% variables even for population-level predictions. Supply one valid
% observed value for each grouping variable; Conditional=false prevents it
% from contributing a fitted random effect.
requiredNames = string(m.Formula.VariableNames);
auxiliaryNames = setdiff(requiredNames, ...
    [string(m.ResponseName) predictorNames],'stable');
for name = auxiliaryNames
    if ismember(name,string(conditions.Properties.VariableNames))
        continue
    end
    observed = m.Variables.(name);
    conditions.(name) = repmat(observed(1),nRows,1);
end

labelParts = strings(nRows,numel(predictorNames));
for predictorIx = 1:numel(predictorNames)
    valuesForLabel = string(conditions.(predictorNames(predictorIx)));
    if ~isCategorical(predictorIx)
        valuesForLabel = predictorNames(predictorIx) + "=" + valuesForLabel;
    end
    labelParts(:,predictorIx) = valuesForLabel;
end
labels = join(labelParts,":",2);
end

function levels = variableLevels(m,variableIx)
levels = m.VariableInfo.Range(variableIx);
if iscell(levels)
    levels = levels{1};
end
levels = string(levels);
end

function values = matchingClass(values,observed,levels)
% Preserve the predictor class expected by the fitted model.
if iscell(observed)
    values = cellstr(values);
elseif isstring(observed)
    values = string(values);
elseif iscategorical(observed)
    values = categorical(values,cellstr(levels));
else
    values = cast(values,'like',observed(1));
end
end

function [low,high] = interceptInterval(m)
interceptIx = string(m.CoefficientNames) == "(Intercept)";
low = m.Coefficients.Lower(interceptIx);
high = m.Coefficients.Upper(interceptIx);
if isa(m,'GeneralizedLinearMixedModel')
    low = m.Link.Inverse(low);
    high = m.Link.Inverse(high);
end
end
