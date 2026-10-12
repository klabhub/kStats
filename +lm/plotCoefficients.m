function [estimate,low,high] =plotCoefficients(m,fe,x,scale,pv)
% Visualize the effects of a generalized linear model by pulling the betas
% back trough the inverse of the link function and displaying the effects
% as bars with error bars.
%
% INPUT
% lm - Linear model
% fe - Names of the fixed effects to show. 
%           Empty (default) will show all except the intercept
%           "*" will show all
% x - Horizontal axis. [1:nrFixedEffects].
% scale - Scale the effects to
%                      'Intercept' ; scale to the fixed intercept effect
%                       'RANDOMSTD' : scale to the standard deviation of
%                       the random effects.
%                       'none' - do not scale (default)
% 'addIntercept' - Set this to true to add the intercept to the
%                   coefficients to show the overall effect at each level.
%                   
% 'sigOnly' - Set to true to show only coefficients that are significantly
%               nonzero (based on their (non-simultaneous) confidence
%               intervals)
% OUTPUT
% estimate - Fxied effect estimates
% lower  - lower edge of CI
% upper - upper edge of CI
%
% BK - Feb 2020.
arguments
    m (1,1) % GeneralizedLinearMixedModel or LinearMixedModel
    fe (1,:) string = string.empty
    x  = [] 
    scale (1,1) string = "none"
    pv.addIntercept (1,1) logical = false
    pv.sigOnly (1,1) logical = false 
end 
allFeNames = string(m.CoefficientNames);
if isscalar(fe) && fe == "*"
    fe = allFeNames(:).';
elseif isempty(fe)
    % Show all except intercept
    fe = setdiff(allFeNames,"(Intercept)",'stable');
else
    % Plot explicitly requested coefficients.
end
if pv.addIntercept
    % The intercept is included in each total estimate, not plotted as a
    % separate effect.
    fe = fe(fe ~= "(Intercept)");
end
nrFe = numel(fe);

if isa(m,'GeneralizedLinearMixedModel')
    link  = m.Link.Inverse;
else
    link = @(x)(x);
end

coeffs = m.Coefficients;
interceptIx = allFeNames == "(Intercept)";


low = NaN(nrFe,1);
high = NaN(nrFe,1);
estimate = NaN(nrFe,1);
for feIx=1:nrFe
    ix = coefficientIndex(allFeNames,fe(feIx));
    if ~isempty(ix)
        if pv.addIntercept && ~interceptIx(ix) && any(interceptIx)
            % Add the intercept and propagate its covariance with this
            % coefficient to obtain the total estimate for this level.
            contrast = zeros(1,numel(allFeNames));
            contrast(ix) = 1;
            contrast(interceptIx) = 1;
            beta = contrast * coeffs.Estimate;
            se = sqrt(max(0,contrast * m.CoefficientCovariance * contrast'));
            critical = coefficientCriticalValue(coeffs,ix);
            low(feIx) = link(beta-critical*se);
            high(feIx) = link(beta+critical*se);
            estimate(feIx) = link(beta);
        else
            low(feIx) = link(coeffs.Lower(ix));
            high(feIx) = link(coeffs.Upper(ix));
            estimate(feIx) = link(coeffs.Estimate(ix));
        end
    else
        error('lm:plotCoefficients:UnknownFixedEffect', ...
            'Fixed effect %s does not exist in this model',fe(feIx));
    end
end

switch upper(scale)
    case 'NONE'
        scale =1;
        units ='';
    case 'INTERCEPT'
        % Scale to the fixed effect intercept
        % How big is the effect relative to the grand mean.
        scale  = 0.01*m.Coefficients.Estimate(strcmpi(m.CoefficientNames,'(Intercept)'));
        units = '(% grand mean)';
    case 'RANDOMSTD'
        % Scale to the standard deviation of the random effects
        % "How big is the effect relative to the variation
        % across subjects?'
        scale = 0.01*std(m.randomEffects);
        units = '(% Random Effects Stdev)';
    otherwise
        error ('Unknown scale %s ', scale);
end

estimate=estimate/scale;
low = low/scale;
high= high/scale;

out = false(size(fe));
if pv.sigOnly
 out = sign(low)~=sign(high);
 estimate(out)= [];
 low(out) = [];
 high(out) =[];
 nrFe = numel(estimate);
 fe(out) = [];
end
if isempty(x) 
    x = 1:nrFe;    
    xLabel = fe;
elseif isstring(x)
    xLabel =x(~out); 
    x = 1:nrFe;    
else
    xLabel =x(~out);    
end


set(gcf, 'defaultTextInterpreter', 'none')
set(gcf, 'defaultLegendInterpreter', 'none')
set(gcf, 'defaultAxesTickLabelInterpreter', 'none')
bar(x,estimate,'FaceColor','w','EdgeColor','k');
hold on
errorbar(x,estimate,estimate-low,high-estimate,'.','CapSize',0,'Color','k');
ylabel([m.ResponseName ' ' units])
set(gca,'XTick',x,'XTickLabel',xLabel,'XtickLabelRotation',45);

end

function critical = coefficientCriticalValue(coeffs,ix)
if any(ix)
    critical = median((coeffs.Upper(ix)-coeffs.Estimate(ix))./coeffs.SE(ix), ...
        'omitnan');
else
    critical = norminv(.975);
end
if isempty(critical) || ~isfinite(critical)
    critical = norminv(.975);
end
end

function ix = coefficientIndex(allFeNames,requested)
% Find a coefficient, treating the factors in an interaction as unordered.
ix = coefficientIndexForTerm(allFeNames,requested);
end

function ix = coefficientIndexForTerm(allFeNames,requested)
% Find a coefficient, treating the factors in an interaction as unordered.
ix = find(allFeNames == requested,1);
if isempty(ix) && contains(requested,":")
    requestedFactors = sort(split(requested,":"));
    for i = 1:numel(allFeNames)
        if contains(allFeNames(i),":") && ...
                isequal(sort(split(allFeNames(i),":")),requestedFactors)
            ix = i;
            return
        end
    end
end
end
