function [p,stat,df,delta,CI,str,c,debugStr] = posthoc(m,A,B,predictedDelta,tail,alpha,scaleMode)
% Perform a posthoc comparison of condition A and B in a Linear Mixed Model,
% using coefTest, but using cell arrays of parm/value pairs to specify conditions A and B
% instead of the contrast
%
%
% INPUT
% lm                = The linear model
% A                 =  Cell array specifying condition A
% B                 = Cell array specifying condition B
% predictedDelta    = Predicted difference between A and B. Defaults to 0.
% tail              = Specify the tail of the distribution. 'left','right' (use one sided T-tests)
%                       or 'both' (use F-test).
% alpha             = Significance level for the confidence interval.
% scaleMode             = Scaling mode for the delta (RAW,INTERCEPT, RANDOMSTD;
%                       see lm.scaleFactor)
%
% OUTPUT
% p                 = The p-value associated with the test. (see coefTest)
% stat              =  F for two-tailed tests,  T for one sided
% df                = Error degrees of freedom (see coefTest)
% delta             = The difference.
% ci                = The 1-alpha confidence interval
% str                = A char that gives the full stats in a publication  ready format
% contrast          = The contrast used for this test.
% debugStr          = Contrast shown together with coefficient names to
%                       help understand why the contrast is the way it is...
% EXAMPLE
% Fit a model with two fixed-effect predictors and a random
%             effect. Test for the significance of the Cylinders term. The
%             p-value is the same as shown in the anova table.
% load carsmall
% T = table(MPG,Weight,Model_Year,Cylinders);
% T.Cylinders = nominal(T.Cylinders);
% glme = fitglme(T,'MPG ~ Weight + Cylinders + (1|Model_Year)','Distribution','Normal','DummyVarCoding','Effects')
% H0 8 and 6 cylinders are the same.
% p = lm.posthoc(glme,{'Cylinders',8},{'Cylinders',6})
%
% EXAMPLE
% Use cell arrays to select specific cells in a multifactorial design. E.g. in a design
% with stim and valid factors, each with multiple levels
% {'stim','VIS','valid',1} represents the cell in which the stim factor has level VIS and the
% valid factor has level 1.
%
%
% BK -  Jan 2021
% Mar 2021- rewrote to use lm.contrast


nrContrast = size(A,1);



nin =nargin;
if nin<7
    scaleMode = 'RAW';
    if nin <6
        alpha = 0.05;
        if nin <5
            tail = 'both';
            if nin <4
                predictedDelta = zeros(nrContrast,1);
            end
        end
    end
end

if nrContrast>1
    % Recursively call this function for each row in the contrast.
    p = nan(nrContrast,1);
    stat =nan(nrContrast,1);
    df =nan(nrContrast,1);
    delta =nan(nrContrast,1);
    CI =nan(nrContrast,2);
    str =cell(nrContrast,1);
    debugStr =cell(nrContrast,1);
    for i=1:nrContrast
        [p(i),stat(i),df(i),delta(i),CI(i,:),str{i},c(i,:),debugStr{i}] = lm.posthoc(m,A(i,:),B(i,:),predictedDelta(i),tail,alpha,scaleMode); %#ok<AGROW>
    end
    return;
end


%% Determine the estimate using the linear model FE
if isnumeric(A)
    if nin < 3
        B =zeros(size(A));
    end
    c =A-B;
else
    c  = lm.contrast(m,A,B); % the linear contrast
end


if isa(m,'LinearModel')
    delta  =c*m.Coefficients.Estimate;
elseif isa(m,'LinearMixedModel') || isa(m,'GeneralizedLinearMixedModel')
    delta = c*m.fixedEffects;
else
    error('Unknown model type')
end


debugStr = strcat(m.CoefficientNames', ' : ' , cellstr(num2str(c')));
predictedDelta = predictedDelta(:).*ones(size(delta));
if ~all(size(delta)==size(predictedDelta))
    error('Predicted delta [nrRows nrCols] has to match the expected delta');
end

if size(delta,1)>1 && ~strcmpi(tail,'both')
    error('One-sided posthoc tests are only defined for conditions, not for factors');
end

%% Assess statistical significance using F or T
switch upper(tail)
    case 'BOTH'
        % Two-sided tests - use the built-in F test.
        [p,stat,~,df] = coefTest(m,c,predictedDelta);
        statName= 'F';
    case {'LEFT','RIGHT'}
        cov  = c*m.CoefficientCovariance*c';
        stat = (delta-predictedDelta)/sqrtm(cov);% Use sqrtm for more stability (less speed)
        df  =  m.DFE;
        if strcmpi(tail,'LEFT')
            p = tcdf(stat,m.DFE);
        else
            p = 1-tcdf(stat,m.DFE);
        end
        statName= 'T';
    otherwise
        error('%s is not a valid tail specification',tail);
end

%% Scalnig
[scale,units] = lm.scaleFactor(m,scaleMode);
delta = delta/scale;

%% Alpha CI
if nargout > 4
    % Compute 1-alpha point estmate confidence intervals, only if requested .
    cov  = c*m.CoefficientCovariance*c';
    criterion = tinv(1-alpha/2,m.DFE);
    updown= criterion.*sqrtm(cov)/scale; % Use sqrtm for more stability (less speed)
    CI = [delta - updown, delta + updown];
end

if nargout > 5
    % Report string
    str = sprintf(['%s(%d) = %.3g, p= %.3g, delta= %.3g' units ' (%d%% CI [%.3g %.3g])'],statName,df,stat,p,delta,100*(1-alpha),CI);
end

end
