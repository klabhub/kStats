function [p,stat,df,delta,CI,str,c,debugStr] = posthoc(m,A,B,predictedDelta,tail,alpha,scaleMode,dfMethod)
% Perform a posthoc comparison of condition A and B in a Linear Mixed Model,
% using cell arrays of parm/value pairs to specify conditions A and B
%
% INPUT
% lm                = The linear model
% A                 =  Cell array specifying condition A
% B                 = Cell array specifying condition B
% predictedDelta    = Predicted difference between A and B. Defaults to 0.
% tail              = Specify the tail of the distribution. 'left','right' or 'both'.
% alpha             = Significance level for the confidence interval.
% scaleMode             = Scaling mode for the delta (RAW,INTERCEPT, RANDOMSTD;
%                       see lm.scaleFactor)
% dfMethod          = Method to determine the error degrees of freedom.
%                   Defaults to 'residual', but can be set to
%                   'satterthwaite'.  Note that 'satterthwaite' will only
%                   work if A and B are conditions in the design (and not
%                   an arbitrary linear combination of the fixed effects).
%
% OUTPUT
% p                 = The p-value associated with the test.
% stat              =  T statistic
% df                = Error degrees of freedom for the T statistics
% delta             = The difference.
% ci                = The 1-alpha confidence interval
% str                = A char that gives the full stats in a publication ready format
% contrast          = The numeric contrast used for this test (derived from A and B specification)
% debugStr          = Contrast shown together with coefficient names to help understand why the contrast is the way it is...
%
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
% Jul 2025- added satterthwaite dof approximation
arguments
    m (1,1)
    A (1,:)
    B (1,:)
    predictedDelta (:,1) double  = zeros(size(A,1))
    tail (1,1) string = "both"
    alpha (1,1) double = 0.05
    scaleMode (1,1) string = "RAW"
    dfMethod (1,1) string  {mustBeMember(dfMethod,["residual","satterthwaite"])} = "residual"
end
nrContrast = size(A,1);
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
        [p(i),stat(i),df(i),delta(i),CI(i,:),str{i},c(i,:),debugStr{i},df(i)] = lm.posthoc(m,A(i,:),B(i,:),predictedDelta(i),tail,alpha,scaleMode,dfMethod);
    end
    return;
end

%% Translate the two conditions A/B into a numeric contrast vector
[c,TA,TB]  = lm.contrast(m,A,B); % the linear contrast
assert((istable(TA) && istable(TB)) || dfMethod=="residual","Satterthwaite residuals can only be computed for conditions in the model. Use cell arrays to define the contrast. ")
assert(~(dfMethod=="satterthwaite" && isa(m,'LinearModel')),"Satterthwaite residuals only apply to mixed effects models" )

%% Estimate the value using the linear model FE
if isa(m,'LinearModel')
    delta  =c*m.Coefficients.Estimate;
    df = m.DFE;
elseif isa(m,'LinearMixedModel') || isa(m,'GeneralizedLinearMixedModel')
    switch dfMethod
        case 'residual'
            delta = c*m.fixedEffects;
            df  =  m.DFE;
        case 'satterthwaite'
            [a,~,dfa] = predict(m,TA,'DFMethod','satterthwaite','Conditional',false);
            [b,~,dfb] = predict(m,TB,'DFMethod','satterthwaite','Conditional',false);
            assert((dfa-dfb)<0.01,"Satterthwaite dof differ between conditions.")
            delta = a-b;
            df = dfa;
    end
else
    error('Unknown model type')
end

debugStr = strcat(m.CoefficientNames', ' : ' , cellstr(num2str(c')));
assert(all(size(delta)==size(predictedDelta)),'Predicted delta [nrRows nrCols] has to match the expected delta');

%% Assess statistical significance using T
cov  = c*m.CoefficientCovariance*c';
se = sqrt(cov);  % scalar standard error
stat = (delta-predictedDelta)/se; % Scalar contrasts only.
pLeft = tcdf(stat,df);
pRight = 1-tcdf(stat,df);
statName= 'T';
switch upper(tail)
    case 'BOTH'
        p = 2*min(pLeft,pRight);
    case 'LEFT'
        p = pLeft;
    case 'RIGHT'
        p= pRight;
end

%% Scalnig
[scale,units] = lm.scaleFactor(m,scaleMode);
delta = delta/scale;

%% Alpha CI
if nargout > 4
    % Compute 1-alpha point estmate confidence intervals, only if requested .
    criterion = tinv(1-alpha/2,df);
    updown= criterion.*se/scale;
    CI = [delta - updown, delta + updown];
end

if nargout > 5
    % Report string
    str = sprintf(['%s(%d) = %.3g, p= %.3g, delta= %.3g' units ' (%d%% CI [%.3g %.3g])'],statName,df,stat,p,delta,100*(1-alpha),CI);
end

end
