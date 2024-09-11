function anova =anova(meta,poolT2)
% Compare summary effects in groups of studies. 
% This ANOVA always uses a random effects analysis (i.e. the true effect sizes of the groups
% are expected to be different).

arguments
    meta (:,1) struct {mustBeVector}  % A vector of meta structs (see  ma.analyze)
    poolT2 (1,1) logical  = false     % Set to true to pool the across-study variance across groups (Borenstein p 171)
end

nrGroups= numel(meta)-1;
ts = [meta.T];
if poolT2
    % Use the pooled T2
    % page 171 in Borenstein et al.
    T2 =  [ts.pooled].^2;
else
    % Use the T2 as already estimated per group
    % page 164 in Borenstein et al.
    T2 =  [ts.value].^2;
end

% We have to recompute the Q's specifically for the ANOVA to include
% random effects (true effect sizes of groups can be different)
% page 164 in Borenstein et al.
qFun =@(w,gg) (w'*gg.^2- (w'*gg).^2./sum(w));
W =cell(nrGroups,1);
g =cell(nrGroups,1);
% Collect effects and weights
for grp=1:nrGroups
    g{grp} =meta(grp).effect.value;
    gVar = meta(grp).effect.variance;
    W{grp} =  1./(gVar +T2(grp)); % Weigh by variance plus between study variance T2
end
% Determine Qstar for each group (page 169 in Borenstein)
Qstar = cellfun(qFun,W,g);
qWithin = sum(Qstar);
% Determine grand mean Qstart
W  = cat(1,W{:});
g =  cat(1,g{:});
qTotal = qFun(W,g);

% Standard anova.
anova.qBetween = qTotal - qWithin;
anova.df = nrGroups-1;
anova.p = 1-chi2cdf(anova.qBetween,anova.df);
% R2 = 1-T2wiithin/T2total
T2Within = T2; % pooled variance (same in all ts)
TTotal = meta(end).T.value.^2; % Variance of the full set
anova.R2 = 1-T2Within/TTotal;

% Also determine the pairwise difference between groups and  of the
% summary effect and the corresponding confidence interval. (Borenstein 19.17)
ms = [meta(1:nrGroups).summary];
M  = [ms.value];
varM = [ms.variance];

anova.diff.value = M -M';
anova.diff.se  = sqrt(varM+varM');
anova.diff.lower = anova.diff.value-1.96*anova.diff.se;
anova.diff.upper= anova.diff.value+1.96*anova.diff.se;
anova.qStar =Qstar;
end