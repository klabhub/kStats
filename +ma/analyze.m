function [meta,anova] = analyze(varargin,pv)
% Given means, standard deviations, and the number of samples for the experimental
% and control groups in each study, do a meta-analysis following Borenstein et al 2009.
%
% 1. Borenstein, M., Hedges, L. V., Higgins, J. P. T. & Rothstein, H. R.
% Introduction to Meta‐Analysis. Principles and Practice of Clinical Trials (Wiley, 2009).
% doi:10.1002/9780470743386.
%
% Currently this can only perform a meta analysis of continuuous data,
% based on Hedges g effect size.
%
% All analyses use a random effects model (i.e. true effect sizes are
% expected to be different across studies).
%
% INPUT
%  6 vectors specifying the number of samples, mean, and standard deviation
%  for group 1 and group 2, in this order n1, m1, sd1, n2, m2, sd2 
%   and then study an indicator variable identifying which study each
%   mean comes from.
% 
% Alternatively, you can pass a single table with columns
% m1,sd1,n1,m2,sd2,n2, and study.
%
% Use the grp input variable (a vector of strings) to specify group membership for a meta
% analysis that includes subgroups. Use the pv.poolT2 input argument to
% specify that the pooled T2 variance should be used (instead of the
% subgroup T2 variance).
% 
% OUTPUT
% meta.nrStudies - number of studies included in the analysis
% meta.name  - The name of the effect size used
% ID - The name of the group of studies (for multi-group analyses only)
% meta.effect:
%       .value - the effect size per study
%       .variance - variance per study
%       .n      -  number of subjects per study
% meta.Q
%       .value - the value of the Q statistic (dispersion across studies)
%       .p     - the p-value associated with the null hypothesis that the dispersion is zero
%       .C      - The C value (See Borenstein)
% meta.T
%       .value -  The value of the T statistic (sqrt of the variance between studies)
%       .ci     - 95% confidence interval for the T statistic
%       .pooled - The pooled variance (only defined for analyses with multiple groups)
% I2
%       .value - The value of the I2 statistic (the fraction of true= heterogeneity across studies relative to the total variance in the observed effects).
%       .ci     - 95% confidence interval
% summary
%       .value - The summary effect size (M)
%       .variance - The varaince of the summary effect size
%       .p      - The p-value associated with the null hypothesis that the summary effect size is zero (by default based on a one-tailed test see pv.tail)
%       .ci - 95% confidence interval of the summary effect size
%       .pi  - 95% prediction interval of the summary effect size (likely range of true effect sizes in a new study)
%
%
% For a group analysis:
%   The first ouput argument is a struct array of meta structs, one for
% each group that has at least minGroupSize observations, plus one (the last item) that
% ignores the subgroups and analyzes all studies together.
%   The second output argument gives the results of an ANOVA over subgroups
% testing the hypothesis that the effects are different across groups:
% anova
%   .qBetween - Sum of squares of Q between groups
%   .df          - Associated degrees of freedom
%   .p   - p-value associated with the null hypothesis that subgroups have the same true effect size. (Chi-squared test)
%   .R2    - R squared ~ explained variance of the anova model
%   .diff
%       .value - Pairwise difference between the effect sizes of all subgroups
%       .se     - Standard error of the pairwise differences
%       .lower - Lower limit of the 95% confidence interval of the pairwise difference
%       .upper  - Upper limit of the 95% confidence interval of the pairwise difference
%
% BK - March 2024
arguments (Repeating)
    varargin
end
arguments
    pv.effectSize (1,1) string {mustBeMember(pv.effectSize,"hedgesg")} = "hedgesg"  % Measure of effect size to use
    pv.tail (1,1) double = 1; % Set to 2 to use a two-tailed test for p-values of the summary effect

    % Options for group analysis only
    pv.group   string = ""  % Vector of string identifiers associating each study with a group.
    pv.minGroupSize (1,1) double = 2 % Only include a group in the ANOVA if it has at least this many studies
    pv.poolT2 (1,1) logical = true; % Use a pooled estimate of the between study variance (T2)
end

% Initialize the output structure
meta = struct('nrStudies',0,'effect',struct,'Q',struct,'T',struct,'I2',struct,'summary',struct,'anova',struct,'study',"");
assert(ismember(numel(varargin),[1 7]),'metaAnalysis requires either a single table or 7 variables specifying the data. See hellp metaAnalysis.')

%% Construct the data table from the input 
if isequal(length(varargin),7) % input is raw statistics, put them into a table
    T = table;
    newname = ["n1", "m1", "sd1", "n2", "m2", "sd2", "study"];
    for i = 1:length(varargin)
        T = addvars(T,varargin{1,i},'NewVariableNames', newname(i));
    end
    assert (all(round(T.n1)==T.n1 & round(T.n2)==T.n2),"n1 and n2 should be integer (number of samples in the study)")
    assert(isstring(T.study),"study should be a string");
elseif isscalar(varargin)
    T = varargin{1,1};
end

%% Determine grouping 
if ~all(pv.group=="")
    totalStudiesBefore = length(unique(T.study));
    assert(height(pv.group)==height(T),"The group input must match the input data (identifying which group each observation belongs to)");
    [G,ID] = findgroups(pv.group);
    nrStudiesInGroup = arrayfun(@(x)length(find(G == x)), unique(G,'stable'), 'Uniform', true);
    enough = nrStudiesInGroup>=pv.minGroupSize;
    if ~all(enough)
        fprintf(2,'%d subgroups have fewer than %d studies to be analyzed and have been removed.\n',sum(~enough),pv.minGroupSize)
    end
    % determine which groups to keep
    validGroups = ID(enough);
    % filter table to keep only valid groups
    stayGroup = ismember(pv.group, validGroups);
    T = T(stayGroup,:);
    pv.group = pv.group(stayGroup);
    % unique studies after pruning
    totalStudiesAfter = length(unique(T.study));
    [G, ID] = findgroups(pv.group);
    nrGroups = numel(ID);
    grpNr = 1:numel(ID);
    % display total number of studies before and after pruning
    fprintf('Total number of studies before pruning: %d\n', totalStudiesBefore);
    fprintf('Total number of studies after pruning: %d\n', totalStudiesAfter);
else
    nrGroups =1;
    grpNr = 1;
    G = ones(height(T),1);
    ID = "not grouped"; % name of the group
end



%% Determine effect size per group and assess heterogeneity
for id = 1:nrGroups
    % Per group
    stay = G==grpNr(id);
    switch upper(pv.effectSize)
        case "HEDGESG"
            [meta(id).effect.value,meta(id).effect.variance,meta(id).effect.n,meta(id).name] = hedgesg(T.m1(stay),T.m2(stay),T.sd1(stay),T.sd2(stay),T.n1(stay),T.n2(stay));
    end
    meta(id).ID =  ID(id);
    meta(id) = ma.heterogeneity(meta(id));
    meta(id).nrStudies = sum(stay);
    meta(id).study = [T.study(stay)];
end

% Add effect sizes and heterogeneity for the ungrouped data as another "group"
if nrGroups>1
    % Also analyze the groups as one ( needed for anova )
    id =nrGroups+1;
    totalStudies = unique(T.study);
    meta(id).nrStudies = length(totalStudies);
    [meta(id).effect.value,meta(id).effect.variance,meta(id).effect.n,meta(id).name] = hedgesg(T.m1,T.m2,T.sd1,T.sd2,T.n1,T.n2);
    meta(id).ID =  "All Studies";
    meta(id).study = "";
    meta(id) = ma.heterogeneity(meta(id));
    % Determine pooled between study variance T2 (option in anova)
    qs = [meta(1:nrGroups).Q];
    Q = sum([qs.value]);
    C = sum([qs.C]);
    df = sum([meta(1:nrGroups).nrStudies]-1);
    pooledT2 = (Q-df)/C;
    for i=1:nrGroups+1
        meta(i).T.pooled =sqrt(pooledT2);
    end
    nrGroups = nrGroups+1;
end


%%  Based on the effect size and heterogeneity, compute the summary effects
for id = 1:nrGroups
    meta(id) = ma.summaryEffect(meta(id),pv.tail,pv.poolT2);
end
%% If requested as an ouput, perform anova of subgroup differences
if nrGroups>1 && nargout >1
    anova = ma.anova(meta,pv.poolT2);
else
    anova =struct; % empty struct
end
end





