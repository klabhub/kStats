%% Data from Borenstein et al. Table 14.1 and used in the worked examples
% of chapter 14 and 18.
% Run this an dcompare with page 94+
m1 = [94 98 98 94 98 96]';
m2 = [92 92 88 82 88 92]';
sd1 = [22 21 28 19 21 21]';
sd2 = [20 22 26 17 22 22]';
n1 = [ 60 65 40 200 50 85]';
n2 = [60 65 40 200 45 85]';
study = ["Carrol" "Grant" "Peck" "Donat" "Stewart" "Young"];
[meta,anova] = ma.analyze(n1,m1,sd1,n2,m2,sd2,study',tail=1, poolT2=false);
ma.plotForest(meta)
title('Figure 14.2 (page 94 Borenstein et al.)')
fprintf('\n *** p 92 Borenstein et al. *** \n  In words, using random-effect weights, the standardized mean difference (Hedges'' g) is %.2f \n with a 95%% confidence interval of %.2f to %.2f. The Z-value is %.2f and the p-value is %f (one-tailed). \n *** \n',meta.summary.value,meta.summary.ci(1),meta.summary.ci(2),meta.summary.value./sqrt(meta.summary.variance),meta.summary.p)
