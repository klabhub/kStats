%% Example 2 - Borenstein et al Chapter 19.
%
% Compare these results with the values and graphs in that chapter.

%% Data From Table 19.1 - testing grouped study comparisons (only g/variance provided)
meta = struct('nrStudies',0,'effect',struct,'Q',struct,'T',struct,'I2',struct,'summary',struct,'anova',struct,'study',"");
% Group A
meta.effect.value = [0.11 0.224 0.338 0.451 0.480]';
meta.effect.variance = [0.01 0.03 0.02 0.015 0.01]';
meta.nrStudies = 5;
meta.study = ["Thornhill" "Kendall" "Vandamm" "Leonard" "Professor"]';
meta.name = "hedgesg";
meta.ID  ="ID";
% Group B
meta(2).effect.value = [0.44 0.492 0.651 0.710 0.740]';
meta(2).effect.variance = [0.015 0.02 0.015 0.025 0.012]';
meta(2).nrStudies = 5;
meta(2).study = ["Jeffries" "Fremont"  "Doyle" "Stella" "Thorwald"];

% Combined A+B
meta(3).effect.value = [meta(1).effect.value;meta(2).effect.value];
meta(3).effect.variance = [meta(1).effect.variance;meta(2).effect.variance];
meta(3).nrStudies = 10;
meta(3).study = ["Thornhill" "Kendall" "Vandamm" "Leonard" "Professor" "Jeffries" "Fremont"  "Doyle" "Stella" "Thorwald"];
%% Compute

% Settings
% To pool the variance across groups set this to true (Page 171+).
% With poolT2 = false, each group uses its own estimate of across-studies
% variance (as done on page 164+)
poolT2 =false;
% Use 1-tailed tests
tail =1;

for i=1:3
    % Determine variance/heterogeneity within each group 
    meta(i)= ma.heterogeneity(meta(i));
    % Determine summar effect (using a random effects model - i.e. including across study variance T2) 
    meta(i) = ma.summaryEffect(meta(i),tail,poolT2);
end

% Show Forest plots 
figure(1);clf;
layout = tiledlayout(2,1);
fprintf('\n\n\t g \t variance \t Q \n');
for i=1:2
    nexttile;
    ma.plotForest(meta(i));
    fprintf('\t %.3f \t %.3f  \t %.3f\n',meta(i).summary.value,meta(i).summary.variance,meta(i).Q.value)
end
linkaxes(gcf().Children().Children())
title(layout,'Figure 19.1');

% Compare groups in an anova
anova = ma.anova(meta,poolT2);
fprintf('\n\n\t Q*: %.3f %.3f \n',anova.qStar(1),anova.qStar(2));
fprintf('\n\n\t Qbetween: %.3f  p = %.3f\n',anova.qBetween,anova.p);
