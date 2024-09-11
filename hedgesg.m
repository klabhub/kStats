function [g,varG,n,name] = hedgesg(m1,m2,sd1,sd2,n1,n2)
% Given mean, standard devs, and n for two groups in each study,
% estimate the effect size as a Hedges' g and its variance.
arguments
    m1 (:,1) double % Mean values of the experimental group, per study
    m2 (:,1) double % Mean values of the control group, per study
    sd1(:,1) double {mustBePositive}% Standard deviations of the experimental group, per study
    sd2(:,1) double {mustBePositive} % Standard deviations of the control group, per study
    n1 (:,1) double {mustBePositive,mustBeInteger} % Number of subjects in the experimental group, per study
    n2(:,1)  double {mustBePositive,mustBeInteger} % Number of subjects in the control group, per study
end

% within study variance
sdWithin = sqrt((((n1-1).*sd1.^2 + (n2-1).*sd2.^2))./(sum([n1 n2],2)-2));
% Standardized mean difference
d = (m1-m2)./sdWithin;
% variance of d
Vd = (n1+n2)./(n1.*n2) + d.^2./(2.*(n1+n2));
% Hedges correction factor
J = 1-3./(4.*(n1+n2-2)-1);
% Hedges g and its variance
g = J.*d;
varG = J.^2.*Vd;
name = 'Hedges'' g';
n =n1+n2;
end
