function plotForest(meta,pv)
% Generate a Forest plot based on a meta analysis.
%  
% BK  - Sep 2024.
arguments
    meta struct             % The results of a meta analysis  (see ma.analyze)
    pv.maxHeight = 0.35     % Set to <0.5 to fill more or less of the distance between studies
    pv.aspectRatio = 2    %
    pv.sort = false         % Sort by N
    pv.yTickLabel cell = cellstr(meta.study);
end

% Sort by number of particpants in the study.
if pv.sort
    [meta.effect.n,ix]  =sort(meta.effect.n,'descend');
    meta.effect.value = meta.effect.value(ix);
    meta.effect.study= meta.effect.study(ix);
    meta.effect.variance = meta.effect.variance(ix);
end

%% Generate the plot
hold on
% Determine scaling and aspect ratio so that studies are shown as
% squares with the surface proportional to their weight.
nrStudies = numel(meta.effect.value);
widths =sqrt(1./meta.effect.variance);    % Scale by 1/var
yRange = nrStudies+1;
% Determine the range of g we will show
se = sqrt(meta.effect.variance);
mn = min(meta.effect.value-1.96*se);
mx = max(meta.effect.value+1.96*se);
lims = [floor(10*mn)/10 ceil(mx*10)/10];
xRange = -diff(lims);
xlim(lims)
ylim([0.5 nrStudies+1.5]) % One more than nrStudies to show summary effect
if ~isempty(pv.aspectRatio)
    pbaspect([1 pv.aspectRatio 1]);
end
aspectRatio = get(gca,'PlotBoxAspectRatio');
widths = pv.maxHeight *widths/max(widths); % Scale to fit in the space per study
for s= 1:nrStudies
    wY = widths(s);
    wX = aspectRatio(2)*xRange/yRange*wY;
    cX = meta.effect.value(s);
    cY = s;
    ci = se(s)*1.96;
    % Patch centered on the effect size, with surface proportional to
    % weight (1./variance)
    patch([cX-wX cX-wX cX+wX cX+wX],[cY-wY cY+wY cY+wY cY-wY],'k')
    % Line showing confidence limits
    line([cX-ci cX+ci],[cY cY],'Color','k')
end

% Add the summary effect as a diamond with tips showing confidence limits
patch([meta.summary.ci(1) meta.summary.value meta.summary.ci(2) meta.summary.value], nrStudies+ [1 1+pv.maxHeight 1 1-pv.maxHeight ],'k');
% And a line showing the prediction interval
line([meta.summary.p(1) meta.summary.pi(2)],(nrStudies+1)*[1 1],'Color','k')
plot([0 0],ylim,'k')
set(gca,'YDir','reverse','YTick',1:nrStudies,'YTickLabel',pv.yTickLabel);
xlabel(sprintf('Effect Size (%s) \n %s',meta.name,meta.ID))
title(sprintf('M = %.2g [%.2g; %.2g] (p=%.3g).\n Q = %.2g (p=%.3g). \n T = %.2g, I^2=%.0f%%' ,meta.summary.value,meta.summary.ci(1),meta.summary.ci(2),meta.summary.p,meta.Q.value,meta.Q.p,meta.T.value,meta.I2.value))
end