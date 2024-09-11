function meta = summaryEffect(meta,tail,poolT2)
% Compute summary effect using random effects
% Weigh by 1/variance weighted average, plus between study variance T2
% INPUT
% meta - meta data struct
% tail - 1 or 2 tails for the p-value of the summary effect
% poolT2 - Set to true to pool variance across groups
% OUPUT
% meta - an updated meta struct
% BK - Sept 2024
arguments
    meta (1,1) struct
    tail (1,1) double = 1
    poolT2 (1,1) logical = false
end
if poolT2 && isfield(meta.T,'pooled')
    T2 = meta.T.pooled^2;% Use the T2 as estimated for all groups in this meta analysis
else
    T2 = meta.T.value^2; % Use the T2 as estimated for these data
    if poolT2
        fprintf('No pooled T2 computed - using T2 for single group instead\n');
    end
end


% Determine summary effects based on a random effects model (include T2)
% Borenstein page 73.
W = 1./(meta.effect.variance+T2);
M = (W'*meta.effect.value)./sum(W);
varM = 1./sum(W);
seM = sqrt(varM);
ciM =  M + norminv(0.975)*[-1 1].*seM;
z = M./seM;
if tail==1
    % One-tailed test
    pM = 1-normcdf(abs(z));
else
    pM = 2*(1-normcdf(abs(z)));
end
% Prediction interval
df = numel(meta.effect.value)-1;
t  = tinv(0.975,df-1);
piM = M + [-1 1]*t*sqrt(varM+poolT2);


meta.summary.value = M;
meta.summary.variance = varM;
meta.summary.p  = pM;
meta.summary.ci  = ciM;
meta.summary.pi  = piM;


end