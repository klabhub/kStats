function m = synthetic(m,opts)
% Given a meta analysis struct, convert the effects to synthetic effects
% based on the options struct
arguments 
    m struct
    opts struct
end

switch upper(opts.mode)
    case "AVERAGE"
        [grpIx,m.study] = findgroups(m.study);
        m.nrStudies= numel(m.study);
        m.effect.value =  splitapply(@(x,w) sum(x.*w)./sum(w),m.effect.value,1./m.effect.variance,grpIx);
        m.effect.variance = splitapply(@(x) weightedVariance(x,opts.r),m.effect.variance,grpIx); 
        m.effect.n = [];%splitapply(@min,m.effect.n,grpIx); 
    otherwise 
        error("Unknown synthetic mode %s",opts.mode)
end

end


function vBar = weightedVariance(v,r)

m = length(v);          % Number of measurements
% We create an m x m matrix where diagonals are 1 and off-diagonals are r
R = ones(m, m) * r;
R(logical(eye(m))) = 1; 

sigmas = sqrt(v);
% 3. Compute the Covariance Matrix (Cov)
cov = diag(sigmas) * R * diag(sigmas);

% 4. Compute the Variance of the Mean (V_bar)
% The sum of all elements in the Covariance Matrix divided by m^2
% gives you (Sum of Vi + Sum of Covariances) / m^2
vBar = sum(cov(:)) / (m^2);

end
