function meta =  heterogeneity(meta)
% Estimate T^2 - between studies variance using the DerSimonian and Laird
% method.  Borenstein et al. page 72
% Used for random effects analaysis and for assessment of heterogeneity
arguments
    meta (1,1) struct  % A meta struct; output of ma.analyze
end
W =  1./meta.effect.variance;
Q = W'*meta.effect.value.^2- (W'*meta.effect.value).^2./sum(W);
df = numel(meta.effect.value)-1;
C  = sum(W) - sum(W.^2)/sum(W);
T2 = (Q-df)/C;
T = sqrt(T2);
I2 = 100*(Q-df)/Q;
pQ = 1- chi2cdf(Q,df);

% Confidence intervals for T2,T & I2. Page 122 of Borenstein et al.
sw1 = sum(W);
sw2 = sum(W.^2);
sw3 = sum(W.^3);
A =df + 2*(sw1-sw2/sw1)*T^2 + (sw2-2*(sw3/sw1)+ sw2^2/sw1^2)*T^4;
varT2 = 2* (A/C^2);
seT2  =sqrt(varT2);
if Q > df+1
    B = 0.5*(log(Q)-log(df))/(sqrt(2*Q)-sqrt(2*df-1));
else
    B = sqrt(1./(2 * (df-1) * (1-(1/(3*(df-1)^2))) ));
end
L = exp(0.5*log(Q/df)-1.96*B);
U = exp(0.5*log(Q/df)+1.96*B);
ciT2 = df/C *[L^2-1 U^2-1];
ciT2 = max(0,ciT2); % Clamp >=0
ciT = sqrt(ciT2);
ciI2 = 100 *[(L^2-1)/L^2 (U^2-1)/U^2];
ciI2 = max(0,ciI2);

% Combine in struct to return
meta.Q.value = Q;
meta.Q.p  =pQ;
meta.Q.C = C;
meta.T.value = T;
meta.T.ci  =ciT;
meta.I2.value = I2;
meta.I2.ci  =ciI2;

end