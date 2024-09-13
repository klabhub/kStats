function m = recodeDummyVars(m,dummyvar)
% Recode dummy vars to effects or reference
arguments
    m
    dummyvar (1,1) string {mustBeMember(dummyvar,["effects" "reference"])}
end

m = fitlme(m.Variables,char(m.Formula), 'exclude',m.ObservationInfo.Excluded,'dummyvarcoding',dummyvar);
