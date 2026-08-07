function pairN = get_closest_integer_dividers(num)

assert(isscalar(num) & isnumeric(num) & rem(num,1)==0, "num must be an integer!");

pairs = arrayfun(@(x) [x, num/x],1:sqrt(num),'UniformOutput',false);
pairs = vertcat(pairs{:});
pairs = sort(pairs,2);
isInt = ~any(rem(pairs,1),2);
pairs = pairs(isInt,:);
pairN = pairs(gen.absargmin(diff(pairs,[],2)),:);

end