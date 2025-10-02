function r = range(arr, dim)

if nargin == 1, dim = 'all'; end

if isnumeric(dim)
    r = cat(dim,min(arr,[],dim), max(arr,[],dim));
elseif strcmp(dim,"all")

    r = [min(arr,[],dim), max(arr,[],dim)];

else
    error("dimension needs to be a numeric value or 'all'.");
end

end