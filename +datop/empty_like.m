function emp = empty_like(x)
% Creates an empty instance of the class of x

emp = feval(sprintf("%s.empty",class(x)));

end