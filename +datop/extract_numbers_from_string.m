function varargout = extract_numbers_from_string(x)
% Extracts all sequences of digits from a string.
%   numericVars = extract_numbers_from_string(x) finds all sequences of one or more
%   digits in the input string x and returns them as a numeric array.
%
%   [numericVars, isNumeric] = extract_numbers_from_string(x) also returns a logical
%   array 'isNumeric' the same size as x, where true indicates that the
%   corresponding character in x is a digit.
%
%   Args:
%     x: The input string.
%
%   Returns:
%     numericVars: A numeric array containing the numbers found in the string.
%                  Returns an empty array if no numbers are found.
%     isNumeric:   (Optional) A logical array indicating the position of digits.

varargout = cell(1, nargout);

varargout{1} = str2double(regexp(x, '\d+', 'match'));

if nargout > 1
    varargout{2} = isstrprop(x,'digit');
end
end