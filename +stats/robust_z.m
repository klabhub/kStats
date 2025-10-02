function varargout = robust_z(arr, dim, var_method, varargin)

%{
Robust Z-score:
Calculated using the median and median absolute deviation (MAD).   
Formula: RobustZ= (x−median)/MAD
​
The median and MAD are less affected by outliers, making the robust Z-score more reliable when dealing with data that may contain extreme values.   

Median: The middle value in a sorted dataset.   
Median Absolute Deviation (MAD): The median of the absolute differences between each data point and the median.
   
%}


if nargin == 1 || isempty(dim)

    dim = 'all';

end


m = median(arr, dim, 'omitnan');

if nargin < 3

    var_method = 'mad';

end


switch var_method

    case 'mad'

        var = median(abs(arr - m), dim, 'omitnan');

    case 'iqr'

        var = iqr(arr, dim);

    case 'std'

        if isempty(varargin)
            scalar = 1/1.349; % iqr calculated from normal standard deviation
        else
            scalar = varargin{1};
        end

        var = scalar * iqr(arr, dim);

    otherwise
        error("Unknown method to compute variance. Pick one of the following: 'mad', 'robust_std', 'iqr'")
end

varargout = cell(1, nargout);
varargout{1} = (arr - m) ./ var;

if nargout > 1, varargout{2} = m; end
if nargout > 2, varargout{3} = var; end

end

