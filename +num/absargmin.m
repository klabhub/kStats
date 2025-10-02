function i = absargmin(varargin)
if nargin == 1
%     varargin = {varargin{1},1,varargin{2:end}};
% else
    varargin{2} = 1;
end
% [~, i] = min(abs(varargin{1}),varargin{2:end});
[~, i] = mink(abs(varargin{1}),varargin{2:end});

end