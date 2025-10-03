function varargout = apply_func_along_dimension(arr, iDim, func, varargin)
dim_size = size(arr);
num_dims = ndims(arr);

% --- Input Validation ---
if iDim > num_dims || iDim <= 0 || ~isnumeric(iDim) || fix(iDim) ~= iDim
    error('iDim must be a positive integer less than or equal to the number of dimensions of arr.');
end
if ~isa(func, 'function_handle')
    error('func must be a function handle.');
end

% --- Determine Output Size (and preallocate) ---
% Apply the function to the *first* slice to determine the output size.
first_slice_indices = repmat({1}, 1, num_dims);
first_slice_indices{iDim} = ':';
first_slice = squeeze(arr(first_slice_indices{:}));
outputN = cell([1,nargout]);
[outputN{:}] = func(first_slice, varargin{:});
output_size_iDim = cellfun(@(x) numel(x), outputN);
% Calculate the size of the output array
op_size = repmat(dim_size,nargout,1);
op_size(:,iDim) = output_size_iDim;
% Preallocate the output array.
varargout = arrayfun(@(i) zeros(op_size(i,:), 'like', outputN{i}), 1:nargout, 'UniformOutput',false);

% --- Prepare for Looping ---
% Create a cell array to index all dimensions *except* iDim.
slice_idx = arrayfun(@(x) 1:x, dim_size, 'UniformOutput',false);
slice_idx{iDim} = ':';
slice_idx = table2cell(combinations(slice_idx{:}));

for iSlice = 1:size(slice_idx,1)

    [outputN{:}] = func(squeeze(arr(slice_idx{iSlice,:})), varargin{:});

    for iOp = 1:nargout

        varargout{iOp}(slice_idx{iSlice,:}) = outputN{iOp};

    end
end

end