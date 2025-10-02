function op = struct_to_varargin(s)

assert(isstruct(s) && isscalar(s), "Input must be a struct array of length 1.")

op = arrayfun(@(f_name, f_val) {f_name(:), handle_empty_val(f_val{:})}, string(fieldnames(s)), struct2cell(s), 'UniformOutput', false);

op = horzcat(op{:})';
end

function op = handle_empty_val(f_val)

if isempty(f_val)

    op = gen.empty_like(f_val);

else

    op = f_val;

end

end