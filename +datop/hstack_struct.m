function s = force_mergestruct(varargin)

filler = missing;
isMerge = cellfun(@(x) isstruct(x), varargin);

idx1_opts = find(~isMerge, 1, 'first');

if ~isempty(idx1_opts)
    opts = varargin(idx1_opts:end);

    ii = 1;
    while ii <= length(opts)

        optN = opts{ii};
        ii = ii + 1;
        switch optN

            case {'fill_with'}

                filler = opts{ii};
                ii = ii + 1;
        end

    end


    s = varargin(1:idx1_opts-1);
else
    
    s = varargin;
end


try

    s = [s{:}];

catch ME

    switch ME.identifier

        case 'MATLAB:catenate:structFieldBad'

            

            fld_names = cellfun(@(x) string(fieldnames(x))', s, 'UniformOutput', false);
            fld_names = unique([fld_names{:}]);

            for ii = 1:length(s)

                nan_fields = setdiff(fld_names, fieldnames(s{ii}));
                for fldN = nan_fields

                    s{ii}.(fldN) = filler;

                end

            end

            s = [s{:}];
        otherwise
            rethrow(ME);
    end


end