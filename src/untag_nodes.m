function [names, tp] = untag_nodes(tagged)
%UNTAG_NODES  Split "NAME_<timepoint index>" back into name and timepoint.
%
%   [names, tp] = UNTAG_NODES(tagged)
%
%   The split is made at the LAST underscore, so protein identifiers that
%   themselves contain underscores (e.g. "RV0001_C") are handled correctly.
%
%   See also TAG_NODES.

tagged = string(tagged(:));
n      = numel(tagged);
names  = strings(n, 1);
tp     = nan(n, 1);

for i = 1:n
    if ismissing(tagged(i)) || tagged(i) == ""
        continue
    end
    idx = strfind(tagged(i), "_");
    if isempty(idx)
        error('untag_nodes:badNode', ...
            'Node id "%s" does not carry a _<timepoint> suffix.', tagged(i));
    end
    p        = idx(end);
    names(i) = extractBefore(tagged(i), p);
    tp(i)    = str2double(extractAfter(tagged(i), p));
end

if any(isnan(tp))
    error('untag_nodes:badNode', 'One or more node ids have a non-numeric suffix.');
end
end
