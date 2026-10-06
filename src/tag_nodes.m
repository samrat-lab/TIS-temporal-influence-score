function tagged = tag_nodes(names, tp)
%TAG_NODES  Build temporal node identifiers "NAME_<timepoint index>".
%
%   tagged = TAG_NODES(names, tp)
%
%   The temporal network is a layered graph: the same protein at two
%   different timepoints is two distinct nodes. Node identity is therefore
%   the (protein, timepoint-index) pair, encoded as NAME_t.
%
%   See also UNTAG_NODES.

names  = normalize_ids(names);
tagged = names + "_" + string(tp);
tagged(names == "") = "";
end
