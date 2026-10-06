function out = normalize_ids(x)
%NORMALIZE_IDS  Canonical protein / TF identifier form: uppercase, trimmed, column.
%
%   out = NORMALIZE_IDS(x)
%
%   Accepts string, char or cellstr input of any shape; always returns an
%   [n x 1] string with missing values replaced by "". All set operations in
%   this pipeline are performed on normalized identifiers, so a protein that
%   appears with inconsistent case across timepoints is treated as one entity.

out = string(x);
out = out(:);
out(ismissing(out)) = "";
out = upper(strtrim(out));
end
