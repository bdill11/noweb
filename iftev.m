function val = iftev(condition,affirmative,negative)
% function val = iftev(condition,affirmative,negative)
%
% Returns affirmative(i) if condition(i) is true (not zero)
% and negative(i) otherwise.
% This is useful for anonymous functions.
aff_idx = find(condition);
neg_idx = find(~condition);
val = zeros(size(condition));
val(aff_idx) = affirmative(aff_idx);
val(neg_idx) = negative(neg_idx);
