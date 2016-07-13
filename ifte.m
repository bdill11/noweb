function val = ifte(condition,affirmative,negative)
% function val = ifte(condition,affirmative,negative)
%
% Returns affirmative if condition is true (not zero)
% and negative otherwise.
% This is useful for anonymous functions.
if condition
    val = affirmative;
else
    val = negative;
end
