function [filtered_pop,idx_filtered_pop] = filter_cv(pop,cv)
% FILTER_CV  Filter solutions by overall constraint violation (CV).
%   Selects solutions in 'pop' whose overall constraint violation is
%   strictly smaller than the given threshold 'cv'.
%
%   Inputs:
%     pop : SOLUTION array (1×N) to be filtered
%     cv  : scalar threshold of overall constraint violation
%
%   Outputs:
%     filtered_pop     : SOLUTION array consisting of selected solutions
%     idx_filtered_pop : indices of selected solutions in the original 'pop'

pop_cv = overall_cv(pop.cons);
idx_filtered_pop = find(pop_cv < cv);
filtered_pop = pop(idx_filtered_pop);
end

% The Overall Constraint Violation
function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end
