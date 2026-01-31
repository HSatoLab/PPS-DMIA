function [filtered_Pop,filtered_Pop_DMcount] = filter_dominated_Pop(Pop_A,Pop_B,Pop_A_DMcount)
% FILTER_DOMINATED_POP  Extract solutions in Pop_A dominated by at least one
% solution in Pop_B (objective space only).
%
%   This function returns the subset of Pop_A that is dominated by Pop_B,
%   where domination is checked using only objective values.
%
%   Inputs:
%     Pop_A          : SOLUTION array (1×N), candidates to be tested
%     Pop_B          : SOLUTION array (1×L), reference set (solutions in Pop_B
%                      are NOT included in the return set)
%     Pop_A_DMcount  : 1×N array, counts how many times each solution in Pop_A
%                      has been selected as a parent for guided mating
%
%   Outputs:
%     filtered_Pop         : SOLUTION array containing Pop_A solutions that
%                            are dominated by at least one solution in Pop_B
%     filtered_Pop_DMcount : corresponding DMcount values for filtered_Pop

    A_objs = Pop_A.objs;
    B_objs = Pop_B.objs;
    [N,M] = size(A_objs);
    [L,~] = size(B_objs);

    % Convert to L*M*N tensors to compare each A(i) against all B solutions
    A_objs_mat = zeros(L, M, N);
    for i = 1:N
        A_objs_mat(:, :, i) = repmat(A_objs(i, :), L, 1);
    end
    B_objs_mat = repmat(B_objs, 1, 1, N);

    R =  A_objs_mat > B_objs_mat;
    S = all(R, 2);
    T = any(S, 1);
    T = squeeze(T);
    filtered_Pop = Pop_A(T);
    filtered_Pop_DMcount = Pop_A_DMcount(T);
end