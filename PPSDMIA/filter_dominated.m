function [M,M_idx] = filter_dominated(P,x)
% FILTER_DOMINATED  Extract solutions in P that are dominated by x.
%   Returns the subset M of solutions in P that are dominated by the given
%   solution x (objective space only).
%
%   Inputs:
%     P : SOLUTION array (1×N), candidate set
%     x : SOLUTION, reference solution
%
%   Outputs:
%     M     : SOLUTION array consisting of solutions in P dominated by x
%     M_idx : indices of solutions in P that are dominated by x
%
%   Note:
%     This implementation assumes a minimization problem and checks
%     "x dominates p" by verifying x_obj(m) <= p_obj(m) for all objectives m.

    p_objs = P.objs;
    x_obj = x.obj;

    [N,M] = size(p_objs);
    M_idx = zeros(N,1);
    set_idx = 1;
    for i = 1 : N
        m = 1;
        while m <= M && x_obj(m) <= p_objs(i,m)
            m = m + 1;
        end
        Dominated = m > M;
        if Dominated
            M_idx(set_idx) = i;
            set_idx = set_idx + 1;
        end
    end
    M_idx = M_idx(M_idx ~= 0);
    M = P(M_idx);
end
