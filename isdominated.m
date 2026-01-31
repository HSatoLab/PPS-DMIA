function dominatedFlag = isdominated(targetSol_obj,Pop_obj)
% ISDOMINATED  Check whether a target solution is dominated by a population.
%   dominatedFlag = ISDOMINATED(targetSol_obj, Pop_obj) returns true if there
%   exists at least one solution in Pop_obj that dominates the target
%   solution in the objective space.
%
%   Inputs:
%     targetSol_obj : 1×M vector of objective values of the target solution
%     Pop_obj       : N×M matrix of objective values of candidate dominators
%
%   Output:
%     dominatedFlag : logical scalar, true if targetSol is dominated by Pop
%
%   Note:
%     This implementation assumes minimization and checks domination using
%     a strict comparison (targetSol_obj > Pop_obj for all objectives).
    
    [N,~] = size(Pop_obj);
    tSol_objmat = repmat(targetSol_obj, N, 1);
    isbigmat = zeros(size(Pop_obj));
    dominatedFlag = false;
    
    if N > 0
        isbigmat(tSol_objmat > Pop_obj) = 1;
        isdominated = all(isbigmat, 2);
        dominatedFlag = any(isdominated == 1);
    end
end
