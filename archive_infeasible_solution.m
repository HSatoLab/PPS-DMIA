function [Population,arch_DMcount] = archive_infeasible_solution(Population,N,arch_DMcount)
% ARCHIVE_INFEASIBLE_SOLUTION  Maintain the infeasible archive A' (size N).
% This function returns a set of *infeasible* solutions selected as follows:
%   1) keep infeasible solutions only (any constraint > 0),
%   2) perform non-dominated sorting in the extended space [objectives, CV],
%      where CV is the overall constraint violation,
%   3) keep only the first front,
%   4) if the archive size exceeds N, truncate it by crowding distance
%      computed in the objective space.
%
% arch_DMcount is filtered/truncated consistently and returned:
%   - arch_DMcount(i) counts how many times the corresponding solution
%     has been selected as a parent for guided mating.
%
% If there is no infeasible solution, this function returns a dummy archive
% of length N (Population filled with SOLUTION()) and arch_DMcount = NaN(1,N).

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is originally written by Wenji Li
% Minor modifications by Ryo Takamiya

    %% Select infeasible solutions
    fIndex           = any(Population.cons > 0,2);
    Population       = Population(fIndex);
    arch_DMcount     = arch_DMcount(fIndex);

    if isempty(Population)
        dummy = SOLUTION();
        Population = repmat(dummy, 1, N);
        arch_DMcount = NaN(1, N);
        return
    else
        % Add overall CV as an additional objective value
        cv_all = overall_cv(Population.cons);
        PopObj_addCon = [Population.objs,cv_all];

        %% Non-dominated sorting
        [FrontNo,~] = NDSort(PopObj_addCon,1);
        Next = (FrontNo == 1);    
        Population = Population(Next);
        arch_DMcount = arch_DMcount(Next);
        if sum(Next) > N
            %% Calculate the crowding distance of each solution
            CrowdDis = CrowdingDistance(Population.objs);
            [~,Rank] = sort(CrowdDis,'descend');
            Population = Population(Rank(1:N));
            arch_DMcount = arch_DMcount(Rank(1:N));
        end
    end
end

% The Overall Constraint Violation
function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end