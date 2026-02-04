classdef PPSDMIA < ALGORITHM
% <2025> <multi/many> <real/integer> <constrained>
% Push and pull search with Directed Mating and Infeasible Archive algorithm
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% R. Takamiya, M. Miyakawa, K. Takadama and H. Sato, PPS-DMIA: Directed
% Mating With Infeasible Solution Archive in Push and Pull Search for
% Constrained Multi-Objective Optimization, in IEEE Access, vol. 13,
% pp. 211302-211316, 2025, doi: 10.1109/ACCESS.2025.3643342.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ryo Takamiya

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [delta,nr] = Algorithm.ParameterSet(0.9,2);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population and initialize archive
            Population = Problem.Initialization();
            arch       = archive_DM(Population,Problem.N,zeros(1,size(Population,2)));
            arch_Con   = archive_infeasible_solution(Population,Problem.N,zeros(1,size(Population,2)));
            Z          = min(Population.objs,[],1);
            arch_DMcount     = zeros(1,size(arch,2));
            arch_Con_DMcount = zeros(1,size(arch_Con,2));
            cntLimit         = 1; % Maximum number of times each archive solution can be used for guided mating

            %% Evaluate the Population
            Tc               = 0.8 * ceil(Problem.maxFE/Problem.N);
            last_gen         = 120;
            change_threshold = 1e-3;
            search_stage     = 1; % 1 for push stage,otherwise,it is in pull stage.
            max_change       = 1;
            epsilon_k        = 0;
            epsilon_0        = 0;
            cp               = 2;
            alpha            = 0.95;
            tao              = 0.1;
            ideal_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            nadir_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen        = ceil(Problem.FE/Problem.N);
                pop_cons   = Population.cons;
                cv         = overall_cv(pop_cons);
                population = [Population.decs,Population.objs,cv];
                rf         = sum(cv <= 1e-6) / Problem.N;
                ideal_points(gen,:) = Z;
                nadir_points(gen,:) = max(population(:,Problem.D + 1 : Problem.D + Problem.M),[],1);

                % The maximumrate of change of ideal and nadir points rk is calculated.
                if gen >= last_gen
                    max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen);
                end

                % The value of e(k) and the search strategy are set.
                if gen < Tc
                    if max_change <= change_threshold && search_stage == 1
                        search_stage = -1;
                        epsilon_0 = max(population(:,end),[],1);
                        epsilon_k = epsilon_0;
                    end
                    if search_stage == -1
                        epsilon_k =  update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp);
                    end
                else
                    epsilon_k = 0;
                end

                % For each solution
                Q = SOLUTION.empty;
                rp = randperm(Problem.N);
                for i = 1 : Problem.N
                    j = rp(i);
                    if rand < delta
                        P = B(j,randperm(size(B,2)));
                    else
                        P = randperm(Problem.N);
                    end
                    p_1 = Population(P(1));
                    p_2 = Population(P(2));

                    % Generate an offspring
                    if search_stage == -1
                        % Extract candidate set M from A ∪ A' for the current subproblem j
                            % From feasible archive A
                            [~,idx_M_fromA] = filter_dominated(arch,Population(j));
                            idx_M_fromA = idx_M_fromA(arch_DMcount(idx_M_fromA) < cntLimit);
                            [M_fromA,idx_tmp] = filter_cv(arch(idx_M_fromA),overall_cv(Population(j).con));
                            idx_M_fromA = idx_M_fromA(idx_tmp);
                            cv_M_fromA = overall_cv(M_fromA.cons);
                            % From infeasible archive A'
                            [~,idx_M_fromA_prime] = filter_dominated(arch_Con,Population(j));
                            idx_M_fromA_prime = idx_M_fromA_prime(arch_Con_DMcount(idx_M_fromA_prime) < cntLimit);
                            [M_fromA_prime,idx_tmp] = filter_cv(arch_Con(idx_M_fromA_prime),overall_cv(Population(j).con));
                            idx_M_fromA_prime = idx_M_fromA_prime(idx_tmp);
                            cv_M_fromA_prime = overall_cv(M_fromA_prime.cons);
                        M = [M_fromA, M_fromA_prime];
                        cv_M = overall_cv(M.cons);
                        if length(M) >= 1
                            % Choose candidates with the minimum constraint violation as M'
                            % (performed separately for A and A', then merged)
                            idx_M_prime_list_fromA = idx_M_fromA(cv_M_fromA == min(cv_M));
                            idx_M_prime_list_fromA_prime = idx_M_fromA_prime(cv_M_fromA_prime == min(cv_M));
                            M_prime = [arch(idx_M_prime_list_fromA), arch_Con(idx_M_prime_list_fromA_prime)];
                            % From M', select the one with the smallest g
                            g_M_prime = max(abs(M_prime.objs-repmat(Z,length(M_prime),1)).*W(j,:),[],2);
                            [~, idx_p_1] = min(g_M_prime);
                            if idx_p_1 <= length(idx_M_prime_list_fromA)
                                idx_p_1 = idx_M_prime_list_fromA(idx_p_1);
                                p_1 = arch(idx_p_1);
                                arch_DMcount(:,idx_p_1) = arch_DMcount(:,idx_p_1) + 1;
                            else
                                idx_p_1 = idx_M_prime_list_fromA_prime(idx_p_1-length(idx_M_prime_list_fromA));
                                p_1 = arch_Con(idx_p_1);
                                arch_Con_DMcount(:,idx_p_1) = arch_Con_DMcount(:,idx_p_1) + 1;
                            end
                            Offspring = OperatorDE(Problem,Population(j),p_1,Population(j));
                        else
                            Offspring = OperatorDE(Problem,Population(j),p_1,p_2);
                        end
                    else
                        Offspring = OperatorDE(Problem,Population(j),p_1,p_2);
                    end
                    Q(1,j) = Offspring;

                    % Update the ideal point
                    Z = min(Z,Offspring.obj);

                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),length(P),1).*W(P,:),[],2);
                    cv_old = overall_cv(Population(P).cons);
                    cv_new = overall_cv(Offspring.con) * ones(length(P),1);
                    % Check whether the offspring is dominated by the feasible archive
                    if search_stage == -1
                        isdominated_offspring = isdominated(Offspring.obj,arch.objs);
                    end

                    if search_stage == 1 % Push Stage
                        Population(P(find(g_old>=g_new,nr))) = Offspring;
                    else  % Pull Stage  &&  An improved epsilon constraint-handling is employed to deal with constraints
                        Population(P(find(((g_old >= g_new) & (((cv_old <= epsilon_k) & (cv_new <= epsilon_k)) | (cv_old == cv_new)) | ((cv_new < cv_old) & (~isdominated_offspring))), nr))) = Offspring;
                    end
                end

                % Output the non-dominated and feasible solutions.
                [arch,arch_DMcount] = archive_DM([arch,Q],Problem.N,[arch_DMcount,zeros(1,size(Q,2))]);
                [R,R_DMcount] = filter_dominated_Pop([Q,arch_Con],Population,[zeros(1,size(Q,2)),arch_Con_DMcount]);
                [arch_Con,arch_Con_DMcount] = archive_infeasible_solution(R,Problem.N,R_DMcount);

                if Problem.FE >= Problem.maxFE
                    Population = arch;
                end
            end
        end
    end
end

% The Overall Constraint Violation
function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end

% Calculate the Maximum Rate of Change
function max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen)
    delta_value = 1e-6 * ones(1,size(ideal_points,2));
    rz = abs((ideal_points(gen,:) - ideal_points(gen - last_gen + 1,:)) ./ max(ideal_points(gen - last_gen + 1,:),delta_value));
    nrz = abs((nadir_points(gen,:) - nadir_points(gen - last_gen + 1,:)) ./ max(nadir_points(gen - last_gen + 1,:),delta_value));
    max_change = max([rz, nrz]);
end

function result = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp)
    if rf < alpha
        result = (1 - tao) * epsilon_k;
    else
        result = epsilon_0 * ((1 - (gen / Tc)) ^ cp);
    end
end
