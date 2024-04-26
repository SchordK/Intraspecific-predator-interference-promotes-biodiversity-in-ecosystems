% This is the Parametre-Control File for Chasing_Intra_Inter Interference Pair Model
% Paramtres:
%           encounter(i, l) = a_il
%           escape(i. l) = d_il
%           eat_resource(i, l) = k_il
%           pair_meet(i, j) = a'_ij ;
%                   when i==j , pair_meet(i, i) = a'_i
%           pair_dessociate(i, j) = d'_ij;
%                   when i==j , pair_dessociate(i, i) = d'_i
%           resource_supply = R_a

%  ------------------------------------------------------------------------------------
%  ------------------------------------------------------------------------------------

%         INTERFERENCE_PAIR
%                /|\
%                 |
%                 |  pair_meet
%                 |   a'_ij
%                 |
%                / \
%      _____----/   \----____                       resource_supply , R_a
%     |   pair_dessociate    |                             |
%     |        d'_ij         |                             |
%     V                      V                             V
%   Cj_FREE               Ci_FREE                   RESOURCE_FREE
%                           /|\                            /|\
%                            |_            escape           _|
%                              -----\      d_il      /-----
%                                    \_            _/
%                                      -____   ___-
%                                           \ /
%                                            |
%                                            |  encounter
%                                            |    a_il
%                                            V
%                                      CHASING_PAIR
%                                            |
%                                            |  eat_resource
%                                            |      k_il
%                                            V
%                                     RESOURCE_EATEN
classdef InterIntraParametre
    properties
        encounter
        escape
        pair_meet
        pair_dessociate
        eat_resource
        resource_supply
        resource_capacity
        resource_dilution
        mortality
        growth_resource_required
    end
    
    methods
        function obj = InterIntraParametre(consumer_number)
            a=0.1;
            obj.encounter = ones(1, consumer_number) * a;                   % a_il
            obj.escape = ones(1, consumer_number) * 2.8;                    % d_il

            obj.pair_meet = zeros(consumer_number, consumer_number);
            obj.pair_meet(triu(true(consumer_number), 0)) = ...
                ones(1, consumer_number * (consumer_number + 1) / 2) ...
                * a*1.4;                                                    % a'_ij
                                                                            % [a'_1, a'_12, a'_13, ... , a'_1n ;
                                                                            %   0,    a'_2, a'_13, ... , a'_2n ;
                                                                            %   0,     0,    a'_3, ... , a'_3n ;
                                                                            %   :      :      :     :      :
                                                                            %   0      0      0     0    a'_n ]

            obj.pair_dessociate = zeros(consumer_number, consumer_number);
            obj.pair_dessociate(triu(true(consumer_number), 1))...
                = ones(1, consumer_number * (consumer_number - 1) / 2) ...
                * 0.15;                                                     % d'_ij
            obj.pair_dessociate(logical(eye(consumer_number))) = 0.02;
                                                                            % d'_i
                                                                            % [d'_1, d'_12, d'_13, ... , d'_1n ;
                                                                            %   0,    d'_2, d'_13, ... , d'_2n ;
                                                                            %   0,     0,    d'_3, ... , d'_3n ;
                                                                            %   :      :      :     :      :
                                                                            %   0      0      0     0    d'_n ]

            obj.eat_resource = ones(1, consumer_number) * 0.15;             % k_il

            obj.resource_supply = ones(1) * 100;                            % R_a
            obj.resource_capacity = ones(1) * 600;                          % K_0
            obj.resource_dilution = ...
                obj.resource_supply / obj.resource_capacity;
            %pd = makedist('Uniform', 'lower', 0.95, 'upper', 1.05);
            %obj.mortality = random(pd, [1, consumer_number])* 0.0125;
            
            rng(100)
            obj.mortality= normrnd(1,0.1,1,consumer_number)*0.009;          % D_i
            obj.growth_resource_required = ones(1, consumer_number) / 0.18; % 1 / w_il
            % growth_resource_required is the inverse of the mass_conversion_ratio w_il
        end
    end
end
