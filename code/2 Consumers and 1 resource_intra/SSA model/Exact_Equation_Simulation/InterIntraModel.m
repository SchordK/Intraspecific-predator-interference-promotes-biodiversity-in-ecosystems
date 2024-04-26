% This is the File for class InterIntraModel,
% Called by Main.m, to construct an InterIntraModel object in Main.m;
% Here, we define the Gillespie's Algorithm for SSA,
%                   and the ODE Simulation process.
% >>> The initialization of InterIntraModel depends on 2 other classes to work properly :
%                               * class InterIntraParametre   and     * class ReactionChannel,
%       * InterIntraParametre.m are used to define model's parametres
%       * ReactionChannel.m has one important method, and one important property :
%           |--> 1, ReactionChannel.reaction_channel_rates():
%                       the INPUT is a given sysytem phase (system state),
%                       the OUTPUT is an array of reaction rates for every reaction channel,
%                   it calculates every reaction rates of all reaction channels, at a given state.
%           |--> 2, ReactionChannel.channel_table;
%                   records every delta(phase_point) for each reaction channel.
%                   the i'th row of ReactionChannel.channel_table:
%                               ReactionChannel.channel_table(i,:) , is the  delta(phase_point) when i'th reaction occurs.

% >>> the method SSA_Simulation, as its name suggests, takes in two arguments : * begin, and  * step_number;
%       * begin, the initial state from where Stochastic Simulation begins.
%       * step_number, the Stochastic Simulation Step Number, during each step, call the method GillespieAlgorithm.

% >>> the method GillespieAlgorithm, the algorithm developed by Daniel T. Gillespie,
%     determines the reaction which will occur after some time, tau.

% >>> the method ODE_Simulation, calls the Matlab function ode45, and the obj.ode_equation

% >>> the method ode_equation defines the dynamic equation to be simulated by ode45.


% |==========================|
%        InterIntraModel    ||
% |=========================||      parametres       |============================|
%            obj.init       ||  <-------<----------  ||  InterIntraParametre     ||
%                |          ||                       |============================|
%               /|          ||                                      |
%              / |          ||                                      |
%           __/  |          ||                                      |
%          /     |          ||                                      V
%    _____/      |          ||                                      |    parametres
%   /            V          ||                                      |
%   |      SSA_Simulation   ||                                      |
%   |            |          ||                                      V
%   |            |          ||                       |===========================================|
%   V            V          ||                       ||               ReactionChannel           ||---->----->_____>______
%   |  if iterration < step number                   |===========================================|                        \
%   |        |    ^             |                    ||                                         ||                         V
%   V    FALSE    |           TRUE                   ||        reaction_channel_rates()         ||                         |
%   |        |    |             |                    ||      ^       |                          ||                         |
%   V        V    ^             V                    ||      |       |          Reaction        ||                         |
%   |       END   |=============================|    || SystemPhase  |            Rates         ||                         |
%   |             |   Gillespie's Algorithm    ||    ||              |          for EVERY       ||                         V
%   |             |============================||    |===============|           Channel        ||                         |  Build
%   V             ^             |              ||           ^        |===========================|                         |  CHANNEL_TABLE
%   |             |             |              ||           |        |<-----  Channel Width ---->|                         |  according to the structure of
%   |             |             V              ||           |              (Number of Channels)                            V                       phase_point and reaction_channel_rates
%   |             |             |              ||           ^                       |                                      |
%   |             |             V              ||    INPUT  |                       | OUTPUT                               |
%   V             ^             \___>_______________>______/                        V                                      |
%   |             |                    phase_point                                  |                                      V
%   |             |                                                                /                                       |
%   V             ^               __<_________<______-________<______________<____/                                        |
%   |             |              /              channel_rates                                                              V
%   |             |              |                                                                     _____       |==================|
%   V             ^              V                                                                       ^  |--->  ||     ~~~~~      ||   <----- delta(phase_point) when channel_index == 1
%   |             |         According to channel_rates,                                          Channel |  |--->  ||    CHANNEL     ||   <----- delta(phase_point) when channel_index == 2
%   |             |              |                  |                                             _Width |  |--->  ||     TABLE      ||   <----- delta(phase_point) when channel_index == 3
%   |             |   Calculate tau    and   channel_index of which the Reaction occur                   |  |--->  ||     ~~~~~      ||   <----- delta(phase_point) when channel_index == 4
%   V             ^              |                  |                                                    |  |--->  ||     ~~~~~      ||   <----- delta(phase_point) when channel_index == 5
%   |             |              |                  |                                                    :  |--->  ||     ~~~~~      ||   <----- delta(phase_point) when channel_index == 6
%   |             |              |                   \_____________  channel_index ------------------------>| :    ||        :       ||
%   |             |              V                                                                       |  | :    ||        :       ||
%   V             ^      Add tau into time_list                                                          v  |--->  ||     ~~~~~      ||   <----- Δ(phase_point) when channel_index == Channel Width
%   |             |         |                                                                          -----       |==================|
%   |             |         V                                                                                              |
%   |             |         |                                                                                              V
%   |             |         |                               check if consumer will grow                                    |
%   V             |         |                   New_phase  <______ delta(phase_point) ____________________________________/
%   |             |         |                       |        check if death happens
%   |             ^         V                       V
%   |             |         |                       |
%   |             |         |           Add New_phase into time_list
%   V             |         |-------<---------<-----|
%   |             |         |                   ||
%   V             |         V                   ||
%    \             \--<---< iterration + 1      ||
%     \                                         ||
%      \                                        ||
%       \========================================|
%        |====================|
%        |_> ODE_Simulation  ||
%               |            ||
%               V            ||
%               |            ||
%             ode45          ||
%               |            ||
%               V            ||
%               |            ||
%              END           ||
% |===========================|






classdef InterIntraModel
    properties
        consumer_number % the number of consumer species, Sc.
        resource_number % the number of resource kinds, SR.
        parameter       % parameters neccessary for model. See Code Comments in InterIntraParametre.m
        channel         % the reaction channel information needed by Gillespie's Algorithm.
                        % See Code Comments in ReactionChannel.m

        time_array_ssa
        evolution_array_ssa
        time_array_ode
        evolution_array_ode
        % time_array_ssa, evolution_array_ssa, time_array_ode, evolution_array_ode
        % record the simulation result data of SSA_Simulation and ODE_Simulation
    end
    
    methods
        % Initialization of an instance of the obj : InterIntraModel
        function obj = InterIntraModel(consumer_number)
            obj.consumer_number = consumer_number;
            obj.resource_number = 1;
            obj.parameter = InterIntraParametre(consumer_number);
            obj.channel = ReactionChannel(consumer_number, obj.parameter);
        end
        
        function [result_time, result_evoltuion] =...
                SSA_Simulation(obj, begin, step_number)

            time_list = zeros(1, step_number);
            evolution_list = zeros(step_number, numel(begin));
            evolution_list(1, 1:numel(begin)) = begin;
            for time_index = 2:step_number
                [new_evolution, tau] = obj.GillespieAlgorithm(...
                    evolution_list(time_index - 1, ...
                    1:numel(begin)));          % Perform the Gillespie's Algorithm
                
                time_list(time_index) = tau;   % tau (i.e. delta_t), the time needed to the occurrence of the next reaction.
                evolution_list(time_index,...
                    1:numel(begin)) = new_evolution;
            end
            result_time = cumsum(time_list);
            result_evoltuion = evolution_list;
        end
        
        function [new_evolution, tau] = GillespieAlgorithm(...
                obj, phase_point)
                                            % Gillespie's Algorithm

            channel_rates =...
                obj.channel.reaction_channel_rates(phase_point);
                                            % reaction_rates for each reaction_channel

            total_rate = sum(channel_rates);
            r = rand(1, 2);
            tau = log(1 / r(1)) / total_rate;
                                            % tau determined by Gillespie's Algorithm

            
            cumulated_rates = cumsum(channel_rates);
            j = find(cumulated_rates > r(2) * total_rate, 1);
                                            % the next reaction index j determined by Gillespie's Algorithm
            
            channel_output = obj.channel.channel_table(j,:);
                                            % the change of state (i.e. delta_phase)
            
            delta_eaten = ...
                channel_output(end - obj.consumer_number: end - 1);
            if any(delta_eaten < 0)         % if the death of consumer occurs.
                                            % When death occurs, the dead consumer individual takes the eaten resource with it,
                                            % by assumming the eaten resource is unifornally distributed among the consumer induviduals,
                                            % the eaten resource will diminish by the amount of : Eaten_Mass / Ci
                eaten_old = phase_point(end - obj.consumer_number:...
                                            end - 1);

                consumer_free_old = phase_point(1:obj.consumer_number);             % Ci_Free
                pair_old = phase_point(obj.consumer_number + 1: ...
                                end - 2 * obj.consumer_number - 1);                 % Intra_Inter Pairs
                square_old = zeros(obj.consumer_number, ...
                                        obj.consumer_number);
                square_old(triu(true(obj.consumer_number), 0)) = pair_old';
                square_row = sum(square_old, 1);
                square_column = sum(square_old, 2)';                                % square_row + square_column = 2 * y_i + sum_j(z_ij)

                chase_old = phase_point(end - 2 * obj.consumer_number:...
                    end - obj.consumer_number - 1);                                 % Chasing Pair: x_i
                
                consumer_old = consumer_free_old + ...
                                chase_old + ...
                                square_row + square_column;                         % Ci = Ci_Free + x_i + 2 * y_i + sum_j(z_ij)
                            
                delta_eaten(delta_eaten<0) = -1 .* ...
                                eaten_old(delta_eaten<0)...
                                        ./consumer_old(delta_eaten<0) ;             % the eaten resource will diminish by the amount of : eaten_old / Ci
            
                channel_output(end - obj.consumer_number: end - 1) =...
                                delta_eaten;
            end
            
            phase_point_new = phase_point + channel_output;
                                            % phase_point_new = phase_point_old + delta_phase_point

            eat_new = phase_point_new(end - obj.consumer_number: end - 1);
            growth = floor(eat_new ./ ...
                obj.parameter.growth_resource_required);
                                            % Caculate if the consumer will grow, determined by the mass conversion ratio from resource to consumer
                                            % THE Consumer will grow, when eaten_resource reach the amount : growth_resource_required.
                                            % growth_resource_required = 1 / resource_consumer_mass_conversion_ratio

            remained_eat = eat_new - ...
                growth .* obj.parameter.growth_resource_required;
                                            % After consumer growth, the amount of eaten_resource is : remained_eat,
                                            % i.e. mod(eat_new ./ obj.parameter.growth_resource_required)

            phase_point_new(1: obj.consumer_number) = ...
                phase_point_new(1: obj.consumer_number) +  growth;
            phase_point_new(end - obj.consumer_number: end - 1) =...
                remained_eat;
            
            new_evolution = phase_point_new;
                                            % the New Phase Point, after dealing with death and growth.

        end
        
        function [result_time, result_evoltuion] = ...
                ODE_Simulation(obj, begin, time_max)
                                            % the ODE simulation using ode45 in Matlab.
            tspan = [0:0.01:time_max];
            [result_time, result_evoltuion] = ...
                ode45(@obj.ode_equation, tspan, begin);
        end

        function dydt = ode_equation(obj, t, phase_point)

            consumer = phase_point(1: obj.consumer_number);                     % Ci

            intra_inter_pair = phase_point(...
                obj.consumer_number + 1: end - obj.consumer_number - 1);        % y_i and z_ij
                                                                                % According to matlab conventions:
                                                                                % [y_1, z_12, y_2, z_13, z_13, ... y_n]^T

            chase = phase_point(end - obj.consumer_number: end - 1);            % x_i

            resource = phase_point(end);                                        % R

            intra_inter_pair_square = zeros(obj.consumer_number, ...
                                                obj.consumer_number);
            intra_inter_pair_square(...
                triu(true(obj.consumer_number), 0)) = intra_inter_pair;         % reshape y_i and z_ij into a square form :
                                                                                % using n to replace Sc for simplicity:
                                                      % intra_inter_pair_square = [y_1, z_12, z_13, ... , z_1n ;
                                                                                %   0,   y_2, z_13, ... , z_2n ;
                                                                                %   0,    0,   y_3, ... , z_3n ;
                                                                                %   :     :     :    :     :
                                                                                %   0     0     0    0    y_n ]

            consumer_free = ...
                consumer ...
                - sum(intra_inter_pair_square, 1)' ...
                - sum(intra_inter_pair_square, 2) ...
                - chase;                                                        % Ci_Free = Ci -  x_i - 2 * y_i - sum_j(z_ij)

            possible_pair_square = consumer_free * consumer_free';
            intra_inter_pair_square_rate = ...
                obj.parameter.pair_meet .* possible_pair_square ...             % diagonal intra-pair : a'_i * (Ci_Free) ^ 2 ;
                ...                                                             % off-diagonal inter-pair : a'_ij * Ci_Free * Cj_Free,
                - obj.parameter.pair_dessociate .* intra_inter_pair_square;     % diagonal intra-pair : - d'_i * y_i ;
                                                                                % off-diagonal : - d'_ij * z_ij
                                                                                % dy_i/dt = a'_i * (Ci_Free) ^ 2 - d'_i * y_i
                                                                                % dz_ij/dt = a'_ij * Ci_Free * Cj_Free - d'_ij * z_ij


            resource_free = resource - sum(chase);                              % R_Free = R - sum(x_i)

            consumer_rate = (obj.parameter.eat_resource' ...
                ./ obj.parameter.growth_resource_required').* chase ...         % w_i * k_i * x_i
                - obj.parameter.mortality' .* consumer;                         % - D_i * Ci
                                                                                % dCi/dt = w_i * k_i * x_i - D_i * C_i

            intra_inter_pair_rate = intra_inter_pair_square_rate(...
                triu(true(obj.consumer_number), 0));                            % Reshape the upper_triangle part of intra_inter_pair_square_rate,
                                                                                % into a column vector, according to Matlab convention.

            chase_rate = ...
                obj.parameter.encounter' .* consumer_free .* resource_free...   % a_i * Ci_Free * R_Free
                - obj.parameter.escape' .* chase ...                            % - d_i * x_i
                - obj.parameter.eat_resource' .* chase;                         % - k_i * x_i
                                                                                % dx_i/dt = a_i * Ci_Free * R_Free - (d_i + k_i) * x_i


            % Abiotic Resource
            resource_rate = ...
                obj.parameter.resource_supply .* (1 ...
                - resource ./ obj.parameter.resource_capacity) ...              % R_a * (1 - R / K0)
                - sum(obj.parameter.eat_resource' .* chase);                    % - sum_i(k_i * x_i)
                                                                                % dR/dt = R_a * (1 - R / K0) - sum_i(k_i * x_i)


            dydt = ...
                [consumer_rate; intra_inter_pair_rate; ...
                chase_rate; resource_rate];
        end

        function dydt = ode_equation_death_in_pairs(obj, t, phase_point)
            %  Reconsider the death events occured in Chasing Pair (x_i), Intra/Inter Interference Pair (y_i, z_ij)
            %  and the dilution event of resource occured in Chasing Pair (x_i),
            %  compared to  obj.ode_equation, there are modifications on dx_i/dt, dy_i/dt, and dz_ij/dt

            consumer = phase_point(1: obj.consumer_number);                     % Ci

            intra_inter_pair = phase_point(...
                obj.consumer_number + 1: end - obj.consumer_number - 1);        % y_i and z_ij
                                                                                % According to matlab conventions:
                                                                                % [y_1, z_12, y_2, z_13, z_13, ... y_n]

            chase = phase_point(end - obj.consumer_number: end - 1);            % x_i
            resource = phase_point(end);                                        % R
            
            intra_inter_pair_square = zeros(obj.consumer_number, ...
                                                obj.consumer_number);
            intra_inter_pair_square(...
                triu(true(obj.consumer_number), 0)) = intra_inter_pair;
                                                                                % reshape y_i and z_ij into a square form :
                                                                                % using n to replace Sc for simplicity:
                                                      % intra_inter_pair_square = [y_1, z_12, z_13, ... , z_1n ;
                                                                                %   0,   y_2, z_13, ... , z_2n ;
                                                                                %   0,    0,   y_3, ... , z_3n ;
                                                                                %   :     :     :    :     :
                                                                                %   0     0     0    0    y_n ]
            consumer_free = ...
                consumer ...
                - sum(intra_inter_pair_square, 1)' ...
                - sum(intra_inter_pair_square, 2) ...
                - chase;                                                       % Ci = Ci_Free + x_i + 2 * y_i + sum_j(z_ij)

            possible_pair_square = consumer_free * consumer_free';
            intra_inter_pair_square_rate = ...
                obj.parameter.pair_meet .* possible_pair_square ...            % a'_i * (Ci_Free) ^ 2 ;  a'_ij * Ci_Free * Cj_Free
                - obj.parameter.pair_dessociate .* intra_inter_pair_square ... % - d'_i * y_i         ;  - d'_ij * z_ij
                - intra_inter_pair_square .* repmat( ...
                                            obj.parameter.mortality, ...
                                            obj.consumer_number, 1) ...        % - D_i * y_i          ;  - D_i * z_ij
                - intra_inter_pair_square .* repmat( ...
                                            obj.parameter.mortality',...
                                            1, obj.consumer_number);            % - D_i * y_i          ;  - D_j * z_ij      death in interference pair
                                                                                % dy_i/dt = a'_i * (Ci_Free) ^ 2 - d'_i * y_i - 2 * D_i * y_i
                                                                                % dz_ij/dt = a'_ij * Ci_Free * Cj_Free - d'_ij * z_ij

            resource_free = resource - sum(chase);                              % R_Free = R - sum(x_i)
            
            consumer_rate = (obj.parameter.eat_resource' ...
                ./ obj.parameter.growth_resource_required').* chase ...         % w_i * k_i * x_i
                - obj.parameter.mortality' .* consumer;                         % - D_i * Ci
                                                                                % dCi/dt = w_i * k_i * x_i - D_i * C_i


            intra_inter_pair_rate = intra_inter_pair_square_rate(...
                triu(true(obj.consumer_number), 0));                            % Reshape the upper_triangle part of intra_inter_pair_square_rate,
                                                                                % into a column vector, according to Matlab convention.

            chase_rate = ...
                obj.parameter.encounter' .* consumer_free .* resource_free...   % a_i * Ci_Free * R_Free
                - obj.parameter.escape' .* chase ...                            % - d_i * x_i
                - obj.parameter.eat_resource' .* chase ...                      % - k_i * x_i
                - obj.parameter.mortality' .* chase ...                         % - D_i * x_i                       mortality in chasing pair
                - obj.parameter.resource_dilution' .* chase;                    % - R_a / K0 * x_i                  dilution in chasing pair
                                                                                % dx_i/dt = a_i * Ci_Free * R_Free - (d_i + k_i + D_i + R_a/K0) * x_i

            resource_rate = ...
                obj.parameter.resource_supply .* (1 ...
                - resource ./ obj.parameter.resource_capacity) ...              % R_a * (1 - R / K0)
                - sum(obj.parameter.eat_resource' .* chase);                    % - sum_i(k_i * x_i)
                                                                                % dR/dt = R_a * (1 - R / K0) - sum_i(k_i * x_i)
            dydt = ...
                [consumer_rate; intra_inter_pair_rate; ...
                chase_rate; resource_rate];
        end
        
    end
end
