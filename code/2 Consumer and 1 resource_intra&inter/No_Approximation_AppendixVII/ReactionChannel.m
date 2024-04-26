% This is the File for class ReactionChannel,
% called by InterIntraModel, provide detailed infomation on different reaction channels.

% >>> the method reaction_channel_rates(phase_point), takes in the current system phase,
%                   and outputs an array of reaction rates for each channel:
%                   [rate_channel_1, rate_channel_2, ... , rate_channel_end]

% >>> the method channel_table_make() produce the channel_table of the sysytem
%                                        _____       |==================|
%                                          ^  |--->  ||     ~~~~~      ||   <----- delta(phase_point) when channel_index == 1
%                                  Channel |  |--->  ||    CHANNEL     ||   <----- delta(phase_point) when channel_index == 2
%                                   _Width |  |--->  ||     TABLE      ||   <----- delta(phase_point) when channel_index == 3
%                                          |  |--->  ||     ~~~~~      ||   <----- delta(phase_point) when channel_index == 4
%                                          |  |--->  ||     ~~~~~      ||   <----- delta(phase_point) when channel_index == 5
%                                          :  |--->  ||     ~~~~~      ||   <----- delta(phase_point) when channel_index == 6
%                                             | :    ||        :       ||
%                                          |  | :    ||        :       ||
%                                          v  |--->  ||     ~~~~~      ||   <----- Δ(phase_point) when channel_index == Channel Width
%                                        -----       |==================|

classdef ReactionChannel
    properties
        parameter
        consumer_number
        channel_width
        channel_table
    end
    
    methods
        function obj = ReactionChannel(consumer_number, parametre)
            obj.consumer_number = consumer_number;
            obj.parameter = parametre;
            obj.channel_width = ...
                8 * consumer_number + 2 * consumer_number^2 + 2;
            obj.channel_table = obj.channel_table_make();
            
            disp('Reaction Channel Ready');
            disp(size(obj.channel_table));
        end
        
        function CHANNEL_BUS = reaction_channel_rates(obj, phase_point)
            % decode the array phase_point:
            % consumer_free         ==== N
            % intra_inter      ==== N * ( N + 1) / 2
            % chase                 ==== N
            % eaten                   ==== N
            % resource_free         ==== 1
            consumer_free = phase_point(1: obj.consumer_number);        % Ci_Free

            intra_inter = phase_point(obj.consumer_number + 1 : ...
                end - 2 * obj.consumer_number - 1);                     % y_i and z_ij
                                                                        % According to matlab conventions:
                                                                        % [y_1, z_12, y_2, z_13, z_13, ... y_n]^T
 
            chase = phase_point(...
                end - 2 * obj.consumer_number: ...
                end - obj.consumer_number - 1);                         % x_i

            resource_free = phase_point(end);                           % R

            square = zeros(obj.consumer_number, obj.consumer_number);
            square(...
                triu(true(obj.consumer_number), 0)) = intra_inter';     % reshape y_i and z_ij into a square form :
                                                                        % using n to replace Sc for simplicity:
                                                      % intra_inter_pair_square = [y_1, z_12, z_13, ... , z_1n ;
                                                                                %   0,   y_2, z_13, ... , z_2n ;
                                                                                %   0,    0,   y_3, ... , z_3n ;
                                                                                %   :     :     :    :     :
                                                                                %   0     0     0    0    y_n ]

            encounter_channel = ...
                obj.parameter.encounter .* consumer_free * resource_free;       % Encounter : Ci_Free----R_Free
                                                                                %                     |
                                                                                %                    \|/
                                                                                %                     V
                                                                                %                Chasing Pair

            escape_channel = obj.parameter.escape .* chase;                     % Escape    : Ci_Free<--------->R _Free
                                                                                %                       /|\
                                                                                %                        |
                                                                                %                        |
                                                                                %                   Chasing Pair

            possible_pairs = consumer_free' * consumer_free;

            possible_pairs(logical(eye(obj.consumer_number))) = ...
                consumer_free .* (consumer_free - 1);
            interference_pair_channel = ...
                possible_pairs .* obj.parameter.pair_meet;                      % Interference Meet: Ci_Free--------Cj_Free
                                                                                %                               |
                                                                                %                              \|/
                                                                                %                               V
                                                                                %           Intra/Inter Interference Pair

            interference_dessociate_channel = ...
                obj.parameter.pair_dessociate .* square;                        % Interference Dessociate: Ci_Free <--------> Cj_Free
                                                                                %                                      /|\
                                                                                %                                       |
                                                                                %                                       |
                                                                                %                       Intra/Inter Interference Pair

            resource_free_supply_channel = obj.parameter.resource_supply;       % Supply Resource  :
                                                                                %         Outside Environment
                                                                                %                   |
                                                                                %                  \|/
                                                                                %                   V
                                                                                %                R_Free

            resource_free_dilution_channel = ...
                obj.parameter.resource_dilution .* resource_free;               % Dilution R_Free

            resource_chase_dilution_channel = ...
                obj.parameter.resource_dilution .* chase;                       % Dilution Resource in Chasin Pair
                                                                                %             Chasing Pair
                                                                                %                 |
                                                                                %                 |    deterioration/dilution
                                                                                %                 |-------------\
                                                                                %                 |              |
                                                                                %                \|/            \|/
                                                                                %                 V              V
                                                                                %              Ci_Free          Nothing

            mortality_free_channel = ...
                obj.parameter.mortality .* consumer_free;                       % Mortality Ci_Free:   ---X Ci_Free X---

            mortality_chase_channel = obj.parameter.mortality .* chase;         % Mortality Consumer in Chasing Pair:
                                                                                %                   Chasing Pair
                                                                                %                       |
                                                                                %                       |
                                                                                %                       V
                                                                                %    Dead-Ci-Dead <------------>  R_Free

            mortality_row_channel = ...
                square .* repmat( ...
                                obj.parameter.mortality', ...
                                1, obj.consumer_number);                        % Mortality Consumer in Interference Pair:
                                                                                %                   Intra/InterInterference Pair
                                                                                %                             Ci-Cj
                                                                                %                               |
                                                                                %                               |
                                                                                %                               V
                                                                                %            Dead-Ci-Dead <------------>  Cj_Free
            mortality_column_channel = ...
                square .* repmat(...
                                obj.parameter.mortality, ...
                                obj.consumer_number, 1);                        % Mortality Consumer in Interference Pair:
                                                                                %                   Intra/InterInterference Pair
                                                                                %                             Ci-Cj
                                                                                %                               |
                                                                                %                               |
                                                                                %                               V
                                                                                %                Ci_Free <------------> Dead-Cj-Dead

            eat_channel = obj.parameter.eat_resource .* chase;                  % Eat    :      Chasing Pair
                                                                                %                    |
                                                                                %                    |
                                                                                %                   \|/
                                                                                %                    V
                                                                                %            Ci_Free & Eaten resource
            
            CHANNEL_BUS = [...
                        encounter_channel, escape_channel,...
                        eat_channel,...
                        interference_pair_channel(triu(true(obj.consumer_number)))',...
                        interference_dessociate_channel(triu(true(obj.consumer_number)))',...
                        mortality_free_channel,...
                        mortality_chase_channel,...
                        mortality_row_channel(triu(true(obj.consumer_number)))',...
                        mortality_column_channel(triu(true(obj.consumer_number)))',...
                        resource_free_supply_channel,...
                        resource_free_dilution_channel,...
                        resource_chase_dilution_channel];
        end
        
        function table = channel_table_make(obj)
            table = [];
            for channel_index = 1:obj.channel_width
                delta_consumer_free = zeros(1, obj.consumer_number);
                delta_square = zeros(...
                                    obj.consumer_number, ...
                                    obj.consumer_number);
                delta_chase = zeros(1, obj.consumer_number);
                delta_eat = zeros(1, obj.consumer_number);
                delta_resource_free = zeros(1);
                
                switch channel_index
                    case num2cell(1: obj.consumer_number)
                        % Encounter : Ci_Free----R_Free
                        %                     |
                        %                    \|/
                        %                     V
                        %                Chasing Pair
                        delta_chase(channel_index) = delta_chase(channel_index) + 1;
                        delta_resource_free(end) = delta_resource_free(end) - 1;
                        delta_consumer_free(channel_index) = delta_consumer_free(channel_index) - 1;
                    
                    case num2cell(obj.consumer_number + 1: 2 * obj.consumer_number)
                        % Escape    : Ci_Free<--------->R _Free       
                        %                       /|\
                        %                        |
                        %                        |
                        %                   Chasing Pair
                        relative_index = channel_index - obj.consumer_number;
                        delta_chase(relative_index) = delta_chase(relative_index) -1;
                        delta_resource_free(end) = delta_resource_free(end) + 1;
                        delta_consumer_free(relative_index) = delta_consumer_free(relative_index) + 1;

                    case num2cell(2 * obj.consumer_number + 1 : 3 * obj.consumer_number)
                        % Eat    :      Chasing Pair     
                        %                    |
                        %                    |
                        %                   \|/
                        %                    V
                        %            Ci_Free & Eaten resource
                        relative_index = channel_index - 2 * obj.consumer_number;
                        delta_chase(relative_index) = delta_chase(relative_index) - 1;
                        delta_eat(relative_index) = delta_eat(relative_index) + 1;
                        delta_consumer_free(relative_index) = delta_consumer_free(relative_index) + 1;

                    case num2cell(3 * obj.consumer_number + 1 : ...
                                   3 * obj.consumer_number ...
                                   + obj.consumer_number * (obj.consumer_number + 1) / 2)
                        % Interference Meet: Ci_Free--------Cj_Free
                        %                               |
                        %                              \|/
                        %                               V
                        %           Intra/Inter Interference Pair
                        relative_index = channel_index - 3 * obj.consumer_number;

                        index_square = -ones(obj.consumer_number, obj.consumer_number);
                        index_square(triu(true(obj.consumer_number), 0)) = ...
                                                    (1 : obj.consumer_number * ( ...
                                                                                obj.consumer_number + 1) / 2);
                        [channel_position_row, channel_position_column] = ...
                                                        find(index_square == relative_index);

                        delta_square(channel_position_row, channel_position_column) = ...
                                                delta_square(channel_position_row, channel_position_column) + 1;
                        delta_consumer_free(channel_position_row) = delta_consumer_free(channel_position_row) - 1;
                        delta_consumer_free(channel_position_column) = delta_consumer_free(channel_position_column) - 1;

                    case num2cell(3 * obj.consumer_number ...
                                   + obj.consumer_number * (obj.consumer_number + 1) / 2 + 1: ...
                                   3 * obj.consumer_number ...
                                   + obj.consumer_number * (obj.consumer_number + 1))
                        % Interference Dessociate: Ci_Free <--------> Cj_Free
                        %                                      /|\
                        %                                       |
                        %                                       |
                        %                       Intra/Inter Interference Pair
                        relative_index = (channel_index ...
                                  - 3 * obj.consumer_number ...
                                  - obj.consumer_number * (obj.consumer_number + 1) / 2);
                        index_square = -ones(obj.consumer_number, obj.consumer_number);
                        index_square(triu(true(obj.consumer_number), 0)) = ...
                                                    (1 : obj.consumer_number * ( ...
                                                                                obj.consumer_number + 1) / 2);
                        [channel_position_row, channel_position_column] = ...
                                                        find(index_square == relative_index);

                        delta_square(channel_position_row, channel_position_column) = ...
                                                delta_square(channel_position_row, channel_position_column) - 1;
                        delta_consumer_free(channel_position_row) = delta_consumer_free(channel_position_row) + 1;
                        delta_consumer_free(channel_position_column) = delta_consumer_free(channel_position_column) + 1;

                    case num2cell(3 * obj.consumer_number ...
                                   + obj.consumer_number * (obj.consumer_number + 1) + 1 : ...
                                   4 * obj.consumer_number ...
                                    + obj.consumer_number * (obj.consumer_number + 1))
                        % Mortality Ci_Free:   ---X Ci_Free X---
                        relative_index = channel_index ...
                                            - 3 * obj.consumer_number ...
                                            - obj.consumer_number * (obj.consumer_number + 1);
                        delta_consumer_free(relative_index) = delta_consumer_free(relative_index) - 1;
                        delta_eat(relative_index) = delta_eat(relative_index)-1;

                    case num2cell(4 * obj.consumer_number ...
                                    + obj.consumer_number * (obj.consumer_number + 1) + 1 : ...
                                    5 * obj.consumer_number ...
                                    + obj.consumer_number * (obj.consumer_number + 1))
                        % Mortality Consumer in Chasing Pair:   
                        %                   Chasing Pair
                        %                       |
                        %                       |
                        %                       V
                        %       X-Ci-X <------------>  R_Free
                
                        relative_index = channel_index ...
                                            - 4 * obj.consumer_number ...
                                            - obj.consumer_number * (obj.consumer_number + 1);
                        delta_chase(relative_index) = delta_chase(relative_index) - 1;
                        delta_resource_free(end) = delta_resource_free(end) + 1;
                        delta_eat(relative_index) = delta_eat(relative_index) - 1;

                    case num2cell(5 * obj.consumer_number ...
                                    + obj.consumer_number * (obj.consumer_number + 1) + 1 : ...
                                     5 * obj.consumer_number ...
                                        + obj.consumer_number * (obj.consumer_number + 1) * 3 / 2)
                        % Mortality Consumer in Interference Pair:   
                        %                   Intra/InterInterference Pair
                        %                             Ci-Cj
                        %                               |
                        %                               |
                        %                               V
                        %                  X-Ci-X <------------>  Cj_Free
                        relative_index = channel_index ...
                                            - 5 * obj.consumer_number ...
                                            - obj.consumer_number * (obj.consumer_number + 1);

                        index_square = -ones(obj.consumer_number, obj.consumer_number);
                        index_square(triu(true(obj.consumer_number), 0)) = ...
                                                    (1 : obj.consumer_number * ( ...
                                                                                obj.consumer_number + 1) / 2);
                        [channel_position_row, channel_position_column] = ...
                                                        find(index_square == relative_index);

                        delta_square(channel_position_row, channel_position_column) = ...
                                    delta_square(channel_position_row, channel_position_column) - 1;
                        delta_consumer_free(channel_position_column) = delta_consumer_free(channel_position_column) + 1;
                        delta_eat(channel_position_row) = delta_eat(channel_position_row) - 1;
                    case num2cell(5 * obj.consumer_number ...
                                    + obj.consumer_number * (obj.consumer_number + 1) * 3 / 2 + 1 : ...
                                    5 * obj.consumer_number ...
                                    + obj.consumer_number * (obj.consumer_number + 1) * 2)
                        % Mortality Consumer in Interference Pair:   
                        %                   Intra/InterInterference Pair
                        %                             Ci-Cj
                        %                               |
                        %                               |
                        %                               V
                        %                Ci_Free <------------> X-Cj-X
                        relative_index = channel_index ...
                                            - 5 * obj.consumer_number ...
                                            - obj.consumer_number * (obj.consumer_number + 1) * 3 / 2;
                        index_square = -ones(obj.consumer_number, obj.consumer_number);
                        index_square(triu(true(obj.consumer_number), 0)) = ...
                                                    (1 : obj.consumer_number * ( ...
                                                                                obj.consumer_number + 1) / 2);
                        [channel_position_row, channel_position_column] = ...
                                                        find(index_square == relative_index);

                        delta_square(channel_position_row, channel_position_column) = ...
                                    delta_square(channel_position_row, channel_position_column) - 1;
                        delta_consumer_free(channel_position_row) = delta_consumer_free(channel_position_row) + 1;
                        delta_eat(channel_position_column) = delta_eat(channel_position_column) - 1;

                    otherwise
                        relative_index = channel_index ...
                                            - 5 * obj.consumer_number ...
                                            - obj.consumer_number * (obj.consumer_number + 1) * 2;
                        switch relative_index
                            case 1
                                % Supply Resource  :
                                %         Outside Environment
                                %                   |
                                %                  \|/
                                %                   V
                                %                R_Free
                                delta_resource_free(end) = delta_resource_free(end) + 1;
                            case 2
                                % Dilution R_Free :
                                %                  R_Free - 1
                                delta_resource_free(end) = delta_resource_free(end) - 1;
                            otherwise
                                % Dilution Resource in Chasin Pair
                                %             Chasing Pair
                                %                 |
                                %                 |    deterioration/dilution                
                                %                 |-------------\
                                %                 |              |
                                %                \|/            \|/
                                %                 V              V
                                %              Ci_Free          Nothing
                                relative_index = relative_index - 2;
                                delta_chase(relative_index) = delta_chase(relative_index) - 1;
                                delta_consumer_free(relative_index) = delta_consumer_free(relative_index) + 1;
                        end
                end

                delta_phase = [delta_consumer_free, ...
                                delta_square(...
                                        triu(true(obj.consumer_number), 0))', ...
                                delta_chase, ...
                                delta_eat, ...
                                delta_resource_free];

                table(end+1,:) = delta_phase;
            end
        end

    end
end
