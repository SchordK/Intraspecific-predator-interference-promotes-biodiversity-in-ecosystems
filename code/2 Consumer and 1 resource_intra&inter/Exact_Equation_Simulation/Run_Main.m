% This is the Main Programme for Chasing_Intra_Inter Interference Pair Model :
%       1, SSA Simulation for Chasing_Intra_Inter Interference Pair Model
%       2, ODE Simulation for Chasing_Intra_Inter Interference Pair Model
%       3, Plot Simulation data

% For bith SSA and ODE Simulation, the Main.m calls the respective methods defined in InterIntraModel.m
% NOTE : the Model Parametres are controlled in InterIntraParametre.m
%
% >>> SSA SIMULATION PROCEDURE ---->
%       1, construct an instance of the class: InterIntraModel
%       2, call the method : SSA_Simulation(begin, simulation_step_number)
%       3, obtain the return value : [time_ssa, evolution_ssa]
%           time_ssa : with size (1, simulation_step_number)
%           evolution_ssa : with size (simulation_step_number, length_of_begin)

% >>> ODE SIMULATION PROCEDURE ----->
%       call the method : ODE_Simulation, consistant with Matlab function ode45()

% >>> PLOT SSA PROCUDURE ---->
%      Collecte the Total Consumer Abundance in evolution_ssa:
%               contains the data of C_Free, Interference_Pair, Chasing_Pair and R_Free
%           1, the data of C_Free, Chasing_Pair and R_Free will be easily handled.
%           2, the data of Interference_Pair will first be reshaped into a 3D array form: square_evolution
%              Suppose, there exists n kinds of Consumer Species.
%              And square_evolution(t) = [C1C1, C1C2. C1C3, ...... C1Cn;
%                                           0,  C2C2, C2C3, ...... C2Cn;
%                                           0,    0,  C3C3, ...... C3Cn;
%                                           :     :     : ,   :      :
%                                           0     0     0,  ...... CnCn]
%                   where we use CiCj denote the interfernce pair between Ci and Cj.
%                                   Naturally, if i==j, CiCi is the Intra-Interference Pair.
%                   Apparently, to calculate the Abundance of Consumer involved in Interference, at t'th moment :
%                                   just simply cuculate the ROW_SUM + COLUMN_SUM of square_evolution(t)
%           3, Collected Data stored in consumer_ssa and resoure_ssa
%clear
tic
consumer_diversity = 2;
simulation_step_number = 1e8;           % simulation steps for SSA

eco = InterIntraModel(consumer_diversity);

% begin_ssa = [consumer_free, intra_inter, chase, eaten, resource_free]
begin_ssa = zeros(1, 3 * consumer_diversity + ...
    consumer_diversity * (consumer_diversity + 1) / 2 + 1);
begin_ssa(1:consumer_diversity) = 20;       % initial state of consumer free, C_Free
begin_ssa(end) =50;                         % initial state of resource free, R_Free
                                         % SSA : evolution_ssa
                                        %  t0   [C_free, CiCj, Chasing_pair, Eaten-Resource, Resource;
                                        %  t1    C_free, CiCj, Chasing_pair, Eaten-Resource, Resource;
                                        %   :
                                        %   :
                                        % t_end  C_free, CiCj, Chasing_pair, Eaten-Resource, Resource]
[time_ssa, evolution_ssa] = ...
    eco.SSA_Simulation(begin_ssa, simulation_step_number);
                                        % SSA Simulation

% begin_ode = [consumer_free, intra_inter, chase, resource_free]
begin_ode = zeros(1, 2 * consumer_diversity + ...
    consumer_diversity * (consumer_diversity + 1) / 2 + 1);
begin_ode(1:consumer_diversity) = 20;       % initial state of consumer, C
begin_ode(end) = 50;                        % initial state of resource, R
options=odeset('MaxStep',0.01);
[time_ode, evolution_ode] = ...
    eco.ODE_Simulation(begin_ode, time_ssa(end));
                                        % ODE Simulation

consumer_free_evolution = evolution_ssa(:, 1:consumer_diversity);
                                 % Consuemr_Free in SSA, Ci_Free


resource_free_evolution = evolution_ssa(:, end);
                                 % Resource_Free in SSA, R_Free


pair_evolution = evolution_ssa(:, ...
            consumer_diversity + 1: end - 2 * consumer_diversity - 1);
                                 % Intra/ Inter - Interference Pair in SSA
        % Reshape the Interference Pair into a 3D Array : square_evolution,
        %                        with shape (length_of_time_ssa, consumer_diversity, consumer_diversity)
square_evolution = zeros(length(time_ssa), consumer_diversity, consumer_diversity);
[i, j, k] = meshgrid(1:consumer_diversity,...
    1:length(time_ssa),...
    1:consumer_diversity);
condition = i <= k;
square_evolution(condition) = reshape(pair_evolution, [1,numel(pair_evolution)]);
        % Fill the simulation data of all Intra_Inter Interference Pairs into square_evolution.


inter_pair_collection = squeeze(sum(square_evolution, 2)) +...
    sum(square_evolution, 3);
        %  Calculate the Abundance of Consumer involved in Interference Pair,  2*y_i + sum_j(z_ij)


chase_evolution = evolution_ssa(:, ...
    end - consumer_diversity*2 :end - consumer_diversity - 1);
        %  Calculate the Abundance of Consumer involved in Chasing Pair, x_i


consumer_ssa = consumer_free_evolution + chase_evolution + inter_pair_collection;
        %  Calculate the Total Abundance for each Consumer: Ci = Ci_Free + x_i + 2*y_i + sum_j(z_ij)


resource_ssa = resource_free_evolution + sum(chase_evolution, 2);
        %  Calculate the Total Abundance of Resource : R = R_Free + x_i



figure
hold on
plot(time_ssa', consumer_ssa)
plot(time_ssa', resource_ssa,'r','linewidth',3)
plot(time_ode, evolution_ode(:,1:consumer_diversity))
plot(time_ode, evolution_ode(:,end),'r--','linewidth',3)
xlabel('time')
ylabel('Consumers and Resource')
set(gca,'XScale','log')
set(gca,'YScale','log')

% --------------------------------------------------------------------------------------
% we can also save the SSA_simulation data, and ODE_data as .mat file in Folder ./data :

% save('./data/time_ssa.mat','time_ssa', '-v7.3')
% save('./data/abundance_ssa.mat','evolution_ssa','-v7.3')
% save('./data/consumer_ssa.mat','consumer_ssa', '-v7.3')
% save('./data/resource_ssa.mat','resource_ssa','-v7.3')

% save('./data/time_ode.mat','time_ode', '-v7.3')
% save('./data/abundance_ode.mat','evolution_ode', '-v7.3')
% --------------------------------------------------------------------------------------

toc
