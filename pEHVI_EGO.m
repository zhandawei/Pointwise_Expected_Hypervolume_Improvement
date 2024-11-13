clearvars;close all;
% objective function
fun_name = 'DTLZ2';
% number of objectives
num_obj = 2;
% number of variables
num_vari = 6;
% lower and upper bounds
lower_bound = zeros(1,num_vari);
upper_bound = ones(1,num_vari);
% Pareto front of the test problem
pareto_front = Calculate_Pareto_Front(1000,num_obj);
% number of initial design points
num_initial = 65;
% number of maximal evaluations
max_evaluation = 150;
% generate initial design points
sample_x = lower_bound + (upper_bound-lower_bound).*lhsdesign(num_initial, num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name, sample_x, num_obj);
% normalize the objective values
sample_y_scaled = (sample_y - min(sample_y))./(max(sample_y)-min(sample_y));
iteration = 0;
evaluation = size(sample_x,1);
kriging_obj = cell(1,num_obj);
% get current non-dominated front
[non_dominated_front,index] = Paretoset(sample_y);
% normalize non-dominated front
non_dominated_front_scaled = sample_y_scaled(index,:);
fprintf('pEHVI-EGO on m=%d %s function, iteration: %d, evaluation: %d, IGD: %f\n',num_obj,fun_name,iteration,evaluation,mean(min(pdist2(pareto_front,non_dominated_front),[],2)));
while evaluation < max_evaluation
    for ii = 1:num_obj
        kriging_obj{ii} = Kriging_Train(sample_x,sample_y_scaled(:,ii),lower_bound,upper_bound,1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    end
    % get a new point by maximizing pEHVI function using GA
    [infill_x,~] = Optimizer_GA(@(x)-Infill_pEHVI(x,kriging_obj,non_dominated_front_scaled),num_vari,lower_bound,upper_bound,40,50);
    % evaluate the new point using the real objective function
    infill_y = feval(fun_name,infill_x,num_obj);
    sample_x = [sample_x;infill_x];
    sample_y = [sample_y;infill_y];
    sample_y_scaled = (sample_y - min(sample_y))./(max(sample_y)-min(sample_y));
    evaluation = evaluation + size(infill_x,1);
    iteration = iteration + 1;
    [non_dominated_front,index] = Paretoset(sample_y);
    non_dominated_front_scaled = sample_y_scaled(index,:);
    fprintf('pEHVI-EGO on m=%d %s function, iteration: %d, evaluation: %d, IGD: %f\n',num_obj,fun_name,iteration,evaluation,mean(min(pdist2(pareto_front,non_dominated_front),[],2)));
end

