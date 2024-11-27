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
% number of initial design points
num_initial = 11*num_vari-1;
% number of maximal evaluations
max_evaluation = 200;
% generate initial design points
sample_x = lower_bound + (upper_bound-lower_bound).*lhsdesign(num_initial, num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name, sample_x, num_obj);
% normalize the objective values
iteration = 0;
evaluation = size(sample_x,1);
kriging_obj = cell(1,num_obj);
% get current non-dominated front
[non_dominated_front,index] = Paretoset(sample_y);
fprintf('pEHVI-EGO on m=%d %s function, iteration: %d, evaluation: %d\n',num_obj,fun_name,iteration,evaluation);
while evaluation < max_evaluation
    for ii = 1:num_obj
        kriging_obj{ii} = Kriging_Train(sample_x,sample_y(:,ii),lower_bound,upper_bound,1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    end
    % get a new point by maximizing pEHVI function using GA
    [infill_x,~] = Optimizer_GA(@(x)-Infill_pEHVI(x,kriging_obj,non_dominated_front),num_vari,lower_bound,upper_bound,40,50);
    % evaluate the new point using the real objective function
    infill_y = feval(fun_name,infill_x,num_obj);
    sample_x = [sample_x;infill_x];
    sample_y = [sample_y;infill_y];
    evaluation = evaluation + size(infill_x,1);
    iteration = iteration + 1;
    [non_dominated_front,index] = Paretoset(sample_y);
    fprintf('pEHVI-EGO on m=%d %s function, iteration: %d, evaluation: %d\n',num_obj,fun_name,iteration,evaluation);
end
if num_obj == 2
    scatter(non_dominated_front(:,1),non_dominated_front(:,2),'ro','filled');
elseif num_obj == 3
    scatter3(non_dominated_front(:,1),non_dominated_front(:,2),non_dominated_front(:,3),'ro','filled');
end




