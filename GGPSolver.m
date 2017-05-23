function [Obj_Optimal, Var_Optimal] = GGPSolver(NumOfVar, Obj, Constraint, LowerBound, UpperBound)
%% General Geometric Programming
% apply the algorithm in
% http://maranas.che.psu.edu/pub/1997/
%Maranas_and_Floudas,_Computers_and_Chem._Eng.,_1997.pdf
% Duo, Hong Kong University of Science and Technology
% dhanaa@connect.ust.hk 
% 2017/05/22


%% Input: 
% NumOfVar: number of variables to be optimized over
% Obj: objective function
% Constraint: a number of constaints
% LowerBound: global lower bounds for all variables
% UpperBound: global upper bounds for all variables

tic
cvx_quiet(true); % shut down cvx console
itermax = 1000; % set a limit for iteration cycles
NumofConstraint = size(Constraint,2); % number of constraints

%% Step 1 initialization
t_l = LowerBound; % global lower bound for the underlying variable
t_u = UpperBound; % globalupper bound
z_l = log(t_l); % transformed bounds
z_u = log(t_u);
z = (z_u-z_l) * rand(1) + z_l; % a random init value
% compute G0_u in DC in page 3
[G0_u,~] = DCproblem(Obj,Constraint,z);
G0_l = 0;
G0_data = G0_u;
z_data = z;
z_l_data = z_l;
z_u_data = z_u;

z_l_temp1 = z_l;
z_l_temp2 = z_l;

z_u_temp1 = z_u;
z_u_temp2 = z_u;

% optimization parameters
T_f = 0.000001;
T_c = 0.0001;
iter = 0; % iteration counter
formatSpec = ' %3$s %2$s %1$s';

% start the iteration
while ((G0_u-G0_l)>T_c)&&(iter<itermax)
    iter = iter + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% step 2 feasibility check and Update upper bound of obj G0_u
    % compute the obj and  constraint based on last iteration z
    [G0_u_LastStep,ConLastStep] = DCproblem(Obj,Constraint,z);
    % check and update G0_u
    if (all(ConLastStep<T_f))
        G0_u = min( G0_u_LastStep, G0_u );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% step 3 partition the current rectangle
    [~,long_index] = max(z_u-z_l); % choose the longest side
    z_l_temp1 = z_l; % box 1
    z_u_temp1 = z_u;
    z_u_temp1(long_index) = (z_u(long_index)+z_l(long_index))/2;
    
    z_l_temp2 = z_l; % box 2
    z_u_temp2 = z_u;
    z_l_temp2(long_index) = (z_u(long_index) + z_l(long_index))/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% step 4 feasibility check of convex relaxazation
    % partition 1, lower rectangle
    [G0_L_temp,ineq_constr] = Rproblem(Obj,Constraint, z_l_temp1, z_u_temp1);
    % check feasibility
    if (G0_L_temp >= G0_u || ...
            any(ineq_constr >= 0))
        ind_1 = 0;
    else
        %we can run cvx on parition 1
        % step 5 cvx inside rectangle (rectangle 1)
        [ObjOpt, zOpt] = CVX_R_problem(NumOfVar, Obj, Constraint, z_l_temp1, z_u_temp1);
        if ObjOpt<G0_u
            G0_data=[G0_data, ObjOpt];
            z_data=[z_data; zOpt];
            z_l_data=[z_l_data ; z_l_temp1];
            z_u_data=[z_u_data; z_u_temp1];
        end
        ind_1 = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % partition 2, upper rectangle
    [G0_L_temp2,ineq_constr2] = Rproblem(Obj,Constraint, z_l_temp2, z_u_temp2);
    % check feasibility
    if (G0_L_temp2 >= G0_u || ...
            any(ineq_constr2 >= 0))
        ind_2 = 0;
    else
        %we can run cvx on parition 2
        % step 5 cvx inside rectangle (rectangle 1)
        [ObjOpt, zOpt] = CVX_R_problem(NumOfVar, Obj, Constraint, z_l_temp2, z_u_temp2);
        if ObjOpt<G0_u
            G0_data=[G0_data, ObjOpt];
            z_data=[z_data ; zOpt];
            z_l_data=[z_l_data ; z_l_temp2];
            z_u_data=[z_u_data ; z_u_temp2];
        end
        ind_2 = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% setp 6 Update lower bound G0_l and opt point
    [G0_l,indx] = min(G0_data);
    z = z_data(indx, :);
    z_l = z_l_data(indx, :);
    z_u = z_u_data(indx, :);
    %G0_u
    %X_temp
    G0_data(indx)=[];
    z_data(indx, :)=[];
    z_l_data(indx, :)=[];
    z_u_data(indx, :)=[];
    
    % compute true obj from last iteration
    [trueObj, trueCon] = DCproblem(Obj, Constraint, z);
    
    %% display iteration information
    if iter==1
        fprintf('Starting the GGP algorithm...\n')
        fprintf('Note: when the algorithm ends, Max(con) should be almost negative 0 \n')
        fprintf('%4s %12s %12s %12s \r\n', 'Iter','Objective',  'Max(Con)', 'Gap_Between_Iteration')
    end
    fprintf('%4d %12.8f %12.8f %12.8f \n', (iter),(trueObj),(max(trueCon)),((G0_u-G0_l)));
end
% return the optimal value needed
Obj_Optimal = trueObj;
Var_Optimal = exp(z);
toc
% display optimization summary
ObjVal_opt = trueObj;
Variables_opt = Var_Optimal;
Constraints = trueCon;
Iter_used = iter;
Summary = table(ObjVal_opt,Variables_opt,Constraints,Iter_used)











%% Convex R problem usage
function [ObjOpt, zOpt] = CVX_R_problem(NumOfVar, Obj, Constraint, z_l, z_u)
cvx_begin
variable x(1, NumOfVar);
ObjFunction = CVX_R_problem_ForOneLine(Obj, z_l, z_u, x);
%compute ConVal given x
NumofConstraint = size(Constraint,2);
expression ConVal(NumofConstraint);
for i = 1: NumofConstraint
    ConVal(i) = CVX_R_problem_ForOneLine(Constraint(i), z_l, z_u, x);
end
minimize ObjFunction;
subject to
x<=z_u;
x>=z_l;
for i = 1: NumofConstraint
    ConVal(i) <= 0;
end
cvx_end
ObjOpt = cvx_optval;
zOpt = x;

function [ObjVal] = CVX_R_problem_ForOneLine(Obj, z_l, z_u, z)
% (R) problem in page 4 and 16
% for cvx optimization
% compute obj first
ObjVal = 0;
NumOfPositiveTermsinObj = size(Obj.positive, 1);
Y_0_L_pos = zeros(1,NumOfPositiveTermsinObj);
Y_0_U_pos = zeros(1,NumOfPositiveTermsinObj);
for row = 1:NumOfPositiveTermsinObj % iterate for all positive terms in Obj
    if (Obj.positive(row,1) > 0) % if not zero
        temp = Obj.positive(row,2:end).* z; % sum(alpha*z_i)
        ObjVal = ObjVal + Obj.positive(row,1) * exp(sum(temp));
    end
end
NumOfNegativeTermsinObj = size(Obj.negative, 1);
Y_0_L_neg = zeros(1,NumOfNegativeTermsinObj);
Y_0_U_neg = zeros(1,NumOfNegativeTermsinObj);
for row = 1:NumOfNegativeTermsinObj % iterate for all positive terms in Obj
    if (Obj.negative(row,1) > 0) % if not zero
        temp_lower = Obj.negative(row,2:end).* z_l; % used for find min and max for Y
        temp_upper = Obj.negative(row,2:end).* z_u;
        if (all(Obj.negative(row,2:end) == 0))
            temp_z = 0;
            A_0 = 1; % when all alpha = 0, we take limit, can have A=B=1
            B_0 = 1;
        else
            temp_z =  Obj.negative(row,2:end).* z;
            temp_min = min(temp_lower,temp_upper);
            temp_max = max(temp_lower,temp_upper);
            Y_0_L_neg(row) = sum(temp_min);
            Y_0_U_neg(row) = sum(temp_max);
            A_0 = (Y_0_U_neg(row) * exp(Y_0_L_neg(row)) - Y_0_L_neg(row)...
                * exp(Y_0_U_neg(row))) / (Y_0_U_neg(row) - Y_0_L_neg(row));
            B_0 = (exp(Y_0_U_neg(row)) - exp(Y_0_L_neg(row)))...
                / (Y_0_U_neg(row) - Y_0_L_neg(row));
        end
        % TYPO in page 16, for G_j^{conv,L}, there should not exist a
        % summation operator before Y_{jk}^L
        ObjVal = ObjVal - Obj.negative(row,1) * (A_0+B_0* sum(temp_z));
    end
end







%% R problem functions
function [ObjVal,ConVal] = Rproblem(Obj, Constraint, z_l, z_u)
% (R) problem in page 4 and 16
ObjVal = RproblemForOneLine(Obj, z_l, z_u);
%compute ConVal given z
NumofConstraint = size(Constraint,2);
ConVal = zeros(1, NumofConstraint);
for i = 1: NumofConstraint
    ConVal(i) = RproblemForOneLine(Constraint(i), z_l, z_u);
end


function [ObjVal] = RproblemForOneLine(Obj, z_l, z_u)
% (R) problem in page 4 and 16
% compute obj first
ObjVal = 0;
NumOfPositiveTermsinObj = size(Obj.positive, 1);
Y_0_L_pos = zeros(1,NumOfPositiveTermsinObj);
Y_0_U_pos = zeros(1,NumOfPositiveTermsinObj);
for row = 1:NumOfPositiveTermsinObj % iterate for all positive terms in Obj
    if (Obj.positive(row,1) > 0) % if not zero
        temp_lower = Obj.positive(row,2:end).* z_l; % used for find min and max for Y
        temp_upper = Obj.positive(row,2:end).* z_u;
        temp_min = min(temp_lower,temp_upper);
        Y_0_L_pos(row) = sum(temp_min);
        ObjVal = ObjVal + Obj.positive(row,1) * exp(Y_0_L_pos(row));
    end
end
NumOfNegativeTermsinObj = size(Obj.negative, 1);
Y_0_L_neg = zeros(1,NumOfNegativeTermsinObj);
Y_0_U_neg = zeros(1,NumOfNegativeTermsinObj);
for row = 1:NumOfNegativeTermsinObj % iterate for all positive terms in Obj
    if (Obj.negative(row,1) > 0) % if not zero
        temp_lower = Obj.negative(row,2:end).* z_l; % used for find min and max for Y
        temp_upper = Obj.negative(row,2:end).* z_u;
        temp_min = min(temp_lower,temp_upper);
        temp_max = max(temp_lower,temp_upper);
        Y_0_L_neg(row) = sum(temp_min);
        Y_0_U_neg(row) = sum(temp_max);
        if (all(Obj.negative(row,2:end) == 0))
            A_0 = 1; % when all alpha = 0, we take limit, can have A=B=1
            B_0 = 1;
        else
            A_0 = (Y_0_U_neg(row) * exp(Y_0_L_neg(row)) - Y_0_L_neg(row)...
                * exp(Y_0_U_neg(row))) / (Y_0_U_neg(row) - Y_0_L_neg(row));
            B_0 = (exp(Y_0_U_neg(row)) - exp(Y_0_L_neg(row)))...
                / (Y_0_U_neg(row) - Y_0_L_neg(row));
        end
        % TYPO in page 16, for G_j^{conv,L}, there should not exist a
        % summation operator before Y_{jk}^L
        ObjVal = ObjVal - Obj.negative(row,1) * (A_0+B_0*Y_0_U_neg(row));
    end
end





%% DC problem
function [ObjVal, ConVal] = DCproblem(Obj, Constraint, z)
% (DC) problem in page 3
% compute objVal given z
ObjVal = 0;
for row = 1:size(Obj.positive, 1) % iterate for all positive terms in Obj
    if (Obj.positive(row,1) > 0) % if not zero
        ObjVal = ObjVal + Obj.positive(row,1) * exp(dot(Obj.positive(row,2:end) , z) );
    end
end
for row = 1:size(Obj.negative, 1)
    if (Obj.negative(row,1) > 0) % if not zero
        ObjVal = ObjVal - Obj.negative(row,1) * exp(dot(Obj.negative(row,2:end) , z) );
    end
end

%compute ConVal given z
NumofConstraint = size(Constraint,2);
ConVal = zeros(1,NumofConstraint);
for i = 1: NumofConstraint
    for row = 1:size(Constraint(i).positive, 1) % iterate for all positive terms in Obj
        if (Constraint(i).positive(row,1) > 0) % if not zero
            ConVal(i)  = ConVal(i)  + Constraint(i).positive(row,1) * exp(dot(Constraint(i).positive(row,2:end), z));
        end
    end
    for row = 1:size(Constraint(i).negative, 1)
        if (Constraint(i).negative(row,1) > 0) % if not zero
            ConVal(i)  = ConVal(i)  - Constraint(i).negative(row,1) * exp(dot(Constraint(i).negative(row,2:end),z));
        end
    end
end

