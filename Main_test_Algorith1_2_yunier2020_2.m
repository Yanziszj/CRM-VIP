%% This file is to test three algorithms for the variation inequlaity problems 
%  Algorithms: 1. Algorithm 1 in our paper (for paramonotone)
%                      2. Algorithm 2 in our paper (for monotone)
%                      3. Extragradient method
%                  wmmm    4. Algoirthm 2 in paper by Yunier and Alfredo (2010)
%

clc;clear;
%% Initialianzation
N =  [5 10 20];                 % the dimension
Q = [2 5 10];                 % the number of ellipsoids      

type_problem = 2;                % the type of problem (mainly means different monotone operators)
                                             % 1: paramonotone with sysmetric matrix
                                             % 2: paramonotone with asysmentric matrix 
                                             % 3: monotone 
% e =sqrt(eps);                         % stopping tolerance
e = 10^(-6);

Max_iteration1 = 30000;          % maximum of the iteration
Max_iteration2 = 30000;
Max_iteration3 = 30000; 
Max_iteration4 = 30000; 
Max_experient = 20;               % Max times to do the experients
Max_run = 20;                      % for each experiment the times we ran

Time = [];
Iteration = [];
Error = [];

Times = zeros(Max_experient,4);
Errors = zeros(Max_experient,4);
Iterations = zeros(Max_experient,4);

Mid_time_ACRM = zeros(3,3);
Mid_time_ECM = zeros(3,3);
Mid_time_Yunier_2010_2 = zeros(3,3);
Mid_time_Yunier_2012 = zeros(3,3);

Mean_error_ACRM = zeros(3,3);
Mean_error_ECM = zeros(3,3);
Mean_error_Yunier_2010_2 = zeros(3,3);
Mean_error_Yunier_2012 = zeros(3,3);

Mid_iteration_ACRM =zeros(3,3);
Mid_iteration_ECM =zeros(3,3);
Mid_iteration_Yunier_2010_2 =zeros(3,3);
Mid_iteration_Yunier_2012 =zeros(3,3);


for k=1:3
    n=N(k);
    for j=1:3
        m = Q(j);
        for i=1:Max_experient
            seed = floor(1e8 * mod(now, 1));
            rng(seed);
            [C, Slater_point]= GeneratingEllipsoid1(n,m);    % Generate m ellipsoid for the feasibility part of the VIP
            x0 = Generatingx0Ellipsoid(C,n,m)*2;   
            M = GeneratingMatrix(n,type_problem);  % Generate the matrix for the operator F in the VIP

            [xk_ACRM,time_ACRM,error_ACRM,iteration_ACRM,index_stopping_ACRM] =...
                ACRM_VIP( C,M,e,Max_iteration1,x0,type_problem);

            [xk_ECM,time_ECM,error_ECM,iteration_ECM,G_yk] =...
                ECM_VIP( C,M,e,Max_iteration2,x0,Slater_point,type_problem);

        %     [xk_Extra,time_Extra,error_Extra,iteration_Extra] =...
        %         Extragradient( C,M,e,Max_iteration3,x0,type_problem);

            [xk_Yunier2010_2,time_Yunier2010_2,error_Yunier2010_2,iteration_Yunier2010_2,index_stopping_Yunier2010_2] =...
                Yunier2010_2( C,M,e,Max_iteration3,x0,type_problem);
            
            [xk_Yunier_2012,time_Yunier_2012,error_Yunier_2012,iteration_Yunier_2012,G_yk] =...
                Yunier_2012( C,M,e,Max_iteration4,x0,Slater_point,type_problem);

            Times(i,1) = time_ACRM;
            Times(i,2) = time_ECM;
        %     Times(i,3) = time_Extra;
            Times(i,3) = time_Yunier2010_2;
            Times(i,4) = time_Yunier_2012;

            Errors(i,1) = error_ACRM;
            Errors(i,2) = error_ECM;
        %     Errors(i,3) = error_Extra;
            Errors(i,3) = error_Yunier2010_2;
            Errors(i,4) = error_Yunier_2012;

            Iterations(i,1) = iteration_ACRM;
            Iterations(i,2) = iteration_ECM;
        %     Iterations(i,3) = iteration_Extra;
            Iterations(i,3) = iteration_Yunier2010_2;
            Iterations(i,4) = iteration_Yunier_2012;
            
            fprintf('Dimension:%d',n);
            fprintf('  number of ellipsoid:%d',m);
            fprintf('  Experiment:%d\n',i);
        end
        Mid_time_ACRM(k,j) = median(Times(:,1));
        Mid_time_ECM(k,j)  = median(Times(:,2));
        Mid_time_Yunier_2010_2(k,j)  =median(Times(:,3));
        Mid_time_Yunier_2012(k,j)  =median(Times(:,4));


        Mean_error_ACRM(k,j) = mean(Errors(:,1));
        Mean_error_ECM(k,j) = mean(Errors(:,2));
        Mean_error_Yunier_2010_2(k,j)= mean(Errors(:,3));
        Mean_error_Yunier_2012(k,j)= mean(Errors(:,4));
        
        Mid_iteration_ACRM(k,j) = median(Iterations(:,1));
        Mid_iteration_ECM(k,j)  = median(Iterations(:,2));
        Mid_iteration_Yunier_2010_2(k,j)  = median(Iterations(:,3));
        Mid_iteration_Yunier_2012(k,j)  = median(Iterations(:,4));
        
        Time = [Time; Times];
        Iteration = [Iteration;Iterations];
        Error = [Error;Errors];
        
    end
end


% Times = log(Times);
PerformanceProfile_time(Time,1);
PerformanceProfile_iteration(Iteration,1);
PerformanceProfile_error(Error,1);





