%% This file is to test three algorithms for the variation inequlaity problems 
%  Algorithms: 1. Algorithm 1 in our paper (for paramonotone)
%                      2. Algorithm 2 in our paper (for monotone)
%                      3. Extragradient method
%                      4. Algoirthm 2 in paper by Yunier and Alfredo (2010)
%

clc;clear;
%% Initialianzation
n =  5;                 % the dimension
m = 5;                 % the number of ellipsoids      

type_problem =2;                % the type of problem (mainly means different monotone operators)
                                             % 1: paramonotone with sysmetric matrix
                                             % 2: paramonotone with asysmentric matrix 
                                             % 3: monotone 

Max_experient = 1;              % Max times to do the experients
e =sqrt(eps);                         % stopping tolerance

Max_iteration1 = 50000;          % maximum of the iteration
Max_iteration2 = 50000;
Max_iteration3 = 50000; 
Max_experient = 1;               % Max times to do the experients
Max_run = 20;                      % for each experiment the times we ran

Times = zeros(Max_experient,3);
Errors = zeros(Max_experient,3);
Iterations = zeros(Max_experient,3);

for i=1:Max_experient

    [C, Slater_point]= GeneratingEllipsoid1(n,m);    % Generate m ellipsoid for the feasibility part of the VIP
    x0 = Generatingx0Ellipsoid(C,n,m)*1.1;   
    M = GeneratingMatrix(n,type_problem);  % Generate the matrix for the operator F in the VIP

    [xk_ACRM,time_ACRM,error_ACRM,iteration_ACRM,index_stopping_ACRM] =...
        ACRM_VIP( C,M,e,Max_iteration1,x0,type_problem);
    
    [xk_ECM,time_ECM,error_ECM,iteration_ECM,G_yk] =...
        ECM_VIP( C,M,e,Max_iteration2,x0,Slater_point,type_problem);
    
    [xk_Extra,time_Extra,error_Extra,iteration_Extra] =...
        Extragradient( C,M,e,Max_iteration3,x0,type_problem);
    
    [xk_Yunier2010_2,time_Yunier2010_2,error_Yunier2010_2,iteration_Yunier2010_2,index_stopping_Yunier2010_2] =...
        Yunier2010_2( C,M,e,Max_iteration1,x0,type_problem);

    Times(i,1) = time_ACRM;
    Times(i,2) = time_ECM;
    Times(i,3) = time_Extra;
    Times(i,4) = time_Yunier2010_2;
    
    Errors(i,1) = error_ACRM;
    Errors(i,2) = error_ECM;
    Errors(i,3) = error_Extra;
    Errors(i,4) = error_Yunier2010_2;
    
    Iterations(i,1) = iteration_ACRM;
    Iterations(i,2) = iteration_ECM;
    Iterations(i,3) = iteration_Extra;
    Iterations(i,4) = iteration_Yunier2010_2;
    
end
% Times = log(Times);
% PerformanceProfile(Times);

% x=[1:1:Max_experient];
% 
% figure(1);
% plot(x,log(Times(:,1)),'r--');
% hold on;
% plot(x,log(Times(:,2)),'g--');
% hold on;
% plot(x,log(Times(:,3)),'k--');
% hold on;
% plot(x,log(Times(:,4)),'c-^');
% hold on;
% 
% legend('Algorithm1','Algorithm2','Extragradient','Yunier_2010');
% xlabel('experiment number');
% ylabel('log(CPU time)');
% title('Performance of times');
% 
% figure(2);
% plot(x,log(Errors(:,1)),'r--');
% hold on;
% plot(x,log(Errors(:,2)),'g--');
% hold on;
% plot(x,log(Errors(:,3)),'k--');
% hold on;
% plot(x,log(Errors(:,4)),'c-^');
% hold on;
% 
% 
% legend('Algorithm1','Algorithm2','Extragradient','Yunier2010');
% xlabel('experiment number');
% ylabel('log(Errors)');
% title('Performance of Errors');



