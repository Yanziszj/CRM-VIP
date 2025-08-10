function [ ] = PerformanceProfile_error(T,type_test)
%PERFORMANCEPROFILE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��


% ģ�����ݣ�ÿ����һ�����⣬ÿ����һ���㷨��CPUʱ�䣨�룩
% T = [
%     1.0, 1.2, 3.5;
%     0.8, 0.9, 2.0;
%     2.0, 2.0, 2.1;
%     1.5, 1.4, 3.0;
%     0.6, 1.1, 1.2;
%     3.0, 2.0, 2.5;
%     1.8, 1.7, 2.0;
%     0.9, 1.0, 1.8;
%     2.2, 2.2, 1.0;
%     1.3, 1.2, 1.1
% ];

[numProblems, numSolvers] = size(T);

% �������ܱ� r_p,s
R = zeros(size(T));
for p = 1:numProblems
    min_time = min(T(p,:));
    R(p,:) = T(p,:) / min_time;
end

% �趨 tau �ķ�Χ
tau = linspace(1, 10, 100);  % ����� 1 �� 5 ��

% ����ÿ���㷨�����ܸſ� rho_s(tau)
rho = zeros(numSolvers, length(tau));
for s = 1:numSolvers
    for i = 1:length(tau)
        rho(s,i) = sum(R(:,s) <= tau(i)) / numProblems;
    end
end

% ��ͼ
figure; hold on;
colors = {'r','b','g','k','m','c'};
methods = {'Alg1','Alg2','Extra'};
linestyles = {'-','--',':','-.'};
switch type_test
    case 1
        for s = 1:numSolvers
        plot(tau, rho(s,:),linestyles{s}, 'Color', colors{s}, 'LineWidth', 2, 'DisplayName', ['Alg ' num2str(s)]);
        end
    case 2        
        for s = 1:numSolvers
            plot(tau, rho(s,:), linestyles{s}, 'Color', colors{s}, 'LineWidth', 2, 'DisplayName', methods{s});
        end
end
xlabel('\tau','FontSize',12);
ylabel('Performance Profile \rho_s(\tau)','FontSize',12);
title('Performance Profile (Error)','FontSize',14);
legend('show','Location', 'southeast');
grid on;

end

