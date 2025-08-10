function [ P ] = CalculatePho_time(Max_experient,r_m,time,l)
%CALCULATEPHO 
%   
% 
% l = 50; % length of the graph of x-axis
step = r_m/l;

% r_ps = zeros(1,Max_experient);
r_ps = times;
% for i=1:Max_experient
%     r_ps(i) = R{8,i};
% end

P = zeros(1,r_m);

for i=1:l
    s=0;
    for j=1:Max_experient
        if r_ps(j)<=1+(i-1)*step
            s=s+1;
        end
    end
    P(i) = s/Max_experient;
end

end
