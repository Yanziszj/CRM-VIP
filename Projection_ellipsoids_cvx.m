function [ x ] = Projection_ellipsoids_cvx(C, n, m, x0)
%PROJECTION_ELLIPSOIDS_CVX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

cvx_begin quiet
    variable x(n)
    minimize( sum_square(x - x0) )
    subject to
        for i = 1:m
            quad_form(x, C{1,i}) + 2 * C{2,i}' * x - C{3,i} <= 0;
        end
cvx_end


end

