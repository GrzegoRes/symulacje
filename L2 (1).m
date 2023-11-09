x_begin = [5, 0, 0, 0, 0, 0];
y_begin = [-5, 0, 0, 0, 0, 0];
z_begin = [200, 0, 0, 0, 0, 0];
vx_begin = [1, 0, 0, 0, 0, 0];
vy_begin = [3, 0, 0, 0, 0, 0];
vz_begin = [0, 0, 0, 0, 0, 0];
g=10;
k = 1;
m = 1;
t_vec = [0, 0, 0, 0, 0];
Ek = [0,0,0,0,0];
Ep = [0,0,0,0,0];
Et = [0,0,0,0,0];

% f(x,y) = x^2 + y^2 + 6*x + 6*y

%{
syms x0 y0 z0 vx0 vy0 vz0 g t;
eqn = z0 + vz0*t - (1/2)*g*t^2 == (x0 + vx0*t)^2 + (y0 + vy0*t)^2 + 6*(x0 + vx0*t) + 6*(y0 + vy0*t);
S = solve(eqn, t);
%}

for i = 1:5
    x0 = x_begin(i);
    y0 = y_begin(i);
    z0 = z_begin(i);
    vx0 = vx_begin(i);
    vy0 = vy_begin(i);
    vz0 = vz_begin(i);
    %{
    t1 =  (vz0 - 2*vx0*x0 - 2*vy0*y0 + (- 4*vx0^2*y0^2 + 4*z0*vx0^2 + 8*vx0*vy0*x0*y0 - 4*vx0*vz0*x0 - 4*vy0^2*x0^2 + 4*z0*vy0^2 - 4*vy0*vz0*y0 + vz0^2 - 2*g*x0^2 - 2*g*y0^2 + 2*g*z0)^(1/2))/(2*vx0^2 + 2*vy0^2 + g);
    t2 = -(2*vx0*x0 - vz0 + 2*vy0*y0 + (- 4*vx0^2*y0^2 + 4*z0*vx0^2 + 8*vx0*vy0*x0*y0 - 4*vx0*vz0*x0 - 4*vy0^2*x0^2 + 4*z0*vy0^2 - 4*vy0*vz0*y0 + vz0^2 - 2*g*x0^2 - 2*g*y0^2 + 2*g*z0)^(1/2))/(2*vx0^2 + 2*vy0^2 + g);
    %}
    
    % rozwiązania analitycznie
%     t1 = -(6*vx0 + 6*vy0 - vz0 + (2*g*z0 - 4*vx0^2*y0^2 - 12*g*x0 - 12*g*y0 - 4*vy0^2*x0^2 + 72*vx0*vy0 - 12*vx0*vz0 - 12*vy0*vz0 - 2*g*x0^2 - 2*g*y0^2 - 24*vy0^2*x0 - 24*vx0^2*y0 + 4*vx0^2*z0 + 4*vy0^2*z0 + 36*vx0^2 + 36*vy0^2 + vz0^2 + 24*vx0*vy0*x0 - 4*vx0*vz0*x0 + 24*vx0*vy0*y0 - 4*vy0*vz0*y0 + 8*vx0*vy0*x0*y0)^(1/2) + 2*vx0*x0 + 2*vy0*y0)/(2*vx0^2 + 2*vy0^2 + g)
%     t2 = -(6*vx0 + 6*vy0 - vz0 - (2*g*z0 - 4*vx0^2*y0^2 - 12*g*x0 - 12*g*y0 - 4*vy0^2*x0^2 + 72*vx0*vy0 - 12*vx0*vz0 - 12*vy0*vz0 - 2*g*x0^2 - 2*g*y0^2 - 24*vy0^2*x0 - 24*vx0^2*y0 + 4*vx0^2*z0 + 4*vy0^2*z0 + 36*vx0^2 + 36*vy0^2 + vz0^2 + 24*vx0*vy0*x0 - 4*vx0*vz0*x0 + 24*vx0*vy0*y0 - 4*vy0*vz0*y0 + 8*vx0*vy0*x0*y0)^(1/2) + 2*vx0*x0 + 2*vy0*y0)/(2*vx0^2 + 2*vy0^2 + g)
    
    
    % rozwiązanie własne
    % f: (vx0^2+vy0^2+g/2)*t^2 + (2*x0*vx0+2*y0*vy0+6*vx0+6*vy0-vz0)*t + (x0^2+y0^2+6*x0+6*y0-z0)
    t0 = 0;
    iterator = 0.1;
    while t0 <= 1e-10
        fun = @(t)(vx0^2+vy0^2+g/2)*t^2 + (2*x0*vx0+2*y0*vy0+6*vx0+6*vy0-vz0)*t + (x0^2+y0^2+6*x0+6*y0-z0); % function
        t0 = fzero(fun, iterator); % initial point
        iterator = iterator + 0.1
    end
    
%     t= [];
%     if isreal(t1)
%         t(length(t)+1) = t1;
%     end
%     if isreal(t2)
%         t(length(t)+1) = t2;
%     end
%     
%     if isempty(t)
%         exit
%     end
% 
%     if i == 1
%         t = min(t(t>=0));
%         
%     else
%         t = max(t);
%     end
%     t_vec(i) = t;


    t_vec(i) = t0;

    x_begin(i + 1) = x_begin(i) + vx_begin(i) * t_vec(i);
    y_begin(i + 1) = y_begin(i) + vy_begin(i) * t_vec(i);
    z_begin(i + 1) = z_begin(i) + vz_begin(i) * t_vec(i) - (1/2)*g*t_vec(i)^2;

    %N = [-2*x_begin(i + 1), -2*y_begin(i + 1), 1];
    x_1 = x_begin(i + 1)-0.01;
    x_2 = x_begin(i + 1)+0.01;
    y_1 = y_begin(i + 1)-0.01;
    y_2 = y_begin(i + 1)+0.01;
    f_prim_po_x = ((x_2^2 + y_begin(i + 1)^2 + 6*x_2 + 6*y_begin(i + 1)) - (x_1^2 + y_begin(i + 1)^2 + 6*x_1 + 6*y_begin(i + 1)))/(x_2-x_1);
    f_prim_po_y = ((x_begin(i + 1)^2 + y_2^2 + 6*x_begin(i + 1) + 6*y_2) - (x_begin(i + 1)^2 + y_1^2 + 6*x_begin(i + 1) + 6*y_1))/(y_2-y_1);
    N = [ -f_prim_po_x, -f_prim_po_y, 1];
    n = N/norm(N);
    v_before = [vx_begin(i), vy_begin(i), vz_begin(i) - g*t_vec(i)];
    v_before_length = sqrt(v_before(1)^2 + v_before(2)^2 + v_before(3)^2);
    v_after = (v_before - 2 * dot(v_before, n) * n)*sqrt(k);
    vx_begin(i + 1) = v_after(1);
    vy_begin(i + 1) = v_after(2);
    vz_begin(i + 1) = v_after(3);
    
    Ek(i) = m * v_before_length^2 / 2;
    Ep(i) = m * g * z_begin(i + 1);
    Et(i) = Ek(i) + Ep(i);
end

x_vec = [];
y_vec = [];
z_vec = [];
for i = 1:5
   for t_iter = 0:0.01:t_vec(i)
       x_vec(length(x_vec)+1) = x_begin(i)+vx_begin(i)*t_iter;
       y_vec(length(y_vec)+1) = y_begin(i)+vy_begin(i)*t_iter;
       z_vec(length(z_vec)+1) = z_begin(i)+vz_begin(i)*t_iter-g*t_iter^2/2;
   end
end

x_plot = -15:0.5:15;
y_plot = x_plot;
[X,Y] = meshgrid(x_plot);
A = X.^2 + Y.^2;
B = 6*X;
C = 6*Y;
F = A + B + C;
mesh(X,Y,F);
hold on
plot3(x_vec, y_vec, z_vec, 'r');
hold off