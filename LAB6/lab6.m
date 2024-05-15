close all;
clear all;
clc;

% Question 4

% Defining the constants of the equation

w = 1.9;
Lx = 4;  Ly = 1;        % dimensions of rectangle
nx = 101; ny = 26;      % no. of grid points
x = linspace(0,Lx,nx);  % x coordinates of nodes
y = linspace(0,Ly,ny);  % y coordinates of nodes
X = (x(:,1:nx-1)+x(:,2:nx))/2;
Y = (y(:,1:ny-1)+y(:,2:ny))/2;

dx = Lx/(nx-1); dy = Ly/(ny-1);
del = dx;

p = zeros(ny+1,nx+1);

Ae = 1; Aw = 1; An = 1; As = 1; Ap = -4;

err = 100;
iter = 0;
while (abs(err) > 1e-10)

    p_prev = p;
    
    for i=2:nx
        for j=2:ny

            Q = (0.1/del)*(sin(x(i)) - sin(x(i-1)) + sin(y(j)) - sin(y(j-1)));
            p(j,i) = (Q - Ae*p(j,i+1) - Aw*p(j,i-1) - As*p(j-1,i) - An*p(j+1,i))/Ap;
            d = p(j,i) - p_prev(j,i);
            p(j,i) = p_prev(j,i) + w*d;

        end
    end

    p(:,nx+1) = -p(:,nx);
    p(:,1) = p(:,2);
    p(1,:) = p(2,:);
    p(ny+1,:) = p(ny,:);

    err = norm(p-p_prev)/sqrt((nx+1)*(ny+1));
    iter = iter + 1;

end

% plotting the isocontours of solution
figure;
contourf(X,Y,p(2:ny,2:nx));
title("Isocontours for p(x,y)");
xlabel("x");
ylabel("y");
colorbar;



% Question 5

u_star = -x + 0.1*sin(x);
v_star = y + 0.1*sin(y);
u = zeros(ny,nx); v = zeros(ny,nx);

for i=1:nx
    for j=1:ny

        u(j,i) = u_star(i) - del*(p(j+1,i+1)-p(j+1,i));
        v(j,i) = v_star(j) - del*(p(j+1,i+1)-p(j,i+1));

    end
end

% plotting the isocontours of solution
figure;
contourf(x,y,u);
title("Isocontours for u(x,y)");
xlabel("x");
ylabel("y");
colorbar;

% plotting the isocontours of solution
figure;
contourf(x,y,v);
title("Isocontours for v(x,y)");
xlabel("x");
ylabel("y");
colorbar;


% Question 6

div_vel = zeros(ny-1,nx-1);

for i=1:nx-1
    for j=1:ny-1
        div_vel(j,i) = (u(j,i+1)-u(j,i)+v(j+1,i)-v(j,i))/del;
    end
end

% plotting the isocontours of solution
figure;
contourf(X,Y,div_vel);
title("Isocontours for \nabla.vel(x,y)");
xlabel("x");
ylabel("y");
colorbar;

div_vel_star = zeros(ny-1, nx-1);
for i=1:nx-1
    for j=1:ny-1
        div_vel_star(j,i) = (u_star(i+1)-u_star(i)+v_star(j+1)-v_star(j))/del;
    end
end

% plotting the isocontours of solution
figure;
contourf(X,Y,div_vel_star);
title("Isocontours for \nabla.vel^*(x,y)");
xlabel("x");
ylabel("y");
colorbar;
