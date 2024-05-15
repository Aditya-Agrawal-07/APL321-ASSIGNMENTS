close all;
clear all;
clc;

Re = 0.01;
u0 = 1; y0 = 0.5; w = 0.3;
Lx = 4;  Ly = 1;        % dimensions of rectangle
nx = 101; ny = 26;      % no. of grid points
x = linspace(0,Lx,nx);  % x coordinates of nodes
y = linspace(0,Ly,ny);  % y coordinates of nodes
X = (x(:,1:nx-1)+x(:,2:nx))/2;
Y = (y(:,1:ny-1)+y(:,2:ny))/2;


dx = Lx/(nx-1); dy = Ly/(ny-1);
del = dx;

u_star = zeros(ny+1,nx+2);
v_star = zeros(ny+2,nx+1);
p_star = zeros(ny+1,nx+1);

for j=2:ny
    u_star(j,2) = u0*exp(-(Y(j-1)-y0)^2/w^2);
end   

vec = u_star(:,2);
u_star = repmat(vec, 1, nx+2);
U = u_star; 
V = v_star;
Residual = 100
ITER = 0;

while ITER < 10
    
    % Prediction step 

    Ae = -0.25; Aw = -0.25; An = -0.25; As = -0.25; Ap = 1;

    % u_star prediction
    omega = 1.53;
    err = 100;
    iter = 0;
    while (abs(err) > 1e-8)
    
        u_prev = u_star;
        
        for i=3:nx+1
            for j=2:ny

                Qp = Re*del*0.25*(p_star(j,i-1)-p_star(j,i));
                u_star(j,i) = (Qp - Ae*u_star(j,i+1) - Aw*u_star(j,i-1) - As*u_star(j-1,i) - An*u_star(j+1,i))/Ap;
                d = u_star(j,i) - u_prev(j,i);
                u_star(j,i) = u_prev(j,i) + omega*d;
    
            end
        end

        u_star(:,2) = vec;  % west
        u_star(:,nx+2) = 2*u_star(:,nx+1) - u_star(:,nx); % east
        u_star(ny+1,:) = -u_star(ny,:);  % north
        u_star(1,:) = -u_star(2,:); % south
        u_star(2:ny,1) = 2*u_star(2:ny,2) - u_star(2:ny,3);

        sz = size(u_prev);

        err = norm(u_star-u_prev)/sqrt(sz(1)*sz(2));
        iter = iter + 1;

    end

    % v_star prediction
    omega = 1.2;
    err = 100;
    iter = 0;
    while (abs(err) > 1e-8)
    
        v_prev = v_star;

        for i=2:nx
            for j=3:ny
                
                Qp = Re*del*0.25*(p_star(j-1,i)-p_star(j,i));
                v_star(j,i) = (Qp - Ae*v_star(j,i+1) - Aw*v_star(j,i-1) - As*v_star(j-1,i) - An*v_star(j+1,i))/Ap;
                d = v_star(j,i) - v_prev(j,i);
                v_star(j,i) = v_prev(j,i) + omega*d;
    
            end
        end

        v_star(:,1) = -v_star(:,2); % west
        v_star(1,:) = 2*v_star(2,:) - v_star(3,:);
        v_star(2,:) = 0;    % south
        v_star(ny+1,:) = 0; % north
        v_star(:,nx+1) = 2*v_star(:,nx) - v_star(:,nx-1);  % east
        

        sz = size(v_prev);

        err = norm(v_star-v_prev)/sqrt(sz(1)*sz(2));
        iter = iter + 1;
    
    end



    % pressure correction step

    Ae = 1; Aw = 1; An = 1; As = 1; Ap = -4; Qd = Re*del*0.25;
    p = p_star;
    omega = 1.947;
    err = 100;
    iter = 0;
    while (abs(err) > 1e-8)
    
        p_prev = p;
        
        for i=2:nx
            for j=2:ny
    
                Q = (1/Qd)*(u_star(j,i+1) - u_star(j,i) + v_star(j+1,i) - v_star(j,i));
                p(j,i) = (Q - Ae*p(j,i+1) - Aw*p(j,i-1) - As*p(j-1,i) - An*p(j+1,i))/Ap;
                d = p(j,i) - p_prev(j,i);
                p(j,i) = p_prev(j,i) + omega*d;
    
            end
        end
    
        p(:,nx+1) = -p(:,nx);
        p(:,1) = p(:,2);
        p(1,:) = p(2,:);
        p(ny+1,:) = p(ny,:);

        sz = size(p_star);
    
        err = norm(p-p_prev)/sqrt(sz(1)*sz(2));
        iter = iter + 1;

    end


    % correct pressure and velocity

    u = zeros(ny+1,nx+2); v = zeros(ny+2,nx+1);

    for i=2:nx+1
        for j=2:ny
            u(j,i) = u_star(j,i) - del*Re*0.25*(p(j,i)-p(j,i-1));    
        end
    end

    for i=2:nx
        for j=2:ny
            v(j,i) = v_star(j,i) - del*Re*0.25*(p(j,i)-p(j-1,i));
        end
    end

    alpha_P = 0.001; alpha_u = 0.001; alpha_v = 0.001;

    U_prev = U; V_prev = V;
    P = p_star + alpha_P*p;
    U = alpha_u*u + (1-alpha_u)*U;
    V = alpha_v*v + (1-alpha_v)*V;

    Residual = norm(U-U_prev)/sqrt((ny-1)*(nx)) + norm(V-V_prev)/sqrt((nx-1)*ny)
    ITER = ITER + 1

    p_star = P;
    u_star = U;
    v_star = V;

end

% Plotting the relevant graphs

figure;
contourf(x,Y,U(2:ny,2:nx+1))
title("Isocontours for u(x,y)");
xlabel("x");
ylabel("y");
colorbar;

figure;
contourf(X,y,V(2:ny+1,2:nx));
title("Isocontours for v(x,y)");
xlabel("x");
ylabel("y");
colorbar;

figure;
contourf(X,Y,p(2:ny,2:nx));
title("Isocontours for p(x,y)");
xlabel("x");
ylabel("y");
colorbar;

u = (U(2:ny,2:nx) + U(2:ny,3:nx+1))/2;
v = (V(2:ny,2:nx) + V(3:ny+1,2:nx))/2;
figure;
quiver(X,Y,u,v);
xlabel('X-axis');
ylabel('Y-axis');
title('Vector Field Plot');
xlim([0,4]);

div_vel = zeros(ny-1,nx-1);

for i=2:nx
    for j=2:ny
        div_vel(j-1,i-1) = (U(j,i+1)-U(j,i)+V(j+1,i)-V(j,i))/del;
    end
end

% plotting the isocontours of solution
figure;
contourf(X,Y,div_vel);
title("Isocontours for \nabla.vel(x,y)");
xlabel("x");
ylabel("y");
colorbar;