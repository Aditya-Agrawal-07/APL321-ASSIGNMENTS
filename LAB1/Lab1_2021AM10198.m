% Question 1

n = 10^(-3);                        % step size
eta_0 = 0;                          % initial eta
eta_f = 10;                         % final_eta
N = (eta_f - eta_0) / n;            % Number of iterations

% Boundary Conditions
f_0 = 0;
g_0 = 0;
h0_0 = 0.1;                         % 1st guess for f''(0)

f0 = zeros(N+1, 1);                 
g0 = zeros(N+1, 1);
h0 = zeros(N+1, 1);

f0(1) = f_0;
g0(1) = g_0;
h0(1) = h0_0;

for i = 1:N                         % Euler's Method to solve system of ODEs
    f0(i+1) = f0(i) + n*g0(i);
    g0(i+1) = g0(i) + n*h0(i);
    h0(i+1) = h0(i) + n*(-0.5)*f0(i)*h0(i);
end

h1_0 = 0.2;                         % 2nd guess for f''(0)

f1 = zeros(N+1, 1);
g1 = zeros(N+1, 1);
h1 = zeros(N+1, 1);

f1(1) = f_0;
g1(1) = g_0;
h1(1) = h1_0;

for i = 1:N                         % Euler's Method to solve system of ODEs
    f1(i+1) = f1(i) + n*g1(i);
    g1(i+1) = g1(i) + n*h1(i);
    h1(i+1) = h1(i) + n*(-0.5)*f1(i)*h1(i);
end

% Newton-Raphson to optimize h_0

hk_1 = h0_0;
hk = h1_0;
gk_1 = g0(N+1);
gk = g1(N+1);
tol = 1e-5;                         % tolerance value for convergence
err = abs(gk - 1);                  % error after each iteration

while (err > tol)

    h_next = hk - ((gk - 1)*(hk - hk_1))/(gk - gk_1);
    
    f = zeros(N+1, 1);
    g = zeros(N+1, 1);
    h = zeros(N+1, 1);

    f(1) = f_0;
    g(1) = g_0;
    h(1) = h_next;

    for i = 1:N                     % Euler's Method to solve system of ODEs
        f(i+1) = f(i) + n*g(i);
        g(i+1) = g(i) + n*h(i);
        h(i+1) = h(i) + n*(-0.5)*f(i)*h(i);
    end

    hk_1 = hk;
    hk = h_next;
    gk_1 = gk;
    gk = g(N+1);
    err = abs(gk - 1);

end

h_0 = hk;                           % Optimized guess value of h_0                           

f = zeros(N+1, 1);
g = zeros(N+1, 1);
h = zeros(N+1, 1);
eta = zeros(N+1, 1);
eta(1) = 0;

f(1) = f_0;
g(1) = g_0;
h(1) = h_0;

for i = 1:N                         % Euler's Method to solve system of ODEs
    eta(i+1) = i*n;
    f(i+1) = f(i) + n*g(i);
    g(i+1) = g(i) + n*h(i);
    h(i+1) = h(i) + n*(-0.5)*f(i)*h(i);
end

solution = zeros(4, 1);              % values of f'=u/U for η = 1,2,3,4

for i = 1:4
    solution(i) = g(i/n + 1);
end

% Plot of u/U_∞ Vs η

figure('Position',[10 10 1600 600]);
subplot(1,2,1);
plot(eta, g);
title('Plot of u/U_∞ Vs η')
ylabel("u/U_∞");
xlabel("η");
set(gca,'FontSize',16);

Expected = [0.32979; 0.62977; 0.84605; 0.95552];
eta_ = [1; 2; 3; 4];

% Expected Vs Computed plots of u/U_∞ Vs η
subplot(1,2,2);
plot(eta_, Expected)                               
hold on
plot(eta_, solution)
legend('Expected', 'Computed');
customTicks = [1, 2, 3, 4];
xticks(customTicks);
ylabel("u/U_∞");
xlabel("η");
title('Expected Vs Computed Solution')
set(gca,'FontSize',16);



% Question 2

U = 10; nu = 1e-5;                  % constants U_∞ and nu

lx = 1;                             % lower limit of x coordinate
Lx = 10;                            % upper limit of x coordinate
Ly = 1e-2;                          % upper limit of y coordinate

% number of points along x and y
nx = 19; ny = 19;

% distance between two consecutive abscissa and ordinate resp.
dx = (Lx-lx)/nx; dy = Ly/ny;

% coordinates along x and y directions
xcoor = (0:nx)*dx + lx;
ycoor = (0:ny)*dy;

% coordinates in a ny-by-nx mesh 
[xm,ym] = meshgrid(xcoor,ycoor);

% eta values in the meshgrid
eta_mesh = fxy(xm,ym,U,nu);

% finding values of u/U_∞ on the mesh using interpolation
fmesh = interp1(eta,g,eta_mesh);

% new figure 
figure('Position',[10 10 1600 600]);

% plot isocontour of 2D function
contourf(xm,ym,fmesh);
xlabel('x');
ylabel('y');
colorbar;
title('Isocontour of u/U_∞');
set(gca,'FontSize',16);
xlim([1.00 10.00])
ylim([0.0000 0.0100])

function f = fxy(x,y,U,nu)              % function to evaluate eta at (x,y)
   f = y .* sqrt(U ./ (nu .* x));
end