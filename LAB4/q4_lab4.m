close all;
clear all;
clc;

% Defining the constants of the equation

rho = 1;                    % density = 1
Gamma = 0.01;               % Gamma = 1
N = [20; 80; 320];          % No. of grids
beta = [1.1; 1.024; 1.006]; % Stretching ratios
ITER = zeros(3,1); iter = zeros(3,1);       % No. of iterations
phi = cell(3,1); PHI = cell(3,1);           % Soln matrices
xcenter = cell(3,1); ycenter = cell(3,1);   % center of CVs
w = [1.51; 1.815; 1.9299]; W = [0.6; 1.657; 1.913]; % relaxation parameter


for k=1:3
    
    n = N(k);           % no. of grids
    r = beta(k);        % Stretching ratio
    nx = n+1; ny = n+1; % no. of coordinate points

    ycoor = linspace(0,1,ny)';  % y coordinates
    xcoor = f(nx,r);            % x coordinates

    xcen = (xcoor(1:nx-1,:) + xcoor(2:nx,:))/2; % x coordinate of CVs 
    ycen = (ycoor(1:ny-1,:) + ycoor(2:ny,:))/2; % y coordinate of CVs

    xcenter{k,1} = xcen;
    ycenter{k,1} = ycen;

    xe = xcoor(2:nx,:);     % x coordinate of mid-point of CV's east face
    xw = xcoor(1:nx-1,:);   % x coordinate of mid-point of CV's west face
    yn = ycoor(2:ny,:);     % y coordinate of mid-point of CV's north face
    ys = ycoor(1:ny-1,:);   % y coordinate of mid-point of CV's south face

    % Matrix for CDS
    Ap=zeros(ny-1,nx-1);
    An=Ap; As=Ap; Aw=Ap; Ae=Ap;

    % Matrix for UDS
    ap=zeros(ny-1,nx-1);
    an=ap; as=ap; aw=ap; ae=ap;

    % Matrix for Source Term Q
    Q = zeros(ny-1,nx-1);

    % For interior nodes
    for i=2:nx-2
        for j=2:ny-2

            dy = ycoor(j+1)-ycoor(j); dx = xcoor(i+1)-xcoor(i);
            
            me = rho*xe(i)*dy;
            mw = -rho*xw(i)*dy;
            mn = -rho*yn(j)*dx;
            ms = rho*ys(j)*dx;

            le = (xe(i)-xcen(i))/(xcen(i+1)-xcen(i));
            lw = (xcen(i)-xw(i))/(xcen(i)-xcen(i-1));
            ln = (yn(j)-ycen(j))/(ycen(j+1)-ycen(j));
            ls = (ycen(j)-ys(j))/(ycen(j)-ycen(j-1));

            Ae(j,i) = me*le - (Gamma*dy)/(xcen(i+1)-xcen(i));
            Aw(j,i) = mw*lw - (Gamma*dy)/(xcen(i)-xcen(i-1));
            An(j,i) = mn*ln - (Gamma*dx)/(ycen(j+1)-ycen(j));
            As(j,i) = ms*ls - (Gamma*dx)/(ycen(j)-ycen(j-1));
            Ap(j,i) = -(Ae(j,i) + Aw(j,i) + As(j,i) + An(j,i));

            ae(j,i) = min(me,0) - (Gamma*dy)/(xcen(i+1)-xcen(i));
            aw(j,i) = min(mw,0) - (Gamma*dy)/(xcen(i)-xcen(i-1));
            an(j,i) = min(mn,0) - (Gamma*dx)/(ycen(j+1)-ycen(j));
            as(j,i) = min(ms,0) - (Gamma*dx)/(ycen(j)-ycen(j-1));
            ap(j,i) = -(ae(j,i) + aw(j,i) + as(j,i) + an(j,i));
        end
    end
    
    % For boundary nodes
    % West Boundary nodes
    for j=1:ny-1
        i = 1;
        dy = ycoor(j+1)-ycoor(j); dx = xcoor(i+1)-xcoor(i);
            
        me = rho*xe(i)*dy;
        mn = -rho*yn(j)*dx;
        ms = rho*ys(j)*dx;

        le = (xe(i)-xcen(i))/(xcen(i+1)-xcen(i));
        Ae(j,i) = me*le - (Gamma*dy)/(xcen(i+1)-xcen(i));
        ae(j,i) = min(me,0) - (Gamma*dy)/(xcen(i+1)-xcen(i));

        if j~=ny-1
        ln = (yn(j)-ycen(j))/(ycen(j+1)-ycen(j));
        An(j,i) = mn*ln - (Gamma*dx)/(ycen(j+1)-ycen(j));
        an(j,i) = min(mn,0) - (Gamma*dx)/(ycen(j+1)-ycen(j));
        end

        if j~=1
        ls = (ycen(j)-ys(j))/(ycen(j)-ycen(j-1));
        As(j,i) = ms*ls - (Gamma*dx)/(ycen(j)-ycen(j-1));
        as(j,i) = min(ms,0) - (Gamma*dx)/(ycen(j)-ycen(j-1));
        end

        Ap(j,i) = -(Ae(j,i) + As(j,i) + An(j,i)) + (Gamma*dy)/(xcen(i));
        ap(j,i) = -(ae(j,i) + as(j,i) + an(j,i)) + (Gamma*dy)/(xcen(i));

        Q(j,i) = (Gamma*dy*(1-ycen(j)))/xcen(i);

    end

    % North Boundary Nodes
    for i=1:nx-1
        j = ny-1;
        dy = ycoor(j+1)-ycoor(j); dx = xcoor(i+1)-xcoor(i);

        me = rho*xe(i)*dy;
        mw = -rho*xw(i)*dy;
        ms = rho*ys(j)*dx;

        if i~=nx-1
        le = (xe(i)-xcen(i))/(xcen(i+1)-xcen(i));
        Ae(j,i) = me*le - (Gamma*dy)/(xcen(i+1)-xcen(i));
        ae(j,i) = min(me,0) - (Gamma*dy)/(xcen(i+1)-xcen(i));
        end

        if i~=1
        lw = (xcen(i)-xw(i))/(xcen(i)-xcen(i-1));
        Aw(j,i) = mw*lw - (Gamma*dy)/(xcen(i)-xcen(i-1));
        aw(j,i) = min(mw,0) - (Gamma*dy)/(xcen(i)-xcen(i-1));
        end

        ls = (ycen(j)-ys(j))/(ycen(j)-ycen(j-1));
        As(j,i) = ms*ls - (Gamma*dx)/(ycen(j)-ycen(j-1));
        as(j,i) = min(ms,0) - (Gamma*dx)/(ycen(j)-ycen(j-1));

        Ap(j,i) = -(Ae(j,i) + Aw(j,i) + As(j,i)) + (Gamma*dx)/(yn(j)-ycen(j));
        ap(j,i) = -(ae(j,i) + aw(j,i) + as(j,i)) + (Gamma*dx)/(yn(j)-ycen(j));
    end
    dy = ycoor(ny) - ycoor(ny-1);
    Ap(ny-1,1) = Ap(ny-1,1) + (Gamma*dy)/(xcen(1));
    ap(ny-1,1) = ap(ny-1,1) + (Gamma*dy)/(xcen(1));

    % East Boundary Nodes
    for j=1:ny-1
        i = nx-1;
        dy = ycoor(j+1)-ycoor(j); dx = xcoor(i+1)-xcoor(i);
 
        mw = -rho*xw(i)*dy;
        mn = -rho*yn(j)*dx;
        ms = rho*ys(j)*dx;

        lw = (xcen(i)-xw(i))/(xcen(i)-xcen(i-1));
        Aw(j,i) = mw*lw - (Gamma*dy)/(xcen(i)-xcen(i-1));
        aw(j,i) = min(mw,0) - (Gamma*dy)/(xcen(i)-xcen(i-1));

        if j~=ny-1
        ln = (yn(j)-ycen(j))/(ycen(j+1)-ycen(j));
        An(j,i) = mn*ln - (Gamma*dx)/(ycen(j+1)-ycen(j));
        an(j,i) = min(mn,0) - (Gamma*dx)/(ycen(j+1)-ycen(j));
        end

        if j~=1
        ls = (ycen(j)-ys(j))/(ycen(j)-ycen(j-1));
        As(j,i) = ms*ls - (Gamma*dx)/(ycen(j)-ycen(j-1));
        as(j,i) = min(ms,0) - (Gamma*dx)/(ycen(j)-ycen(j-1));
        end

        Ap(j,i) = -(Aw(j,i) + As(j,i) + An(j,i));
        ap(j,i) = -(aw(j,i) + as(j,i) + an(j,i));
    end
    dx = xcoor(nx)-xcoor(nx-1);
    Ap(ny-1,nx-1) = Ap(ny-1,nx-1) + (Gamma*dx)/(yn(ny-1)-ycen(ny-1));
    ap(ny-1,nx-1) = ap(ny-1,nx-1) + (Gamma*dx)/(yn(ny-1)-ycen(ny-1));
 
    % South Boundary Nodes
    for i=1:nx-1
        j = 1;
        dy = ycoor(j+1)-ycoor(j); dx = xcoor(i+1)-xcoor(i);

        me = rho*xe(i)*dy;
        mw = -rho*xw(i)*dy;
        mn = -rho*yn(j)*dx;

        if i~=nx-1
        le = (xe(i)-xcen(i))/(xcen(i+1)-xcen(i));
        Ae(j,i) = me*le - (Gamma*dy)/(xcen(i+1)-xcen(i));
        ae(j,i) = min(me,0) - (Gamma*dy)/(xcen(i+1)-xcen(i));
        end

        if i~=1
        lw = (xcen(i)-xw(i))/(xcen(i)-xcen(i-1));
        Aw(j,i) = mw*lw - (Gamma*dy)/(xcen(i)-xcen(i-1));
        aw(j,i) = min(mw,0) - (Gamma*dy)/(xcen(i)-xcen(i-1));
        end
        
        ln = (yn(j)-ycen(j))/(ycen(j+1)-ycen(j));
        An(j,i) = mn*ln - (Gamma*dx)/(ycen(j+1)-ycen(j));
        an(j,i) = min(mn,0) - (Gamma*dx)/(ycen(j+1)-ycen(j));

        Ap(j,i) = -(Ae(j,i) + Aw(j,i) + An(j,i));
        ap(j,i) = -(ae(j,i) + aw(j,i) + an(j,i));
    end
    dy = ycoor(2) - ycoor(1);
    Ap(1,1) = Ap(1,1) + (Gamma*dy)/(xcen(1));
    ap(1,1) = ap(1,1) + (Gamma*dy)/(xcen(1));


    phi1 = zeros(nx-1,ny-1);
    err = 100;
    iter1 = 0;
    w1 = w(k,1);

    % UDS solution using SOR
    while (err > 1e-10)
        phi10 = phi1;

        for i=2:nx-2
            for j=2:ny-2
                phi1(j,i) = (Q(j,i) - (ae(j,i)*phi1(j,i+1) + an(j,i)*phi1(j+1,i) + as(j,i)*phi1(j-1,i) + aw(j,i)*phi1(j,i-1)))/ap(j,i); 
                d = phi1(j,i) - phi10(j,i);
                phi1(j,i) = phi10(j,i) + w1*d;
            end
        end

        for i=2:nx-2
            j = 1;
            phi1(j,i) = (Q(j,i) - (ae(j,i)*phi1(j,i+1) + an(j,i)*phi1(j+1,i) + aw(j,i)*phi1(j,i-1)))/ap(j,i);
            d = phi1(j,i) - phi10(j,i);
            phi1(j,i) = phi10(j,i) + w1*d;
            j = ny-1;
            phi1(j,i) = (Q(j,i) - (ae(j,i)*phi1(j,i+1) + as(j,i)*phi1(j-1,i) + aw(j,i)*phi1(j,i-1)))/ap(j,i);
            d = phi1(j,i) - phi10(j,i);
            phi1(j,i) = phi10(j,i) + w1*d;
        end

        for j=2:ny-2
            i = 1;
            phi1(j,i) = (Q(j,i) - (ae(j,i)*phi1(j,i+1) + an(j,i)*phi1(j+1,i) + as(j,i)*phi1(j-1,i)))/ap(j,i);
            d = phi1(j,i) - phi10(j,i);
            phi1(j,i) = phi10(j,i) + w1*d;
            i = nx-1;
            phi1(j,i) = (Q(j,i) - (an(j,i)*phi1(j+1,i) + as(j,i)*phi1(j-1,i) + aw(j,i)*phi1(j,i-1)))/ap(j,i); 
            d = phi1(j,i) - phi10(j,i);
            phi1(j,i) = phi10(j,i) + w1*d;
        end

        phi1(1,1) = (Q(1,1) - ae(1,1)*phi1(1,2) - an(1,1)*phi1(2,1))/ap(1,1); i = 1; j = 1;
        d = phi1(j,i) - phi10(j,i);
        phi1(j,i) = phi10(j,i) + w1*d;
        phi1(1,nx-1) = (Q(1,nx-1) - aw(1,nx-1)*phi1(1,nx-2) - an(1,nx-1)*phi1(2,nx-1))/ap(1,nx-1); i = nx-1; j = 1;
        d = phi1(j,i) - phi10(j,i);
        phi1(j,i) = phi10(j,i) + w1*d;
        phi1(ny-1,1) = (Q(ny-1,1) - ae(ny-1,1)*phi1(ny-1,2) - as(ny-1,1)*phi1(ny-2,1))/ap(ny-1,1); i = 1; j = ny-1;
        d = phi1(j,i) - phi10(j,i);
        phi1(j,i) = phi10(j,i) + w1*d;
        phi1(ny-1,nx-1) = (Q(ny-1,nx-1) - aw(ny-1,nx-1)*phi1(ny-1,nx-2) - as(ny-1,nx-1)*phi1(ny-2,nx-1))/ap(ny-1,nx-1); i = nx-1; j = ny-1;
        d = phi1(j,i) - phi10(j,i);
        phi1(j,i) = phi10(j,i) + w1*d;

        err = norm(phi1 - phi10)/sqrt((nx-1)*(ny-1));
        iter1 = iter1 + 1;
    end   
    phi{k,1} = phi1;
    iter(k,1) = iter1


    % CDS solution using SOR
    phi1 = zeros(nx-1,ny-1);
    err = 100;
    iter1 = 0;
    w1 = W(k,1);

    while (err > 1e-10)
        phi10 = phi1;

        for i=2:nx-2
            for j=2:ny-2
                phi1(j,i) = (Q(j,i) - (Ae(j,i)*phi1(j,i+1) + An(j,i)*phi1(j+1,i) + As(j,i)*phi1(j-1,i) + Aw(j,i)*phi1(j,i-1)))/Ap(j,i); 
                d = phi1(j,i) - phi10(j,i);
                phi1(j,i) = phi10(j,i) + w1*d;
            end
        end

        for i=2:nx-2
            j = 1;
            phi1(j,i) = (Q(j,i) - (Ae(j,i)*phi1(j,i+1) + An(j,i)*phi1(j+1,i) + Aw(j,i)*phi1(j,i-1)))/Ap(j,i);
            d = phi1(j,i) - phi10(j,i);
            phi1(j,i) = phi10(j,i) + w1*d;
            j = ny-1;
            phi1(j,i) = (Q(j,i) - (Ae(j,i)*phi1(j,i+1) + As(j,i)*phi1(j-1,i) + Aw(j,i)*phi1(j,i-1)))/Ap(j,i);
            d = phi1(j,i) - phi10(j,i);
            phi1(j,i) = phi10(j,i) + w1*d;
        end

        for j=2:ny-2
            i = 1;
            phi1(j,i) = (Q(j,i) - (Ae(j,i)*phi1(j,i+1) + An(j,i)*phi1(j+1,i) + As(j,i)*phi1(j-1,i)))/Ap(j,i);
            d = phi1(j,i) - phi10(j,i);
            phi1(j,i) = phi10(j,i) + w1*d;
            i = nx-1;
            phi1(j,i) = (Q(j,i) - (An(j,i)*phi1(j+1,i) + As(j,i)*phi1(j-1,i) + Aw(j,i)*phi1(j,i-1)))/Ap(j,i); 
            d = phi1(j,i) - phi10(j,i);
            phi1(j,i) = phi10(j,i) + w1*d;
        end

        phi1(1,1) = (Q(1,1) - Ae(1,1)*phi1(1,2) - An(1,1)*phi1(2,1))/Ap(1,1); i = 1; j = 1;
        d = phi1(j,i) - phi10(j,i);
        phi1(j,i) = phi10(j,i) + w1*d;
        phi1(1,nx-1) = (Q(1,nx-1) - Aw(1,nx-1)*phi1(1,nx-2) - An(1,nx-1)*phi1(2,nx-1))/Ap(1,nx-1); i = nx-1; j = 1;
        d = phi1(j,i) - phi10(j,i);
        phi1(j,i) = phi10(j,i) + w1*d;
        phi1(ny-1,1) = (Q(ny-1,1) - Ae(ny-1,1)*phi1(ny-1,2) - As(ny-1,1)*phi1(ny-2,1))/Ap(ny-1,1); i = 1; j = ny-1;
        d = phi1(j,i) - phi10(j,i);
        phi1(j,i) = phi10(j,i) + w1*d;
        phi1(ny-1,nx-1) = (Q(ny-1,nx-1) - Aw(ny-1,nx-1)*phi1(ny-1,nx-2) - As(ny-1,nx-1)*phi1(ny-2,nx-1))/Ap(ny-1,nx-1); i = nx-1; j = ny-1;
        d = phi1(j,i) - phi10(j,i);
        phi1(j,i) = phi10(j,i) + w1*d;
        
        err = norm(phi1 - phi10)/sqrt((nx-1)*(ny-1));
        iter1 = iter1 + 1;
    end
    PHI{k,1} = phi1;
    ITER(k,1) = iter1

end

% Plotting the graphs for both methods at different grid sizes
figure;

for i=1:3
    subplot(2,3,i);
    contourf(xcenter{i,1},ycenter{i,1},phi{i,1});
    subplot(2,3,i+3);
    contourf(xcenter{i,1},ycenter{i,1},PHI{i,1});
end

subplot(2,3,1); ylabel({'\bf{UDS}', 'y'}); xlabel('x');
subplot(2,3,2); ylabel('y'); xlabel('x');
subplot(2,3,3); ylabel('y'); xlabel('x');
subplot(2,3,4); ylabel({'\bf{CDS}','y'}); xlabel({'x','\bf{N = 20}'});
subplot(2,3,5); ylabel('y'); xlabel({'x','\bf{N = 80}'});
subplot(2,3,6); ylabel('y'); xlabel({'x','\bf{N = 320}'});
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
colorbar('Position', [0.93 0.1 0.02 0.8]);
sgtitle('\bf{Isocontours for ɸ(x,y)}','Fontsize', 16);


phi1 = phi{1,1}; phi2 = phi{2,1}; phi3 = phi{3,1}; 
PHI1 = PHI{1,1}; PHI2 = PHI{2,1}; PHI3 = PHI{3,1};
x1 = xcenter{1,1}; x2 = xcenter{2,1}; x3 = xcenter{3,1};

% Saving the solution matrix for Question 5
save('solution.mat','phi1','phi2','phi3','PHI1','PHI2','PHI3','x1','x2','x3');


function coor = f(n,r)  % function to build coordinate axis with constant stretching ratio
    x = zeros(n,1);
    x(n,1) = 1;
    delx = (r - 1)/(r^(n-1) - 1);
    for i=2:n-1
        x(i) = x(i-1) + delx;
        delx = r*delx;
    end
    coor = x;
end
