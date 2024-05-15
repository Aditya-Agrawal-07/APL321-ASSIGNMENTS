close all;
clear all;
clc;

RE = [5 ; 50];
alpha_P = 0.001; alpha_u = 0.001; alpha_v = 0.001;
u0 = 1; y0 = 0.5; w = 0.1;
Lx = 1;  Ly = 1;                % dimensions of rectangle
nx = 50; ny = 50;               % no. of grid points
dx = Lx/nx; dy = Ly/ny;         % cell width
del = dx;                       % cell width

xp = dx/2.0+(0:nx-1)*dx;        % x coordinate of pressure nodes
yp = dy/2.0+(0:ny-1)*dy;        % y coordinate of pressure nodes

xu = linspace(0,Lx,nx+1);       % x coordinates of u nodes
yu = yp;                        % y coordinates of u nodes

xv = xp;                        % x coordinates of v nodes
yv = linspace(0,Ly,ny+1);       % y coordinates of v nodes

for idx=2:2
    
    Re = RE(idx);
    
    u = zeros(ny+2,nx+3);
    v = zeros(ny+3,nx+2);

    p = zeros(ny+2,nx+2);       % ghost cells for pressure

    uin = u0*exp(-(yu-y0).^2/w^2);  % initial condition of u_inlet
    uin = uin';
    for i=1:nx+3
        u(2:ny+1,i) = uin(:,:);
    end
    
    u_star = u; 
    v_star = v; 
    p_star = p;

    Niter = 400;
    Residual = zeros(Niter,1);
    ITER = 0;

    while ITER < Niter

        uold = u; vold = v; pold = p_star;

        aue = zeros(ny+2,nx+3); auw = aue; aun = auw; aus = aun; qu = aue;
        ave = zeros(ny+3,nx+2); avw = ave; avn = ave; avs = avn; qv = ave;
        avp = ave; aup = aue;

        % u_star prediction

        omega = 1.5;
        err = 100;
        iter = 0;

        for i=2:nx+2
                for j=2:ny+1

                    me = 0.5*(u_star(j,i)+u_star(j,i+1))*dy;
                    mw = -0.5*(u_star(j,i)+u_star(j,i-1))*dy;
                    mn = 0.5*(v_star(j+1,i)+v_star(j+1,i-1))*dx;
                    ms = -0.5*(v_star(j,i)+v_star(j,i-1))*dx;

                    Aue = 1/Re - min(0,me); aue(j,i) = Aue;
                    Auw = 1/Re - min(0,mw); auw(j,i) = Auw;
                    Aun = 1/Re - min(0,mn); aun(j,i) = Aun;
                    Aus = 1/Re - min(0,ms); aus(j,i) = Aus;
                    Aup = -4/Re - (max(0,me) + max(0,ms) + max(0,mn) + max(0,mw));  aup(j,i) = Aup;
                    Qup = (p_star(j,i)-p_star(j,i-1))*dy; qu(j,i) = Qup;

                end
        end

        aup(:,2) = 1; aue(:,2) = 0; aun(:,2) = 0; aus(:,2) = 0; auw(:,2) = 0;
        
%         aup(:,nx+2) = aup(:,nx+2) + 2*aue(:,nx+2);
%         auw(:,nx+2) = auw(:,nx+2) - aue(:,nx+2);
%         aue(:,nx+2) = 0;

        while abs(err) > 1e-10

            u_prev = u_star;

            for i=3:nx+2
                for j=2:ny+1

                    u_star(j,i) = (qu(j,i) - auw(j,i)*u_star(j,i-1) - aus(j,i)*u_star(j-1,i) - aue(j,i)*u_star(j,i+1) - aun(j,i)*u_star(j+1,i))/aup(j,i);

                    d = u_star(j,i) - u_prev(j,i);
                    u_star(j,i) = u_prev(j,i) + omega*d;

                end
            end

            u_star(2:ny+1,2) = uin(:,:);      % west bc
            u_star(2:ny+1,1) = 2*u_star(2:ny+1,2) - u_star(2:ny+1,3);
            u_star(1,:) = -u_star(2,:);
            u_star(ny+2,:) = -u_star(ny+1,:);
            u_star(:,nx+2) = 2*u_star(:,nx+1) - u_star(:,nx);

            sz = size(u_star);
            err = norm(u_star-u_prev)/sqrt(sz(1)*sz(2));
            iter = iter + 1;

        end

        % v_star prediction

        omega = 1.5;
        err = 100;
        iter = 0;

        for i=2:nx+1
            for j=2:ny+2

                me = 0.5*(u_star(j-1,i+1)+u_star(j,i+1))*dy;
                mw = -0.5*(u_star(j-1,i)+u_star(j,i))*dy;
                mn = 0.5*(v_star(j,i)+v_star(j+1,i))*dx;
                ms = -0.5*(v_star(j-1,i)+v_star(j,i))*dx;
        
                Ave = 1/Re - min(0,me); ave(j,i) = Ave;
                Avw = 1/Re - min(0,mw); avw(j,i) = Avw;
                Avn = 1/Re - min(0,mn); avn(j,i) = Avn;
                Avs = 1/Re - min(0,ms); avs(j,i) = Avs;
                Avp = -4/Re - (max(0,me) + max(0,ms) + max(0,mn) + max(0,mw)); avp(j,i) = Avp;
                Qvp = (p_star(j,i)-p_star(j-1,i))*dy;   qv(j,i) = Qvp;
            end
        end

        avp(2,:) = 1; ave(2,:) = 0; avw(2,:) = 0; avn(2,:) = 0; avs(2,:) = 0;

        
        while abs(err) > 1e-10

            v_prev = v_star;

            for i=2:nx+1
                for j=3:ny+1

                    v_star(j,i) = (qv(j,i) - avw(j,i)*v_star(j,i-1) - avs(j,i)*v_star(j-1,i) - ave(j,i)*v_star(j,i+1) - avn(j,i)*v_star(j+1,i))/avp(j,i);

                    d = v_star(j,i) - v_prev(j,i);
                    v_star(j,i) = v_prev(j,i) + omega*d;

                end
            end

            v_star(:,1) = -v_star(:,2); % west
            v_star(2,:) = 0;    % south
            v_star(ny+2,:) = 0; % north
            v_star(:,nx+2) = 2*v_star(:,nx+1) - v_star(:,nx);  % east  
            v_star(1,:) = 2*v_star(2,:) - v_star(3,:);

            sz = size(v_star);

            err = norm(v_star-v_prev)/sqrt(sz(1)*sz(2));
            iter = iter + 1;

        end


        % pressure correction step

        Ape = zeros(ny+2,nx+2);
        Apw = Ape;
        Apn = Ape;
        Aps = Ape;
        App = Ape;

        p_dash = zeros(ny+2,nx+2);
        omega = 1.9;
        err = 100;
        iter = 0;

        for i=2:nx+2
                for j=2:ny+1

                    me = 0.5*(u_star(j,i)+u_star(j,i+1))*dy;
                    mw = -0.5*(u_star(j,i)+u_star(j,i-1))*dy;
                    mn = 0.5*(v_star(j+1,i)+v_star(j+1,i-1))*dx;
                    ms = -0.5*(v_star(j,i)+v_star(j,i-1))*dx;

                    Aup = -4/Re - (max(0,me) + max(0,ms) + max(0,mn) + max(0,mw));  aup(j,i) = Aup;

                end
        end

        aup(:,2) = 1;

        for i=2:nx+1
            for j=2:ny+2

                me = 0.5*(u_star(j-1,i+1)+u_star(j,i+1))*dy;
                mw = -0.5*(u_star(j-1,i)+u_star(j,i))*dy;
                mn = 0.5*(v_star(j,i)+v_star(j+1,i))*dx;
                ms = -0.5*(v_star(j-1,i)+v_star(j,i))*dx;
        
                Avp = -4/Re - (max(0,me) + max(0,ms) + max(0,mn) + max(0,mw)); avp(j,i) = Avp;
            end
        end

        avp(2,:) = 1;

        for i=2:nx+1
                for j=2:ny+1

                    Ape(j,i) = dy/aup(j,i+1);
                    Apw(j,i) = dy/aup(j,i);
                    Apn(j,i) = dx/avp(j+1,i);
                    Aps(j,i) = dx/avp(j,i);
                    App(j,i) = -(Aps(j,i) + Apn(j,i) + Apw(j,i) + Ape(j,i));
                end
        end

        while abs(err) > 1e-10

            p_prev = p_dash;

            for i=2:nx+1
                for j=2:ny+1

                    Qpp = u_star(j,i) - u_star(j,i+1) + v_star(j,i) - v_star(j+1,i);

                    p_dash(j,i) = (Qpp - Ape(j,i)*p_dash(j,i+1) - Apw(j,i)*p_dash(j,i-1) - Apn(j,i)*p_dash(j+1,i) - Aps(j,i)*p_dash(j-1,i))/App(j,i);

                    d = p_dash(j,i) - p_prev(j,i);
                    p_dash(j,i) = p_prev(j,i) + omega*d;

                end
            end

            p_dash(:,nx+2) = -p_dash(:,nx+1);
            p_dash(:,1) = p_dash(:,2);
            p_dash(1,:) = p_dash(2,:);
            p_dash(ny+2,:) = p_dash(ny+1,:);

            sz = size(p_dash);

            err = norm(p_dash-p_prev)/sqrt(sz(1)*sz(2));
            iter = iter + 1;

        end
        
        
        % correct pressure and velocity

        for i=2:nx+2
          for j=2:ny+1
            u(j,i) = u_star(j,i) + dy/aup(j,i)*(p_dash(j,i)-p_dash(j,i-1));
          end
        end
    
        for i=2:nx+1
          for j=2:ny+2
            v(j,i) = v_star(j,i) + dx/avp(j,i)*(p_dash(j,i)-p_dash(j-1,i));
          end
        end

        % pressure and velocity under-relaxation

        p = alpha_P*p_dash + pold;
        p_star = p;

        u_star = alpha_u*u + (1-alpha_u)*uold; 
        v_star = alpha_v*v + (1-alpha_v)*vold;


        % RHS of x-momentum eqn for residual
        Qu = zeros(ny+2,nx+3);
    
        for ix=2:nx+2
            for iy=2:ny+1
                Qu(iy,ix) = dy*(p_star(iy,ix)-p_star(iy,ix-1));
            end
        end

        resuavg = calcres(aup,aue,auw,aus,aun,Qu,u_star,nx+3,ny+2);
        

        % RHS of y-momentum eqn for residual
        Qv = zeros(ny+3,nx+2);
    
        for ix=2:nx+1
            for iy=2:ny+2
                Qv(iy,ix) = dx*(p_star(iy,ix)-p_star(iy-1,ix));
            end
        end

        resvavg = calcres(avp,ave,avw,avs,avn,Qv,v_star,nx+2,ny+3);

        ITER = ITER + 1
        Residual(ITER,1) = resvavg + resuavg;
       
        % Error in u,v and p

        erru = norm(u_star-uold)/sqrt((nx+1)*ny);
        errv = norm(v_star-vold)/sqrt(nx*(ny+1));
        errp = norm(p_star-pold)/sqrt(nx*ny);

    end

    
    % Plotting Residual Vs Iterations for Re = 50
    if idx==2

        figure;
        plot(1:Niter,Residual);
        title("Residual Vs No. of Iterations");
        xlabel("Iterations");
        ylabel("Residual");
    end

    % Plotting isocontours of u(x,y)
    figure;
    contourf(xu,yp,u(2:ny+1,2:nx+2));
    title("Isocontours for u(x,y)");
    xlabel("x");
    ylabel("y");
    colorbar;
    
    % Plotting isocontours of v(x,y)
    figure;
    contourf(xp,yv,v(2:ny+2,2:ny+1));
    title("Isocontours for v(x,y)");
    xlabel("x");
    ylabel("y");
    colorbar;
    
    % Plotting isocontours of p(x,y)
    figure;
    contourf(xp,yp,p(2:ny+1,2:nx+1));
    title("Isocontours for p(x,y)");
    xlabel("x");
    ylabel("y");
    colorbar;

    % Interpolation of velocity field at pressure nodes
    vp = (v(2:ny+1,2:nx+1) + v(3:ny+2,2:nx+1))/2;
    up = (u(2:ny+1,2:nx+1) + u(2:ny+1,3:nx+2))/2;

    % Plotting streamlines of the velocity field
    [x,y] = meshgrid(xp,yp);
    startx = linspace(min(xp),max(xp),20);
    starty = linspace(min(yp),max(yp),20);
    [StartX,StartY] = meshgrid(startx, starty);
    figure;
    streamline(x,y,up,vp,StartX(:),StartY(:));
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Streamline Plot');

end


% function to evaluate residual (taken from code shared for lab7)
function resavg=calcres(ap,ae,aw,as,an,rhs,c,nxp,nyp)

  resavg=0.0; 

       for ix=1:nxp
          for iy=1:nyp
             res=rhs(iy,ix);

             if (ix>1) 
               res=res-aw(iy,ix)*c(iy,ix-1); 
             end

             if (iy>1) 
               res=res-as(iy,ix)*c(iy-1,ix);
             end

             if (iy<nyp) 
               res=res-an(iy,ix)*c(iy+1,ix);
             end

             if (ix<nxp) 
               res=res-ae(iy,ix)*c(iy,ix+1);
             end

             res=res-ap(iy,ix)*c(iy,ix);

             resavg=resavg+res^2;
           end
       end
     
     resavg=sqrt(resavg/(nxp*nyp));  
   
end