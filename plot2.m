clc
clear all
close all 
% Plot for first case where we have the exact solutions

% Initial data set
% Choose the case (1 for this part)
data = 1;
g=1;
switch data
    case 1
        u=0.25;
        hIC =@(x) 1+0.5*sin(pi*x);
        mIC =@(x) u*hIC(x);
        S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
        bc  = 'Periodic';
    case 2
        u=0.25;
        hIC =@(x) 1+0.5*sin(pi*x);
        mIC =@(x) u*hIC(x);
        S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
        bc  = 'Open';
    case 3
        u=0;
        hIC=@(x) 1-0.1*sin(pi*x);
        mIC=@(x) 0.*x;
        S=@(x,t) 0*x;
        bc = 'Periodic';
    case 4
        u=0.25;
        hIC=@(x) 1-0.2*sin(2*pi*x);
        mIC=@(x) 0.5+0*x;
        S=@(x,t) 0*x;
        bc = 'Periodic';
end


% compute for two dx differents the solution at final time
M=[100,500];
for i=1:2
dx=2/M(i);
xf = 0:dx:2;
xc = (0.5*dx):dx:(2-0.5*dx);
N  = length(xc);
Tfinal = 2.0;
CFL    = 0.5;

% Averaging initial conditions
% Cell-center values sufficient for second-order schemes
Uavg = [hIC(xc);mIC(xc)];

% U1 --> solution with the choosen limiter
U1 = Uavg;
%Choose the limiters change based on the limiter we want to use
 limiter1='NONE';
% limiter1='MINMOD';
% limiter1='MUSCL';
%limiter1='MINMODTVB';

%initial time and iteration
time = 0;
iter=0;

% Solve
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(U1(2,:)./U1(1,:))+ sqrt(g*U1(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(U1,k,dx,g,N,bc,limiter1,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = U1 + k*RHS1;
    q2=(3/4)*U1+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter1,time+k,xc));
    q3=(1/3)*U1+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter1,time+k/2,xc));
    U1=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
if (i==1)
    U=U1;
    xc1=xc;
end
end

%plot the solutions-- change the legend depending on the case and the flux
figure(1)
subplot(2,2,1)
plot(xc1,U(1,:),'r');
hold on;
plot(xc1,hIC(xc1-time),'--k');
legend('h obtained by Roe flux', 'true solution h')
title('solution with dx=0.02 and T=2');
subplot(2,2,2)
plot(xc,U1(1,:),'r');
hold on;
plot(xc,hIC(xc-time),'--k');
legend('h obtained by Roe flux', 'true solution h')
title('solution with dx=0.004 and T=2');
subplot(2,2,3)
plot(xc1,U(2,:),'r');
hold on;
plot(xc1,mIC(xc1-time),'--k');
legend('m obtained by Roe flux', 'true solution m')
subplot(2,2,4)
plot(xc,U1(2,:),'r');
hold on;
plot(xc,mIC(xc-time),'--k');
legend('m obtained by Roe flux', 'true solution m')
%print('21bRoeminmodtvb1','-dpdf') 


%% repeat with initial conditions different so we need to compute a reference solution
clc
clear all
close all 

% Initial data set choose the data corresponding to the right case
data = 3;
g=1;
switch data
    case 1
        u=0.25;
        hIC =@(x) 1+0.5*sin(pi*x);
        mIC =@(x) u*hIC(x);
        S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
        bc  = 'Periodic';
    case 2
        u=0.25;
        hIC =@(x) 1+0.5*sin(pi*x);
        mIC =@(x) u*hIC(x);
        S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
        bc  = 'Open';
    case 3
        u=0;
        hIC=@(x) 1-0.1*sin(pi*x);
        mIC=@(x) 0.*x;
        S=@(x,t) 0*x;
        bc = 'Periodic';
    case 4
        u=0.25;
        hIC=@(x) 1-0.2*sin(2*pi*x);
        mIC=@(x) 0.5+0*x;
        S=@(x,t) 0*x;
        bc = 'Periodic';
    case 5
        hIC=@(x) 1+0.*x;
        mIC=@(x) -1.5*(x<=1);
        S=@(x,t) 0*x;
        bc = 'Open';
        
end


%First compute reference solution
dx=0.0005;
xf = 0:dx:2;
xc = (0.5*dx):dx:(2-0.5*dx);
N  = length(xc);
Tfinal = 2.0;
%Tfinal = 0.5;
CFL    = 0.5;

% Averaging initial conditions
% Cell-center values sufficient for second-order schemes
Uavg = [hIC(xc);mIC(xc)];

Uex = Uavg;

%Choose the limiters we use MUSCL for reference solution
%limiter1='NONE';
%limiter1='MINMOD';
limiter1='MUSCL';
% limiter1='MINMODTVB';

%initial time and iteration
time = 0;
iter=0;

% Solve
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(Uex(2,:)./Uex(1,:))+ sqrt(g*Uex(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(Uex,k,dx,g,N,bc,limiter1,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = Uex + k*RHS1;
    q2=(3/4)*Uex+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter1,time+k,xc));
    q3=(1/3)*Uex+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter1,time+k/2,xc));
    Uex=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
% Save the space cells too
xcex=xc;

%compute for different dx
M=[100,500];
for i=1:2
dx=2/M(i);
xf = 0:dx:2;
xc = (0.5*dx):dx:(2-0.5*dx);
N  = length(xc);
Tfinal = 2.0;
%Tfinal = 0.5;
CFL    = 0.5;

% Averaging initial conditions
% Cell-center values sufficient for second-order schemes
Uavg = [hIC(xc);mIC(xc)];

% U1 --> solution with the coosen slope limiter
U1 = Uavg;

%Choose the limiters
% limiter1='NONE';
 limiter1='MINMOD';
% limiter1='MUSCL';
%limiter1='MINMODTVB';

%initial time and iteration
time = 0;
iter=0;

% Solve for NONE
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(U1(2,:)./U1(1,:))+ sqrt(g*U1(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(U1,k,dx,g,N,bc,limiter1,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = U1 + k*RHS1;
    q2=(3/4)*U1+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter1,time+k,xc));
    q3=(1/3)*U1+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter1,time+k/2,xc));
    U1=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
if (i==1)
    U=U1;
    xc1=xc;
end
end

%plot the solution against reference one again change the legend depending
%on the flux used.

figure(1)
subplot(2,2,1)
plot(xc1,U(1,:),'r');
hold on;
plot(xcex,Uex(1,:),'--k');
legend('h obtained by Roe flux', 'true solution h')
title('solution with dx=0.02 and T=0.5');
subplot(2,2,2)
plot(xc,U1(1,:),'r');
hold on;
plot(xcex,Uex(1,:),'--k');
legend('h obtained by Roe flux', 'true solution h')
title('solution with dx=0.004 and T=0.5');
subplot(2,2,3)
plot(xc1,U(2,:),'r');
hold on;
plot(xcex,Uex(2,:),'--k');
legend('m obtained by Roe flux', 'true solution m')
subplot(2,2,4)
plot(xc,U1(2,:),'r');
hold on;
plot(xcex,Uex(2,:),'--k');
legend('m obtained by Roe flux', 'true solution m')
%print('Roeminmodtvb4','-dpdf') 



