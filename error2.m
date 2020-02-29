clc
clear all
close all 

% define the vector where we plug in the errors
error_h1=zeros(1,4);
error_h2=zeros(1,4);
error_h3=zeros(1,4);
error_h4=zeros(1,4);
error_m1=zeros(1,4);
error_m2=zeros(1,4);
error_m3=zeros(1,4);
error_m4=zeros(1,4);


% choose the case for the initial conditions and boundary conditions
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

% Define for compute the error
M=[100,200,500,1000];
for i=1:4
    dx=2./M(i);
    xf = 0:dx:2;
    xc = (0.5*dx):dx:(2-0.5*dx);
    N  = length(xc);
    Tfinal = 2.0;
    CFL    = 0.5;

    % Averaging initial conditions
    % Cell-center values sufficient for second-order schemes
    Uavg = [hIC(xc);mIC(xc)];

    % U1 --> solution with no slope limiting
    % U2 --> solution with minmod limiting
    % U3 --> solution with MUSCL limiting
    % U4 --> solution with MUSCL TVB limiter
    U1 = Uavg;
    U2 = Uavg;
    U3 = Uavg;
    U4 = Uavg;

    %Choose the limiters
    limiter1='NONE';
    limiter2='MINMOD';
    limiter3='MUSCL';
    limiter4='MINMODTVB';

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
error_h1(i)=abs(max(abs(U1(1,:)-hIC(xc-time))));
error_m1(i)=abs(max(abs(U1(2,:)-u*hIC(xc-time))));

time = 0;
iter=0;

% Solve for MINMOD
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(U2(2,:)./U2(1,:))+ sqrt(g*U2(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(U2,k,dx,g,N,bc,limiter2,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = U2 + k*RHS1;
    q2=(3/4)*U2+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter2,time+k,xc));
    q3=(1/3)*U2+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter2,time+k/2,xc));
    U2=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
error_h2(i)=abs(max(abs(U2(1,:)-hIC(xc-time))));
error_m2(i)=abs(max(abs(U2(2,:)-u*hIC(xc-time))));
%solve for MUSCL
time = 0;
iter=0;

% Solve for NONE
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(U3(2,:)./U3(1,:))+ sqrt(g*U3(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(U3,k,dx,g,N,bc,limiter3,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = U3 + k*RHS1;
    q2=(3/4)*U3+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter3,time+k,xc));
    q3=(1/3)*U3+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter3,time+k/2,xc));
    U3=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
error_h3(i)=abs(max(abs(U3(1,:)-hIC(xc-time))));
error_m3(i)=abs(max(abs(U3(2,:)-u*hIC(xc-time))));
%solve for MINMODTVB
time = 0;
iter=0;
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(U4(2,:)./U4(1,:))+ sqrt(g*U4(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(U4,k,dx,g,N,bc,limiter4,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = U4 + k*RHS1;
    q2=(3/4)*U4+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter4,time+k,xc));
    q3=(1/3)*U4+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter4,time+k/2,xc));
    U4=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
error_h4(i)=abs(max(abs(U4(1,:)-hIC(xc-time))));
error_m4(i)=abs(max(abs(U4(2,:)-u*hIC(xc-time))));
end


%plot the errors change the number of the error to change the error for the limiter we plot
figure();
subplot(121)
loglog(1./M,error_h1);
hold on;
loglog(1./M,(1./M).^1)
legend('error','O(Ne-1)')
title('Error of numerical solution h(x,t) at T=2')
subplot(122)
loglog(1./M,error_m1);
hold on;
loglog(1./M,(1./M).^1)
legend('error','O(Ne-1)')
title('Error of numerical solution m(x,t) at T=2')
%print('21bLF0err','-dpdf') 



%% error case 2

%Here we do kind of the same thing as before but we need to compute a
%reference solution 
clc
clear all
close all 
error_h1=zeros(1,4);
error_h2=zeros(1,4);
error_h3=zeros(1,4);
error_h4=zeros(1,4);
error_m1=zeros(1,4);
error_m2=zeros(1,4);
error_m3=zeros(1,4);
error_m4=zeros(1,4);


%compute reference solution
data = 4;
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


%First compute reference solution using muscl

dx=0.0005;
xf = 0:dx:2;
xc = (0.5*dx):dx:(2-0.5*dx);
N  = length(xc);
Tfinal = 2.0;
CFL    = 0.5;

% Averaging initial conditions
% Cell-center values sufficient for second-order schemes
Uavg = [hIC(xc);mIC(xc)];

Uex = Uavg;


%Choose the limiters
limiter3='NONE';
limiter2='MINMOD';
limiter1='MUSCL';
limiter4='MINMODTVB';

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
xcex=xc;

% now compute the error
M=[100,200,500,1000];
for i=1:4
    dx=2./M(i);
    xf = 0:dx:2;
    xc = (0.5*dx):dx:(2-0.5*dx);
    N  = length(xc);
    Tfinal = 2.0;
    CFL    = 0.5;

    % Averaging initial conditions
    % Cell-center values sufficient for second-order schemes
    Uavg = [hIC(xc);mIC(xc)];

    % U1 --> solution with no slope limiting
    % U2 --> solution with minmod limiting
    % U3 --> solution with MUSCL limiting
    % U4 --> solution with MUSCL TVB limiter
    U1 = Uavg;
    U2 = Uavg;
    U3 = Uavg;
    U4 = Uavg;

    %Choose the limiters
    limiter1='NONE';
    limiter2='MINMOD';
    limiter3='MUSCL';
    limiter4='MINMODTVB';

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
error_h1(i)=abs(max(abs(U1(1,:)-interp1(xcex,Uex(1,:),xc,'linear'))));
error_m1(i)=abs(max(abs(U1(2,:)-interp1(xcex,Uex(2,:),xc,'linear'))));

time = 0;
iter=0;

% Solve for MINMOD
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(U2(2,:)./U2(1,:))+ sqrt(g*U2(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(U2,k,dx,g,N,bc,limiter2,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = U2 + k*RHS1;
    q2=(3/4)*U2+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter2,time+k,xc));
    q3=(1/3)*U2+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter2,time+k/2,xc));
    U2=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
error_h2(i)=abs(max(abs(U2(1,:)-interp1(xcex,Uex(1,:),xc,'linear'))));
error_m2(i)=abs(max(abs(U2(2,:)-interp1(xcex,Uex(2,:),xc,'linear'))));
%solve for MUSCL
time = 0;
iter=0;

% Solve for muscl
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(U3(2,:)./U3(1,:))+ sqrt(g*U3(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(U3,k,dx,g,N,bc,limiter3,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = U3 + k*RHS1;
    q2=(3/4)*U3+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter3,time+k,xc));
    q3=(1/3)*U3+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter3,time+k/2,xc));
    U3=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
error_h3(i)=abs(max(abs(U3(1,:)-interp1(xcex,Uex(1,:),xc,'linear'))));
error_m3(i)=abs(max(abs(U3(2,:)-interp1(xcex,Uex(2,:),xc,'linear'))));
%solve for MINMODTVB
time = 0;
iter=0;
while time < Tfinal
    % compute timestep
    k = CFL*dx/max(abs(U4(2,:)./U4(1,:))+ sqrt(g*U4(1,:)));
    if(time + k > Tfinal)
        k = Tfinal - time;
    end
    %compute RHS
    RHS1 = evalRHS(U4,k,dx,g,N,bc,limiter4,time,xc);	
    % Update solution to k by SSP-RK3
    q1 = U4 + k*RHS1;
    q2=(3/4)*U4+(1/4)*(q1+k*evalRHS(q1,k,dx,g,N,bc,limiter4,time+k,xc));
    q3=(1/3)*U4+(2/3)*(q2+k*evalRHS(q2,k,dx,g,N,bc,limiter4,time+k/2,xc));
    U4=q3;
    %update time and iteration
    time = time + k;
    iter=iter+1;
end
error_h4(i)=abs(max(abs(U4(1,:)-interp1(xcex,Uex(1,:),xc,'linear'))));
error_m4(i)=abs(max(abs(U4(2,:)-interp1(xcex,Uex(2,:),xc,'linear'))));
end

%plot the errors changing the number of the vector to change limiter
figure();
subplot(121)
loglog(1./M,error_h3);
hold on;
loglog(1./M,(1./M).^2)
legend('error','O(Ne-2)')
title('Error of numerical solution h(x,t) at T=2')
subplot(122)
loglog(1./M,error_m3);
hold on;
loglog(1./M,(1./M).^2)
legend('error','O(Ne-2)')
title('Error of numerical solution m(x,t) at T=2')
print('31bLFminmodtvberr3','-dpdf') 




