% Function evaluates the RHS corresponding to a second-order MUSCL scheme
function RHS = evalRHS(U,k,h,g,N,bc,limiter,time,xc)


% Need to extend by 2 ghost cells on either side
U_ext = apply_bc(U,bc,2);
u=0.25;
hIC =@(x) 1+0.5*sin(pi*x);
S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
% Obtain limited slope for N+2 cells
dU      = zeros(2,N+2);
dU(1,:) = SlopeLimiter(U_ext(1,:),limiter,N,h);
dU(2,:) = SlopeLimiter(U_ext(2,:),limiter,N,h);

UL=U_ext(:,2:N+3)-dU/2;
UR=U_ext(:,2:N+3)+dU/2;

%change the flux we want to use
% LF Flux
flux1=LaxFriedFlux(UR(:,2:N+1),UL(:,3:N+2),k,g,h);
flux2=LaxFriedFlux(UR(:,1:N),UL(:,2:N+1),k,g,h);
%Roe Flux
%  flux1=RoeFlux(UR(:,2:N+1),UL(:,3:N+2),g);
%  flux2=RoeFlux(UR(:,1:N),UL(:,2:N+1),g);

% change RHS if we use initial condition with zero source term or not
RHS=-(1/h)*(flux1- flux2)+S(xc,time);
%RHS=-(1/h)*(flux1- flux2);

return

