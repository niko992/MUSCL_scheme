function dU = SlopeLimiter(U,limiter,N,h)
M=5;


dUL = U(:,2:N+3)-U(:,1:N+2);
dUR = U(:,3:N+4)-U(:,2:N+3);
switch limiter
    case 'NONE'
        dU = 0*dUL;
    case 'MINMOD'
        dU = minmod([dUL',dUR'])';
    case 'MUSCL'
        dU = minmod([0.5*(dUL'+dUR'),2*dUL',2*dUR'])';
    case 'MINMODTVB'
        dU = minmodTVB([dUL',dUR'],M,h);
    otherwise
        error('Unknown limiter function requested!')
        
end
        

return