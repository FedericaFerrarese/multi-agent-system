function [dz]=modelAR_2D(t,z,par,P)
% INPUT
% t = time
% z = initial position
% par = parameters 
% P = number of predators 
% OUTPUT
% dz = preys and predators forces 
D = 1; 

N  = length(z(1:end-P*2))/(2*D); %number of agents

% preys and predators positions 
x  = z(1:N);
y  = z(N+1:2*N);
zx = z(2*N+1:2*N+P);
zy = z(2*N+P+1:2*N+2*P);


x   = x(:); zx = zx(:); % column vectors
y   = y(:); zy = zy(:); % column vectors

% compute the forces 
[vx,vy,vzx,vzy] = forceAR_2D_vec(t,x,y,zx,zy,par,P,N);

dzx  = vzx;
dzy  = vzy;
dx   = vx;
dy   = vy;

% preys and predators forces 
dz = [dx;dy;dzx;dzy];
end