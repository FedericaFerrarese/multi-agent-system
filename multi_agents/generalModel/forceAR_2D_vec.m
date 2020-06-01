function [vx,vy,vzx,vzy] = forceAR_2D_vec(t,x,y,zx,zy,par,P,N)
% INPUT 
% t = time
% x,y = preys positions
% zx,zy = predators positions
% par = parameters 
% P = number of predators 
% N = number of preys 
% OUTPUT
% vx,vy = preys forces
% vzx,vzy = predators forces
B   = @(r) par.b./(r.^2);        % preys predators repulsion
A  = @(r) r.^(-2)-par.a;         % preys preys repulsion and attraction
C= @(r)r.^(-par.p);              % attraction between predators and preys


Rx   = sqrt((x-x').^2 + (y-y').^2+eye(N));    % matrix of distances
Rz   = sqrt((x-zx').^2 + (y-zy').^2);         % matrix of distances


FA   = A(Rx);  % repulsion and attraction force between preys
FB = B(Rz);    % repulsion force between preys and predators 
FC = C(Rz);    % attraction force between predators and preys

if P==1
    vx = x; vy = y;
    % preys forces 
for i=1:N
    vx (i) =  (sum(FA(i,:)'.*(x(i)-x)))./N + FB(i).*(x(i)-zx);
    vy (i) =  (sum(FA(i,:)'.*(y(i)-y)))./N + FB(i).*(y(i)-zy);
end
    % predator force 
    vzx = sum(par.c*FC.*(x-zx))/N;
    vzy = sum(par.c*FC.*(y-zy))/N;   
end


if P>1
    vx = x; vy = y; vzx=zx; vzy =zy;
    
    % preys forces 
    for i=1:N
        vx (i) =  (sum(FA(i,:)'.*(x(i)-x)))./N + sum(FB(i,:)'.*(x(i)-zx));
        vy (i) =  (sum(FA(i,:)'.*(y(i)-y)))./N + sum(FB(i,:)'.*(y(i)-zy));
    end
    
    % predators force 
     for k=1:P
         vzx(k) =sum(par.c(k)*FC(:,k).*(x-zx(k)))/N+par.R*sum((zx(k)-zx));
         vzy(k) =sum(par.c(k)*FC(:,k).*(y-zy(k)))/N+par.R*sum((zy(k)-zy));
    end
end   
end