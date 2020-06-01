% main file 
clc
close all
clear all

% Fix parameters
Nstar = 3;                  % number of nearest preys 


N     = 500;                % number of preys
dt    = 0.1;                % time step
tf    = 10; t0 = 0;         % time interval 
steps = (tf-t0)/dt + 1;     % number of steps


par.p      = 3;         % distance influence
par.a      = 1;         % attraction between preys
par.b      =1.5;        % repulsion between preys and predators 
par.c    = 1;           % attraction between predators and preys 

par.g = 1e-4;           % penalization term for control
    
% selection of the control
disp ('[1] Control 1')
disp('[2] Control 2')
disp('[3] Control 3')
P = input('initial data 1-3: ');


% initial data
% preys positions
x           =  -0.5+rand(N,1);        
y           =  -1.5+3*rand(N,1);

% predator position
zx          = -1+2*rand;     
zy          =  -1+2*rand;     

% plot of initial data
figure
plot(x,y,'bo','LineWidth',2.5)
hold on
plot(zx,zy,'ro','LineWidth',2.5)
axis([-1,1,-1,1])
axis('equal')
xlabel('x','FontSize',12,'FontWeight','bold')
ylabel('y','FontSize',12,'FontWeight','bold')
 title([' Preys-predator model for t=' num2str(1),',for a=' num2str(par.a),...
    ',b =' num2str(par.b), ',c =' num2str(par.c),',p =' num2str(par.p)])
set(gca,'FontSize',12,'FontWeight','bold')


B   = @(r) par.b./(r.^2);        % preys predators repulsion
A  = @(r) r.^(-2)-par.a;         % preys preys repulsion and attraction
C= @(r)r.^(-par.p);              % attraction between predators and preys


figure
for i = 1:steps-1
     
    Rx   = sqrt((x-x').^2 + (y-y').^2+eye(N)); % matrix of distances
    Rz   = sqrt((x-zx').^2 + (y-zy').^2);      % matrix of distances
    FA   = A(Rx);  % repulsion and attraction force between preys
    FB = B(Rz);    % repulsion force between preys and predators 
    FC = C(Rz);    % attraction force between predators and preys
    
    % control 1 
    if P==1
        % centre of mass 
        xbar =sum(x)/N;
        ybar=sum(y)/N;
    end
    
    % control 2 
    if P==2
        R = (max(Rz)+min(Rz))/2; % radius
        V=Rz<=R; % truth vector
        % center of mass of the ball
        xbar=sum(V.*x)/sum(V);
        ybar=sum(V.*y)/sum(V);
    end
    
    % control 3
    if P==3 
        R = (max(Rz)+min(Rz))/2; % radius
        [s,in]=sort(Rz); % sort the distances 
        Vx=0;
        Vy=0;
        % Choose the N* nearest values 
        for k=1:Nstar
            Vx=Vx+x(in(k));
            Vy=Vy+y(in(k));
        end
        % centre of mass of the N* nearest preys
        xbar = Vx/Nstar;
        ybar = Vy/Nstar;
    end
    
    % preys forces
    vx = x; vy = y;
    for j=1:N
        vx (j) =  (sum(FA(j,:)'.*(x(j)-x)))./N + FB(j).*(x(j)-zx);
        vy (j) =  (sum(FA(j,:)'.*(y(j)-y)))./N + FB(j).*(y(j)-zy);
    end
    
    % Euler method 
    x=x+dt*vx;
    y=y+dt*vy;

    % compute the control 
    ux= (-zx-dt*sum(FC.*(x-zx))/N+xbar)/(dt+par.g);
    uy= (-zy-dt*sum(FC.*(y-zy))/N+ybar)/(dt+par.g);
    
    % Euler method 
    zx = zx +dt/N*sum(FC.*(x-zx))+dt*ux;
    zy= zy +dt/N*sum(FC.*(y-zy))+dt*uy;
    
    % plot of the time evolution
    plot(x,y,'bo','LineWidth',2.5)
    hold on
    plot(zx,zy,'ro','LineWidth',2.5)
    axis([-5,5,-5,5])
    title([' Preys-predator model for t=' num2str(i+1),',for a=' num2str(par.a),...
    ',b =' num2str(par.b), ',c =' num2str(par.c),',p =' num2str(par.p)])
    pause(0.1)
    hold off
    drawnow
end
xlabel('x','FontSize',12,'FontWeight','bold')
ylabel('y','FontSize',12,'FontWeight','bold')
 
set(gca,'FontSize',12,'FontWeight','bold')
