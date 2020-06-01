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
par.c    = [4 4.1];     % attraction between predators and preys 
par.R = 0.1;            % repulsion between predators 

par.g =1e-10;               % penalization term 

% selection of the control
disp ('[1] Control 1')
disp('[2] Control 2')
disp('[3] Control 3')
P1 = input('initial data 1-3: ');

% initial data
% preys positions
x = -0.5+rand(N,1);       
y = -1.5+3*rand(N,1);

% predators position
zx = [-1 ; 1];
zy = [0;0];

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
     
    Rx   = sqrt((x-x').^2 + (y-y').^2+eye(N));  % matrix of distances
    Rz   = sqrt((x-zx').^2 + (y-zy').^2);       % matrix of distances
    Rzz = sqrt((zx-zx').^2+(zy-zy').^2);        % matrix of distances
    
    FA   = A(Rx);  % repulsion and attraction force between preys
    FB = B(Rz);    % repulsion force between preys and predators 
    FC = C(Rz);    % attraction force between predators and preys
    
    % control 1 
    if P1==1
        xbar =sum(x)/N;
        ybar=sum(y)/N;
    end
    
    % control 2
    if P1==2 
        V1 = ((sum(zx)/2-x).^2+(sum(zy)/2-y).^2).^(1/2);
        V=V1<=(min(V1)+max(V1))/2;
        xbar=sum(V.*x)/sum(V);
        ybar=sum(V.*y)/sum(V);
    end
    
    % control 3
    if P1==3
        V1 = ((sum(zx)/2-x).^2+(sum(zy)/2-y).^2).^(1/2);

        [s,in]=sort(V1);
        Vx=0;
        Vy=0;
        for k=1:Nstar
            Vx=Vx+x(in(k));
            Vy=Vy+y(in(k));
        end
        xbar = Vx/Nstar;
        ybar = Vy/Nstar;
    end

    % preys forces
    vx = x; vy = y;
    for j=1:N
        vx (j) =  (sum(FA(j,:)'.*(x(j)-x)))./N + sum(FB(j,:)'.*(x(j)-zx));
        vy (j) =  (sum(FA(j,:)'.*(y(j)-y)))./N + sum(FB(j,:)'.*(y(j)-zy));
    end
    
    % Euler method 
    x=x+dt*vx;
    y=y+dt*vy;

    % predators force 
    for k=1:2
        vzx(k) =sum(par.c(k)*FC(:,k).*(x-zx(k)))/N-par.R*sum(zx(k)-zx);
        vzy(k) =sum(par.c(k)*FC(:,k).*(y-zy(k)))/N-par.R*sum(zy(k)-zy);

    end
    
    % control
    ux= (-sum(zx)-dt*sum(vzx)+2*xbar)/(2*dt+par.g);
    uy= (-sum(zy)-dt*sum(vzy)+2*ybar)/(2*dt+par.g);
    
    % Euler method 
    zx=zx+dt*sum(vzx)+dt*ux;
    zy=zy+dt*sum(vzy)+dt*uy;
    
    
    % time evolution
    plot(x,y,'bo','LineWidth',2.5)
    hold on
    plot(zx,zy,'ro','LineWidth',2.5)
    axis([-6,6,-6,6])
    %axis('equal')
    title([' Preys-predator model for t=' num2str(i+1),',for a=' num2str(par.a),...
    ',b =' num2str(par.b), ',c =' num2str(par.c),',p =' num2str(par.p)])
    pause(0.1)
    hold off
    drawnow
end
xlabel('x','FontSize',12,'FontWeight','bold')
ylabel('y','FontSize',12,'FontWeight','bold')
 
set(gca,'FontSize',12,'FontWeight','bold')
