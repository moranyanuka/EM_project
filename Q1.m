clear ;close all

% Question 1
%----------------- A1  -----------------

%initialize variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nvariables for Q1
p1.a = 1 + 0.5*ceil((0.25*(0+7+1))); %[m]
p1.b = (3*p1.a)/16; %[m]
p1.sigma = 1; %[1/ohm*m]
p1.v0 = 1;%[v]
p1.N = 140;
p1.h = p1.a/(p1.N-1);%[m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%variables for the RMS error
pt.a = 1 + 0.5*ceil((0.25*(0+7+1))); %[m]
pt.b = (3*p1.a)/16; %[m]
pt.sigma = 1; %[1/ohm*m]
pt.v0 = 1;%[v]
pt.N = 70;
pt.h = pt.a/(pt.N-1);%[m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pot = setPlatePot(3,2,p1); %return the potential Matrix when the electrode is at (n=2,m=3)
potTest = setPlatePot(3,2,pt);%potential matrix with half the number of points
E_rms = Erms(pot,potTest);%check the relative error
figure(1)
imagesc(0:p1.h:p1.a,p1.a:-p1.h:0,pot) %plot the potential
hold on;
set(gca, 'YDir','reverse');%flip the y axis
axis xy;
colorbar;
title('Potential Distribution');
xlabel('x[m]');
ylabel('y[m]');
hold off;

[Jx,Jy] = gradient(-1.*pot,p1.h); %sigma = 1 => sigma*E = j => E=J
figure(2)
hold on
streamslice(0:p1.h:p1.a,0:p1.h:p1.a,Jx,Jy,2) %plot current lines
set(gca, 'YDir','reverse');%flip the y axis
set(gca,'YtickLabel',p1.a:-0.5:0)
xlabel('x[m]')
ylabel('y[m]')
title('Current Flow')
%b
[Z,currConserv] = calcImp(p1);% Z - impedance matrix, currConserv - current conservation (%)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions



function pot = setPlatePot(m,n,p)

xn = ((6*n-2)*p.a)/32; %ground center position
ym = ((6*m-2)*p.a)/32; %electrode center position
sol_v = zeros(1,p.N^2); % the solutions vector
M = zeros(p.N^2,p.N^2); %the potential coefficients

%electrode and ground positions
elecEnd = ym + p.b/2; 
elecStart = ym - p.b/2;
groundStart = xn - p.b/2;
groundEnd = xn + p.b/2;

% initiate equations matrix estimation
i = 1;
for x = 0:p.N-1
    for y = 0:p.N-1
        if x==0 && ((y * p.h <= elecEnd) && (y * p.h >= elecStart)) %the electrode zone - dirichle condition
            sol_v(i) = p.v0;
            M(i,i) = 1;
        elseif  y == p.N-1 && ((x * p.h <= groundEnd) && (x * p.h >= groundStart)) %the ground zone- dirichle condition
            sol_v(i) = 0;
            M(i,i) = 1;
        elseif x==0 %left edge - noyman condition
            M(i,i) = 1;
            M(i,i+p.N) = -1;
        elseif x == p.N-1 %right edge - noyman condition
            M(i,i) = 1;
            M(i,i-p.N) = -1;
        elseif y == 0 %upper edge - noyman condition
            M(i,i+1) = -1;
            M(i,i) = 1;
        elseif y == p.N-1 %lower edge - noyman condition
            M(i,i-1) = -1;
            M(i,i) = 1;
        else %inside the plate - laplace estimation
            M(i,i) = -4;
            M(i,i+1) = 1;
            M(i,i+p.N) = 1;
            M(i,i-1) = 1;
            M(i,i-p.N) = 1;
        end
        i = i +1;
    end
end

M = sparse(M);
phi = M\sol_v'; % calc the potential vector (solve system of equations)  
pot = zeros(p.N,p.N); 

%reorganize the potentials back in the plate Matrix
i = 1;
for y = 1:p.N
    for x = 1:p.N
        pot(x,y) = phi(i);
        i = i+1;
    end
end           
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z,current] = calcImp(p)
current = zeros(5);
Z = zeros(5);
elec = 0; 
ground = 0;
for m = 1:5
    for n = 1:5
        curr_pot = setPlatePot(m,n,p);
        [Jx,Jy] = gradient(-1.*curr_pot,p.h);       
        for i = 1:p.N
            elec = elec + Jx(i,1)*p.h; %sum over the current out of the electrode
            ground = ground + Jy(p.N,i)*p.h; %sum over the current that goes in to the ground
        end
        
        current(m,n) = abs(elec-ground)/((elec+ground)/2); %deviation from conservation of current
        Z(m,n) = p.v0/((elec+ground)/2);%find the impedance
        elec = 0;
        ground =0;
    end
end
Z = flipud(Z);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = Erms(pot,potTest)
newPot = zeros(length(potTest));%build a new matrix with the size of potTest
x = 1;
y = 1;
for i = 1:2:length(pot)-1
    for j = 1:2:length(pot)-1        
        newPot(x,y) = (pot(i,j)+pot(i+1,j)+pot(i,j+1)+pot(i+1,j+1))/4;
        y=y+1;
    end
    y=1;
    x=x+1;
end
Eror = sqrt(sum((potTest-newPot).^2));
dividor = sqrt(sum(potTest.^2));
E = Eror/dividor;
end
