
clear ;close all
% Question 2
%-----------------   -----------------

p2.a = 1 + 0.5*ceil((0.25*(0+7+1)));%[m]
p2.b = (3*p2.a)/16;%[m]
p2.sigma = 1;%[1/ohm]
p2.v0 = 1;%[v]
p2.N = 141;
p2.h = p2.a/(p2.N-1);%[m]

%the foreign body position:
dx = 0.25*(ceil(0.25*(0+1)));%[m]
dy = 0.25*(ceil(0.25*(7+1)));%[m]
xc = p2.a/2; 
yc = p2.a/2; 
p2.sigma1 = (1+0.8)*1; %z3-z7 = 6-2 = 4 => sign(4) = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p1.a = 1 + 0.5*ceil((0.25*(0+7+1)));
p1.b = (3*p1.a)/16;
p1.sigma = 1;
p1.v0 = 1;
p1.N = 141;
p1.h = p1.a/(p1.N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2.1
%plot the potential

pot = setPlatePot2(3,2,p2,xc,yc,dx,dy);
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

%plot the current

[Jx,Jy] = gradient(-1.*pot,p1.h); %sigma = 1 => sigma*E = j => E=J
xelec = floor((xc-dx/2)/p2.h):ceil((xc+dx/2)/p2.h);
yelec = floor((yc-dy/2)/p2.h):ceil((yc+dy/2)/p2.h);
Jx(xelec,yelec) = p2.sigma1.*Jx(xelec,yelec);
Jy(xelec,yelec) = p2.sigma1.*Jy(xelec,yelec);
figure(2)
hold on
streamslice(0:p1.h:p1.a,0:p1.h:p1.a,Jx,Jy,2) %plot current lines
set(gca, 'YDir','reverse');%flip the y axis
set(gca,'YtickLabel',p1.a:-0.5:0)
title('Current Flow') 
xlabel('x[m]');
ylabel('y[m]');
hold off


%calculate  the impedance matrix
[Z1,current1] = calcImp(p1);
[Z2,current2] = calcImp2(p2,xc,yc,dx,dy);

D = calcD(Z1,Z2);

%2.2

xc = 0.35*p2.a; 
yc = 0.4*p2.a; 

pot = setPlatePot2(3,2,p2,xc,yc,dx,dy);
figure(3)
imagesc(0:p1.h:p1.a,p1.a:-p1.h:0,pot) %plot the potential
hold on;
set(gca, 'YDir','reverse');%flip the y axis
axis xy;
colorbar;
title('Potential Distribution');
xlabel('x[m]');
ylabel('y[m]');
hold off;

%plot the current

[Jx,Jy] = gradient(-1*pot,p2.h); %sigma = 1 => sigma*E = j => E=J
xelec = floor((xc-dx/2)/p2.h):ceil((xc+dx/2)/p2.h);
yelec = floor((yc-dy/2)/p2.h):ceil((yc+dy/2)/p2.h);
Jx(xelec,yelec) = p2.sigma1.*Jx(xelec,yelec);
Jy(xelec,yelec) = p2.sigma1.*Jy(xelec,yelec);
figure(4)
hold on
streamslice(0:p1.h:p1.a,0:p1.h:p1.a,Jx,Jy,2) %plot current lines
set(gca, 'YDir','reverse');%flip the y axis
set(gca,'YtickLabel',p1.a:-0.5:0)
title('Current Flow')
xlabel('x[m]');
ylabel('y[m]');
hold off

[Z3,current3] = calcImp2(p2,xc,yc,dx,dy);

D2 = calcD(Z1,Z3);


%2.3

xc = p2.a/2; 
yc = p2.a/2; 
dxx = 0.3*dx; %[m]
dyy = 0.3*dy; %[m]
p2.sigma2 = 0.5*p2.sigma1; %[1/ohm]


pot = setPlatePot3(3,2,p2,xc,yc,dx,dy,dxx,dyy);
figure(5)
imagesc(0:p1.h:p1.a,p1.a:-p1.h:0,pot) %plot the potential
hold on;
set(gca, 'YDir','reverse');%flip the y axis
axis xy;
colorbar;
title('Potential Distribution');
xlabel('x[m]');
ylabel('y[m]');
hold off;

%plot the current

[Jx,Jy] = gradient(-1.*pot,p1.h); %sigma = 1 => sigma*E = j => E=J
xelec = floor((xc-dx/2)/p2.h):ceil((xc+dx/2)/p2.h);
yelec = floor((yc-dy/2)/p2.h):ceil((yc+dy/2)/p2.h);
xcover = floor((xc-dx/2-dxx)/p2.h):ceil((xc+dx/2+dxx)/p2.h);
ycover = floor((yc-dy/2-dyy)/p2.h):ceil((yc+dy/2+dyy)/p2.h);
Jx(xcover,ycover) = p2.sigma2.*Jx(xcover,ycover);
Jy(xcover,ycover) = p2.sigma2.*Jy(xcover,ycover);
Jx(xelec,yelec) = p2.sigma1.*Jx(xelec,yelec);
Jy(xelec,yelec) = p2.sigma1.*Jy(xelec,yelec);
figure(6)
hold on
streamslice(0:p1.h:p1.a,0:p1.h:p1.a,Jx,Jy,2) %plot current lines
set(gca, 'YDir','reverse');%flip the y axis
set(gca,'YtickLabel',p1.a:-0.5:0)
title('Current Flow')
xlabel('x[m]');
ylabel('y[m]');
hold off

[Z4,current4] = calcImp3(p2,xc,yc,dx,dy,dxx,dyy);

D3 = calcD(Z1,Z4);



%|--------------functions-------------------|

function pot = setPlatePot(m,n,p)

xn = ((6*n-2)*p.a)/32; %ground center position
ym = ((6*m-2)*p.a)/32; %electrode center position
sol_v = zeros(1,p.N^2); % the solutions vector
M = zeros(p.N^2,p.N^2); %the potential coefficients

% electrode and ground positions
elecEnd = ym + p.b/2; 
elecStart = ym - p.b/2;
groundStart = xn - p.b/2;
groundEnd = xn + p.b/2;

% initiante equations matrix estimation
i = 1;
for x = 0:p.N-1
    for y = 0:p.N-1
        if x==0 && ((y * p.h <= elecEnd) && (y * p.h >= elecStart)) %the electrode zone - dirichle condition
            sol_v(i) = p.v0;
            M(i,i) = 1;
        elseif  y == p.N-1 && ((x * p.h <= groundEnd) && (x * p.h >= groundStart)) %the ground zone- dirichle condition
            sol_v(i) = 0;
            M(i,i) = 1;
        elseif x == 0 %left edge - noyman condition
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



function pot = setPlatePot2(m,n,p,xc,yc,dx,dy)

xn = ((6*n-2)*p.a)/32; %grond center position
ym = ((6*m-2)*p.a)/32; %electrode center position
sol_v = zeros(1,p.N^2); % the solutions vector
M = zeros(p.N^2,p.N^2); %the potential coefficients

% electrode's and ground's starting and ending positions
elecStart = ym - p.b/2;
elecEnd = ym + p.b/2;
groundStart = xn - p.b/2;
groundEnd = xn + p.b/2;

% the foreign body's starting and ending positions
xstart = xc - (dx/2);
xend = xc + (dx/2);
ystart = yc -(dy/2);
yend = yc +(dy/2);


% initiante equations matrix estimation
i = 1;
for x = 0:p.N-1
    for y = 0:p.N-1
        if ( (abs(x * p.h - xstart) <= p.h/2) && ( (y*p.h < yend) && (y*p.h > ystart)) )%the left edge of the foreign body
            M(i,i) = -1*p.sigma1 - p.sigma;
            M(i,i-p.N) = p.sigma;
            M(i,i+p.N) = p.sigma1;
        elseif ((abs(x * p.h - xend) <= p.h/2) && ( (y*p.h < yend) && (y*p.h > ystart)) )%the right edge of the foreign body
            M(i,i) = -1*p.sigma1 - p.sigma;
            M(i,i-p.N) = p.sigma1;
            M(i,i+p.N) = p.sigma; 
        elseif ((abs(y * p.h - ystart) <= p.h/2) && ( (x*p.h < xend) && (x*p.h > xstart)) )%the upper edge of the foreign body
            M(i,i) = -1*p.sigma1 - p.sigma;
            M(i,i-1) = p.sigma;
            M(i,i+1) = p.sigma1;
        elseif ((abs(y * p.h - yend) <= p.h/2) && ( (x*p.h < xend) && (x*p.h > xstart)) )%the lower edge of the foreign body
            M(i,i) = -1*p.sigma1 - p.sigma;
            M(i,i-1) = p.sigma1;
            M(i,i+1) = p.sigma;
        elseif x == 0 && ((y * p.h < elecEnd) && (y * p.h > elecStart)) %the electrode zone - dirichle condition
            sol_v(i) = p.v0;
            M(i,i) = 1;
        elseif  y == p.N-1 && ((x * p.h < groundEnd) && (x * p.h > groundStart)) %the ground zone- dirichle condition
            sol_v(i) = 0;
            M(i,i) = 1;
        elseif x == 0 %left edge - noyman condition
            M(i,i) = 1;
            M(i,i+p.N) = -1;
        elseif x == p.N-1 %right edge - noyman condition
            M(i,i-p.N) = 1;
            M(i,i) = -1;
        elseif y == 0 %upper edge - noyman condition
            M(i,i+1) = 1;
            M(i,i) = -1;
        elseif y == p.N-1 %lower edge - noyman condition
            M(i,i-1) = 1;
            M(i,i) = -1;
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

function pot = setPlatePot3(m,n,p,xc,yc,dx,dy,dxx,dyy)

xn = ((6*n-2)*p.a)/32; %grond center position
ym = ((6*m-2)*p.a)/32; %electrode center position
sol_v = zeros(1,p.N^2); % the solutions vector
M = zeros(p.N^2,p.N^2); %the potential coefficients

% electrode's and ground's starting and ending positions
elecStart = ym - p.b/2;
elecEnd = ym + p.b/2;
groundStart = xn - p.b/2;
groundEnd = xn + p.b/2;

% the foreign body's starting and ending positions
xstart = xc - (dx/2);
xend = xc + (dx/2);
ystart = yc -(dy/2);
yend = yc +(dy/2);

%the cover's starting and ending positions
xcstart = xstart - dxx;
xcend = xend + dxx;
ycstart = ystart - dyy;
ycend = yend + dyy;

% initiante equations matrix estimation
i = 1;
for x = 0:p.N-1
    for y = 0:p.N-1
        if ( (abs(x * p.h - xcstart) <= p.h/2) && ( (y*p.h < ycend) && (y*p.h > ycstart)) )%the left edge of the cover
            M(i,i) = -1*p.sigma2 - p.sigma1;
            M(i,i-p.N) = p.sigma1;
            M(i,i+p.N) = p.sigma2;
        elseif ((abs(x * p.h - xcend) <= p.h/2) && ( (y*p.h < ycend) && (y*p.h > ycstart)) )%the right edge of the cover
            M(i,i) = -1*p.sigma2 - p.sigma1;
            M(i,i-p.N) = p.sigma2;
            M(i,i+p.N) = p.sigma1; 
        elseif ((abs(y * p.h - ycstart) <= p.h/2) && ( (x*p.h < xcend) && (x*p.h > xcstart)) )%the upper edge of the cover
            M(i,i) = -1*p.sigma2 - p.sigma1;
            M(i,i-1) = p.sigma1;
            M(i,i+1) = p.sigma2;
        elseif ((abs(y * p.h - ycend) <= p.h/2) && ( (x*p.h < xcend) && (x*p.h > xcstart)) )%the lower edge of the foreign body
            M(i,i) = -1*p.sigma2 - p.sigma1;
            M(i,i-1) = p.sigma2;
            M(i,i+1) = p.sigma1; 
        elseif ( (abs(x * p.h - xstart) <= p.h/2) && ( (y*p.h < yend) && (y*p.h > ystart)) )%the left edge of the foreign body
            M(i,i) = -1*p.sigma1 - p.sigma;
            M(i,i-p.N) = p.sigma;
            M(i,i+p.N) = p.sigma1;
        elseif ((abs(x * p.h - xend) <= p.h/2) && ( (y*p.h < yend) && (y*p.h > ystart)) )%the right edge of the foreign body
            M(i,i) = -1*p.sigma1 - p.sigma;
            M(i,i-p.N) = p.sigma1;
            M(i,i+p.N) = p.sigma; 
        elseif ((abs(y * p.h - ystart) <= p.h/2) && ( (x*p.h < xend) && (x*p.h > xstart)) )%the upper edge of the foreign body
            M(i,i) = -1*p.sigma1 - p.sigma;
            M(i,i-1) = p.sigma;
            M(i,i+1) = p.sigma1;
        elseif ((abs(y * p.h - yend) <= p.h/2) && ( (x*p.h < xend) && (x*p.h > xstart)) )%the lower edge of the foreign body
            M(i,i) = -1*p.sigma1 - p.sigma;
            M(i,i-1) = p.sigma1;
            M(i,i+1) = p.sigma;
        elseif x == 0 && ((y * p.h < elecEnd) && (y * p.h > elecStart)) %the electrode zone - dirichle condition
            sol_v(i) = p.v0;
            M(i,i) = 1;
        elseif  y == p.N-1 && ((x * p.h < groundEnd) && (x * p.h > groundStart)) %the ground zone- dirichle condition
            sol_v(i) = 0;
            M(i,i) = 1;
        elseif x == 0 %left edge - noyman condition
            M(i,i) = 1;
            M(i,i+p.N) = -1;
        elseif x == p.N-1 %right edge - noyman condition
            M(i,i-p.N) = 1;
            M(i,i) = -1;
        elseif y == 0 %upper edge - noyman condition
            M(i,i+1) = 1;
            M(i,i) = -1;
        elseif y == p.N-1 %lower edge - noyman condition
            M(i,i-1) = 1;
            M(i,i) = -1;
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

function [Z,current] = calcImp2(p,xc,yc,dx,dy)

current = zeros(5);
Z = zeros(5);
elec = 0; 
ground = 0;
for m = 1:5
    for n = 1:5
        curr_pot = setPlatePot2(m,n,p,xc,yc,dx,dy);
        [Jx,Jy] = gradient(-1*curr_pot,p.h,p.h); %actually, it is Ex,Ey
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

function [Z,current] = calcImp3(p,xc,yc,dx,dy,dxx,dyy)

current = zeros(5);
Z = zeros(5);
elec = 0; 
ground = 0;
for m = 1:5
    for n = 1:5
        curr_pot = setPlatePot3(m,n,p,xc,yc,dx,dy,dxx,dyy);
        [Jx,Jy] = gradient(-1*curr_pot,p.h); %actually, it is Ex,Ey
        for i = 1:p.N
            elec = elec + Jx(i,1)*p.h; %sum over the current out of the electrode
            ground = ground + Jy(p.N,i)*p.h; %sum over the current that goes in to the ground
        end
        current(m,n) = abs(elec-ground)/((elec+ground)/2); %deviation from conservation of current
        Z(m,n) = p.v0/((elec+ground)/2); %find the impedance
        elec = 0;
        ground =0;
    end
end
Z = flipud(Z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = calcD(Z1,Z2)

D = 0;
for m=1:5
    for n=1:5
        D = D+(Z1(m,n)-Z2(m,n))^2;
    end
end
D = sqrt(D);
end

