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
      
