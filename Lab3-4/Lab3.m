clear all
close all
clc

% mol_dyn1.m
% mol_dyn27.m - 27febr NEW:  mas0centr=sum(mas2r);
% modified  mo26dyn //: min(r) min(v) ?
% mol dynamics ( nstep=...) mas(1)>=<mas(2)for N = 2 
% functions used:     forc2,  sep2,  
% funct-ns tested in:  tsforc2, ..f=forc2(r); f(1,2)-LJpot 
% Vervet tested at Lx=Ly=0.4;{0.3 0.3}; dt=0.01; vmax=0.35; nstep=25  

N = 2;  % Number of particles
mas = [1.0 1.0];  % Masses of the molecules
%mas(1, N) = 1.0;
Lx = 1;  % Length of the box in x-direction
Ly = 1;  % Length of the box in y-direction
s2xy = Lx * Lx + Ly * Ly;  % Square distance for normalization
L = sqrt(s2xy);  % Total length
vmax = 0.12;  % Maximum velocity
dt = 0.015;  % Time step for integration
dt2 = dt * dt;  % Square of the time step for calculations
dr = 0;  % Change in distance
nstep = 500;  % Total number of steps
rn(1:nstep) = 0;  % Initialize distance array

   % Initial state

% x0 = [-0.3,0.3]; y0 = [-0.3,0.3]; x
v0 = [vmax, -vmax];  % Initial velocities of the molecules
x0 = [-0.1, 0.5];  % Initial x positions
y0 = [-0.1, 0.5];  % Initial y positions
r2xy = [x0 .* x0 + y0 .* y0];  % Square distances
rxy = sqrt(r2xy);  % Distances
mas2r = mas .* rxy;  % Weighted distances by mass
%mas0centr = sum(mas2r);  % Center of mass
acold = [0 0];  % Initial accelerations
xt = x0;  % Current x positions
yt = y0;  % Current y positions
r0 = sep2(xt, yt);  % Initial separation between molecules

xx = -Lx:Lx;  % X-axis for plotting
yy = -Ly:Ly;  % Y-axis for plotting

% Initial state plot
figure()
plot(xx, yy)  % Plot the box
hold on
plot(x0, y0, 'g+')  % Plot initial positions of molecules
title('+ - molecules start position') 
hold off

% f0 = forc1(r0);
f0 = forc2(r0);  % Calculate initial forces
ac(1) = acold(1) + f0(1,1) * dr ./ mas(1);  % Acceleration for molecule 1
ac(2) = acold(2) - f0(1,1) * dr ./ mas(2);  % Acceleration for molecule 2 
             %! minus from 3-rd Newton's Law    
x = x0;  % Set current positions
y = y0;  % Set current positions
vt = v0;  % Set current velocities

for k = 1:nstep
    % Update positions using Verlet integration
    xt = x + vt.*dt + ac.*dt.*dt/2;                      
    yt = y + vt.*dt + ac.*dt.*dt/2;

    % temporary  state xt, yt
    figure()  % Screen ++
    plot(xx,yy)
    hold on
    plot(x0, y0, 'g+')
    plot(xt(1), yt(1), 'bo')
    plot(xt(2), yt(2), 'ro')
    title('Moleculas start +G and current *R positions') 
    hold off 

    step = k;  % Current step number
    acold = ac;  % Save previous accelerations
    x = xt;  % Update positions
    y = yt;  % Update positions
    r = sep2(x, y);  % Calculate new separations
    rn(k) = r;  % Store separation for this step
    dr = r - r0;  % Change in separation

   % Periodic boundary conditions
    if dr > L
        r = r - L;  % Adjust separation
        xt = x - Lx;  % Adjust x position
        yt = y - Ly;  % Adjust y position
    end

    r0 = r;  % Update previous separation
    f = forc2(r);  % Calculate new forces
    ac(1) = acold(1) + f(1,1) * dr ./ mas(1);  % Update acceleration for molecule 1
    ac(2) = acold(2) - f(1,1) * dr ./ mas(2);  % Update acceleration for molecule 2
    
    v = vt;  % Store previous velocities
    vt = v + (acold + ac) .* dt / 2;  % Update velocities using average acceleration
end

% Save distances to a file
save('distances1.mat', 'rn');

mas0centr = sum(mas2r);  % Final center of mass
new_r2xy = [x .* x + y .* y];  % Square distances for final positions
rxynew = sqrt(new_r2xy);  % Final distances
mas2rnew = mas .* rxynew;  % Weighted final distances
mas2centr = sum(mas2rnew);  % Final center of mass

figure()
plot(rn)
title('Distance r(k) [dt = 0.015, nstep = 500] ')
