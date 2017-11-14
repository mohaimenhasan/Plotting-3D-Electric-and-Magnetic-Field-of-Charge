function [Etot, Ex, Ey, Ez] = ringofcharge(a, rhos, x, y, z, N)
epsilon = 8.854e-12;
dtheta = 2*pi/N;
thetaPrime= linspace(0, 2*pi, N);
for k=1:length(thetaPrime)
    integrand=dtheta/((x-a*cos(k))^2 + (y-a*sin(k))^2 + z^2)^(3/2);   
    dEx(k) = integrand*(x-a*cos(k));
    dEy(k) = integrand*(y-a*sin(k));
    dEz(k) = integrand*z;
end


Ex = ((rhos*a)/(4*pi*epsilon))*sum(dEx);
Ey = ((rhos*a)/(4*pi*epsilon))*sum(dEy);
Ez = ((rhos*a)/(4*pi*epsilon))*sum(dEz);

Etot = (Ex^2+Ey^2+Ez^2)^0.5;
