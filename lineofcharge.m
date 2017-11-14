function [Etot, Ex, Ey, Ez] = lineofcharge(h, rhol, x, y, z, N)
epsilon = 8.854e-12;
dz = 2*h/N;
zprime = linspace(-h, h, N)
for k=1:length(zprime)
    integrand=dz/((x^2 + y^2 + (z-zprime(k))^2)^(3/2));
    dEx(k)=integrand;
    dEy(k)=integrand;
    dEz(k)=(z- zprime(k))*integrand;
end
Ex = ((rhol*x)/(4*pi*epsilon))*sum(dEx);
Ey = ((rhol*x)/(4*pi*epsilon))*sum(dEy);
Ez = (rhol/(4*pi*epsilon))*sum(dEz);
Etot = (Ex^2+Ey^2+Ez^2)^0.5;