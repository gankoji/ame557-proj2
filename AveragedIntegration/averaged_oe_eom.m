function xDot = averaged_oe_eom(~,x,mu,BC,J2,R,i0)

a = x(1);
e = x(2);

%Error Handling
if e < 0 
    xDot = zeros(4,1);
    return
end

[rho_a,~] = atmosphere(a - R);

rho_a = rho_a*1000^3;   %kg/km^3

H = a*(1-e) - R;

n = sqrt(mu/(a^3));
nu = a*e/H;
p = a*(1-e^2);

aDot = -(1/BC)*n*a^2*rho_a*(besseli(0,nu)+2*e*besseli(1,nu)+0.75*e^2*(besseli(0,nu)+besseli(2,nu)));
eDot = -(1/BC)*n*a*(1-e^2)*rho_a*(besseli(1,nu)+0.5*e*(besseli(0,nu)+besseli(2,nu))+0.125*e^2*(3*besseli(1,nu)+besseli(3,nu)));
OmegaDot = -3*J2*R^2*n/(2*p^2)*cosd(i0);
omegaDot = 3*J2*R^2*n/(4*p^2)*(5*cosd(i0)^2-1);

xDot = [aDot; eDot; OmegaDot; omegaDot];