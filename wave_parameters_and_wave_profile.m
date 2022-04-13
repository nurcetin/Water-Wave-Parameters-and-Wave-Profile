clear
T = input(' Enter wave period T (sec)');
H = input(' Enter wave height H (m)');
% d : the distance from the still water level (SWL) to sea bottom
d = input(' Enter water depth d (m)');
g = input(' Enter acceleration of gravity g (m/s^2)');
rho = input(' Enter water density rho (kg/m^3)');
% z is taken as positive for above SWL and taken as negative for below SWL
z = input(' Enter z coordinate (m)');
% x axes is the SWL axes, x=0 is at crest, x=L/2 is at through, x=L is at crest
t = input(' Enter time t (sec)');
purpose = input('Do you want to wave profile or wave parameters?','s');
% L0 : deep water wave length (m)
L0=g/2/pi*T*T;
% C0 : deep water celerity (m/s)
C0=g/2/pi*T;
% sgm : angular wave frequency (rad/s)
sgm=2*pi/T;
% s : wave steepness
s=d/L0;
% a : wave amplitude (m)
a=H/2;
if s>=0.5
    % L : wave length at the depth taken into consideration
    L=L0;
    % C : wave celerity (m/s)
    C=C0;
    % k : wave number (rad/m)
    k=2*pi/L;
    n=0.5;
    % A : particle orbital major (horizontal) semi axis (m)
    A = a*exp(k*z);
    % B : particle orbital minor (vertical) semi axis (m)
    B = a*exp(k*z);
    if isequal(purpose, 'wave profile')
        x=input('Enter the x range (m), ex. 0:L');
    else 
    x = input(' Enter x coordinate (m), it can be in terms of L');
    % u : horizontal particle velocity (m/s)
    u=a*sgm*exp(k*z)*cos(k*x-sgm*t);
    % w : vertical particle velocity (m/s)
    w=a*sgm*exp(k*z)*sin(k*x-sgm*t);
    end
elseif (s<0.5) && (s>0.0157)
    % L : wave length at the depth taken into consideration
    syms L
    eqn=g/2/pi*T*T*tanh(2*pi*d/L)==L;
    L = solve(eqn,L);
    L=double(-L);
    % k : wave number (rad/m)
    k=2*pi/L;
    n=0.5*(1+2*k*d/sinh(2*k*d));
    % A : particle orbital major (horizontal) semi axis (m)
    A = a*cosh(k*(z+d))/sinh(k*d);
    % B : particle orbital minor (vertical) semi axis (m)
    B = a*sinh(k*(z+d))/sinh(k*d);
    % C : wave celerity (m/s)
    C=g/2/pi*T*tanh(2*pi*d/L);
    if isequal(purpose, 'wave profile')
      x=input('Enter the x range (m), ex. 0:L');  
    else
    x = input(' Enter x coordinate (m), it can be in terms of L');
    % u : horizontal particle velocity (m/s)
    u=a*sgm*cosh(k*(z+d))/sinh(k*d)*cos(k*x-sgm*t);
    % w : vertical particle velocity (m/s)
    w=a*sgm*sinh(k*(z+d))/sinh(k*d)*sin(k*x-sgm*t);
    end
else
    % L : wave length at the depth taken into consideration
    L=T*sqrt(g*d);
    % k : wave number (rad/m)
    k=2*pi/L;
    % C : wave celerity (m/s)
    C=sqrt(g*d);
    n=1;
    % A : particle orbital major (horizontal) semi axis (m)
    A = a/k/d;
    % B : particle orbital minor (vertical) semi axis (m)
    B = a*(1+z/d);
    if isequal(purpose, 'wave profile')
        x=input('Enter the x range (m), ex. 0:L');
    else
    x = input(' Enter x coordinate (m), it can be in terms of L');
    % u : horizontal particle velocity (m/s)
    u=a*sgm/k/d*cos(k*x-sgm*t);
    % w : vertical particle velocity (m/s)
    w=a*sgm*(1+z/d)*sin(k*x-sgm*t);
    end
end
% Cg : group velocity (m/s)
Cg=n*C;
% Et = total wave energy obtained for one wave length joule/m
Et=rho*g*H*H*L/8;
% Ea = energy density (specific energy) per unit surface area joule/m^2
Ea=rho*g*H*H/8;
% Pd : dynamic presssure (Pa)
Pd=rho*g*cosh(k*(z+d))/cosh(k*d)*a*cos(k*x-sgm*t);
% Ph : hydrostatic pressure (Pa)
Ph=-rho*g*z;
% P = power for a unit crest width watt/m
P=rho*g*H*H/8*Cg;
% eta(x,t) : instantaneous vertical displacement of water surface above SWL
% (m)
eta=a*cos(k*x-sgm*t);
% Pt : total pressure (Pa)
Pt=Pd+Ph;
if isequal(purpose, 'wave profile')
    plot(x,eta)
else
value =[L0;C0;a;s;L;C;sgm;k;u;w;n;Cg;Et;Ea;P;eta;Pd;Ph;Pt;A;B];
variablenames = {'L0 (m)';'C0 (m/s)';'a (m)';'s';'L (m)';'C (m/s)';'sgm (rad/s)';'k (rad/m)';'u (m/s)';'w(m/s)';'n';'Cg (m/s)';'Et (J/m)';'Ea (J/m^2)';'P (Watt/m)';'eta (m)';'Pd (Pa)';'Ph (Pa)';'Pt (Pa)';'A (m)';'B (m)'};
results=table(value,'RowNames',variablenames)
end
