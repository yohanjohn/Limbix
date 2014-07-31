%function dydt = hodgkinHuxley7(t,y,tInp,Inp, Wex, Win, Xtau, Itau,N)
function dydt = hodgkinHuxley15(t,y,tInp,Inp, pset)
%N = 20;      % Number of cells
% Traub parameters

taur = pset.taur;
tauf = pset.tauf;
N = pset.N;

vT = pset.vT;
C   = pset.C;
gNa = pset.gNa;
gK  = pset.gK;
gL  = pset.gL;
vNa = pset.vNa;
vK =  pset.vK;
vL = pset.vL;
Vsyn = pset.Vsyn;

gM   = pset.gM;
tmax = pset.tmax;

gE = pset.gE;
gI = pset.gI;

Vx = 2; gT = 0.4; vCa = 120;
%Inp = interp1(tInp, Inp,t);

% The stimulating current (from the synaptic input). Varying this

Vindex = 1:N;
mindex = (N+1):2*N;
hindex = (2*N+1):3*N;
nindex = (3*N+1):4*N;
Xindex = (4*N+1):5*N;
Iindex = (5*N+1):6*N;
pindex = (6*N+1):7*N;
Gindex = (7*N+1):8*N;
uindex = (8*N+1):9*N;
Inhindex = (9*N+1):10*N;
Excindex = (10*N+1):11*N;
%Eindex = (5*N+1):6*N;


V = y(Vindex);
Spik = heaviside(V-20);

aM = -0.32.*(V-vT-13)./(-1+exp(-(V-vT-13)./4));
bM = 0.28.*(V-vT-40)./(exp((V-vT-40)./5)-1);
aH = 0.128.*exp(-(V-vT-17)./18);
bH = 4./(1+exp(-(V-vT-40)./5));
aN = -0.032.*(V-vT-15)./(-1+exp(-(V-vT-15)./5));
bN = 0.5.*exp(-(V-vT-10)./40);

pinf = 1./(1 + exp(-(V+35)./10));
tp = tmax./(3.3.*exp((V+35)./20) + exp(-(V+35)./20));

sinf = 1./(1+exp(-(V + Vx +57)./6.2));
uinf = 1./(1+exp((V + Vx +81)./4));
%tu = (30.8 + (211.4 + exp((V + Vx + 113.2)./5)))./(3.7.*(1+exp((V+Vx+84)./3.2)));
tu = 30.8 + ((211.4 + exp((V + Vx + 113.2)./5)))./(3.7.*(1+exp((V+Vx+84)./3.2)));

IT = gT*sinf.*sinf.*y(uindex).*(V-vCa);
%I   = Inp' + Wex*yI.*(50-V);
%gE = Wex*yI;

I   = Inp' + gE.*(50-V);

%inhcurr = - (Win*y(Iindex)).*(V - Vsyn);
inhcurr = -gI.*(V - Vsyn);

dydt(Inhindex,1) = inhcurr;
dydt(Excindex,1) = I;

dydt(Gindex,1) = sinf.*sinf.*y(uindex);

dydt(Vindex,1) = (1./C).*(I ...
    + inhcurr ...
    - gNa.*y(mindex).^3.*y(hindex).*(y(Vindex)-vNa) ...
    - gK.*y(nindex).^4.*(y(Vindex)-vK) ...
    - gL.*(y(Vindex)-vL) ...
    - gM.*y(pindex).*(y(Vindex) - vK)  ...
    - IT); 

dydt(mindex,1) = (1-y(mindex)).*aM - y(mindex).*bM;
dydt(hindex,1) = (1-y(hindex)).*aH - y(hindex).*bH;
dydt(nindex,1) = (1-y(nindex)).*aN - y(nindex).*bN;

dydt(pindex,1) = (pinf - y(pindex))./tp;
dydt(uindex,1) = (uinf - y(uindex))./(tu);

dydt(Xindex,1) = (1-y(Xindex)).*Spik - y(Xindex)./taur;
%dydt(Eindex,1) = Etau*(((1-y(Eindex)).*(Ws*Is(:,t) + (Wex*X))) - 2*y(Eindex));
dydt(Iindex,1) = ((tauf+taur)./tauf).*((1-y(Iindex)).*(y(Xindex)).*(2./taur) -y(Iindex)./tauf);