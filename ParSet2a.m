function pset = ParSet2a(initflag,N,pset,indx,neurontype)
Vsyn = -80;
switch neurontype
    case 'other'
        %len = 61.4e-4; diam = len; %Units: cm
        %Cm = pi*len*diam*1e3; %Units: nF
        C   = 1;
        gNa = 50;
        gK  = 5;
        gL  =0.01;
        vNa = 50;
        vK = - 90;
        vL = - 70.3;
        vT = -56.2; % mV
        gM = 0.02;
        tmax =  608;
        taur = 0.1;
        tauf = 2;
    case 'gLtest'
        %len = 61.4e-4; diam = len; %Units: cm
        %Cm = pi*len*diam*1e3; %Units: nF
        C   = 1;
        gNa = 50;
        gK  = 5;
        gL  =0.05;
        vNa = 50;
        vK = - 90;
        vL = - 70.3;
        vT = -56.2; % mV
        gM = 0.02;
        tmax =  608;
    
    case 'default'
        %len = 61.4e-4; diam = len; %Units: cm
        %Cm = pi*len*diam*1e3; %Units: nF
        C   = 1;
        gNa = 50;
        gK  = 5;
        gL  =0.01;
        vNa = 50;
        vK = - 90;
        vL = - 70.3;
        vT = -56.2; % mV
        gM = 0;
        tmax =  608;
        taur = 0.1;
        tauf = 2;
   
    
    case 'RS'
        len = 61.4e-4; diam = len; %Units: cm
        Cm = pi*len*diam*1e3; %Units: nF
        C   = Cm;
        gNa = 56;
        gK  = 6;
        gL  =0.0205;
        vNa = 50;
        vK = - 100;
        vL = - 70.3;
        vT = -56.2; % mV
        gM = 0.075;
        tmax =  608;
        taur = 0.1;
        tauf = 2;
   
    case 'RSIN'
        len = 61.8e-4; diam = len; %Units: cm
        Cm = pi*len*diam*1e3; %Units: nF
        C   = Cm;
        gNa = 21;
        gK  = 10;
        gL  =0.133;
        vNa = 50;
        vK = - 100;
        vL = - 56.2;
        vT = -56.2; % mV
        gM = 0.075;
        tmax =  608;
        taur = 0.1;
        tauf = 2;
   
        
    case 'FSIN'
        len = 67e-4; diam = len; %Units: cm
        Cm = pi*len*diam*1e3; %Units: nF
        C   = Cm;
        gNa = 56;
        gK  = 10;
        gL  =0.15;
        vNa = 50;
        vK = - 100;
        vL = - 70.3;
        vT = -56.2; % mV
        gM = 0.0;
        tmax =  608;
        taur = 0.1;
        tauf = 2;
       
    case 'IN'
        len = 67e-4; diam = len; %Units: cm
        Cm = pi*len*diam*1e3; %Units: nF
        C   = Cm;
        gNa = 56;
        gK  = 10;
        gL  =0.15;
        vNa = 50;
        vK = - 100;
        vL = - 70.3;
        vT = -56.2; % mV
        gM = 0.0;
        tmax =  608;
        taur = 0.2;
        tauf = 4;
        
     case 'IN2'
        len = 67e-4; diam = len; %Units: cm
        Cm = pi*len*diam*1e3; %Units: nF
        C   = Cm;
        gNa = 56;
        gK  = 10;
        gL  =0.15;
        vNa = 50;
        vK = - 100;
        vL = - 70.3;
        vT = -56.2; % mV
        gM = 0.0;
        tmax =  608;
        taur = 0.2;
        tauf = 8;
   
   
end


if initflag == 1
    pset.vT   = vT.*ones(N,1);
    pset.C    = C.*ones(N,1);
    pset.gNa  = gNa.*ones(N,1);
    pset.gK   = gK.*ones(N,1);
    pset.gL   = gL.*ones(N,1);
    pset.vNa  = vNa.*ones(N,1);
    pset.vK   = vK.*ones(N,1);
    pset.vL   = vL.*ones(N,1);
    pset.Vsyn = Vsyn.*ones(N,1);
    pset.gM   = gM.*ones(N,1);
    pset.tmax = tmax.*ones(N,1);
    pset.taur = taur.*ones(N,1);
    pset.tauf = tauf.*ones(N,1);
else
    
    Ng = length(indx);
    
    pset.vT(indx)   = vT.*ones(Ng,1);
    pset.C(indx)    = C.*ones(Ng,1);
    pset.gNa(indx)  = gNa.*ones(Ng,1);
    pset.gK(indx)   = gK.*ones(Ng,1);
    pset.gL(indx)   = gL.*ones(Ng,1);
    pset.vNa(indx)  = vNa.*ones(Ng,1);
    pset.vK(indx)   = vK.*ones(Ng,1);
    pset.vL(indx)   = vL.*ones(Ng,1);
    pset.Vsyn(indx) = Vsyn.*ones(Ng,1);
    pset.gM(indx)   = gM.*ones(Ng,1);
    pset.tmax(indx) = tmax.*ones(Ng,1);
    pset.taur = taur.*ones(N,1);
    pset.tauf = tauf.*ones(N,1);
end
end