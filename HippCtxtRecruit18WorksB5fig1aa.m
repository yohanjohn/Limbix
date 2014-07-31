% Based on HippConjunctive34b.m
clear all; close all;

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

tic;
init1 = [0 .5 .5 .5];

spik = 0; rat = 1;

Nregions = 4;
Nreg = zeros(Nregions,1);
Nreg(1) = 30;
Nreg(2) = 10;
Nreg(3) = Nreg(1);
Nreg(4) = Nreg(2);
Nreg(5) = Nreg(2); %BLA
Nreg(6) = Nreg(5); %BLA

N = sum(Nreg); % Number of neurons
region1 = 1:Nreg(1);
region2 = Nreg(1)+1:Nreg(1)+Nreg(2);
region1IN = region2(end)+1:region2(end)+Nreg(3);
region2IN = region1IN(end)+1:region1IN(end)+Nreg(4);
region3 = region2IN(end)+1:region2IN(end)+Nreg(5);
region3IN = region3(end)+1:region3(end)+Nreg(5);

Nseconds = 120;
Nmsec = Nseconds*1000; % milliseconds

tlen = Nmsec*500; %number of time steps for HH
ll = linspace(0,Nmsec,tlen);
h = ll(2);


downsamp1 = 200; subsamp = 1;
y = zeros(tlen/downsamp1,N);

Wex = (zeros(N)); Win = (zeros(N));

Wins = ones(Nreg(1)) - 0*eye(Nreg(1));
Wexs = 0.0001*rand(Nreg(1));
W12  = 0.001*rand(Nreg(2),Nreg(1));
W21  = 0.001*rand(Nreg(1),Nreg(2));

W23 = zeros(Nreg(4),1);

toc
if spik == 1
    tic
    initialvalues = [zeros(1,N) 0.5*zeros(1,N) 0.5*zeros(1,N) 0.5*zeros(1,N) zeros(1,N) zeros(1,N) zeros(1,N) zeros(1,N)];
    z = initialvalues';
    tInp = linspace(0,Nmsec,tlen);
    Inp = zeros(tlen,N);
    randInp = ll(2)*50/1000 > rand(tlen,N);
    randInp = double(randInp);
    toc
    tic
    for ii = 1:N
        randInp(:,ii) =conv(randInp(:,ii),ones(round(1000),1),'same'); %Inp = f(tInp)
    end
    randInp = randInp ./max(max(randInp));
    toc
    tic
    Vsyn = -20;
    
    %     taur = 0.1;
    %     tauf = 2;
    %     taur = 0.1;
    %     tauf = 5;
    %
    pset = [];
    pset = ParSet2a(1,N,pset,1:N,'RS');
    pset = ParSet2a(0,N,pset,region3,'IN');
    
    pset.Wex = Wex;
    pset.Win = Win;
    %     pset.taur = taur;
    %     pset.tauf = tauf;
    pset.gE = zeros(N,1);
    pset.gI = zeros(N,1);
    
    pset.N = N;
    
    Vindex = 1:N;
    mindex = (N+1):2*N;
    hindex = (2*N+1):3*N;
    nindex = (3*N+1):4*N;
    Xindex = (4*N+1):5*N;
    Iindex = (5*N+1):6*N;
    pindex = (6*N+1):7*N;
    Gindex = (7*N+1):8*N;
    
    y(1,Vindex) = z(Vindex);
    
    %%
    
    
    
    r1 = 3:13;
    r2 = 17:27;
    r12 = intersect(region1(r1),region1(r2));
    
    Imax = 5;
    fr = 32;
    Inp(:,region1(r1)) = Imax*0.5*(1+square(fr*pi*tInp/Nmsec - pi))'*ones(1,length(r1));
    Inp(:,region1(r2)) = Imax*0.5*(1+square(fr*pi*tInp/Nmsec ))'*ones(1,length(r1));
    
    Inp(:,region2(1)) = 5*(1+square(fr*pi*tInp/Nmsec - pi))';
    Inp(:,region2(2)) = 5*(1+square(fr*pi*tInp/Nmsec ))';
    
    Inp(:,r12) = Imax;
    Inp(tlen*0.5:tlen,:) = 0;
    
    
    
    theta = 0.5*(1+square(5*pi*tInp/1000 + pi));
    DA = 1;
    %Inp = (theta'*ones(1,N)).*(Inp + randInp);
    Inp = (Inp + 3*randInp);
    Inp(1:1000,:) = 0;
    toc
    tic
    DA = 1; phas = 0;
    for tt = 1:tlen
        phas = phas + 1/tlen;
        %         pset.gI(region1) = Wins*z(Iindex(region1));
        %         pset.gI(region2) = ones(Nreg(2))*z(Iindex(region2));
        %         pset.gE(region1) = Wexs*z(Iindex(region1));
        
        pset.gI(region1) = 2*Wins*z(Iindex(region3));
        pset.gI(region2) = ones(Nreg(2))*z(Iindex(region4));
        %
        pset.gE(region1) = 0.7*Wexs*z(Iindex(region1));
        
        pset.gE(region3) = 0.3*eye(Nreg(3))*z(Iindex(region1));
        pset.gE(region4) = 0.3*eye(Nreg(4))*z(Iindex(region2));
        %
        if tt> tlen*0.55
            Inp(tt,region1) = Imax*0.5*(1+square(pi*(1:Nreg(1))/Nreg(1) -15*phas +pi/2,10)');
            DA = 0;
        end
        
        if tt>0.7*tlen && tt<0.71*tlen
            phas = 0;
        end
        
        
        qq = hodgkinHuxley14b1(ll(tt),z,tInp,Inp(tt,:),pset);
        z = z + h*qq;
        
        mx = z(Iindex(region1));
        
        %Wexs  = Wexs + DA*0.05*h.*(ones(Nreg(1),1)*mx'-Wexs).*(mx*ones(1,Nreg(1)));
        
        
        
        y(subsamp,:) = z(Vindex);
        subsamp = subsamp + 1-sign(mod(tt,downsamp1));
    end
    toc
end
%%
if rat == 1
    tic
    secslice = 2000;
    tlenx = Nseconds*secslice; tInp2 = linspace(0,Nseconds,tlenx);
    downsamp = 40; st = 1;
    ttlen = tlenx/downsamp;
    
    llx = linspace(0,1,secslice); hx = llx(2);
    
    x = zeros(N,1); Inpx = x; Inhx = x;
    yy = zeros(N,1);
    
    InhW = zeros(N,1);
    InhW(region1) = 40;    InhW(region2) = 20; InhW(region3) = 20;
    PassDec = ones(N,1);  PassDec(region1) = 1;  PassDec(region2) = 10;
    PassDec(region1IN) = 18;  PassDec(region2IN) = 18; PassDec(region3) = 5;
    
    gate = zeros(N,1); gates = zeros(ttlen,N);
    gfac = zeros(N,1); gfac(region2) = 1;
    %PassDec = PassDec.*(1 + 0.85*(rand(N,1)-0.5));
    
    
    xs = zeros(ttlen,N); yys = zeros(ttlen,N);
    
    Inp2 = zeros(tlenx,N);
    
    XintRate = 0.5*ones(N,1);
    
    %XintRate(region1IN) = 0.5*(1 + 0.5*(rand(Nreg(3),1)-0.5));
    %XintRate(region2IN) = 2.1;
   
    r1 = (3:14) +2;
    r2 = (16:27)-2;
    r12 = intersect(region1(r1),region1(r2));
    fr1 = 32; fr2 = 2;
    
    LR = 3;    fracc = 20;    drct = -1; fracadd = 0.002;
    INtonic = 0;
    INtonic1 =0;
    
    Inmem = zeros(Nreg(1),1);    Inmem(r12) = 1;
        
    randInpx = ll(2)*350/1000 > rand(tlenx,Nreg(2));
    randInpx = double(randInpx);
    for ii = 1:Nreg(2)
        randInpx(:,ii) =conv(randInpx(:,ii),ones(1000,1),'same'); %Inp = f(tInp)
    end
    randInpx = randInpx ./max(max(randInpx));
    
    %Inp2(:,region2) = randInpx;
    
    %    Inp2(:,r12) = 1;
    %Inp2(1:tlenx*0.23,r1) = 1;    Inp2(tlenx*0.27:tlenx*0.5, r2) = 1;
    Inp2(1:round(tlenx*0.120),r1) = 1;    Inp2(round(tlenx*0.130):tlenx*0.25, r2) = 1;
    theta = 0.5*(1+square(0.35*pi*tInp2 + pi));

    DA = 10*theta;
    DA(round(tlenx*0.120):tlenx) = 0;
    thr = 0.1;
    phas = 0;
    
    Wext = zeros(ttlen, length(Wexs(:)));    W12t = zeros(ttlen, length(W12(:)));    W21t = zeros(ttlen, length(W21(:)));
    
    conflict = zeros(ttlen,1);
    %conflict(tlenx*0.25:tlenx*0.3) = 1;
    conf = 0; ACh = 0.1;
    AChs = zeros(ttlen,1);
    
    Lelig = zeros(Nreg(2),1);
    Leligs = zeros(ttlen,Nreg(2));
    reg2ton = zeros(Nreg(2),1);
    eta = 2;
    sigpE = 20; sigpI = 20;
    sfun = @(ecks) sigmerf(ecks,1.1,5);
    %sfun = @(ecks) sigmerf(ecks,0.8,10);
    
    pattcomp = 40;
    %pattcomp = zeros(Nreg(2),1);
    %pattcomp(Nreg(2)*0.5+1:Nreg(2),1) = 20;
        
    
    Inp2(tlenx*0.75:tlenx,r2(end):region1(end)) = 0;
    boost = 0;
    
    %Lrate = linspace(0.01,5,Nreg(2))';
    amyf = 2;
    
    Rboost = 1;
    Lboost = 1;
    patt = r12;
    
    for tt = 1:tlenx
        phas = phas + 1/tlenx;
        mx = max(x,0);
% 
%         if tt>0.124*tlenx && tt <0.126*tlenx
%             x(:) = 0;
%         end
        
        if tt> 0.25*tlenx && tt<0.27*tlenx
            %x(:) = 0;
            phas = 0;
            %x(:) = 0;
        end
        
        if tt> tlenx*0.27
            Inp2(tt,region1) = 0.5*(1+square(pi*(1:Nreg(1))/Nreg(1) +drct*12*phas ,fracc)');
             drct = -1;
         end
        
   
        if tt>tlenx*0.5 && tt<tlenx*0.51
             Inp2(tt,region1) = 0;
             %boost = 0.2;
             x(:) = 0;
        end
 
        if tt>tlenx*0.51 && tt< tlenx *0.72
             Inp2(tt,region1) = 0;
             Inp2(tt,patt) = 1;
        end
        
        if tt>tlenx*0.70 && tt< tlenx *0.73
             Inp2(tt,region1) = 0;
             phas = 0;
        end
        
        
        if tt>tlenx*0.73
            Inp2(tt,region1) = 0.5*(1+square(pi*(1:Nreg(1))/Nreg(1) +drct*12*phas ,fracc)');
            drct = -1;
        end
        
        if tt>tlenx*0.88 && tt< tlenx *0.9
             Inp2(tt,region1) = 0;
        end
        
        
        if tt>tlenx*0.9
             Inp2(tt,region1) = 0;
             Inp2(tt,patt) = 1;
        end
        
        
        compfac = (max(0.4-max(mx(region2)),0));
        
        %conf = conf + 0.01*hx.*((1-conf).*(compfac*5 + 0*sum(mx(region3))) - 5*conf - 10*(conf+1).*max(mx(region2)));
        conf = conf + 0.05*hx.*((1-conf).*(compfac*10)  - 15*conf);
        conf = max(conf,0);
        conflict(st) = conf;
        
        ACh = ACh + 0.25*hx.*((1-ACh).*max(max(mx(region2)-0.4,0))*10  - 15*ACh);
        AChs(st) = ACh;
        
        %Inpx(region1) = 8.*Inp2(tt,region1)' + max(conf,0).*pattcomp*Wexs*mx(region1);
        Inpx(region1) = 8.*Inp2(tt,region1)' + max(0.1-ACh,0).*pattcomp*Wexs*mx(region1);
        Inpx(region1IN) = 0.4*Wins*mx(region1) +INtonic1;
 
        Inpx(region2) = 30*sfun(W12*mx(region1)) + 10*Inp2(tt,region2)' + Rboost*randInpx(tt,:)' + amyf*mx(region3) +0*sigmerf(mx(region2),0.4,sigpE);
        Inpx(region2IN) = 20*(ones(Nreg(2)) -eye(Nreg(2)))*sigmerf(mx(region2),0.4,sigpI) +INtonic; 
        
%         Inpx(region2) = 5*W12*mx(region1) + 10*Inp2(tt,region2)' + 0*mx(region3) +0*sigmerf(mx(region2),0.4,sigpE);
%         Inpx(region2IN) = 5*(ones(Nreg(2)))*mx(region2); 
        
        
        Inpx(region3) =50*sigmerf(W23.*mx(region2),0.25,20); 
        Inpx(region3IN) = 0.4*ones(Nreg(6))*mx(region3);
        
        yy(region1) = mx(region1IN);
        yy(region2) = mx(region2IN);
        yy(region3) = mx(region3IN);
        
        
        x = x + XintRate.*hx.*( (1-x).*(Inpx) - PassDec.*x - (x+1).*(InhW.*yy + 0*gate.*gfac ));
        
        Lelig = Lelig + 0.25*hx.*((1-Lelig).*sigmerf(mx(region2),0.4,20));
        %Lelig = Lelig + 0.15*hx.*((mx(region2)-Lelig).*sigmerf(mx(region2),0.4,20));
        
        %Wexs  = Wexs + 0.01*LR*hx.*(ones(Nreg(1),1)*mx(region1)'-Wexs).*(mx(region1)*ones(1,Nreg(1)));
        Wexs  = Wexs + max(ACh - 0.1,0)*LR*hx.*(ones(Nreg(1),1)*mx(region1)'-Wexs).*(mx(region1)*ones(1,Nreg(1)));
        
        W12   = W12  +Lboost*LR*hx.*((ones(Nreg(2),1)*mx(region1)')-W12).*((max(1-Lelig,0).*mx(region2))*ones(1,Nreg(1)));
       %W12   = W12  + LR*hx.*((ones(Nreg(2),1)*mx(region1)')-W12).*DA(tt).*(max(mx(region2)-Lelig,0)*ones(1,Nreg(1)));
        
        W21 = W21 + LR*hx.*(mx(region1)*ones(1,Nreg(2)) - W21).*(ones(Nreg(1),1)*(mx(region2))'); % Outstar with presynaptic gating
        
       
        %W23   = W23  +10*LR*hx.*(mx(region2)-W23).*mx(region3);
        W23 = W23 + DA(tt).*20.*hx.*((1-W23).*max(mx(region2) -0.35,0));
        
%         W12 = W12+ DA(tt).*LR.*hx.*((1-W12).*(mx(region2)*mx(region1)')*5 - 0.1*(W12+1).*(mx(region2)*max(thr-mx(region1),0)' + max(thr - mx(region2),0)*mx(region1)'));
%         W21 = W21+ DA(tt).*LR*hx.*((1-W21).*(mx(region1)*mx(region2)')*5 -  0.1*(W21+1).*(mx(region1)*max(thr-mx(region2),0)' + max(thr - mx(region1),0)*mx(region2)'));
%         
%         Wexs = Wexs + DA(tt).*LR*hx.*((1-Wexs).*(mx(region1)*mx(region1)') -  0.5*(Wexs+1).*(mx(region1)*max(thr-mx(region1),0)' + max(thr - mx(region1),0)*mx(region1)'));
%         
%         Wexs = max(Wexs,0);
%         W12 = max(W12,0);
%         W21 = max(W21,0);
        
        Wext(st,:) = Wexs(:);
        W12t(st,:) = W12(:);
        W21t(st,:) = W21(:);
        
        xs(st,:) = x;
        Leligs(st,:) = Lelig;
        
        yys(st,region2) =W12*mx(region1);
        
        st = st + 1-sign(mod(tt,downsamp));
    end
    toc
end
%%
tic
t = 1:tlen;

if spik == 1
    figure
    NFig = 4;
    subplot(NFig,1,1)
    %imagesc(Inp')
    ld = downsample(ll,downsamp1);
    plot(ld,y(:,1),ld,y(:,region2(1)),'r')
    
    subplot(NFig,1,2)
    imagesc(y')
    %imagesc(y(:,(4*N+1):5*N)');
    %imagesc(heaviside(y(:,Vindex)')); colormap(gray); colormap(flipud(colormap)); %cbfreeze
    
    subplot(NFig,1,3)
    %imagesc(y(:,1:N)');
    imagesc(y(:,Vindex(region1))'); colormap(bone); %colormap(flipud(colormap)); cbfreeze
    %imagesc(Inp')
    subplot(NFig,1,NFig)
    imagesc(Inp'); %colormap(jet)
end
%%
if rat ==1
    figure('Position',[ 1100 100 700 850])
    Nfig = 5;
    subplot(Nfig,1,1)
    imagesc(Inp2(:,region1)') 
    ylabel('Neuron #')
    title('Input')
    subplot(Nfig,1,2)
    imagesc(max(xs(:,region1)',0))
    title('CA3 activity')
    ylabel('Neuron #')
    subplot(Nfig,1,4)
    imagesc(max(xs(:,region2)',0))
    title('CA1 activity')
    ylabel('Neuron #')
    
    
     subplot(Nfig,1,3)
     %plot(1:ttlen,conflict,1:ttlen,AChs)
     plot(1:ttlen,AChs,'LineWidth',2)
     axis([0 ttlen 0 1.1*max(AChs)])
     title('Modulation of CA3 recurrence strength (ACh)')
     ylabel('Activity')
     subplot(Nfig,1,Nfig)
     %plot(Leligs)
     %axis([0 ttlen 0 1.1*max(max(Leligs))])
     imagesc(max(xs(:,region3)',0))
     title('BLA activity')
     xlabel('Time (arbitrary units)')
     ylabel('Neuron #')
      print -depsc Fig1.eps
     
     figure('Position',[ 1300 500 500 500])
    subplot(221)
    imagesc(Wexs)
    subplot(222)
    imagesc(W12)
    subplot(223)
    imagesc(W21)
    
%     figure
%     subplot(311)
%     plot(Wext)
%     subplot(312)
%     plot(W12t)
%     subplot(313)
%     plot(W21t)
%     
    
    
end
toc