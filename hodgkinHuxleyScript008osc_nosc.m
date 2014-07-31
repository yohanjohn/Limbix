clear all; close all;
tic;
init1 = [0 .5 .5 .5];
N = 60; % Number of neurons

Ng = N/3;

%initialvalues = [zeros(1,N) 0.5*zeros(1,N) 0.5*zeros(1,N) 0.5*zeros(1,N) zeros(1,N) zeros(1,N) zeros(1,N) zeros(1,N)];
initialvalues = [-55*ones(1,N) zeros(1,N) zeros(1,N) zeros(1,N) ...
    zeros(1,N) zeros(1,N) zeros(1,N) zeros(1,N) zeros(1,N) zeros(1,N) zeros(1,N)];
durat = 2000; tlen = durat*500;
ll = linspace(0,durat,tlen);
h = ll(2);
z = initialvalues';

y = zeros(tlen,N);
Inhs = y; Excs = y;


region1 = 1:round(N/3);
region2 = round(N/3)+1:round(2*N/3);
region3 = round(2*N/3)+1:N;

tInp = linspace(0,durat,tlen);

Inp = zeros(tlen,N);
randInp = ll(2)*150/1000 > rand(tlen,N);
randInp = double(randInp);


for ii = 1:N
    randInp(:,ii) =conv(randInp(:,ii),ones(round(.005*tlen),1),'same'); %Inp = f(tInp)
end
randInp = randInp./max(max(randInp));
% 
% for tt = 1:tlen
%     Inp(tt,:) =2*((1+sin(10*pi*(1:N)/N - 0.050*pi*tInp(tt))).*max(sin(0.01*pi*tInp(tt)' +(1:N)/3),0));
% end

Inp(tlen*0.4:tlen*0.6,region1(3:6)) = 5;
Inp(tlen*0.5:tlen*0.8,region1(10:13)) = 3;
Inp(0.72*tlen:tlen*0.9,region1(17:20)) = 12;

  Inp(tlen*0.1:tlen*0.2,region3) = -4;
  Inp(:,region2) = 0;
%  Inp(:,region1) = 0;


Inp = Inp+randInp;

% Inp(:,region3) = 0;
% Inp(tlen*0.1:tlen*0.2,region3) = -4;

Wex = (zeros(N)); Win = (zeros(N));


Wins = rand(Ng);
Wexs = rand(Ng);
WEI = rand(Ng);

for ii = 1:Ng
    
    Wins(ii,region1) = (1-exp(-(((region1)-region1(ii))/(0.30*Ng)).^2));
    Wins(ii,region1) = Wins(ii,region1)./sum(Wins(ii,region1));
    
    Wexs(ii,:) = 0.1*(exp(-(((region2)-region2(ii))/(0.2*Ng)).^2));
end

WEI = Wexs;
Wins = Wins.*(ones(Ng) - eye(Ng));


Vsyn = -20;

taur = 0.1;
tauf = 2;

pset = [];
pset = ParSet2a(1,N,pset,1:N,'RS');
pset = ParSet2a(0,N,pset,region3,'IN2');

pset.Wex = Wex;
pset.Win = Win;
% pset.taur = taur;
% pset.tauf = tauf;
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
uindex = (8*N+1):9*N;
Inhindex = (9*N+1):10*N;
Excindex = (10*N+1):11*N;

y(1,Vindex) = z(Vindex);
toc
%%
WI = 25;
spik = 1; rat = 0;
if spik == 1
    tic
    for tt = 1:tlen
        
%         if tt == tlen*0.3
%             WI = 1;
%         end
%         
%         if tt == tlen*0.6
%             WI = 0;
%         end

        pset.gI(region2) = WI*Wins*z(Iindex(region3));
        %pset.gI(region3) = 0.02*WI*Wins*z(Iindex(region3));
        
        pset.gE(region3) = 0.1*WEI*z(Iindex(region2));% + 0.5*z(Iindex(region1));
        
        pset.gE(region2) = 0.1*Wexs*z(Iindex(region2)) + 0.04*z(Iindex(region1));
        %pset.gE(region2) = 0.05*Wexs*z(Iindex(region2)) + 0.005*z(Iindex(region1));
        
        qq = hodgkinHuxley15(ll(tt),z,tInp,Inp(tt,:),pset);
        z = z + h*qq;
        y(tt,:) = z(Vindex);
        Inhs(tt,:) = qq(Inhindex);
        Excs(tt,:) = qq(Excindex);
    end
    toc
end
if rat == 1
    tic
    tlenx = durat*400; tInp2 = linspace(0,durat,tlenx);
    llx = linspace(0,durat,tlenx); hx = llx(2);
    
    x = zeros(N,1); Inpx = x; Inhx = x;
    yy = zeros(N,1); 
    xs = zeros(tlenx,N);
    yys = zeros(tlenx,N);
    
    Inp2 = zeros(tlenx,N);
    for tt = 1:tlenx
        Inp2(tt,:) =2*((1+sin(30*pi*(1:N)/N - 0.010*pi*tInp2(tt))).*max(sin(0.02*pi*tInp2(tt)' +(1:N)/5),0));
    end
    Inp2(:,region3) = 0;
    Inp2(:,region2) = 0;
    
    Inp2(round(tlenx*0.5):tlenx,region2(1):region2(Ng)) = 5;
    Inp2(round(tlenx*0.7):tlenx,region1(1:4)) = 7;
    Inp2(round(tlenx*0.9):tlenx,region1(Ng-4:Ng)) = 8;
    
    
    for tt = 1:tlenx
        
        mx = max(x,0);
        
        Inpx(region1) = 0.10*Wexs*mx(region3);
        Inpx(region3) = 0.10*Wexs*mx(region1);
        Inpx(region2) = 0.01*Wexs*mx(region3);
        Inpx(region2) = 0.01*Wexs*mx(region1);
        
        Inhx(region1) = 2*Wins*mx(region2);
        Inhx(region2) = 2*Wins*mx(region2);
        
        x = x + 0.1*hx.*( (10-x).*(30*Inpx + 8*Inp2(tt,:)') - 50*x - 10*(x+10).*(Inhx +0.7*yy ));
        
        yy = yy + 0.01*hx.*((10-yy).*mx*0.1 - 10*yy);
        
        xs(tt,:) = x;
        yys(tt,:) = yy;
    end
    toc
end
%%
t = 1:tlen;
if spik == 1
    figure
    NFig = 4;
    subplot(NFig,1,1)
    %imagesc(Inp')
    plot(ll,y(:,1),ll,y(:,region2(1)),'r')
    
    subplot(NFig,1,2)
    %imagesc(xs')
    %imagesc(y(:,(4*N+1):5*N)');
    imagesc(heaviside(y(:,Vindex)')); colormap(gray); colormap(flipud(colormap)); %cbfreeze
    
    subplot(NFig,1,3)
    %imagesc(y(:,1:N)');
    imagesc(y(:,Vindex)'); %colormap(gray); colormap(flipud(colormap)); cbfreeze
    %imagesc(Inp')
    subplot(NFig,1,NFig)
    imagesc(Inp'); %colormap(jet)
    
%     figure
%     subplot(311)
%     imagesc(Excs');
%     subplot(312)
%     imagesc(Inhs');
%     subplot(313)
%     imagesc(Excs' + Inhs');
    
% figure
% subplot(211)
% imagesc(Wexs)
% subplot(212)
% imagesc(Wins)
   
end
%%
if rat ==1
    figure
    subplot(311)
    imagesc(max(xs',0))
    subplot(312)
    imagesc(yys')
    subplot(313)
    imagesc(Inp2')
end
