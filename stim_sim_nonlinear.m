% this one randomly sets a certain fraction of the weights to zero
% AND makes the first part a simulation without perturbations

%% first: run ratemodel_4
% Q's:  use sqrt or squared to measure weight differences?
% averaging:  simple or need to change to use a gaussian smoothing?
%  main thing for connections:  use all, or only excitatory?

%clear all; clf

% it doesn't actually start here - the code currently does random eye
% positions at each iteration
Estartpos=-20;


duration=2500;
delt=1;
n=100;
t=0:delt:duration;
samplenum=250; % this is the number of timesteps of the linear network simulation to average in order to use for saving r and int

numperturb=20;

%% set from the rate model!
w=Wij(1:n,1:n)';
zeromatrix=(w~=0);
length(reshape(find(zeromatrix==0),1,[]))/n^2
tonic=Ti(1:n);

%%



dub=w; % use "dub" as a variable

nums=10; % number of cycles
initial=1; % number of initial ones where you don't perturb


%init=zeros(n,M);
r_na=[];int_na=[];
rr0_na=[]; int0_na=[];

u=[]; u0=[];
s=[]; s0=[];
D=[]; Dc=[]; D2=[]; Dc2=[]; D0=[]; Dc0=[]; D20=[]; Dc20=[];
init=[];
beta=zeros(n,n);
betac=zeros(n,n);
beta0=zeros(n,n);
betac0=zeros(n,n);
for k=1:nums
    if k<=initial
        init(:,k)=HalfWaveRectify((r0(1:n)+(Estartpos+40*rand)*slopes(1:n)));
        inits=sinfE(init(:,k));
        [R,INT,U,S]=nonlinear_network_simulation(n,dub,tonic,slopes,init(:,k),init(:,k),inits,t); 
        %[R,INT,U,S]=nonlinear_network_simulation(n,dub,tonic,slopes,R(end,:)',U(end,:)',S(end,:)',0:delt:(duration*initial));
        R0=R; INT0=INT; U0=U; S0=S;
    else
        % simulation without peturbations: change eye position
        init(:,k)=HalfWaveRectify((r0(1:n)+(Estartpos+40*rand)*slopes(1:n)));
        inits=sinfE(init(:,k));
        [R0,INT0,U0,S0]=nonlinear_network_simulation(n,dub,tonic,slopes,init(:,k),init(:,k),inits,t); 
        %[R0,INT0,U0,S0]=nonlinear_network_simulation(n,dub,tonic,slopes,rr0_na(end,:)',u0(end,:)',s0(end,:)',t);
        % simulation with perturbations
        toperturb=ceil(rand(1,numperturb)*n);
        %init(:,k)=r_na(end,:)';
        %inits=s(end,:)';
        %perturbation
        init(toperturb,k)=zeros(1,numperturb);
        inits(toperturb)=zeros(1,numperturb);
        [R,INT,U,S]=nonlinear_network_simulation(n,dub,tonic,slopes,init(:,k),init(:,k),inits,t);
    end
    r_na=[r_na;R];    int_na=[int_na;INT];    u=[u;U];  s=[s;S];
    rr0_na=[rr0_na;R0];   int0_na=[int0_na;INT0];   u0=[u0;U0]; s0=[s0;S0];
    
    r_na=R;
    int_na=INT;
    rr0_na=R0;
    int0_na=INT0;
    
    
    %%%%%%%%%%%%%%% do averaging so that only part of the rate data is used:
    if samplenum>1
         tempi=[];
         tempr=[];
         tempi0=[];
         tempr0=[];
         q=samplenum;
         for i=1:(floor(size(r_na,1)/q))
             tempi=[tempi; mean(int_na(((i-1)*q+1):(q*i),:))]; %Important!  to subtract tonic current
             tempr=[tempr; mean(r_na(((i-1)*q+1):(q*i),:))];
             tempi0=[tempi0; mean(int0_na(((i-1)*q+1):(q*i),:))];
             tempr0=[tempr0; mean(rr0_na(((i-1)*q+1):(q*i),:))];
         end
        r=tempr;
        rr0=tempr0;
        int=FInverse(tempi)-repmat(tonic',size(tempi,1),1);
        int0=FInverse(tempi0)-repmat(tonic',size(tempi,1),1);
    else
        r=r_na;
        rr0=rr0_na;
        int=FInverse(int_na)-repmat(tonic',size(int_na,1),1);
        int0=FInverse(int0_na)-repmat(tonic',size(int0_na,1),1);
    end
    
    % linear regression ish fit;
     mmax=max(max(abs(w)));
    parfor i=1:n
       % betac linear regression
        Amat=zeros(n);
        bvector=ones(n,1);
        Aeq=1.0*diag(~zeromatrix(:,i));
        beq=zeros(n,1);
        options = optimset('LargeScale','off');
        if n==25
            lb=zeros(n,1);
            ub=mmax*ones(n,1);
        else
            lb=-mmax*(w(:,i)<0);
            ub=mmax*(w(:,i)>0);
        end
        f=1/(initial+k-1);
        betac(:,i)=(1-f)*betac(:,i)+f*lsqlin(sinfE(r),int(:,i),Amat,bvector,Aeq,beq,lb,ub);
        betac0(:,i)=(1-f)*betac0(:,i)+f*lsqlin(sinfE(rr0),int0(:,i),Amat,bvector,Aeq,beq,lb,ub);
        % beta linear regression
        Aeq=zeros(n);
        lb=-mmax*ones(n,1);
        ub=mmax*ones(n,1);
        beta(:,i)=(1-f)*beta(:,i)+f*lsqlin(sinfE(r),int(:,i),Amat,bvector,Aeq,beq,lb,ub);
        beta0(:,i)=(1-f)*beta0(:,i)+f*lsqlin(sinfE(rr0),int0(:,i),Amat,bvector,Aeq,beq,lb,ub);
        
        
        tosolve=find(zeromatrix(:,i)>0);
%         temp=zeros(n,1);
%         temp(tosolve)=pinv(r(:,tosolve))*(int(:,i));
%         betac(:,i)=temp;    
%         temp=zeros(n,1);
%         temp(tosolve)=pinv(rr0(:,tosolve))*(int0(:,i));
%         betac0(:,i)=temp;    
    end
%     beta0=pinv(r)*int;  %connectivity derived from activity
%     beta0=pinv(rr0)*int0;  %connectivity derived from activity
    a=mean(mean(dub.^2));
    d=((dub-beta).^2/a);
    dc=((dub-betac).^2/a);    
    d0=((dub-beta0).^2/a);
    dc0=((dub-betac0).^2/a);

    
    % (find(zeromatrix>0))
    sqrtmean=@(x) sqrt(mean(reshape(x,1,[])));
    sqrtmean2=@(x) sqrt(mean(reshape(x(find(zeromatrix>0)),1,[])));
    D(k)=sqrtmean(d);  %real vs derived difference
    Dc(k)=sqrtmean(dc);
    D2(k)=sqrtmean2(d);
    Dc2(k)=sqrtmean2(dc);
    D0(k)=sqrtmean(d0);
    Dc0(k)=sqrtmean(dc0);
    D20(k)=sqrtmean2(d0);
    Dc20(k)=sqrtmean2(dc0);


     figure(1); 
     clf(1)
     colormap redbluecmap
     lim=max(max(abs(w)))*0.8;
    subplot(1,5,1);
     imagesc((beta0),[-lim lim]);
     title('W*')
     subplot(1,5,2);
     imagesc((beta),[-lim lim]);
     title('W*p')
     subplot(1,5,3);
     imagesc((betac0),[-lim lim]);
     title('W*c')
	 subplot(1,5,4);
     imagesc((betac),[-lim lim]);
     title('W*cp')
	 subplot(1,5,5);
     imagesc((dub),[-lim lim]);
     title('W (actual network connectivity)')
     if k==1
         pause(0.5)
     end
     figure(2);clf(2)
     subplot(1,2,1)
     plot(r)
     subplot(1,2,2)
     plot(rr0)
end

% visualize the vectors of activity after each random perturbation
% x:  neuron number
% y:  activity measure
figure(2);clf(2)
plot(r)

figure(3)
clf(3)
% this plots the decrease of D, the median of the difference between the "true" and "estimated" connectivity matrices, as a function of the number of perturbations.  It probably decreases very sharply at first and then not so much around k=10, because that is where Emre said the fit starts to level off
%axis square
%subplot(1,2,1);
hold on
plot(0:nums-1,D0,'kd-','linewidth',2,'markersize',7);
plot(0:nums-1,D,'ko-','linewidth',2,'markersize',7);
plot(0:nums-1,Dc0,'kx-','linewidth',2,'markersize',7);
plot(0:nums-1,Dc,'k.-','linewidth',2,'markersize',20);
%set(gca,'box','off');axis square
title('Connection weight matrix difference')
xlabel('Number of perturbations (or time)')
ylabel('Average difference from true weights')
legend('Unconstrained, regular activity','Unconstrained, perturbations','Zeros known, regular activity','Zeros known, perturbations')
figure(4)
clf(4)
%subplot(1,2,2);
hold on
plot(0:nums-1,D20,'kd-','linewidth',2,'markersize',7);
plot(0:nums-1,D2,'ko-','linewidth',2,'markersize',7);
plot(0:nums-1,Dc20,'kx-','linewidth',2,'markersize',7);
plot(0:nums-1,Dc2,'k.-','linewidth',2,'markersize',20);
%set(gca,'box','off');axis square
title('Connection weight matrix difference:  only nonzero weights')
legend('Unconstrained, regular activity','Unconstrained, perturbations','Zeros known, regular activity','Zeros known, perturbations')
ylabel('Average difference from true weights')
xlabel('Number of perturbations (or time)')
hold off
D
Dc

