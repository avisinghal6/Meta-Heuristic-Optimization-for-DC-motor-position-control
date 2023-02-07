rng('shuffle');
clear;
clc;
Max_iter=2;            % maximum generations
N=20;                   %BAT numbers
dim=5;
lb=[18 0 0 0.001 0.001];
ub=[20 1 1 0.98 0.98];
Fmax=1;                 %maximum frequency
Fmin=0;                 %minimum frequency
A=rand(N,1);            %loudness for each BAT
r=rand(N,1);            %pulse emission rate for each BAT
alpha=0.7;              %constant for loudness update
gamma=0.7;              %constant for emission rate update
ro=0.002; 
wmax=0.3;
wmin=0.1;%initial pulse emission rate
% Initializing arrays
F=zeros(N,1);           % Frequency
v=zeros(N,dim);           % Velocities
 for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        x(:,i)=rand(N,1).*(ub_i-lb_i)+lb_i;
    end
   
for ii=1:N
    Kp=x(ii,1);
    Ki=x(ii,2);
    Kd=x(ii,3);
    mun=x(ii,4);
    lam=x(ii,5);
    sim('fpid');                                                     % proses simülasyonu
    load('error1.mat');                                                      % prosesten hata deðerleri çekiliyor
    load('error.mat');                                                      % prosesten hata deðerleri çekiliyor
    e = ITAE(2,:);                                                          % error vector
    %SSE = (0.75*ST)+(0.18*e(end))+(0.07*PO); 
    SSE=e(end);
    fitness(ii)=SSE;
end
[fmin,index]=min(fitness);          %find the initial best fitness value,
bestsol=x(index,:);                 %find the initial best solution for best fitness value
%%

iter=1;             % start the loop counter
while iter<=Max_iter 
          iter
          fmin
%            w=wmax-(wmax-wmin)*iter/Max_iter;
          %start the loop for iterations
    for ii=1:N
        
        F(ii)=Fmin+(Fmax-Fmin)*rand;   %randomly chose the frequency
        
        v(ii,:)=v(ii,:)+(x(ii,:)-bestsol)*F(ii);  %update the velocity
        x(ii,:)=x(ii,:)+v(ii,:);    
         %update the BAT position
        %         x(ii,:)=round(x(ii,:));
        % Apply simple bounds/limit
        Flag4up=x(ii,:)>ub;
        Flag4low=x(ii,:)<lb;
        x(ii,:)=(x(ii,:).*(~(Flag4up+Flag4low)))+ub.*Flag4up+lb.*Flag4low;
        
        %check the condition with r
        if rand>r(ii)
            % The factor 0.001 limits the step sizes of random walks
            %               x(ii,:)=bestsol+0.001*randn(1,dim);
             epsi=-1+(1-(-1))*rand;
             x(ii,:)=bestsol+epsi*mean(A);
               xnew=bestsol+epsi*mean(A);
%            xnew=bestsol+0.01*randn(1,dim);
            
        else 
            for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        xnew(1,i)=rand(1,1).*(ub_i-lb_i)+lb_i;
            end
        
        
        end
        Flag4up=xnew>ub;
        Flag4low=xnew<lb;
        xnew=(xnew.*(~(Flag4up+Flag4low)))+ub.*Flag4up+lb.*Flag4low;
        
        Kp=xnew(1);
        Ki=xnew(2);
        Kd=xnew(3);
        mun=xnew(4);
        lam=xnew(5);
        sim('fpid');                                                     % proses simülasyonu
    load('error1.mat');                                                      % prosesten hata deðerleri çekiliyor
    load('error.mat');                                                      % prosesten hata deðerleri çekiliyor
    e = ITAE(2,:);                                                          % error vector
    
    SSE=e(end);
        fitnessnew=SSE;  % calculate the objective function
        % Update if the solution improves, or not too loud
        if (fitnessnew<=fmin) && (rand<A(ii)) ,
            
            fitness(ii)=fitnessnew;
            A(ii)=alpha*A(ii);
            r(ii)=ro*(1-exp(-gamma*iter));
            x(ii,:)=xnew;
            bestsol=x(ii,:);
            fmin=fitnessnew;
           
        end
%         if fitnessnew<=fmin,
%             bestsol=x(ii,:);
%             fmin=fitnessnew;
%         end
%         
    end
   
    
    iter=iter+1;
    fmin;          % update the while loop counter
end
%
[bestfit]=(fmin);
BestPositions=bestsol;
Kp=BestPositions(1)
Ki=BestPositions(2)
Kd=BestPositions(3)
mun=BestPositions(4)
lam=BestPositions(5)
sim('fpid');
 load('error1.mat');
  load('error.mat');
  load('error2.mat');
e = ITAE(2,:);                                                          % error vector
    time = ITAE(1,:);                                                       % time vector
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Peak;
    e2 = fark2(2,:);
    time2 = fark2(1,:);
    STI2 = stepinfo(e2,time2,1);
    ST2=STI2.SettlingTime
    PO2=STI2.Overshoot
    RT2 = STI2.RiseTime
