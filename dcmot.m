%**************************************************************************
% Global Neighborhood Algorithm (by A.Alazzam and H.W.Lewis, 2013)
%**************************************************************************
%clear all 0.01980  0.147 4

%model_parameters_14012017;
%**************************************************************************

N = 4;                                                                      % number of parameters to optimize 
M = 20;                                                                     % number of particles
T = 10 ;                                                                      % iteration number
x_min = [3 0.1 3 0];                                                            % parameter minimum limit
x_max = [4 0.25 4 200];                                                         % parameter maximum limit
delta_min = 0;                                                              % Neighborhood minimum limit
delta_max = 0.3;                                                            % Neighborhood maximum limit
t = 1;   
Ni=100;% iteration counter start
%**************************************************************************
for m = 1:M,
    for n = 1:N,
        x(m,n) = x_min(1,n)+(x_max(1,n)-x_min(1,n))*rand(1,1);
    end
end
fx_b = ones(1,M)';                                                          % initial eligibility values
xf = [x fx_b];                                                              % initial candidate matrix (including conformances)
%**************************************************************************
for i = 1:M,
    Kp = xf(i,1);
    Ki = xf(i,2);
    Kd = xf(i,3);
   Ni=xf(i,4);
 % Ni=500;
    sim('pidn');                                                     % proses simülasyonu
    load('error1.mat');                                                      % prosesten hata deðerleri çekiliyor
    load('error.mat');                                                      % prosesten hata deðerleri çekiliyor
    e = ITAE(2,:);                                                          % error vector
    time = ITAE(1,:);                                                       % time vector
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Peak;
    %SSE = (0.75*ST)+(0.18*e(end))+(0.07*PO); 
    SSE=e(end);
% R=1-0.3032*Kp+0.2813*Ki-0.0229*Kd;
% I=0.2435*Kp+0.0742*Ki-0.5252*Kd;
% SSE=abs(R)+abs(I);
         % the starting value of the total error
    f(i) = SSE;                                                             % total error
    xf(i,N+1) = f(i);
    
end
xf_sira = sortrows(xf,[N+1]);                                                 % ranking by relevance
x_best = xf_sira(1:M,1:N+1);                                                  % best order parameter values
f_best_son = xf_sira(1,N+1);                                                  % best fit
g_son = x_best(1,:);                                                        %the solution with the best fit
fbs = 0.2*ones(T-1,1);
%**************************************************************************
while t <= T && SSE>0,                                                               % iteration start
    t=t
    iyi_xf = xf_sira(1:M/2,1:N+1);                                          % best 50%
    for md = 1:M/2,                                                         % plus / minus delta neighborhood
        for nd = 1:N+1,
            delta = delta_min+(delta_max-delta_min)*unifrnd(-1,1);
            komsu_xf(md,nd) = iyi_xf(md,nd)+delta;
            
%              if komsu_xf(md,nd)<x_min(nd)
%                 komsu_xf(md,nd)=x_min(nd);
%              elseif komsu_xf(md,nd)>x_max(nd)
%                    komsu_xf(md,nd)=x_max(nd);
%              end
             
        end
    end
   for lamda =1:M/2
       Flag4up=komsu_xf(lamda,1:N)>x_max;
        Flag4low=komsu_xf(lamda,1:N)<x_min;
        komsu_xf(lamda,1:N)=(komsu_xf(lamda,1:N).*(~(Flag4up+Flag4low)))+x_max.*Flag4up+x_min.*Flag4low;
   end
    %**********************************************************************
    for m = 1:M/2,                                                          % the remaining 50% are reproduced
        for n = 1:N,
            x_kalan(m,n) = x_min(1,n)+(x_max(1,n)-x_min(1,n))*rand(1,1);
        end
    end
    yeni_xf = [x_kalan ones(M/2,1)];                                  %the remaining 50% were randomly produced
    xf_sira_yeni = [komsu_xf;yeni_xf];                               % new candidate matrix (including eligibility)
    %**********************************************************************
    for i = 1:M,
        Kp = xf_sira_yeni(i,1);
        Ki = xf_sira_yeni(i,2);
        Kd = xf_sira_yeni(i,3);
       Ni = xf_sira_yeni(i,4);
      % Ni=500;
        sim('pidn');                                                 % proses simülasyonu
        load('error1.mat');                                                  % prosesten hata deðerleri çekiliyor                                          
         load('error.mat');
        e = ITAE(2,:);                                                          % error vector
    time = ITAE(1,:);                                                       % time vector
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Peak;
   % SSE = (0.75*ST)+(0.18*e(end))+(0.07*PO); 
   SSE=e(end);
 % the starting value of the total error
%  R=1-0.3032*Kp+0.2813*Ki-0.0229*Kd;
% I=0.2435*Kp+0.0742*Ki-0.5252*Kd;
% SSE=abs(R)+abs(I);
    f(i) = SSE;                                                         % total error
        xf_sira_yeni(i,N+1) = f(i);
    end
    xf_sira = sortrows(xf_sira_yeni,[N+1]);                                % uygunluða göre sýralama
    x_best = xf_sira(1:M,1:N);                                              % en iyi sýralamalý parametre deðerleri
    f_best = xf_sira(1,N+1);                                                  % en iyi uygunluk
    g = x_best(1,:);                                                        % en iyi uygunluða sahip çözüm
    %**********************************************************************   
    if f_best < f_best_son,                                                 % en iyi uygunluk karþýlaþtýr
        f_best_son = f_best;                                                % en iyi uygunluðu sakla
        g_son = g;                                                          % en iyi sonucu sakla
    end
    t=t+1;                                                                  % iterasyon sayacýný bir arttýr
    f_best_son
    fbs(t) = f_best_son;
end
sonuc = g_son
Kp = g_son(1);         
Ki = g_son(2);
Kd = g_son(3);
Ni=g_son(4);
sim('pidn');
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
%subplot(2,1,1), plot(tout,freq_1);
%subplot(2,1,2), plot(fbs);