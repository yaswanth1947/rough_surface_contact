% %plotting kappa^{\prime\prime} as a function of contact area evolution.
%The input files are those generated using 'dcfft_g.m'. 

%Cauchy rf
%fractal dim
% als = [1 0.6 0.2];%[1.8 1.4 1 0.6 0.2];
% bets =[1 1 1];%[1 1 1 1 1];
%hurst exp
% als = [1 1 1 1 1];
% bets = [1.8 1.4 1 0.6 0.2];

%Dagum rf
%fractal dim
% als = [1 0.8 0.6 0.4 0.2];
% bets = [1 1 1 1 1];
%hurst exp
als = [0.2 0.2 0.2 0.2 0.2];
bets = [1 0.8 0.6 0.4 0.2];

cols = ['b' 'r' 'y' 'm' 'g'];
% cols = 'r';

for j = 1:length(als)
    al = als(j);
    bet = bets(j);
    size = 2048;
    yaxisall = zeros(10,20);

    for fi =1:10
        ii = fi;
%         load(strcat('C:\Users\Yaswanth\OneDrive - University of Illinois - Urbana\Desktop\Fractal_Mechanics\Contact_Implementation\new_simulations\al',num2str(al),'b',num2str(bet),'\result_',num2str(size),'_al',num2str(al),'b',num2str(bet),'_',num2str(ii),'_ini_af.mat'));
        load(strcat('C:\Users\Yaswanth\OneDrive - University of Illinois - Urbana\Desktop\Fractal_Mechanics\Contact_Implementation\new_simulations\gam',num2str(al),'b',num2str(bet),'\result_',num2str(size),'_gam',num2str(al),'b',num2str(bet),'_',num2str(ii),'_ini_af.mat'));
        yaxisall(fi,:) = yaxis;
    end

    stdm = std(yaxisall);
    stder = stdm/sqrt(10);
    meanv = mean(yaxisall);
    
    slopev = rms_slope*[-3*meanv(1)+4*meanv(2)-meanv(3) meanv(3:end)-meanv(1:end-2) meanv(end-2)-4*meanv(end-1)+3*meanv(end)]/(2*(xaxis(2)-xaxis(1)));
    
    figure(1)
    errorbar(xaxis/rms_slope,slopev,stder,'o-k','MarkerSize',5,'MarkerFaceColor',cols(j),'LineWidth',1)
    hold on
end
figure(1)
% %Cauchy
% leg1 = legend({'D=2.1','D=2.3','D=2.5','D=2.7','D=2.9'},'Location','best');
% leg1 = legend({'D=2.5','D=2.7','D=2.9'},'Location','best');
% leg2 = legend({'H=0.1','H=0.3','H=0.5','H=0.7','H=0.9'},'Location','best');
% title(leg1,'H=0.5')
% title(leg2,'D=2.5')
% % %Dagum
% leg3 = legend({'D=2.5','D=2.6','D=2.7','D=2.8','D=2.9'},'Location','best');
leg4 = legend({'H=0.5','H=0.6','H=0.7','H=0.8','H=0.9'},'Location','best');
% title(leg3,'H=0.5')
title(leg4,'D=2.9')
ax = gca;
ax.FontSize = 12; 
ax.FontName = 'Times';
xlabel('$$W/(E^{\prime}A_0\sqrt{\langle|\nabla h|^2\rangle})$$','Interpreter','latex','FontSize',14)
ylabel('$$E^{\prime}\sqrt{\langle|\nabla h|^2\rangle}{\partial A}/{\partial W}$$','Interpreter','latex','FontSize',14)
% %Cauchy
% title('Cauchy random fields')
% % % title('Cauchy random fields: H=0.5')
% % % title('Cauchy random fields: D=2.5')
% % %Dagum
% title('Dagum random fields')
% % % title('Dagum random fields: H=0.5')
% % % title('Dagum random fields: D=2.9')
