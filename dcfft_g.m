% This code solves the contact problem. All the variables are
% self-explanatory.Comments are mentioned where additional details are
% necessary. The input file is the height surface data generated using the
% R code. The output file contains contact evolution data.

rms_slope = 0.02;

al = 0.6;
bet = 1;
size = 2048;

nu1 = 0;
ESTAR = 1;

Ws = (0.005:0.0025:0.0525)*rms_slope;
wlength = length(Ws);
xaxis = zeros(1,wlength);   
yaxis = zeros(1,wlength);
zaxis = zeros(1,wlength);
waxis = zeros(1,wlength);

M = size;
N = size;

for fi = 1:10
    
    tic
    
    ii = fi;
%     load the height surface file generated using the 'rough_sruface_cefn.r' file 
    load(strcat('C:\Users\yjetti2\OneDrive - University of Illinois - Urbana\Desktop\Fractal_Mechanics\Contact_Implementation\new_simulations\al',num2str(al),'b',num2str(bet),'\hs_al',num2str(al),'b',num2str(bet),'_',num2str(size),'_',num2str(ii),'.mat'));

    c = 1/(M);
    dx = c;
    dy = c;
    delta = (M-1)*c;

    xs = (0:c:delta);
    ys = (0:c:delta);
    xlength = length(xs);
    ylength = length(ys);
    xsn = repmat(xs',length(ys),1);
    ysn = repmat(ys,1,length(xs));

    hxs = (max(max(hs))-(hs));

    hxsnew = hxs;

    hxslength = length(hxs);

    A0 = (M)^2*dx*dy;
    epsilon = 10^(-5);

    parfor q = 1:wlength
        disp(q)
        W = Ws(q);
        xaxis(q) = W/(ESTAR*A0);

        data_p = zeros(2*M,2*N);
        P0 = W/((dx*dy)*(M*N))*ones(M,N);
        pij = P0;
        G0 = 1;
        d = 0;
        Gkm1 = G0;

        error = inf;

        hijkm1 = -ones(M,N);
        tij = zeros(M,N);

        whind = 0;

        while error>epsilon
            whind = whind + 1;

            uzf = dcfft(data_p,M,N,ESTAR,pij);
            uzf = uzf';

            gij = hxs + uzf;
            Nc = -sum(sum(hijkm1));
            gbar = sum(sum(gij(hijkm1<0)))/Nc;
            gij = gij-gbar;

            Gk = sum(sum(gij(hijkm1<0).^2));
            betak = Gk/Gkm1;

            tij(hijkm1<0) = gij(hijkm1<0) + d*betak*tij(hijkm1<0);
            tij(hijkm1==0) = 0;
            Gkm1 = Gk; 

            rij = dcfft(data_p,M,N,ESTAR,tij);
            rij = rij';
            rijbar = sum(sum(rij(hijkm1<0)))/Nc;
            rij = rij-rijbar;
            tau = sum(sum(gij(hijkm1<0).*tij(hijkm1<0)))/sum(sum(rij(hijkm1<0).*tij(hijkm1<0)));

            pijol = pij;
            pij(hijkm1<0) = pij(hijkm1<0) - tau*tij(hijkm1<0);

            d = 1;
            for i = 1:M
                for j = 1:N
                    if pij(i,j)<0
                        pij(i,j) = 0;
                        if gij(i,j)<0
                            d = 0;
                            pij(i,j) = pij(i,j) - tau*gij(i,j);
                            hijkm1(i,j) = -1;
                        else
                            hijkm1(i,j) = 0;
                        end
                    end
                end
            end

            Wt = dx*dy*sum(sum(pij));
            pij = W/Wt*pij;                          

            error = dx*dy/W*sum(sum(abs(pij-pijol)));   
        end
        
        Sd = (sum(sum(abs(diff(hijkm1,1,1))))+sum(sum(abs(diff(hijkm1,1,2)))))*dy;
        At = -sum(sum(hijkm1))*dx*dy;
        Ap = At-(pi-1+log(2))/24*Sd*dx;
        yaxis(1,q) = Ap/A0;
        zaxis(1,q) = At;
        waxis(1,q) = Sd;
        disp(xaxis(q))
        disp(yaxis(1,q))
        disp(At)
        disp('***')
    end
%     saving the contact result as .mat file which will be used for
%     generating the contact evolution figures
    save(strcat('result_2048_al',num2str(al),'b',num2str(bet),'_',num2str(fi),'_ini_af_sd.mat'),'xaxis','yaxis','zaxis','waxis','Ws','rms_slope')    
    toc
end

function uzf = dcfft(data_p,M,N,ESTAR,input)
for i = 1:M
    for j = 1:N
        data_p(i,j) = input((i-1)*N+j);
    end
end
load(strcat('data_emwes',num2str(M),'tm.mat'),'data_emwes')
data_em = data_emwes/ESTAR;
uz = ifft2(fft2(data_em).*fft2(data_p));
uzf = (uz(1:M,1:N));
end