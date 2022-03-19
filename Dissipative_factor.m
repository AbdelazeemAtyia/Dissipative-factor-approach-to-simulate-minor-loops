clc
clear
load H_hr;load B_hr; % the neasured data of the harmonic hystresis
load Rt; % the dissipative factor determined for each synnetric minor loop
load Bt; % the maximum of B at each minor loop
load p35h; % the parameter of the major loop without harmonic
   
pnhr=p35h;
Ms=pnhr(1);k=pnhr(2);c=pnhr(3);alpha=pnhr(4);a=pnhr(5);
 nr=6;% the number of repeated loops of Hm,and Bm
 Hm=H_hr(:,4);Bm=B_hr(:,4);
n_points = length(Hm);
mu0=4*pi*10^-7;
 Hm=repmat(Hm,[nr,1]);Bm=repmat(Bm,[nr,1]);
 
      B=zeros(size(Hm));
  jj=0;nl=length(Hm)/6;
M=zeros(size(Hm));
M(1)=(B(1)/mu0)-(Hm(1));
 for(j=1:size(Hm,2))
 for i=2:size(Hm,1)
  
  He=Hm(i-1,j)+alpha*M(i-1,j);
  Man=Ms*(coth((He)/a)-a/(He));
  dMandHe=Ms*(1-(coth((He)/a)^2)+(a/(He))^2)/a;
  Mirr=(M(i-1,j)-c*Man)/(1-c); 
if i==(nl+1)+jj
         jj=jj+nl;
     end
 dH=Hm(i)-Hm(i-1);dB=Bm(i)-Bm(i-1);
 if(dH>0)
 delta=1;
 else
 delta=-1;
 end
 R=1;
 if i>=2+jj & i<=20+jj       % the first part of the descending branch    
RR=1;  % RR is the disspative factor
     elseif i>=21+jj & i<=50+jj         % the first minor loop on  the descending branch   
 RR=interp1(Bt,Rt,abs(Bm(i,j)),'linear');
      elseif i>=78+jj & i<=102+jj        % the second minor(-ve quarter) loop on  the descending branch   
 RR=interp1(Bt,Rt,abs(Bm(i,j)),'linear'); 
    elseif i>=145+jj & i<=175+jj         % the first minor loop on  the ascending branch  
 RR=interp1(Bt,Rt,abs(Bm(i,j)),'linear');
   elseif i>=203+jj & i<=227+jj         % the first minor loop on  the ascending branch  
 RR=interp1(Bt,Rt,abs(Bm(i,j)),'linear');
 else                                    % the remaing parts of the hystresis loop
  RR=1;
 end
      if (Man-Mirr)*dH>0
      deltaM=1;
   end
  if (Man-Mirr)*dH<=0
      deltaM=1;  
  end
dMirrdHe=deltaM*(Man-RR*Mirr)/(k*delta);
dMdH=(((1-c)*(dMirrdHe))+(c*dMandHe))/(1-(alpha*c*dMandHe)-(alpha*(1-c)*dMirrdHe));
   M(i,j)=M(i-1,j)+dMdH*dH;
  B(i,j)=mu0*(Hm(i,j)+M(i,j));
   end
 end
 Bs=B(1:end,1:j);Hm=Hm(1:end,1:j);
 figure(2)
plot(Hm(end-n_points:end),Bm(end-n_points:end),'r','LineWidth',1.5);
hold on
 plot(Hm(end-n_points:end),Bs(end-n_points:end),'b','LineWidth',1.5);
 hold on
  hold off
legend('Meas','Simu','Location','northwest')
title(['JA Res at B = ', num2str(max(Bm)),'T'],'fontweight','bold','fontsize',15);
xlabel('H[A/m]'),ylabel('B[T]')
set(gca,'FontSize',15,'fontweight','bold') 
