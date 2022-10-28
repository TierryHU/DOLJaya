function [x0,c,sig]=Singleparticlegeneration(c,sig,y0)

%% Inverse transform sampling (Inverse function method) 
% % global Icdf a b d
% %     % definite integral value from PDF
% %     a_s=a;
% %     b_s=b;
% %     d_s=d;
% %      d=eval(d);
% %      a=eval(a);
% %      b=eval(b);
% %         x_all=eval(Icdf);% solution including complex numbers    
% %         x0=real(x_all);
% %         b=b_s;
% %         a=a_s;
% %         d=d_s;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Accompany global variables %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % global Icdf
% % syms x y cdf sig c y0
% %     cdf =-y0+0.5+(1592262918131443*((1764878998669257*(c - x)^2)/(36028797018963968*sig^2) - 2069847464455053/562949953421312)*(c - x))/...
% %         (2251799813685248*sig*((c - x)^2/(2*sig^2) + 1834349511160241/562949953421312)*((2^(1/2)*((1764878998669257*(c - 1)^2)/...
% %         (36028797018963968*sig^2) - 2069847464455053/562949953421312)*(c - 1))/(2*sig*((c - 1)^2/(2*sig^2) + 1834349511160241/562949953421312))...
% %         - (2^(1/2)*((1764878998669257*(c + 1)^2)/(36028797018963968*sig^2) - 2069847464455053/562949953421312)*(c + 1))/...
% %         (2*sig*((c + 1)^2/(2*sig^2) + 1834349511160241/562949953421312))));
% %     % definite integral value from PDF
% %     Icdf_all=solve(cdf,x);
% % Icdf=Icdf_all(2);

%% PDF Sampler (Mont Calro)
%         Lowerbound=-1;
%         Upperbound=1;
%         splits=10;
% %       bins=10;
%         x=linspace(Lowerbound,Upperbound,splits);
%         Deltax=x(2)-x(1);
%         c=c_0*ones(1,splits);%2;%c_0;%PV(i,1);%
%         sig=sig_0*ones(1,splits);%3;%sig_0;%PV(i,2);
%         pdf=(exp(-(x-c).^2./(2*(sig.^2))).*sqrt(2/pi))./(sig.*(erf((c+1)./(sqrt(2).*sig))-erf((c-1)./(sqrt(2).*sig))));
% 
%         %         sampler=PDFsampler(y,9);
% %         x0=sampler.nextRandom;
% %         [pdf_times,pdf_quartiles]=hist(y,bins);
% %         pdf=pdf_times/sum(pdf_times);
%         cdf=zeros(splits,1);
% %         cdf(1)=pdf(1)*Deltax;
%         for i=2:splits
%             cdf(i)=cdf(i-1)+(Deltax*(pdf(i)+pdf(i-1))/2);
%         end
% %         plot(x,cdf);
%         x0=-2;
%         r=rand;%rand;
%         if r<cdf(1)
%             x0 = x(1)+(r/cdf(1))*(x(2)-x(1));            
%         else if r>cdf(end)
%                 x0=x(splits);
%              else
%                 for j=1:splits-1
%                     if(r>=cdf(j) && r<cdf(j+1))
%                         x0 = x(j)+((r-cdf(j))/(cdf(j+1)-cdf(j)))*(x(j+1)-x(j));
%                     end
%                 end      
%              end
%         end
%         if x0==-2
%             x0=0;
%             return; 
%         end
%% Accept and Reject
% global generationflag
% pdfunifyflag=0;
% if sig>0.2
%     if c>=1
%         F_c=(exp(-(1-c)^2/(2*(sig^2)))*sqrt(2/pi)/(sig*(erf((c+1)/(sqrt(2)*sig))-erf((c-1)/(sqrt(2)*sig)))));
%         c=1;
%     else if c<=-1
%             F_c=(exp(-(-1-c)^2/(2*(sig^2)))*sqrt(2/pi)/(sig*(erf((c+1)/(sqrt(2)*sig))-erf((c-1)/(sqrt(2)*sig)))));
%             c=-1;
%         else
%             F_c=(exp(-(c-c)^2/(2*(sig^2)))*sqrt(2/pi)/(sig*(erf((c+1)/(sqrt(2)*sig))-erf((c-1)/(sqrt(2)*sig)))));
%         end
%     end
%     if sig>2%generationflag<10
%         M=0.5/F_c;
%     else
%         M=(normpdf(c,c,sig))/F_c;
%     end
% else
%     pdfunifyflag=1; % the distribution see the same distribution with the normal distribution
% end
% if pdfunifyflag==0
%     while(1)
%         if sig>2%generationflag<10 % The distribution has converged into the -1~1 range
%             x=2*(rand-0.5);%normrnd(c,sig,1,1);
%         else 
%             x=normrnd(c,sig,1,1);
%             if x>1 || x<-1
%                 x=2*(rand-0.5);%normrnd(c,sig,1,1);
%             end
%         end
%             u=rand;
%             F=(exp(-(x-c)^2/(2*(sig^2)))*sqrt(2/pi)/(sig*(erf((c+1)/(sqrt(2)*sig))-erf((c-1)/(sqrt(2)*sig)))));
%         if sig>2%generationflag<10    
%             MG=M*0.5;%(normpdf(x,c(1),sig(1)));
%         else
%             MG=M*(normpdf(x,c(1),sig(1)));
%         end
%         if u<F/MG %&& x>=-1 && x<=1
%             x0=x;
%             break;
%         end
% %         disp(['Something happens'] );   
%         if generationflag<10
%             generationflag=generationflag+1;
%         end
%     end
% else
%     if c>1
%         c=1-0.1*rand;
%     else if c<-1
%             c=-1+0.1*rand;
%         end
%     end
%     x0=normrnd(c,sig,1,1);        
% end
% %% Faster direct generate
%% Approximation
% x=0.8;
% %x<0.5
% erfstand=erf(x);
% erfap=x*(3.6767877+(-0.097970465)*(x^2))/(3.2584593+x^2);
% erfcstand=erfc(x);
% %0.5<=x<=4
% erfcap=exp(-x^2)*(0.73033+(-0.023877)*x)/(0.66211+x);
% %x>4
% erfcap=exp(-x^2)/x*(1/(sqrt(pi))+1/x^2*((-0.124368544)+(-0.0968210364)*x^(-2))/(0.440917061+x^(-2)))
% 
% syms x y cdf sig c y1 y2 y3 a b
% a=(2^(1/2)*(c - 1))/(2*sig);
% b=(2^(1/2)*(c + 1))/(2*sig);
% d=(1/(2*sig^2))^(1/2)*(c - x);
% y1=a*(3.6767877+(-0.097970465)*(a^2))/(3.2584593+a^2);
% y2=b*(3.6767877+(-0.097970465)*(b^2))/(3.2584593+b^2);
% y3=d*(3.6767877+(-0.097970465)*(d^2))/(3.2584593+d^2);
% pdfa=(exp(-(x-c)^2/(2*(sig^2)))*sqrt(2/pi))/(sig*y3);
% %         pdf=(exp(-(x-c).^2./(2*(sig.^2))).*sqrt(2/pi))./(sig.*(erf((c+1)./(sqrt(2).*sig))-erf((c-1)./(sqrt(2).*sig))));
% % pdfa =(7186705221432913*exp(-(c - x)^2/(2*sig^2)))/(9007199254740992*sig*((2^(1/2)*((1764878998669257*(c - 1)^2)/...
% %     (36028797018963968*sig^2) - 2069847464455053/562949953421312)*(c - 1))/(2*sig*((c - 1)^2/(2*sig^2) +...
% %     1834349511160241/562949953421312)) - (2^(1/2)*((1764878998669257*(c + 1)^2)/(36028797018963968*sig^2) -...
% %     2069847464455053/562949953421312)*(c + 1))/(2*sig*((c + 1)^2/(2*sig^2) + 1834349511160241/562949953421312))));
% cdfa=(7186705221432913*pi^(1/2)*(-((1/(2*sig^2))^(1/2)*((1764878998669257*(c - x)^2)/(36028797018963968*sig^2)...
%     - 2069847464455053/562949953421312)*(c - x))/((c - x)^2/(2*sig^2) + 1834349511160241/562949953421312)))/...
%     (18014398509481984*sig*(-(2^(1/2)*((1764878998669257*(c - 1)^2)/(36028797018963968*sig^2) - 2069847464455053/562949953421312)*(c - 1))...
%     /(2*sig*((c - 1)^2/(2*sig^2) + 1834349511160241/562949953421312))-(-(2^(1/2)*((1764878998669257*(c + 1)^2)/(36028797018963968*sig^2)...
%      -2069847464455053/562949953421312)*(c + 1))/(2*sig*((c + 1)^2/(2*sig^2)+ 1834349511160241/562949953421312))))*(1/(2*sig^2))^(1/2));
% 
% cdfa=(1592262918131443*((1764878998669257*(c - x)^2)/(36028797018963968*sig^2) - 2069847464455053/562949953421312)...
%     *(c - x))/(2251799813685248*sig*((c - x)^2/(2*sig^2) + 1834349511160241/562949953421312)*((2^(1/2)*((1764878998669257*(c - 1)^2)/...
%     (36028797018963968*sig^2) - 2069847464455053/562949953421312)*(c - 1))/(2*sig*((c - 1)^2/(2*sig^2) + 1834349511160241/562949953421312))...
%     - (2^(1/2)*((1764878998669257*(c + 1)^2)/(36028797018963968*sig^2) - 2069847464455053/562949953421312)*(c + 1))/(2*sig*((c + 1)^2/(2*sig^2)...
%     + 1834349511160241/562949953421312))))
% 
% Icdfa=finverse(cdfa,x)
%% Novel quick method
        Lowerbound=-1;
        Upperbound=1;
        NotAchieveFlag=1;% The required solution is not achieved
        while(NotAchieveFlag)
            if sig<=Upperbound
                x=normrnd(c,sig,1,1);
                if x>Lowerbound&&x<Upperbound
                    NotAchieveFlag=0; % The required solution is achieved
                end
            else 
                x=2*rand-1;
                NotAchieveFlag=0; % The required solution is achieved
            end
        end
        x0=x;
end