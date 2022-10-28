function [Best_score]=Jaya(Gm_o,D,Np,lb,ub,fobj,func_num)

% Coded by Tianyu Hu (Thierry) and Zhile Yang (Slerk)
% Modified on 07-03-2020
                disp('                    Jaya Algorithm                    ')
% R. Venkata Rao 'Jaya: A simple and new optimization algorithm for solving
% constrained and unconstrained optimization problems'
%***************Initialization**************%
% roundnumber=20;
% best=100;
% Optimal=zeros(1,roundnumber);
% sumall=0;
% worst=0;
%     for q=1:roundnumber
% % Gm=200; 
% % Lowerbound=-5;Upperbound=10;
% % Np=10; 
% % D=30;
% % St=(Upb-Lowb)*rand(Np,D)+Lowb;
Lowerbound=ones(1,D)*lb;
Upperbound=ones(1,D)*ub;
pop=repmat(Lowerbound,Np,1)+rand(Np,D).*(repmat(Upperbound,Np,1)-repmat(Lowerbound,Np,1));
Best_score=inf;
Gm=Gm_o;
G=1; 
ge=zeros(1,Gm);
Best=zeros(Gm,D);
Worst=zeros(Gm,D);
% Tr_new=zeros(Gm,D);
Fit_Tr=zeros(Np,1);
St=pop;
% St_Mean=zeros(1,D);
St_new=zeros(Np,D);
St_selfnew=zeros(Np,D);
Pop_Mean=zeros(Gm,D);
Pop_Variance=zeros(Gm,D);
     while G<=Gm
    %***************Teacher Phase**************%
        for i=1:D
            Sum_Column=0;
            for j=1:Np
              Sum_Column=Sum_Column+St(j,i);
            end
            Pop_Mean(G,i)=Sum_Column/Np;
        end
    %**********Calculate the variance of particles in each generation*********%  
        for i=1:D
            Sum_Column_v=0;
            for j=1:Np
              Sum_Column_v=Sum_Column_v+(St(j,i)-Pop_Mean(G,i))^2;
            end
            Pop_Variance(G,i)=Sum_Column_v/Np;
        end
            for j=1:Np
              Fit_Tr(j)=fobj((St(j,:))',func_num);
            end
         [Yb,Bindex]=min(Fit_Tr); % Best index
         [Yw,Windex]=max(Fit_Tr); % Worst index        
         Best(G,:)=St(Bindex,:);
         Worst(G,:)=St(Windex,:);
         for k=1:Np
            r1=rand;
            r2=rand;
            %St_new(k,:)=St(k,:)+r1*(Bindex(G,:)-((round(1+r2)*Pop_Mean(G,:))));
            St_new(k,:)=St(k,:)+r1*(Best(G,:)-abs(St(k,:)))-r2*(Worst(G,:)-abs(St(k,:)));
            [St_new]=Checkbound(St_new,Lowerbound,Upperbound,Np,D,G);
            if fobj(St_new(k,:)',func_num)<fobj(St(k,:)',func_num)
                St(k,:)=St_new(k,:);
            end
         end
        LocalBest=Fit_Tr(Bindex);
        if G==1
            ge(1)=LocalBest;
        end
        if G>1
            if LocalBest<ge(G-1)
                ge(G)=LocalBest;
             else 
                ge(G)=ge(G-1);   
            end
        end
        Best_score(G)=ge(G);
        G=G+1;
     end         
    %         optimal(segnumber,:)=Tr(Gm,:);
    %         Fit_value(segnumber)=Bat_Second_ord(bestx(Gm,:),y,n,v_T,cur,ocv_s);
% %         ii=linspace(1,Gm,Gm);
% %         plot(ii,ge,'r-')
%         Optimal(1,q)=ge(G-1);
%         if Optimal(1,q)<best
%             best=Optimal(1,q);
%         end            
%         if worst<Optimal(1,q)
%             worst=Optimal(1,q); 
%         end
%         sumall=sumall+Optimal(1,q);
%     mean=sumall/roundnumber;
 end
% function y=f(Z)
%     y= sphere(Z);
% end
