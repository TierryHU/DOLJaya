function [ge]=DOLJaya(Gm,D,Np,lb,ub,fobj,func_num)
% Coded by Tianyu Hu (Thierry)
% Update on 2021
                disp('                    DOLJaya Algorithm                    ')
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
G=1; 
ge=zeros(1,Gm);
Tr=zeros(Gm,D);
Fit_Tr0=zeros(2*Np,1);
StO=zeros(Np,D);
StO_new=zeros(2*Np,D);
Fit_Tr=zeros(Np,1);
Upperbound1=zeros(1,D);
Lowerbound1=zeros(1,D);
St=pop;
op=zeros(Np,D);
Jr=0.3; %Jumping rate
w=15; %Weight
St_new=zeros(Np,D);
Pop_Mean=zeros(Gm,D);
Pop_Variance=zeros(Gm,D);
    %***************Initialization**************%
       for i=1:Np
              for l=1:D
                  Upperbound1(l)=max(St(:,l));
                  Lowerbound1(l)=min(St(:,l)); 
                  op(i,l)=Upperbound1(l)+Lowerbound1(l)-St(i,l);
              end
            StO(i,:)=St(i,:)+w*rand*(rand*op(i,:)-St(i,:));
            StO_new(i,:)=St(i,:);
            StO_new(Np+i,:)=StO(i,:);
            [StO_new]=Checkbound(StO_new,Lowerbound1,Upperbound1,2*Np,D,G);
        end
          for l=1:2*Np
               Fit_Tr0(l)=fobj(StO_new(l,:)',func_num);
               if Fit_Tr0(l)<Best_score
               Best_score=Fit_Tr0(l);
               end
          end
             [Value,Index]=sort(Fit_Tr0);
          for i=1:Np
                 St(i,:)=StO_new(Index(i),:);
          end
     while G<=Gm
    %***************Teacher Phase**************%
        for i=1:D
            Sum_Column=0;
            for l=1:Np
              Sum_Column=Sum_Column+St(l,i);
            end
            Pop_Mean(G,i)=Sum_Column/Np;
        end
    %**********Calculate the variance of particles in each generation*********%  
        for i=1:D
            Sum_Column_v=0;
            for l=1:Np
              Sum_Column_v=Sum_Column_v+(St(l,i)-Pop_Mean(G,i))^2;
            end
            Pop_Variance(G,i)=Sum_Column_v/Np;
        end
            for l=1:Np
              %Fit_Tr(j)=fobj(St(j,:),func_num);
              Fit_Tr(l)=fobj(St(l,:)',func_num);
            end
         [Yb,Bindex]=min(Fit_Tr); % Best index
         [Yw,Windex]=max(Fit_Tr); % Worst index        
         Best(G,:)=St(Bindex,:);
         Worst(G,:)=St(Windex,:);
         for m=1:Np
            r1=rand;
            r2=rand;
            %St_new(k,:)=St(k,:)+r1*(Bindex(G,:)-((round(1+r2)*Pop_Mean(G,:))));
            St_new(m,:)=St(m,:)+r1*(Best(G,:)-abs(St(m,:)))-r2*(Worst(G,:)-abs(St(m,:)));
            [St_new]=Checkbound(St_new,Lowerbound,Upperbound,Np,D,G);
            if fobj(St_new(m,:)',func_num)<fobj(St(m,:)',func_num)
                St(m,:)=St_new(m,:);
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
        G=G+1;
     
                %***************Dynamic opposite phase**************%  
    if rand<Jr
       for i=1:Np
              for l=1:D
                  Upperbound1(l)=max(St(:,l));
                 Lowerbound1(l)=min(St(:,l)); 
                  op(i,l)=Upperbound1(l)+Lowerbound1(l)-St(i,l);
              end
            StO(i,:)=St(i,:)+w*rand*(rand*op(i,:)-St(i,:));
            StO_new(i,:)=St(i,:);
            StO_new(Np+i,:)=StO(i,:);
            [StO_new]=Checkbound(StO_new,Lowerbound1,Upperbound1,2*Np,D,G);
        end
          for l=1:2*Np
               Fit_Tr0(l)=fobj(StO_new(l,:)',func_num);
               if Fit_Tr0(l)<Best_score
               Best_score=Fit_Tr0(l);
               end
               
          end
          
          
             [Value,Index]=sort(Fit_Tr0);
          for i=1:Np
                 St(i,:)=StO_new(Index(i),:);
          end            
        LocalBest=Fit_Tr(Bindex);
        if G==1
            %ge(1)=LocalBest;
            ge(1)=Best_score;
        end
        if G>1
            if Best_score<ge(G-1)
                %LocalBest<ge(G-1)
                ge(G)=Best_score;
                %ge(G)=LocalBest;
             else 
                ge(G)=ge(G-1);   
            end
        end
        %Best_score(G)=ge(G);
        %Best_score(:,Gm+1)=[];
        %Best_score(:,Gm+2)=[];
        %Best_score(:,Gm+3)=[];
        %best_score(:,G)=Best_score;
        [l,m]=size(ge);
if m==Gm+1
    ge(:,1)=[];
else
    if m==Gm+2
        ge(:,1:2)=[];
    end
        %ge(:,Gm+1)=[];
        G=G+1;        

end
    
    end
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
