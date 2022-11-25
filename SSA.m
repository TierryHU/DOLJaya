%_________________________________________________________________________________
%  Salp Swarm Algorithm (SSA) source codes version 1.0
%
%  Developed in MATLAB R2016a
%
%  Author and programmer: Seyedali Mirjalili
%
%         e-Mail: ali.mirjalili@gmail.com
%                 seyedali.mirjalili@griffithuni.edu.au
%
%       Homepage: http://www.alimirjalili.com
%
%   Main paper:
%   S. Mirjalili, A.H. Gandomi, S.Z. Mirjalili, S. Saremi, H. Faris, S.M. Mirjalili,
%   Salp Swarm Algorithm: A bio-inspired optimizer for engineering design problems
%   Advances in Engineering Software
%   DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.07.002
%____________________________________________________________________________________

function [Best_score]=SSA(Gm,D,Np,lb,ub,fobj,func_num)

disp('                    SSA Algorithm                    ')
% if size(ub,1)==1
%     ub=ones(dim,1)*ub;
%     lb=ones(dim,1)*lb;
% end
Upperbound=ones(1,D)*ub;
Lowerbound=ones(1,D)*lb;

% Gm=Gm_o;
Best_score = zeros(1,Gm);

%Initialize the positions of salps
% pop=initialization(N,dim,ub,lb);
pop=repmat(Lowerbound,Np,1)+rand(Np,D).*(repmat(Upperbound,Np,1)-repmat(Lowerbound,Np,1));


FoodPosition=zeros(1,D);
FoodFitness=inf;


%calculate the fitness of initial salps

for i=1:size(pop,1)
    SalpFitness(1,i)=fobj(pop(i,:)',func_num);
end

[sorted_salps_fitness,sorted_indexes]=sort(SalpFitness);

for newindex=1:Np
    Sorted_salps(newindex,:)=pop(sorted_indexes(newindex),:);
end

FoodPosition=Sorted_salps(1,:);
FoodFitness=sorted_salps_fitness(1);

%Main loop
G=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness of salps
while G<=Gm
    
    c1 = 2*exp(-(4*G/Gm)^2); % Eq. (3.2) in the paper
    
    for i=1:size(pop,1)
        
        pop= pop';
        

        if i<=Np/2
            for j=1:1:D
                c2=rand();
                c3=rand();
                %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                if c3<0.5 
                    pop(j,i)=FoodPosition(j)+c1*((Upperbound(j)-Lowerbound(j))*c2+Lowerbound(j));
                else
                    pop(j,i)=FoodPosition(j)-c1*((Upperbound(j)-Lowerbound(j))*c2+Lowerbound(j));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        elseif i>Np/2 && i<Np+1
            point1=pop(:,i-1);
            point2=pop(:,i);
            
            pop(:,i)=(point2+point1)/2; % % Eq. (3.4) in the paper
        end
        
        pop= pop';
    end
    
    for i=1:size(pop,1)
        
        Tp=pop(i,:)>Upperbound;Tm=pop(i,:)<Lowerbound;
        pop(i,:)=(pop(i,:).*(~(Tp+Tm)))+Upperbound.*Tp+Lowerbound.*Tm;
        
        SalpFitness(1,i)=fobj(pop(i,:)',func_num);
        
        if SalpFitness(1,i)<FoodFitness
            FoodPosition=pop(i,:);
            FoodFitness=SalpFitness(1,i);
            
        end
    end
    
    Best_score(G)=FoodFitness;
    G = G + 1;
end



