clear;

phi_2=[1, 1, 1]';
z_2  =[1, 2, 3, 4]';

A    =1.*...
      [ 1,  1,  5;...
        2,  1,  1;...
        1,  2,  5;...
        1,  1,  1];

q_2     =A*phi_2;

[phi_new,q_new] =mem_solve(q_2,A,1,4*[1,1,1]',1,10000);

% phi_new =phi_2;
% 
% q_2     =[1, 1, 1, 1];
% phi_new =ones(size(phi_2,1),1);
% beta    =0.1;
% t       =zeros(size(z_2,1),1);
% c       =zeros(size(z_2,1),1);
% var_q   =1;
% w=rand(size(z_2,1),1);
% w=w./sum(w);
% k=0;
% chi=5;
%  while chi>0.01
% 
%      for i=1:1:size(z_2,1)
%          a=(A(i,:).*phi_new').^-1;
%         t(i,1)=min(a(A(i,:)~=0));
%      end;
%      
%      c=beta.*(1-(A*phi_new)./q_2').*t;
%      
%      for j=1:1:size(phi_2,1)
%          L(j)=0;
%          for i=1:1:size(z_2,1)
%             L(j)=L(j)+w(i).*c(i).*A(i,j);
%          end;
%          phi_new(j)=phi_new(j)./(1-phi_new(j).*L(j)');
%      end;
%      chi=sum((q_2'-A*phi_new).^2)./var_q;
%      k=k+1;
%      chi
%  end