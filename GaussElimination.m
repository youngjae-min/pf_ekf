function A=GaussElimination(A,type)
%parallel computing Gauss elimination
p=size(A,1);
a=ones(size(A));
a(A==0)=0;
for i=1:p
    a(:,i)=a(:,i).*[1:p]';
end
P=zeros(1,p);
k=1;j=1;m=ones(1,p);
while true
    if j==1
        if a(m(j),j)~=0
            P(1,j)=a(m(j),j);
            j=j+1;
        else
            m(j)=m(j)+1;
            continue
        end
    end
    
    k=setdiff(a(:,j),P(1,:));
    k(k==0)='';
    if isempty(k)|m(j)>length(k)
        j=j-1;
        m(j)=m(j)+1;
        m(j+1)=1;
        if j==0
            j=1;
        end
        P(1,j)=0;
        continue
    else
        P(1,j)=k(m(j));
        j=j+1;
    end
    if isempty(find(P==0,1))
        break
    end
end
%%
%elimination in loop
A=A(P,:);
for i=1:p-1
    A1=A;
    parfor j=i+1:p
        A1(j,:)=((A(j,:)-A(i,:)/A(i,i)*A(j,i)));
    end
    A=A1;
end
if strcmp(type,'full')
    for i=1:p-1
        A1=A;
        parfor j=i:p-1
            A1(p-j,:)=A(p-j,:)-A(p-i+1,:)/A(p-i+1,p-i+1)*A(p-j,p-i+1);
        end
        A=A1;
    end
end
end