function [position num]=simulation_model_movie(init,width,tmesh,par,time);
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate stochastic population dynamics
% the 'type' variable determines of species or pair abundance

rng(10)  % 10,60,350

num=zeros(7,tmesh(end)-1);
num(:,1)=init;
t=0;

number_R=init(1);
number_C1=init(2);
number_C2=init(3);
number_C1R=0;
number_C2R=0;
number_C1C1=0;
number_C2C2=0;
newR=0;   %new born resource
newC1=0;   %new born consumer1
deathC1=0;  %death consumer1
newC2=0;   %new born consumer2
deathC2=0;  %death consumer2


sign=zeros(2,1);
position_C1=randi([-width,width],2,number_C1); %initialize positions of Consumer1
position_C1=[position_C1 sign];

position_C2=randi([-width,width],2,number_C2); %initialize positions of Consumer2
position_C2=[position_C2 sign];

position_R=randi([-width,width],2,number_R);%initialize positions of Resource
position_R=[position_R sign];

position_C1C1=zeros(2,1);  %last colomn means nothing
position_C1R=zeros(2,1);  %last colomn means nothing

position_C2C2=zeros(2,1);  %last colomn means nothing
position_C2R=zeros(2,1);  %last colomn means nothing

index = 1;
position{1,1}=position_R;
position{2,1}=position_C1;
position{3,1}=position_C2;
position{4,1}=position_C1R;
position{5,1}=position_C2R;
position{6,1}=position_C1C1;
position{7,1}=position_C2C2;
ppi=2;
index=index+1
while index<length(tmesh)
    
    %consumer1+consumer1 formation
    %find equal position consumer1 &consumer1
    %if equal -->C1C1 pair
    i=1;
    seq=randperm(number_C1);
    k=1;
    while i<=number_C1
        j=seq(k);
        if j==i && k<number_C1
            k=k+1;
            j=seq(k);
        elseif j==i && k==number_C1
            i=i+1;
            if i>number_C1
                break
            end
            seq=randperm(number_C1);
            k=1;
        end
        if (position_C1(1,i)-position_C1(1,j))^2+(position_C1(2,i)-position_C1(2,j))^2<=par.r^2
            temp=position_C1(:,i);
            position_C1C1=[temp position_C1C1];
            position_C1(:,[i,j])=[];
            number_C1=number_C1-2;
            number_C1C1=number_C1C1+1;
            seq=randperm(number_C1);
            k=1;
        else
            k=k+1;
        end
        if k>number_C1
            i=i+1;
            seq=randperm(number_C1);
            k=1;
        end
    end
    
    
    %consumer1+consumer1 devorce
    p=rand(1,number_C1C1);
    pp=p<=par.d1;   %-->C1+C1
    temp=position_C1C1(:,pp);
    position_C1=[temp position_C1];
    theta=rand(size(temp(1,:)))*2*pi;
    temp=round(temp+par.r*[cos(theta);sin(theta)]);
    position_C1=[temp position_C1];
    position_C1C1(:,pp)=[];
    p(:,pp)=[];
    number_C1=number_C1+2*sum(pp);
    number_C1C1=number_C1C1-sum(pp);
    
    pp=p<=par.d1+par.k1 & p>par.d1;   %-->C1
    temp=position_C1C1(:,pp);
    position_C1=[temp position_C1];
    position_C1C1(:,pp)=[];
    p(:,pp)=[];
    number_C1=number_C1+sum(pp);
    number_C1C1=number_C1C1-sum(pp);
    
    
    %consumer2+consumer2 formation
    %find equal position consumer2 &consumer2
    %if equal -->C2C2 pair
    i=1;
    seq=randperm(number_C2);
    k=1;
    while i<=number_C2
        j=seq(k);
        if j==i && k<number_C2
            k=k+1;
            j=seq(k);
        elseif j==i && k==number_C2
            i=i+1;
            if i>number_C2
                break
            end
            seq=randperm(number_C2);
            k=1;
        end
        if (position_C2(1,i)-position_C2(1,j))^2+(position_C2(2,i)-position_C2(2,j))^2<=par.r^2
            temp=position_C2(:,i);
            position_C2C2=[temp position_C2C2];
            position_C2(:,[i,j])=[];
            number_C2=number_C2-2;
            number_C2C2=number_C2C2+1;
            seq=randperm(number_C2);
            k=1;
        else
            k=k+1;
        end
        if k>number_C2
            i=i+1;
            seq=randperm(number_C2);
            k=1;
        end
    end
    
    
    %consumer2+consumer2 devorce
    p=rand(1,number_C2C2);
    pp=p<=par.dd1;   %-->C2+C2
    temp=position_C2C2(:,pp);
    position_C2=[temp position_C2];
    theta=rand(size(temp(1,:)))*2*pi;
    temp=round(temp+par.r*[cos(theta);sin(theta)]);
    position_C2=[temp position_C2];
    position_C2C2(:,pp)=[];
    p(:,pp)=[];
    number_C2=number_C2+2*sum(pp);
    number_C2C2=number_C2C2-sum(pp);
    
    pp=p<=par.dd1+par.kk1 & p>par.dd1;   %-->C2
    temp=position_C2C2(:,pp);
    position_C2=[temp position_C2];
    position_C2C2(:,pp)=[];
    p(:,pp)=[];
    number_C2=number_C2+sum(pp);
    number_C2C2=number_C2C2-sum(pp);
    
    
    %consumer1+resource formation
    %find equal position consumer1 & resource
    %if equal -->C1R pair
    i=1;
    if number_R~=0
    seq=randperm(number_R);
    k=1;
    while i<=number_C1
        j=seq(k);
        if (position_C1(1,i)-position_R(1,j))^2+(position_C1(2,i)-position_R(2,j))^2<=par.r^2
            temp=position_C1(:,i);
            position_C1R=[temp position_C1R];
            position_C1(:,i)=[];
            position_R(:,j)=[];
            number_C1=number_C1-1;
            number_R=number_R-1;
            number_C1R=number_C1R+1;
            if number_R==0
                break
            end
            seq=randperm(number_R);
            k=1;
        else
            k=k+1;
        end
        if k>number_R
            i=i+1;
            seq=randperm(number_R);
            k=1;
        end
    end
    end
    
    
    %consumer1+resource devorce
    p=rand(1,number_C1R);
    pp=p<=par.d;   %-->C1+R
    temp=position_C1R(:,pp);
    position_C1=[temp position_C1];
    theta=rand(size(temp(1,:)))*2*pi;
    temp=round(temp+par.r*[cos(theta);sin(theta)]);
    position_R=[temp position_R];
    position_C1R(:,pp)=[];
    p(:,pp)=[];
    number_C1=number_C1+sum(pp);
    number_R=number_R+sum(pp);
    number_C1R=number_C1R-sum(pp);
    
    pp=p<=par.d+par.k & p>par.d;   %-->C1+
    temp=position_C1R(:,pp);
    position_C1=[temp position_C1];
    position_C1R(:,pp)=[];
    p(:,pp)=[];
    number_C1=number_C1+sum(pp);
    number_C1R=number_C1R-sum(pp);
    newC1=newC1+par.ww*sum(pp);
    
    
    %consumer2+resource formation
    %find equal position consumer2 & resource
    %if equal -->C2R pair
    i=1;
    if number_R~=0
    seq=randperm(number_R);
    k=1;
    while i<=number_C2
        j=seq(k);
        if (position_C2(1,i)-position_R(1,j))^2+(position_C2(2,i)-position_R(2,j))^2<=par.r^2
            temp=position_C2(:,i);
            position_C2R=[temp position_C2R];
            position_C2(:,i)=[];
            position_R(:,j)=[];
            number_C2=number_C2-1;
            number_R=number_R-1;
            number_C2R=number_C2R+1;
            if number_R==0
                break
            end
            seq=randperm(number_R);
            k=1;
        else
            k=k+1;
        end
        if k>number_R
            i=i+1;
            seq=randperm(number_R);
            k=1;
        end
    end
    end
    
    
    %consumer2+resource devorce
    p=rand(1,number_C2R);
    pp=p<=par.dd;   %-->C2+R
    temp=position_C2R(:,pp);
    position_C2=[temp position_C2];
    theta=rand(size(temp(1,:)))*2*pi;
    temp=round(temp+par.r*[cos(theta);sin(theta)]);
    position_R=[temp position_R];
    position_C2R(:,pp)=[];
    p(:,pp)=[];
    number_C2=number_C2+sum(pp);
    number_R=number_R+sum(pp);
    number_C2R=number_C2R-sum(pp);
    
    pp=p<=par.dd+par.kk & p>par.dd;   %-->C2+
    temp=position_C2R(:,pp);
    position_C2=[temp position_C2];
    position_C2R(:,pp)=[];
    p(:,pp)=[];
    number_C2=number_C2+sum(pp);
    number_C2R=number_C2R-sum(pp);
    newC2=newC2+par.ww*sum(pp);
    
    
    %new born consumer1
    if newC1>=1
        n=floor(newC1);
        newC1=newC1-n;
        temp=randi([-width,width],2,n);
        position_C1=[temp position_C1];
        number_C1=number_C1+n;
    end
    
    
    %new born consumer2
    if newC2>=1
        n=floor(newC2);
        newC2=newC2-n;
        temp=randi([-width,width],2,n);
        position_C2=[temp position_C2];
        number_C2=number_C2+n;
    end
    
    
    %new born resource
    newR=newR+par.R00*(1-number_R/par.K0);
    if newR>=1
        n=floor(newR);
        newR=newR-n;
        temp=randi([-width,width],2,n);
        position_R=[temp position_R];
        number_R=number_R+n;
    end
    
    %dead consumer1
    deathC1=deathC1+par.D*number_C1;
    if deathC1>=1
        n=floor(deathC1);
        deathC1=deathC1-n;
        if n<number_C1
            t=randperm(number_C1);
            t=t(1:n);
            position_C1(:,t)=[];
            number_C1=number_C1-n;
        else
            position_C1(:,1:number_C1)=[];
            number_C1=0;
        end
    end
    
    
    %dead consumer2
    deathC2=deathC2+par.DD*number_C2;
    if deathC2>=1
        n=floor(deathC2);
        deathC2=deathC2-n;
        if n<number_C2
            t=randperm(number_C2);
            t=t(1:n);
            position_C2(:,t)=[];
            number_C2=number_C2-n;
        else
            position_C2(:,1:number_C2)=[];
            number_C2=0;
        end
    end
    
    
    %consumer1 move
    x=randi([0,1],1,number_C1);
    dis=randi([-1,1],1,number_C1);
    temp=find(x==0);
    position_C1(1,temp)=position_C1(1,temp)+dis(temp);
    position_C1(1,find(position_C1(1,:)>width))=-width;
    position_C1(1,find(position_C1(1,:)<-width))=width;
    temp=find(x==1);
    position_C1(2,temp)=position_C1(2,temp)+dis(temp);
    position_C1(2,find(position_C1(2,:)>width))=-width;
    position_C1(2,find(position_C1(2,:)<-width))=width;
    
    
    %consumer2 move
    x=randi([0,1],1,number_C2);
    dis=randi([-1,1],1,number_C2);
    temp=find(x==0);
    position_C2(1,temp)=position_C2(1,temp)+dis(temp);
    position_C2(1,find(position_C2(1,:)>width))=-width;
    position_C2(1,find(position_C2(1,:)<-width))=width;
    temp=find(x==1);
    position_C2(2,temp)=position_C2(2,temp)+dis(temp);
    position_C2(2,find(position_C2(2,:)>width))=-width;
    position_C2(2,find(position_C2(2,:)<-width))=width;
    %
    
    %resource move
    x=randi([0,1],1,number_R);
    dis=randi([-1,1],1,number_R);
    temp=find(x==0);
    position_R(1,temp)=position_R(1,temp)+dis(temp);
    position_R(1,find(position_R(1,:)>width))=-width;
    position_R(1,find(position_R(1,:)<-width))=width;
    temp=find(x==1);
    position_R(2,temp)=position_R(2,temp)+dis(temp);
    position_R(2,find(position_R(2,:)>width))=-width;
    position_R(2,find(position_R(2,:)<-width))=width;
    
    
    %C1R move
    if number_C1R>0
        x=randi([0,1],1,number_C1R);
        dis=randi([-1,1],1,number_C1R);
        temp=find(x==0);
        position_C1R(1,temp)=position_C1R(1,temp)+dis(temp);
        position_C1R(1,find(position_C1R(1,:)>width))=-width;
        position_C1R(1,find(position_C1R(1,:)<-width))=width;
        temp=find(x==1);
        position_C1R(2,temp)=position_C1R(2,temp)+dis(temp);
        position_C1R(2,find(position_C1R(2,:)>width))=-width;
        position_C1R(2,find(position_C1R(2,:)<-width))=width;
    end
    
    
    %C2R move
    if number_C2R>0
        x=randi([0,1],1,number_C2R);
        dis=randi([-1,1],1,number_C2R);
        temp=find(x==0);
        position_C2R(1,temp)=position_C2R(1,temp)+dis(temp);
        position_C2R(1,find(position_C2R(1,:)>width))=-width;
        position_C2R(1,find(position_C2R(1,:)<-width))=width;
        temp=find(x==1);
        position_C2R(2,temp)=position_C2R(2,temp)+dis(temp);
        position_C2R(2,find(position_C2R(2,:)>width))=-width;
        position_C2R(2,find(position_C2R(2,:)<-width))=width;
    end
    
    
    %C1C1 move
    if number_C1C1>0
        x=randi([0,1],1,number_C1C1);
        dis=randi([-1,1],1,number_C1C1);
        temp=find(x==0);
        position_C1C1(1,temp)=position_C1C1(1,temp)+dis(temp);
        position_C1C1(1,find(position_C1C1(1,:)>width))=-width;
        position_C1C1(1,find(position_C1C1(1,:)<-width))=width;
        temp=find(x==1);
        position_C1C1(2,temp)=position_C1C1(2,temp)+dis(temp);
        position_C1C1(2,find(position_C1C1(2,:)>width))=-width;
        position_C1C1(2,find(position_C1C1(2,:)<-width))=width;
    end
    
    
    %C2C2 move
    if number_C2C2>0
        x=randi([0,1],1,number_C2C2);
        dis=randi([-1,1],1,number_C2C2);
        temp=find(x==0);
        position_C2C2(1,temp)=position_C2C2(1,temp)+dis(temp);
        position_C2C2(1,find(position_C2C2(1,:)>width))=-width;
        position_C2C2(1,find(position_C2C2(1,:)<-width))=width;
        temp=find(x==1);
        position_C2C2(2,temp)=position_C2C2(2,temp)+dis(temp);
        position_C2C2(2,find(position_C2C2(2,:)>width))=-width;
        position_C2C2(2,find(position_C2C2(2,:)<-width))=width;
    end
    
    if index==time(ppi)
        position{1,ppi}=position_R;
        position{2,ppi}=position_C1;
        position{3,ppi}=position_C2;
        position{4,ppi}=position_C1R;
        position{5,ppi}=position_C2R;
        position{6,ppi}=position_C1C1;
        position{7,ppi}=position_C2C2;
        ppi=ppi+1;
    end
    
    num(:,index)=[number_R;number_C1;number_C2;number_C1R;number_C2R;number_C1C1;number_C2C2];
    
    index=index+1;
end
