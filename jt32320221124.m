%%%%%%% decision-making algorithm
clc
clear
%%%%%%% A1 and A2 represent the membership degrades and nonmembership degrades, repectively.
  A1=[   ];

  A2=[  ];
  
%%%%%%% The values of two thresholds
vv1=0.9;
vv2=0.1;
vv3=0.45;
vv4=0;
[m,n]=size(A1);
aaaa1=[];aaaa2=[];
for i=1:m
    aaa1=[];aaa2=[];
    for j=1:m
        aa1=[];aa2=[];
        for o=1:n
            a=A1(i,o);b=A1(j,o);c=A2(i,o);d=A2(j,o);
            A=[1;1+b-a;1+c-d];
            a1=min(A);
            a2=max(0,d-c);
            aa2=[aa2;a2];aa1=[aa1;a1];
        end
        a11=min(aa1);a22=max(aa2);
        aaa1=[aaa1;a11];aaa2=[aaa2;a22];
    end
    aaaa1=[aaaa1;aaa1'];aaaa2=[aaaa2;aaa2']; 
end
C1=[];C2=[];C3=[];C4=[];
for ii=1:m
    c11=[];c22=[];c33=[];c44=[];
    for jj=1:n
        c1=max(A1(ii,jj),vv1);
        c2=min(A2(ii,jj),vv2);
        c3=min(A1(ii,jj),vv3); 
        c4=max(A2(ii,jj),vv4); 
        c11=[c11;c1];c22=[c22;c2];c33=[c33;c3];c44=[c44;c4];
    end
    C1=[C1;c11'];C2=[C2;c22'];C3=[C3;c33'];C4=[C4;c44'];
end
[nc,mc]=size(C1);
dd=ones(nc,1);ddd=zeros(nc,1);
D1=[];D2=[];D3=[];D4=[];
for oo=1:mc
    d11=[];d22=[];d33=[];d44=[];
    for p=1:m
        b1=(dd+C1(:,oo)-aaaa1(p,:)');
        b2=(dd+aaaa2(p,:)'-C2(:,oo));
        B1=[dd,b1,b2];
        dd1=min(min(B1));
        d11=[d11;dd1];
        b3=(C2(:,oo)-aaaa2(p,:)');
        B2=[ddd,b3];
        dd2=max(max(B2));
        d22=[d22;dd2];
        b4=(C3(:,oo)+aaaa1(p,:)'-dd);
        B3=[ddd,b4];
        dd3=max(max(B3));
        d33=[d33;dd3];
        b5=(aaaa2(p,:)'+C4(:,oo));
        B4=[dd,b5];
        dd4=min(min(B4));
        d44=[d44;dd4];
    end
    D1=[D1,d11];D2=[D2,d22];D3=[D3,d33];D4=[D4,d44]; 
end
D1=D1';D2=D2';D3=D3';D4=D4';
[md,nd]=size(D1);
E1=[];E2=[];E3=[];E4=[];
for e1=1:md
    E11=[];E22=[];E33=[];E44=[];
    for e2=1:nd
        ee1=aaaa1(:,e2)+D1(e1,:)'-dd;
        EE1=[ddd,ee1];
        eee1=max(max(EE1));
        E11=[E11,eee1];
        ee2=aaaa2(:,e2)+D2(e1,:)';
        EE2=[dd,ee2];
        eee2=min(min(EE2));
        E22=[E22,eee2];
        ee4=-aaaa1(:,e2)+D3(e1,:)'+dd;
        ee5=aaaa2(:,e2)-D4(e1,:)'+dd;
        EE3=[dd,ee4,ee5];
        eee3=min(min(EE3));
        E33=[E33,eee3];
        ee6=D4(e1,:)'-aaaa2(:,e2);
        EE4=[ddd,ee6];
        eee6=max(max(EE4));
        E44=[E44,eee6];
    end
    E1=[E1;E11];E2=[E2;E22];E3=[E3;E33];E4=[E4;E44]; 
end
F1=E1+E3-E1.*E3; 
F2=E2.*E4; 
[mf,nf]=size(F1);
G1=[];G2=[];
for g1=1:mf
    G22=[];G11=[];
    for g2=1:nf
        gg1=F1(g1,g2).*ones(1,nf);
        ggg1=gg1.*F2(g1,:);
        G11=[G11;ggg1];
    end
       G1=[G1;G11];  
end
p1=0.5;% The value of preference 
[mg,ng]=size(G1);
H=[];
for h1=1:mg
    H1=[];
    for h2=1:ng
        if G1(h1,h2)<=p1
            h11=G1(h1,h2)/p1;
        else
            h11=1;
        end
        H1=[H1,h11];
    end
    H=[H;H1]; 
end

PD1=H(1:6,1:6);
PD2=H(7:12,1:6);
PD3=H(13:18,1:6);
PD4=H(19:24,1:6);
PD5=H(25:30,1:6);

[nRows,nCols] = size(PD1);
PD1(1:(nRows+1):nRows*nCols)=0.5;

[nRows,nCols] = size(PD2);
PD2(1:(nRows+1):nRows*nCols)=0.5;

[nRows,nCols] = size(PD3);
PD3(1:(nRows+1):nRows*nCols)=0.5;

[nRows,nCols] = size(PD4);
PD4(1:(nRows+1):nRows*nCols)=0.5;

[nRows,nCols] = size(PD5);
PD5(1:(nRows+1):nRows*nCols)=0.5;

pd1=PD1';
pd2=PD2';
pd3=PD3';
pd4=PD4';
pd5=PD5';
W=[0.11 0.24 0.30 0.24 0.11]; % The weight vector of 5 attributes
[mk,nk]=size(PD1);
lg=ones(mk,nk);
KK1=[];KK2=[];

for ll=1:mk
  for  nn=1:nk
      kk1(ll,nn)=1-((1-PD1(ll,nn)).^W(1,1)).*((1-PD2(ll,nn)).^W(1,2)).*((1-PD3(ll,nn)).^W(1,3)).*((1-PD4(ll,nn)).^W(1,4)).*((1-PD5(ll,nn)).^W(1,5));
      kk2(ll,nn)=(pd1(ll,nn).^W(1,1)).*(pd2(ll,nn).^W(1,2)).*(pd3(ll,nn).^W(1,3)).*(pd4(ll,nn).^W(1,4)).*(pd5(ll,nn).^W(1,5));
  end
end
KK1=[KK1;kk1];KK2=[KK2;kk2];

[mk,nk]=size(KK1);
RR=ones(mk,nk);
RR1=RR-KK1;
R1=cumprod(RR1,2);R3=cumprod(RR1);
Rr1=ones(mk,1)-R1(:,end);Rr3=ones(1,nk)-R3(end,:);
R2=cumprod(KK2,2);R4=cumprod(KK2);
Rr2=R2(:,end);Rr4=R4(end,:);
[mr,nr]=size(Rr1);
TT=ones(mr,nr);
T1=TT-(TT-Rr1).^(0.2); 
T2=(Rr2).^(0.2); 
T3=TT-(TT-Rr3').^(0.2); 
T4=(Rr4').^(0.2); 
TT1=max(T1);TT2=min(T2); 
TT3=min(T3);TT4=max(T4); 
Z1=sqrt(0.5*((T1-TT1.*TT).^2+(T2-TT2.*TT).^2+0.5*((T1+T2)-(TT1+TT2).*TT).^2)); 
Z2=sqrt(0.5*((T3-TT3.*TT).^2+(T4-TT4.*TT).^2+0.5*((T3+T4)-(TT3+TT4).*TT).^2)); 
ZZ1=Z2./(Z1+Z2); 
A=ZZ1'
[qn,wn]=sort(A,'descend') % The ranking order of alternatives











