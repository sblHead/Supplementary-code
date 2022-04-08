%%
%-------TMM----Based on a code by S. Rao, adapted by S. Head-------------
%----REF: Sathyanarayan Rao, Transmittance and Reflectance Spectra of---- 
%-------Multilayered Dielectric Stack using Transfer Matrix Method------- 
%--------------------MATLAB Central File Exchange------------------------
%------------------------------------------------------------------------ 
% clear all
% close all
% clc
%------------------------------------------------------------------------ 
%%

%Load / import data
params = importdata('');
% columns of n and k data for each material between 380 and 780 nm
struct = importdata('');
% rows containing material depths for each layer [d1,d2,d3...,dN]

total_runs = length(struct); 

%Data Output
len = length(380:5:780);
Out = zeros(len,total_runs);

% SOURCE PARAMETERS
lam0=(380:5:780)*1e-9;%free space wavelength
k0=(2*pi)./lam0;
pte = 1/sqrt(2); %amplitude of TE polarization
ptm = 1i*pte; %amplitude of TM polarization
ni=1.0; % incident medium refractive index
Ur0 = 1;
Er0 = 1;
Urt = 1;
Ert = 1;


place = 1;
for r=1:total_runs 
    
    N=6;   % number of material layers
    materials = struct(r,:);
    L = struct(r,:)*1e-7;
    Ur=ones(1,N);
    NR=zeros(len,N);
    NI=zeros(len,N);
    
    count = 1;
    for c=1:N
        nR = (params(:,(materials(count)*2)-1));
        NR(:,count)=nR;
        nI = (params(:,(materials(count)*2)));
        NI(:,count)=nI;
        count = count+1; 
    end
    ER=(NR+i*NI).^2;  % gives a vector of complex permitivities
    
    % Free Space Parameters
    Kz0=1;
    Ohm0=1i*Kz0*eye(2);
    Q0 = [0 1; -1 0];
    P0 = [0 1; -1 0];
    V0 = Q0*(Ohm0^-1);
    W0 = eye(2);
    
    
    for l=1:length(lam0)
        % Initialize global scattering matrix
        Sg11=zeros(2,2); Sg12=eye(2); Sg21=eye(2); Sg22=zeros(2,2);
        
        % medium 1
        Kz1=sqrt(Ur0*Er0);
        Ohm1=1i*Kz1*eye(2);
        Q1 = (1/Ur0)*[0 (Ur0*Er0); -(Ur0*Er0) 0];
        V1 = Q1*(Ohm1^-1);
        W1 = eye(2);
        A1=eye(2)+(V0^-1)*V1;
        B1=eye(2)-(V0^-1)*V1;
        
        Sr11 = -(A1^-1)*B1;
        Sr12 = 2*(A1^-1);
        Sr21 = 0.5*(A1-(B1*(A1^-1)*B1));
        Sr22 = B1*(A1^-1);
        
        Sg11=Sr11;
        Sg12=Sr12;
        Sg21=Sr21;
        Sg22=Sr22;
        
        % For all central layers
        for i=1:N
            Kz=sqrt(Ur(i)*ER(l,i));
            Q=(1/Ur(i))*[ 0 Ur(i)*ER(l,i) ;-Ur(i)*ER(l,i) 0];
            Ohm=1i*Kz*eye(2);
            V=Q*(Ohm^-1);
            A=eye(2)+((V^-1)*V0);
            B=eye(2)-((V^-1)*V0);
            X=expm(Ohm*k0(l)*L(i));

            S11=((A-((X*B*(A^-1)*X*B)))^-1)*((X*B*(A^-1)*X*A)-B);
            S22=S11;
            S12=(((A-((X*(B/A)*X*B)))^-1)*X)*(A-B*(A^-1)*B);
            S21=S12;
            
            % updating global scattering matrix
            SAB11=Sg11+(Sg12*((eye(2)-(S11*Sg22))^-1)*S11*Sg21);
            SAB12=Sg12*((eye(2)-(S11*Sg22))^-1)*S12;
            SAB21=S21*((eye(2)-(Sg22*S11))^-1)*Sg21;
            SAB22=S22+(S21*((eye(2)-(Sg22*S11))^-1)*Sg22*S12);
            
            Sg11=SAB11;
            Sg12=SAB12;
            Sg21=SAB21;
            Sg22=SAB22;
            
        end
        
       % medium 2    
        Kz2=sqrt(Urt*Ert);
        Ohm2=1i*Kz2*eye(2);
        Q2 = (1/Urt)*[0 Urt*Ert; -(Urt*Ert) 0];
        V2 = Q2*(Ohm1^-1);
        W2 = eye(2);
        A2=eye(2)+(V0^-1)*V2;
        B2=eye(2)-(V0^-1)*V2;
        
        
        St11 = B2*(A2^-1);
        St12 = 0.5*(A2-(B2*(A2^-1)*B2));
        St21 = 2*(A2^-1);
        St22 = -(A2^-1)*B2;
        
        % updating global scattering matrix
        SAB11=Sg11+(Sg12*((eye(2)-(St11*Sg22))^-1)*St11*Sg21);
        SAB12=Sg12*((eye(2)-(St11*Sg22))^-1)*St12;
        SAB21=St21*((eye(2)-(Sg22*St11))^-1)*Sg21;
        SAB22=St22+(St21*((eye(2)-(Sg22*St11))^-1)*Sg22*St12);
        
        Sg11=SAB11;
        Sg12=SAB12;
        Sg21=SAB21;
        Sg22=SAB22;
        
        erf = Sg11*[pte;ptm];
        etf = Sg21*[pte;ptm];
        erx=erf(1);
        ery=erf(2);
        etx=etf(1);
        ety=etf(2);

        
        R=abs(erx)^2+abs(ery)^2;
        Rx(l)=abs(R);   % Refelctance value
        
    end
    Out(1:len,place)=Rx;
    place = place + 1;
    
end
