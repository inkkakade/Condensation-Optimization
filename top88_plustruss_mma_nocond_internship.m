function x=top88_plustruss_mma_nocond_internship(nelx,nely,volfrac,penal,rmin,ft)
close all
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
P=4;
Sl=1;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
D=E0/(1-nu^2)*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
B=1/2*[-1 0 1 0 1 0 -1 0;0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
DB=D*B;
B=1/2*[-1 0 1 0 1 0 -1 0;0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
DB=D*B;
Cvm=[1 -0.5 0;-0.5 1 0;0 0 3];
Sel=DB'*Cvm*DB;
DE = (1/(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
fixeddofs = [2*(fix(nelx/2)*(nely+1)+1:(nely+1):((nelx)*(nely+1)+1))-1;2*(fix(nelx/2)*(nely+1)+1:(nely+1):((nelx)*(nely+1)+1))];
fixeddofs=fixeddofs(:);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
A=1;                    %Initializing area decision variable
I=100;                  %Initializing inertia desicion variable
x = repmat(1,nely,nelx); %vector of decision variables used for topology optimization
xPhys = x;
loop = 0;
change = 1;
%All vectors and constants needed for MMA Optimization
m = 1;
n = length(xPhys(:));
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = xPhys(:);
xold1   = xval;
xold2   = xval;
xmin    = zeron;
xmax    = eeen;
low     = xmin;
upp     = xmax;
C       = 1000*eeem;
d       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 100000;
kkttol  = 0.001;
%%%% The iterations start:
kktnorm = kkttol+10;
% kktnorm = kkttol;
outit = 0;
change=1;
I_plot = [];
compliance_plot = [];
A_plot = [];
i = 1;  
%% START ITERATION
while kktnorm > kkttol && outit < maxoutit && change>0.001
    outit   = outit+1;
    outeriter = outeriter+1;
    if ft == 1
        xPhys(:)=xval;
    elseif ft == 2
        xPhys(:) = (H*xval(:))./Hs;
    end
    %% Generate the truss stiffness matrix
    [Kcc,Kce,Kee,Fe,Fc,Recovery_matrixc,Recovery_matrixe]=truss_stiffness_no_condensation(nelx,nely,E0,A,I);
%% Generate the Projection matrix, load vector and stiffness matrix ready for assembly
    coupling_nodes = [nely+1:nely+1:(nely+1)*(nelx+1)];
    coupling_dofs=[coupling_nodes*2-1;coupling_nodes*2];
    coupling_dofs=coupling_dofs(:);
    Pr=sparse(coupling_dofs,1:2*(nelx+1),ones(1,2*(nelx+1)),2*(nelx+1)*(nely+1),2*(nelx+1));
    Kt=[Pr*Kcc*Pr' Pr*Kce;Kce'*Pr' Kee];
    Kt=(Kt+Kt')/2;
    Ft=[Pr*Fc;Fe];
    
    Rt=[Recovery_matrixc*Pr' zeros(size(Recovery_matrixc,1),size(Recovery_matrixe,2))
        zeros(size(Recovery_matrixe,1),size(Pr',2)) Recovery_matrixe];
    
    U = zeros(2*(nely+1)*(nelx+1)+length(Fe),1);
    F = U;
    Lambda = [U,U];
    F=F+Ft;
    alldofs = [1:length(F)];
    freedofs = setdiff(alldofs,fixeddofs);  
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K0 = sparse(iK,jK,sK,size(Kt,1),size(Kt,2)); K=K0+Kt;K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    %Compliance of DZ only
    dv = ones(nely,nelx);
    %Total Compliance
    c=F'*U;
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    %% Gather info for MMA
    f0val = c;
    fval=(mean(xPhys(:))-volfrac)/volfrac;%;1*(-0.38-FAN_Ax_disp)/0.38
    df0dx=dc(:);
    dfdx = dv(:)'/length(dv(:))/volfrac;%;-1*dfandisp(:)'/0.38
    innerit=0;
    outvector1 = [outeriter innerit f0val fval'];
    outvector2 = xval;
    %% MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,C,d);
    last = x(nely, :); %Option 1: Reduction of A and I based on the
                       %density of the DZ near interface. 
    A = 1;
    I = 100;
    A = sum(sum(x(nely-5:nely,:)))*A/(5*nelx);
    I = sum(sum(x(nely-5:nely,:)))*I/(5*nelx);   
    compliance_plot(i) = c;
    A_plot(i) = A;
    I_plot(i) = I;
    i = i+1;
    xold2 = xold1;
    xold1 = xval;
    xval  = x(:);
    change = norm(xval-xold1);
    %% Area Optimization %Option 3: Basaed on fmincon 'stress based' optimization.
    
% %     s = (U(edofMat)*(DE*B)').*repmat(E',1,3);
% %     vms = reshape(sqrt(sum(s.^2,2)-s(:,1).*s(:,2)+2.*s(:,3).^2),nely,nelx); %von Mises Stress in each element from top88PTO
% %     
% %     vmsinterface = vms(nely,:);
% %     vmstotal = sum(vmsinterface); %Stress in the elements at the interface is assumed to be stress in the non design zone beam.
% %     
% %     xstart = [sqrt(Area) 0.5*(sqrt(Area)) vmstotal]; %assuming a beam of square cross section
% %     Acons = [-1 0 0; 0 -1 0; 0 0 0]; %A stress based non-linear constraint taking into account Bending and Normal Stress for Aluminium
% %     bcons = [0; 0];
% %     Aeqcons = [0.5 -1 0];
% %     beqcons = [0];
% %     options = optimoptions(@fmincon,'Algorithm','sqp');
% %     [xarea,fval] = fmincon(@objfun,xstart,Acons,bcons,Aeqcons,beqcons,[],[],... 
% %    @confuneq,options);
% % 
% %     Area Update
% %     Area = (xarea(1)*2*xarea(2))*10^9;%Forced the order of magnitude to be correct as I am not sure about the initial system of units
% %     I = (xarea(1)^4)/12*10^22;
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f kktnorm.:%7.3f ch.:%7.3f A.:%7.3f\n',outit,c, ...
        mean(xPhys(:)),kktnorm,change,A);
    figure(2)
    hold on
    plot(outit,c,'bo','MarkerFaceColor','b')
    plot(outit,mean(xPhys(:))*100,'ro','MarkerFaceColor','r')
    plot(outit,A*100,'go','MarkerFaceColor','g')
    title(['Convergence V = ',num2str(mean(xPhys(:))*100),',compliance =',num2str(c),', iter = ', num2str(outit),', Area = ', num2str(A)])
    grid on
    legend('compliance','V %', 'Area')
    xlabel('iter')
    %% %% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,x(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,d);
    outvector1 = [outeriter innerit f0val fval(:)'];
    outvector2 = xval;
    
    %% PLOT DENSITIES
    h = figure(1); set(h,'Color',[1 1 1]);
    [Yy,Xx]=find(nodenrs);
    Yy=nely+1-Yy;
    Xx=Xx-1;
    hold on; patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; drawnow;hold off
    title(['Design zone at iteration ',num2str(outit)])
    colormap(gray);
end
%%
figure
hold on
plot(A_plot*100,compliance_plot, '-r')
xlabel('Area (Scaled)')
ylabel('Compliance')
title('Sensitivty Analysis of A w.r.t. Compliance')
