%keplerian turbulence synthetic model (optimized)

%wave number parameters
N=512;
kmin=5;
kmax=50;

%cell parameters
M=64;
Lx=2;
Ly=2;
Lz=2;

for j=1:M
    x(j)=-Lx/2+(j-1)*Lx/(M-1);
    y(j)=-Ly/2+(j-1)*Ly/(M-1);
    z(j)=-Lz/2+(j-1)*Lz/(M-1);
end

%wave vector amplitudes
kamp(1)=kmin;
kamp(N)=kmax;
for j=2:(N-1)
    kamp(j)=kamp(1)*(kamp(N)/kamp(1))^((j-1)/(N-1));
end

%energy dissipation
eps=(1.5/((kmin)^(-2/3)-(kmax)^(-2/3)))^(3/2);

%Kolmogorov energy spectrum and frequencies
for j=1:N
    E(j)=1.5*(eps)^(2/3)*kamp(j)^(-5/3);
    w(j)=0.5*sqrt(kamp(j)^3*E(j));
end

%random direction of the wave vector
for j=1:N
fi=2*pi*rand(1);
tetta=pi-pi*rand(1);
kedx(j)=sin(tetta)*sin(fi);
kedy(j)=sin(tetta)*cos(fi);
kedz(j)=cos(tetta);
kx(j)=kedx(j)*kamp(j);
ky(j)=kedy(j)*kamp(j);
kz(j)=kedz(j)*kamp(j);
end

%poloidal-toroidal-divergent decomposition (Craya-Herring frame)
for j=1:N
    alfa=2*pi*rand(1);
    fi=2*pi*rand(1);
    
    %velocity components in this decomposition
    u1(j)=sqrt(E(j)/(4*pi*(kamp(j))^2))*cos(alfa)*exp(1i*fi);
    u2(j)=sqrt(E(j)/(4*pi*(kamp(j))^2))*sin(alfa)*exp(1i*fi);
    
    %the unit vectors connected to toroidal and poloidal modes
    e1x(j)=ky(j)/sqrt(kx(j)^2+ky(j)^2);
    e1y(j)=-kx(j)/sqrt(kx(j)^2+ky(j)^2);
    e1z(j)=0;
    
    e2x(j)=(kx(j)*kz(j))/(sqrt(kx(j)^2+ky(j)^2)*kamp(j));
    e2y(j)=(ky(j)*kz(j))/(sqrt(kx(j)^2+ky(j)^2)*kamp(j));
    e2z(j)=-sqrt(kx(j)^2+ky(j)^2)/kamp(j);
    
    %velocity amplitudes in Fourier space
    vampx(j)=(u1(j)*e1x(j)+u2(j)*e2x(j));
   vampy(j)=(u1(j)*e1y(j)+u2(j)*e2y(j));
    vampz(j)=(u1(j)*e1z(j)+u2(j)*e2z(j));
    end

%synthetic velocity field in real space
dt=0.01;
im(1,1,1,501) = 0; 
hF=figure; 

clear kamp kedx kedy kedz u1 u2 e1x e1y e1z e2x e2y e2z;

for m=1:M
    for l=1:M
        for k=1:M
            for j=1:N
vx1(m,l,k,j)=vampx(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddxvx1(m,l,k,j)=1i*kx(j)*vampx(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddyvx1(m,l,k,j)=1i*ky(j)*vampx(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddzvx1(m,l,k,j)=1i*kz(j)*vampx(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));

            end
        end
    end
    m
end

vx=real(sum(vx1,4));
ddxvx=real(sum(ddxvx1,4));
ddyvx=real(sum(ddyvx1,4));
ddzvx=real(sum(ddzvx1,4));
clear vx1 ddxvx1 ddyvx1 ddzvx1;

isocaps(x,y,z,vx,-10)
view(20,20);
title(['Ux(x,y,z,0)']);
colorbar('vert');
caxis([-0.8, 0.8]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;
clear vx1;

for m=1:M
    for l=1:M
        for k=1:M
            for j=1:N

vy1(m,l,k,j)=vampy(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddxvy1(m,l,k,j)=1i*kx(j)*vampy(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddyvy1(m,l,k,j)=1i*ky(j)*vampy(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddzvy1(m,l,k,j)=1i*kz(j)*vampy(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));

            end
        end
    end
    m
end

%vy=2.0*real(sum(vy1,4))+u0;
vy=real(sum(vy1,4));
ddxvy=real(sum(ddxvy1,4));
ddyvy=real(sum(ddyvy1,4));
ddzvy=real(sum(ddzvy1,4));

hF=figure;
isocaps(x,y,z,vy,-10)
view(20,20);
title(['Uy(x,y,z)']);
colorbar('vert');
caxis([-1.5, 1.5]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;
clear vy1 ddxvy1 ddyvy1 ddzvy1;

for m=1:M
   for l=1:M
       for k=1:M
           for j=1:N
vz1(m,l,k,j)=vampz(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddxvz1(m,l,k,j)=1i*kx(j)*vampz(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddyvz1(m,l,k,j)=1i*ky(j)*vampz(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
ddzvz1(m,l,k,j)=1i*kz(j)*vampz(j)*exp(1i*(kx(j)*x(m)+ky(j)*y(l)+kz(j)*z(k)));
            end
        end
    end
    m
end

vz=real(sum(vz1,4));
ddxvz=real(sum(ddxvz1,4));
ddyvz=real(sum(ddyvz1,4));
ddzvz=real(sum(ddzvz1,4));
clear vz1 ddxvz1 ddyvz1 ddzvz1;

hF=figure; 
isocaps(x,y,z,vz,-10)
view(20,20);
title(['Uz(x,y,z,0)']);
colorbar('vert');
caxis([-0.8, 0.8]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;

%kinetic energy distribution at the physical space
for m=1:M
   for l=1:M
        for k=1:M
            EE(m,l,k)=(1/2)*(vx(m,l,k)^2+vy(m,l,k)^2+vz(m,l,k)^2);
        end
    end
end

for j=1:M
    for l=1:M
        for m=1:M
Density(j,l,m)=0.02*(ddxvx(j,l,m)^2+ddyvy(j,l,m)^2+ddzvz(j,l,m)^2+2*(ddxvy(j,l,m)*ddyvx(j,l,m)+ddxvz(j,l,m)*ddzvx(j,l,m)+ddyvz(j,l,m)*ddzvy(j,l,m))-3*ddyvx(j,l,m));
        end
    end
end

hF=figure; 
isocaps(x,y,z,Density,-200)
view(20,20);
title(['Density(x,y,z)']);
colorbar('vert');
caxis([-30, 30]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;

for j=1:M
    z0(j)=0;
end
 
hF=figure; 
isocaps(x,y,z0,Density,-200)
view(0,90);
title(['Density(x,y,z)']);
colorbar('vert');
caxis([-30, 30]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;
