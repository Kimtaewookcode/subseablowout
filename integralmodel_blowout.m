%% Integral plume model-MJ Fridedl
% Authors: Kim Taewook
% Linkedin : https://www.linkedin.com/in/taewook-kim/
% GITHUB : github.com/Kimtaewookcode
% Email : kimtaewook87@gmail.com
% Based on  
% http://www.sciencedirect.com/science/article/pii/S014111879900022X
clear
close all
clc
g = 9.81;%gravity[m/s2]
rho_w=998;%density of water[kg/s]
m_release=0.208;%0.71milgram%0.208;%mass rate[kg/s]
phi=3.14;%phi
p0=101325;%atmopheric pressure[pa]
rho_g0=1.225;%density of releasing gas[kg/m3]
vs=0.35;%slip velocity[m/s]
h=7;%height of the water surface
hp=10.33;
x1 = 0;%centerline
hoff=0;%h offset
z0=1.75;%z value where you want to know about profile at
%coefficients%
alpha=0.1285;%0.13;%0.1285;entrainment coefficient
gam=1.5;%1.5;
beta=0.39;%0.5;theoretical value%0.39;experimental value : for a loss-free rise;1for instantneous
lambda=0.8;

%%calculation in x(radial) direction%%
vs1=vs*((g*m_release*(lambda^2+1))/((h+hp)*phi*gam*rho_g0*2*alpha^2))^(-1/3);
s1=(1+lambda^2)*vs1;
z1=z0/(h+hp);
z2=h/(h+hp);
v1=((25/12)^(1/3))*(z1^(-1/3))*(1+11*z1/39+511/2*(z1/39)^2)-s1*7/22*(1+345/343*z1/13+86175/11662*(z1/13)^2)+(s1^2)*13/121*(12/25)^(1/3)*(z1^(1/2))*(1-59489/1436*z1/39-2825583625/23347324*(z1/39)^2);
v2=((25/12)^(1/3))*(z2^(-1/3))*(1+11*z2/39+511/2*(z2/39)^2)-s1*7/22*(1+345/343*z2/13+86175/11662*(z2/13)^2)+(s1^2)*13/121*(12/25)^(1/3)*(z2^(1/2))*(1-59489/1436*z2/39-2825583625/23347324*(z2/39)^2);
v=v1*((g*m_release*(lambda^2+1))/((h+hp)*phi*gam*rho_g0*2*alpha^2))^(1/3);
vh=v2*((g*m_release*(lambda^2+1))/((h+hp)*phi*gam*rho_g0*2*alpha^2))^(1/3);
b1=3/5*z1*(1-z1/13-7*(z1/13)^2)+s1*3/110*((12/25)^(1/3))*(z1^(4/3))*(1-1046/49*z1/39-227726/833*(z1/39)^2)-(s1^2)*48/15121*((25/12)^(1/3))*(z1^(5/3))*(1-34663/9408*z1+225707803/240143904*z1^2);
b2=3/5*z2*(1-z2/13-7*(z2/13)^2)+s1*3/110*((12/25)^(1/3))*(z2^(4/3))*(1-1046/49*z2/39-227726/833*(z2/39)^2)-(s1^2)*48/15121*((25/12)^(1/3))*(z2^(5/3))*(1-34663/9408*z2+225707803/240143904*z2^2);
b=b1*(2*alpha*(h+hp));%plume width at z0
bh=b2*(2*alpha*(h+hp));%plume width at the surface
x = [0 : 0.01: 200];%x range- need to be set
y=v*exp(-(x.^2)/(b^2));%velocity at z0
yh=vh*exp(-(x1^2)/(b^2));%velocity at the surface

hf=beta*gam*(yh.^2)/g;%the peak of fountain profile
hr=hf*exp(-(x.^2)/(bh^2))-hoff;%fountain profile

void_1=1/(1-z1)/((b1^2)*(v1+s1));
void1=void_1*((((1+lambda^2)^2)*gam*(m_release/rho_g0)^2)/((phi^2)*(lambda^2)*(2^5)*(alpha^4)*((h+hp)^5)*g))^(1/3);
voidr=void1*exp(-(x.^2)/((lambda^2)*(b.^2)));%void fraction profile at z0

void_2=1/(1-z2)/((b2^2)*(v2+s1));
void2=void_2*((((1+lambda^2)^2)*gam*(m_release/rho_g0)^2)/((phi^2)*(lambda^2)*(2^5)*(alpha^4)*((h+hp)^5)*g))^(1/3);
voidh=void2*exp(-(x.^2)/((lambda^2)*(b^2)));%void fraction profile at the surface
mflux_surf=voidh*yh*rho_g0;%mass flux per unit area at the surface

plot(x, y, '.-')
xlim([0 4])
title(['velocity at z=',num2str(z0)]);
xlabel('radial position[m]');
ylabel('velocity[m/s]');

figure;plot(x, hr, '.-')
xlim([0 4])
title('fountain ');
xlabel('radial position[m]');
ylabel('z[m]');

figure;plot(x, voidr, '.-')
xlim([0 4])
title(['void fraction at z=',num2str(z0)]);
xlabel('radial position[m]');
ylabel('voidfraction');

figure;plot(x, mflux_surf, '.-')
xlim([0 4])
title(['mass flux at the surface']);
xlabel('radial position[m]');
ylabel('massflux[kg/s/m2]');

%%%calculation in z direction%%
initial=0;%z start
last=h;%z end
dz=0.01;%delta z
n=(last-initial)/dz;%the number of data
zarr=zeros(n,1);% n x 1 array for z
veloarr=zeros(n,1);%n x 1 array for velocity in z direction
voidarr=zeros(n,1);%n x 1 array for voidfraction in z direction
i=1;

for z=initial:dz:last
zarr(i)=z;
z1=z/(h+hp);
vs1=vs*((g*m_release*(lambda^2+1))/((h+hp)*phi*gam*rho_g0*2*alpha^2))^(-1/3);
s1=(1+lambda^2)*vs1;
v1_z=((25/12)^(1/3))*(z1.^(-1/3))*(1+11*z1/39+511/2*(z1/39).^2)-s1*7/22*(1+345/343*z1/13+86175/11662*(z1/13).^2)+(s1^2)*13/121*(12/25)^(1/3)*(z1.^(1/2))*(1-59489/1436*z1/39-2825583625/23347324*(z1/39).^2);
v_z=v1_z*((g*m_release*(lambda^2+1))/((h+hp)*phi*gam*rho_g0*2*alpha^2))^(1/3);
b1_z=3/5*z1*(1-z1/13-7*(z1/13)^2)+s1*3/110*((12/25)^(1/3))*(z1^(4/3))*(1-1046/49*z1/39-227726/833*(z1/39)^2)-(s1^2)*48/15121*((25/12)^(1/3))*(z1^(5/3))*(1-34663/9408*z1+225707803/240143904*z1^2);
b_z=b1_z*(2*alpha*(h+hp));
bz(i)=b_z;%designate br's 'i' th value
veloarr_z(i)=v_z*exp(-(x1^2)/(b_z^2));%designate veloarr's 'i' th value

void1_z=1/(1-z1)/((b1_z.^2)*(v1_z+s1));
void_z=void1_z*((((1+lambda^2)^2)*gam*(m_release/rho_g0)^2)/((phi^2)*(lambda^2)*(2^5)*(alpha^4)*((h+hp)^5)*g))^(1/3);
voidz(i)=void_z*exp(-(x1^2)/((lambda^2)*(b_z^2)));%designate voidr's 'i' th value


 i=i+1;
end
figure;plot(bz,zarr, '.-')
ylim([0 h])
title('plume width');
ylabel('z position[m]');
xlabel('width[m]');
figure;plot(veloarr_z,zarr, [1 1],'.-','.')
xlim([0 4])
ylim([0 h])
title('velocity at the centerline ');
ylabel('z position[m]');
xlabel('velocity[m/s]');
figure;plot(voidz,zarr, '.-')
xlim([0 0.5])
ylim([0 h])
title('voidfraction at the centerline');
ylabel('z position[m]');
xlabel('voidfraction');
instantneous_maximum_fountain=3.1*hf;
instantneous_average_fountain=2*hf;