%% 2D coninjection with random walk model
% Random walk model's isotropic behaviour is applied in a cone injection.
% The injection describes the cone injection that is used in CFD model
% (ANSYS Fluent )
% Authors: Kim Taewook
% Linkedin : https://www.linkedin.com/in/taewook-kim/
% GITHUB : github.com/Kimtaewookcode
% Email : kimtaewook87@gmail.com

clear
clf
clc
hold on; 
steps = 100;
timestep=0.01;%s
trials=10; 
x=zeros(trials, steps);
y=zeros(trials, steps); 
J=0;
swirl=0;%0.3;%m/s
coneangle=13;%degree
conev=1.87;%m/s
xconev=conev*cosd(coneangle);
yconev=conev*sind(coneangle);
k=0.02;
k2=sqrt(k*2/3);
for t = 1:trials 
 

for i=1:(steps-1)
J=rand;

if J<0.25
%x(t,i+1)=x(t,i)+(-k2*timestep);
%y(t,i+1)=y(t,i)+(k2*timestep);
x(t,i+1)=x(t,i)+(-k2+xconev)*timestep;
y(t,i+1)=y(t,i)+(k2+yconev)*timestep;


elseif J<0.5 
%x(t,i+1)=x(t,i)+(k2*timestep);
%y(t,i+1)=y(t,i)+(k2*timestep);
x(t,i+1)=x(t,i)+(k2+xconev)*timestep;
y(t,i+1)=y(t,i)+(k2+yconev)*timestep;

elseif J<0.75
%x(t,i+1)=x(t,i)+(-k2*timestep);
%y(t,i+1)=y(t,i)+(-k2*timestep);
x(t,i+1)=x(t,i)+(-k2+xconev)*timestep;
y(t,i+1)=y(t,i)+(-k2+yconev)*timestep; 

else
%x(t,i+1)=x(t,i)+(k2*timestep);
%y(t,i+1)=y(t,i)+(-k2*timestep);       
x(t,i+1)=x(t,i)+(k2+xconev)*timestep;
y(t,i+1)=y(t,i)+(-k2+yconev)*timestep;     


end 
end 
coneangle=coneangle-2.6;
xconev=conev*cosd(coneangle);
yconev=conev*sind(coneangle);
plot(x(t,:),y(t,:),'Color', [rand rand rand]) 
plot(x(t,100),y(t,100), 'ko')

%xlim([0 0.05])
%ylim([-0.05 0.05])
xlim([0 2])
ylim([-0.45 0.45])

end
%title('trajectories by instantaneous velocity(random walk) ');
title('trajectories by mean velocity+instantaneous velocity(random walk) ');
xlabel('X Displacement[m]');
ylabel('Y Displacement[m]');
