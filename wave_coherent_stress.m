clear;clc;
%% the wind speed at 10 m  above the surface
U10=6;
%% wind stress computed from Smith (1980)
Cd=0.61+0.061*U10;
u_star=sqrt( Cd*U10^2  *10^-3); %% friction velocity 

zL_part=-10/u_star^3/(0.1); %% atmosphere stability
wind_stress=u_star^2;
%% z
y=(-16):10^-2:log(200);
z=exp(y);
%% wind profile given by MOST under neutral condition
z0=exp( -( U10*0.4/u_star-log(10) ) );

U_log=u_star/0.4*log(exp(y)/z0);

z0_local=find(abs(U_log)==min(abs(U_log)) );


%%  dimensionless wing gradient according to Hogstrom (1988)
zL_profile=z./zL_part/u_star^3;
if zL_part<=0
phi=(1-19*zL_profile).^(-1/4);
else
phi=(1+6*zL_profile);
end
%% wave spectra computed from Donelan et al. (1985)
[w,wp,cp,spectrum_o]=spectrum_Donelan(30,0.4);
k=w.^2/9.8;c=w./k;c(1)=0;
[w,wp,cp,spectrum_wind]=spectrum_Donelan(U10,1.2);
spectrum_swell=spectrum_o.*exp(-w.^3/0.5);
spectrum=spectrum_swell+spectrum_wind*0;

%% solve Eqs. (23) and (24) iterative. The boundary conditon at 10 m above the surfae
%% is U10. This is done by first assuming that the wind profile follows U_log and 
%% the turbulent stress equals to total wind stress
z_10_loca=find( abs(z-10)==min(abs(z-10)) );
u_star_t(1:length(y))=u_star;
U=U_log;
dU_dz=u_star_t/0.4./exp(y);
yita=-41.61*u_star^2;
D=1.74;
%% to generate wave coherent stress

mm=2;nn=600;

for u_star_loop=1:1:5
for z_i=z0_local:1:length(y)

  singula=1-exp(-1*k*z(z_i).*abs(mm*(1-U(z_i)./c)).^nn); %% used to remove any singularity 
   WC_stress_part1= D*k.*exp(-D*k*z(z_i))./(c-U(z_i)).*singula-exp(-D*k*z(z_i))./(c-U(z_i)).^2*dU_dz(z_i).*singula;   
   WC_stress_part2= yita*0.4*z(z_i)/u_star_t(z_i)*w.^2.*spectrum/phi(z_i);
   
   WC_stress_all(z_i,:)=( WC_stress_part1(2:end).*WC_stress_part2(2:end))*diff(w(1:2));
   WC_stress_all(z_i,WC_stress_all(z_i,:)>0)=0; %% used to remove positive component 

   
   WC_stress(z_i)=sum(WC_stress_all(z_i,:));

end

u_star_t=sqrt(wind_stress./(1-WC_stress/u_star^2));

WC_stress=WC_stress/u_star^2.*u_star_t.^2;

dz=diff(z);
%%% based on u_star_t, the new wind profile is computed to determine z0
for z_i=z_10_loca:-1:1
       U_zi(z_i)=sum( u_star_t(z_i:z_10_loca)/0.4./z((z_i:z_10_loca)).*dz(z_i:z_10_loca).*phi(z_i:z_10_loca));
    end
z_0=find( abs(U_zi-U10)==min(abs(U_zi-U10)) );
clearvars U

for z_i=z_0:1:length(z)
    U(z_i)=sum(u_star_t(z_0:z_i)/0.4./z((z_0:z_i)).*dz(z_0-1:z_i-1).*phi(z_0-1:z_i-1));
    end
z0_loca=z_0;

end


%% f ~ z/Lturb

z_3=find( abs(z-3)==min(abs(z-3)) );

zL=-2:0.1:2; %% atmosphere stability

for i=1:1:length(zL)
 %%  dimensionless wing gradient according to Hogstrom (1988)
   
if zL(i)<=0
phi_f=(1-19*zL(i)).^(-1/4);
else
phi_f=(1+6*zL(i));
end

f_zL(i)=1/ (phi_f.^4.*(1-zL(i)./phi_f)).^(1/4); %% without wave effect

beta2_2=-WC_stress(z_3)'/u_star_t(z_3)^2;
f_zL_wave(i)=1/ (phi_f.^4.*abs(1-(zL(i)--abs(beta2_2))./phi_f)).^(1/4); %% containing wave effect

end

plot(zL,f_zL,'*');hold on;

plot(zL,f_zL_wave);hold on;