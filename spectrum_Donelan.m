function [w,wp,cp,spectrum]=spectrum_Donelan(U,wave_age)

f=1;
w=1*f*(0:2:2048*f*10)/1024/f*2*pi;
cp=wave_age*U;
kp=9.8/cp^2;
wp=sqrt(9.8*kp);
yita=0.08*(1+4/(U/cp)^3);
gama=1.7+6*log(U/cp);
alpha=0.006*(U/cp)^0.55;
da_tao=exp(-(w-wp).^2/2/yita^2/wp^2);
spectrum=alpha*9.8^2*w.^-4*wp^-1.*exp(-(wp./w).^4).*gama.^da_tao;

end

