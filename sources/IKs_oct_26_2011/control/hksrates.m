
% this function takes input as starting time t, initial conditions x, and
% value of voltage v. Put in parameters into it now.
function y=hksrates(t,x,v,ploc)

global FoRT;

alpha=ploc(1)*exp(v*FoRT*ploc(2));
beta=ploc(3)*exp(v*FoRT*ploc(4));
gamma=ploc(5)*exp(v*FoRT*ploc(6));
delta=ploc(7)*exp(v*FoRT*ploc(8));
theta=ploc(9);
eta=ploc(10)*exp(v*FoRT*ploc(11));
psi=ploc(12)*exp(v*FoRT*ploc(13));
omega=ploc(14)*exp(v*FoRT*ploc(15));



c1=x(1);
c2=x(2);
c3=x(3);
c4=x(4);
c5=x(5);
c6=x(6);
c7=x(7);
c8=x(8);
c9=x(9);
c10=x(10);
c11=x(11);
c12=x(12);
c13=x(13);
c14=x(14);
a=x(15);
o1=x(16);
o2=x(17);

do2dt = psi*o1-omega*o2;
do1dt = theta*a+omega*o2-(psi+eta)*o1;
dadt = gamma*c14+eta*o1-(4*delta+theta)*a;
dc14dt = alpha*c13+4*delta*a+2*gamma*c12-(beta+3*delta+gamma)*c14;
dc13dt = beta*c14+gamma*c11-(alpha+3*delta)*c13;
dc12dt = alpha*c11+3*delta*c14+3*gamma*c9-(2*beta+2*delta+2*gamma)*c12;
dc11dt = 2*alpha*c10+2*beta*c12+2*gamma*c8+3*delta*c13-(beta+alpha+gamma+2*delta)*c11;
dc10dt = beta*c11+gamma*c7-(2*alpha+2*delta)*c10;
dc9dt = alpha*c8+2*delta*c12+4*gamma*c5-(3*beta+delta+3*gamma)*c9;
dc8dt = 2*alpha*c7+3*beta*c9+3*gamma*c4+2*delta*c11-(2*beta+alpha+2*gamma+delta)*c8;
dc7dt = 3*alpha*c6+2*beta*c8+2*gamma*c3+2*delta*c10-(beta+2*alpha+delta+gamma)*c7;
dc6dt = beta*c7+gamma*c2-(3*alpha+delta)*c6;
dc5dt = alpha*c4+delta*c9-(4*beta+4*gamma)*c5;
dc4dt = 2*alpha*c3+4*beta*c5+delta*c8-(3*beta+alpha+3*gamma)*c4;
dc3dt = 3*alpha*c2+3*beta*c4+delta*c7-(2*beta+2*alpha+2*gamma)*c3;
dc2dt = 4*alpha*c1+2*beta*c3+delta*c6-(beta+3*alpha+gamma)*c2;
dc1dt = beta*c2-4*alpha*c1;

y(1,1)=dc1dt;
y(2,1)=dc2dt;
y(3,1)=dc3dt;
y(4,1)=dc4dt;
y(5,1)=dc5dt;
y(6,1)=dc6dt;
y(7,1)=dc7dt;
y(8,1)=dc8dt;
y(9,1)=dc9dt;
y(10,1)=dc10dt;
y(11,1)=dc11dt;
y(12,1)=dc12dt;
y(13,1)=dc13dt;
y(14,1)=dc14dt;
y(15,1)=dadt;
y(16,1)=do1dt;
y(17,1)=do2dt;

return