clc
close all; clear all;


sigma = 4;% Conductivity in S/m - steel 1.45e6, sea water 4, pinewood 1.754e-3, Medium-density fibreboard 2.9e-3, brick 0.0175, carbon fibre 200
epsr = 80; % relative electric permittivity - steel 1, sea water 80, pinewood 1.699, Medium-density fibreboard 3, brick 4.6, carbon fibre 3
mur = 1; % relative magnetic permeability - steel 1000, 1 for all others
t=0.004; % thickness of layer in meters

eps0 = 8.8541878128e-12; %electric permittivity of freespace in F/m
mu0 = 1.25663706212e-6; %magnetic permeability of freespace in N/A^2
c = 3e8;%m/s (speed of light in vacuum)
imp0=377; % freespace impedence in Ohms
fnum = 100; % number of points in frequency plot
xnum = 1000; % number of points in amplitude decay plot
f = linspace(1e8, 5e9, fnum);
x = linspace(0, 5e3, xnum);
fplot = 100; % frequency in Hz used in amplitude decay plot
Pant = 1; % Transmitting antenna power in watts
Pdet = 0.01; % Minimum detectable power in watts


alpha = attenuation_coef(f,sigma,epsr,mur);%Calculating attenuation coefficient of material with conductivity sigma, relative permittivity epsr and relative permeability mur
beta = phase_coef(f,sigma,epsr,mur);%Calculating phase constant of material with conductivity sigma, relative permittivity epsr and relative permeability mur
impl = impedence(f,sigma,epsr,mur);%Calculating impedence of material with conductivity sigma, relative permittivity epsr and relative permeability mur


slabtrans = slabtransmisivity(transmisivity(impl,imp0),transmisivity(imp0,impl),reflectivity(imp0,impl),reflectivity(imp0,impl),alpha,beta,t);
slabref = slabreflectivity(transmisivity(impl,imp0),transmisivity(imp0,impl),reflectivity(impl,imp0),reflectivity(imp0,impl),reflectivity(imp0,impl),alpha,beta,t);
% Calculating reflection and transmission coefficients for wave incident on
% slab of material with thickness t, impedence impl, and propagation coefficient alpha+j*beta, sandwiched between air

figure(1)%frequency spectrum for slab transmittance and reflectance
plot(f,abs(slabref).^2,'b','LineWidth',4);
hold on
plot(f,abs(slabtrans).^2,'k','LineWidth',4);
legend('Reflectance','Transmittance','Location','east')
xlabel('Frequency (Hz)')

figure(2)%Power decay of fixed frequency electromagnetic wave with distance 
plot(x,Pant*(exp(-attenuation_coef(fplot,sigma,epsr,mur).*x)).^2,'b','LineWidth',4);
hold on
plot(x,Pdet*x./x,'k','LineStyle',':','LineWidth',3);
ylim([0 Pant])
xlim([0 100])
legend('100 Hz','Location','southwest')
xlabel('Distance (m)')

figure(3)%Power decay of fixed frequency electromagnetic wave with distance (y = log scale)
semilogy(x,Pant*(exp(-attenuation_coef(fplot,sigma,epsr,mur).*x)).^2,'b','LineWidth',4);
hold on
semilogy(x,Pdet*x./x,'k','LineStyle',':','LineWidth',3);
ylim([1e-10 Pant])
xlim([0 500])
legend('100 Hz','detection limit','Location','southwest')
xlabel('Distance (m)')

function alpha = attenuation_coef(freq,sigma,epsr,mur)
    c = 3e8;%m/s (speed of light in vacuum)
    mu0 = 1.25663706212e-6;%magnetic permeability of freespace in N/A^2
    eps0 = 8.8541878128e-12;%electric permittivity of freespace in F/m
    alpha = 2*pi*(freq/c).*sqrt((epsr*mur/2).*(sqrt(1+(sigma./(2*pi*freq*epsr*eps0)).^2)-1));
end

function beta = phase_coef(freq,sigma,epsr,mur)
    c = 3e8;%m/s (speed of light in vacuum)
    mu0 = 1.25663706212e-6;%magnetic permeability of freespace in N/A^2
    eps0 = 8.8541878128e-12;%electric permittivity of freespace in F/m
    beta = 2*pi*(freq/c).*sqrt((epsr*mur/2).*(sqrt(1+(sigma./(2*pi*freq*epsr*eps0)).^2)+1));
end

function imp = impedence(freq,sigma,epsr,mur)
    c = 3e8;%m/s (speed of light in vacuum)
    mu0 = 1.25663706212e-6;%magnetic permeability of freespace in N/A^2
    eps0 = 8.8541878128e-12;%electric permittivity of freespace in F/m
    imp = sqrt(1i*2*pi*freq*mur*mu0./(sigma+1i*2*pi*freq*epsr*eps0));
end

function ref = reflectivity(imp2,imp1)%function computing reflection coeff for wave passing from material with impedence imp1 to imp2
    c = 3e8;%m/s (speed of light in vacuum)
    mu0 = 1.25663706212e-6;%magnetic permeability of freespace in N/A^2
    eps0 = 8.8541878128e-12;%electric permittivity of freespace in F/m
    ref = (imp2-imp1)./(imp1+imp2);
end

function trans = transmisivity(imp2,imp1) %function computing transmission coeff for wave passing from material with impedence imp1 to imp2
    c = 3e8;%m/s (speed of light in vacuum)
    mu0 = 1.25663706212e-6;%magnetic permeability of freespace in N/A^2
    eps0 = 8.8541878128e-12;%electric permittivity of freespace in F/m
    trans = 2*imp2./(imp1+imp2);
end

function slabtrans = slabtransmisivity(t1,t2,r2,r3,alph,bet,d) %function computing transmission coeff for wave passing from freespace through a slab of material with propogation constant alph+j*bet and thickness d
%t1 - transmissivity into slab, t2 - transmissivity out of slab, r2 - reflection from within slab off of incident wall, r3 - reflection from within slab off of transmission wall (same as r2 if same material both sides)   
    slabtrans = t1.*t2.*exp(-alph*d-1i*bet*d)./(1-r2.*r3.*exp(-alph*2*d-1i*bet*2*d));
end

function slabref = slabreflectivity(t1,t3,r1,r2,r3,alph,bet,d) %function computing reflection coeff for wave passing from freespace through a slab of material with propogation constant alph+j*bet and thickness d
%t1 - transmissivity into slab, t3 - transmissivity back out of slab, r1 - reflection from outside slab off of incident wall, r2 - reflection from within slab off of incident wall, r3 - reflection from within slab off of transmission wall (same as r2 if same material both sides)    
    slabref = r1+t1.*r2.*t3.*exp(-alph*2*d-1i*bet*2*d)./(1-r2.*r3.*exp(-alph*2*d-1i*bet*2*d));
end


    