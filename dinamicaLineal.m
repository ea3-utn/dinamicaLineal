##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    ANALISIS LINEAL DINAMICO 
## /_____/_/  |_/___/___/___/     |
##                                |    : script principal
##---------CICLO LECTIVO 2020----------------------------------------------------

## CONFIGURACION

clear

pkg load symbolic; # Carga de paquete que me permite hacer operaciones algebraicas

warning ('off','OctSymPy:sym:rationalapprox');

syms x t tau

T=linspace(0,0.2,4000);

## DECLARACIONES

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

## CARACTERISTICAS DEL MATERIAL

densidad=8614.7329956;

young=189978213061.0;

cargaArmonica=0.209426;

alfa=0.000199641900697;

beta=29.8774375097;

###-----------

base=15e-3;

altura=1e-3;

E=young;

A=base*altura;

L=240e-3;

RHO=densidad;

I=(1/12)*base*altura^3;

g=9.81;

##---- CARACTERISTICAS FISICO-GEOMETRICAS

## NODO=[0 0;(L/2) 0;L 0]; % [Xi Yi] en fila i define nodo i

## ELEMENTO=[1 2 E A RHO I;2 3 E A RHO I]; % [NodoInicial NodoFinal E A I] en fila i define ubicacion de la viga "i-esima" y sus propiedades


NODO=[0 0;(L/4) 0;2*(L/4) 0;3*(L/4) 0;L 0]; % [Xi Yi] en fila i define nodo i

ELEMENTO=[1 2 E A RHO I;2 3 E A RHO I;3 4 E A RHO I;4 5 E A RHO I]; % [NodoInicial NodoFinal E A I] en fila i define ubicacion de la viga "i-esima" y sus propiedades

tipoElemento=2;

#########%  CONDICIONES DE CONTORNO Y DESPLAZAMIENTOS ---> G L O B A L E S

f=83;

w=2*pi*f;

Acce=cargaArmonica*sin(w*tau);

#Acce=taylor(acce, tau, 0, 'order', 3);

Vel=int(Acce,tau);

Disp=int(Vel,tau);

## ---- CALCULO DE LAS CONDICIONES INICIALES

velocidad=function_handle(Vel);

UPO=velocidad(0);

desplazamiento=function_handle(Disp);

UO=desplazamiento(0);




CCx=[1 0]; % [Nodo Ux] define condicion de contorno en nodo i

CCy=[1 Disp]; % [Nodo Uy] define condicion de contorno en nodo i 

CCw=[1 0]; % [Nodo Wz] define condicion de contorno en nodo i 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  SCRIPT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

barra=2; # Grados de libertad por nodo 

vigaPortico=3;

matrizElementos=[barra;vigaPortico]; # El orden vertical de los componentes indica su codigo

GL=size(NODO,1)*matrizElementos(tipoElemento); % Cant. de grados de libertad globales

[KG,MG]=rigidezGlobal(ELEMENTO,NODO,GL,tipoElemento);

[Px]=cargaEquivalente(KG,MG,CCx,CCy,CCw,GL,tipoElemento);

[KG,MG]=condContorno(KG,MG,CCx,CCy,CCw,GL,tipoElemento); % Aplicacion de las condiciones de contorno

GLNN=size(KG,1);

########### AUTOVECT/AUTOVALORES -->  W2  ###############

[phi,w2]=eig(inv(MG)*KG);

frecuenciasNaturales=sort(sqrt(w2)/(2*pi)*ones(GLNN,1));

printf("\nTotal %d Frecuencias Naturales:\n\n",GLNN)

printf("%d Hz\n",frecuenciasNaturales)


########## MATRIZ DE AMORTIGUAMIENTO POR METODO DE RAYLEIGH #######

velocidadAngular=sort(sqrt(w2)*ones(GLNN,1));

[C,alfa,beta]=amortiguamientoRayleigh(2,velocidadAngular(1),velocidadAngular(2),alfa,beta,MG,KG);


########### SUPERPOSICION MODAL ##########################

U=zeros(size(frecuenciasNaturales,1),size(T,2));

for i=1:GLNN 

  
	     #----- PROYECCIONES ------------------------------------#
  
  Mr=phi(:,i)'*MG*phi(:,i); # Proyeccion modal de la masa
  
  Kr=phi(:,i)'*KG*phi(:,i); # Proyeccion modal de la rigidez
  
  Pr=phi(:,i)'*Px; # Proyeccion modal de la carga

  Wr=sqrt(Kr/Mr);

  Epn=(beta/(2*Wr))+alfa*(Wr/2);

  printf("\nFrecuencia=%d Hz | Amortiguamiento=%d",(Wr/(2*pi)),Epn)
  
  Wd=Wr*sqrt(1-Epn^2);
  
	     # ---- SOLUCION DE LA INTEGRAL DE DUHAMEL --------------#
  
  #yr=function_handle((1/(Mr*(Wr)))*int(Pr*e^(-Epn*Wr*(t-tau))*sin((Wd)*(t-tau)),tau,0,t));
  
  homogenea=function_handle((e^(-Epn*Wr*t))*(UO*cos(Wd*t)+(UPO/Wd+(Epn*Wr*UO)/Wd)*sin(Wd*t)));

  homogeneaNumerica=homogenea(T);

  homogeneaNumerica(isnan(homogeneaNumerica))=0;
  
  particular=(1/(Mr*(Wd)))*Pr*e^(-Epn*Wr*(t-tau))*sin((Wd)*(t-tau));
  
  yr=function_handle(particular);

  yrIntegrado=integrador(yr,T);

  yrIntegrado(isnan(yrIntegrado))=0;
  
             # ---- PROYECCION DEL DESPLAZAMIENTO TEMPORAL ----------#

  U+=phi(:,i)*(yrIntegrado+homogeneaNumerica); # Reemplazando ya, el dominio temporal y convirtiendolo en numerico.

  #U+=phi(:,i)*(yrIntegrado); # Reemplazando ya, el dominio temporal y convirtiendolo en numerico.
endfor


############## COMPARACION - ASTER ################



Shaker=function_handle(Disp);

VelU5=diff(U(11,:)-Shaker(T))./diff(T);

acceU5=diff(VelU5)./diff(T(1:end-1));

TVEL=T(1:end-1);

TACCE=T(1:end-2);

figure (1);clf;hold on;grid on;

title ('N3 OCTAVE')

plot(T,U(11,:)-Shaker(T),["--" markStyle(4) color(4) ";DISP;"]);

hold off


figure (2);clf;hold on;grid on;

title ('N3 OCTAVE')

plot(TVEL,VelU5,["--" markStyle(3) color(3) ";VEL;"]);

hold off

Simul=dlmread('/home/nico/_org/investigacionDesarrollo/labEstructuras/caia2021/reglaLarga3Elementos/resultados/visualizacionOctave/plotFreqSave.resu',',',0,0);

figure (3);clf;hold on;grid on;

title ('N3 OCTAVE')

plot(TACCE,acceU5,["--" markStyle(5) color(5) ";ACCE;"]);

plot(Simul(:,1),Simul(:,2),["--" markStyle(5) color(4) ";Acce Codeaster;"]);
hold off

keyboard
