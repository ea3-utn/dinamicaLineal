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

T=linspace(0,0.05,500);

## DECLARACIONES

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

## CARACTERISTICAS DEL MATERIAL

b=1;

E=100;

A=3;

L=2;

RHO=7;

I=100;

##---- CARACTERISTICAS FISICO-GEOMETRICAS


NODO=[0 0;L/2 0;L 0]; % [Xi Yi] en fila i define nodo i

ELEMENTO=[1 2 E A RHO I;2 3 E A RHO I]; % [NodoInicial NodoFinal E A I] en fila i define ubicacion de la viga "i-esima" y sus propiedades

tipoElemento=2; 

#########%  CONDICIONES DE CONTORNO Y DESPLAZAMIENTOS ---> G L O B A L E S

CCx=[1 0;3 0]; % [Nodo Ux] define condicion de contorno en nodo i

CCy=[1 (-2e5)*tau^2;3 sym(0)]; % [Nodo Uy] define condicion de contorno en nodo i 

CCw=[1 0;3 0]; % [Nodo Wz] define condicion de contorno en nodo i 

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

########### SUPERPOSICION MODAL ##########################


U=zeros(size(frecuenciasNaturales,1),size(T,2));

for i=1:GLNN 
  
  
	     #----- PROYECCIONES ------------------------------------#

  
  
  #y0=(phi(:,i)'*Cinit)/(phi(:,i)'*phi(:,i)) # Proyeccion de la condicion inicial de la funcion del tiempo del modo

  y0=0;
  
  Mr=phi(:,i)'*MG*phi(:,i); # Proyeccion modal de la masa
  
  Kr=phi(:,i)'*KG*phi(:,i); # Proyeccion modal de la rigidez
  
  Pr=phi(:,i)'*Px; # Proyeccion modal de la carga

	     # ---- SOLUCION DE LA INTEGRAL DE DUHAMEL --------------#
  
  yr=function_handle(y0*cos((sqrt(Kr/Mr))*t)+(1/(Mr*(sqrt(Kr/Mr))))*int(Pr*sin((sqrt(Kr/Mr))*(t-tau)),tau,0,t));


             # ---- PROYECCION DEL DESPLAZAMIENTO TEMPORAL ----------#

  U+=phi(:,i)*yr(T); # Reemplazando ya, el dominio temporal y convirtiendolo en numerico.
  
  
endfor


aster=dlmread('../codeAster/analisisArmonico/despDY_N2.resu',',',0,0);

figure (1);clf;hold on;grid on;

title ('COMPARACION ASTER-OCTAVE')

N1y=function_handle((-2e5)*t^2);

plot(T,N1y(T),["--" markStyle(1) color(1) ";N1 Octave;"]);

plot(T,U(2,:),["--" markStyle(2) color(2) ";N2 Octave;"]);

plot(aster(1:500,1),aster(1:500,2),["--" markStyle(3) color(3) ";N1 Aster;"]);

plot(aster(1:500,1),aster(1:500,3),["--" markStyle(4) color(4) ";N2 Aster;"]);

hold off
