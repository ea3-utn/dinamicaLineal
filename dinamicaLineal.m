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

syms x

## DECLARACIONES

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

## CARACTERISTICAS DEL MATERIAL

rho=7800;

A=0.0079;

I=4.9087e-6;

E=2.1e11;

L=2.25;

##---- CARACTERISTICAS FISICO-GEOMETRICAS


NODO=[0 0;0 2*L;L 2*L;L L;L 0]; % [Xi Yi] en fila i define nodo i

ELEMENTO=[1 2 E A rho I;2 3 E A rho I;3 4 E A rho I;4 5 E A rho I;2 4 E A rho I;1 4 E A rho I]; % [NodoInicial NodoFinal E A Rho I] 

tipoElemento=1; 

#########%  CONDICIONES DE CONTORNO Y CARGAS ---> G L O B A L E S

%-------- COND./CARGAS NULAS SE DEJAN V A C I A S = [] -------------#

CCx=[5 0]; % [Nodo Ux] define condicion de contorno en nodo i

CCy=[1 0;5 0]; % [Nodo Uy] define condicion de contorno en nodo i 

CCw=[0 0]; % [Nodo Wz] define condicion de contorno en nodo i 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  SCRIPT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

barra=2; # Grados de libertad por nodo 

vigaPortico=3;

matrizElementos=[barra;vigaPortico]; # El orden vertical de los componentes indica su codigo

GL=size(NODO,1)*matrizElementos(tipoElemento); % Cant. de grados de libertad globales

[KG,MG]=rigidezGlobal(ELEMENTO,NODO,GL,tipoElemento);

[KG,MG]=condContorno(KG,MG,CCx,CCy,CCw,GL,tipoElemento); %Aplicacion de las condiciones de contorno

########### AUTOVECT/AUTOVALORES -->  W2  ###############

[phi,w2]=eig(inv(MG)*KG);


sqrt(w2)/(2*pi)



