##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    AMORTIGUAMIENTO DE RAYLEIGH
## /_____/_/  |_/___/___/___/     |
##                                |    : crea la matriz C por el metodo de Rayleigh
##---------CICLO LECTIVO 2020----------------------------------------------------

# Clough & Penzien - Dynamics of structures - Pag.256

#  --  USO: FLAG=1 --> Entran Epsilons, salen Alfas | FLAG=2 --> Inversa | FLAG=0 --> Sistema sin amort.

function [C,alfa,beta]=amortiguamientoRayleigh(FLAG,wn,wm,en,em,M,K)

  if (FLAG==1)

    constantes= @(Wn,Wm,En,Em) ((2*Wm*Wn)/(Wn^2-Wm^2))*[Wn -Wm;(-1/Wn) (1/Wm)]*[Em;En];

    numerica=constantes(wn,wm,en,em);alfa=numerica(2);beta=numerica(1);

    C=alfa*[M]+beta*[K];

    printf("\nCoeficientes de la matriz de amortiguamiento de Rayleigh:\n\n",alfa)

    printf("alfa=%d \n",alfa)

    printf("beta=%d \n",beta)
    
  elseif (FLAG==2)

    alfa=en;

    beta=em;
    
    constantes= @(Wn,Wm,alfa,beta) (1/((2*Wm*Wn)/(Wn^2-Wm^2)))*inv([Wn -Wm;(-1/Wn) (1/Wm)])*[beta;alfa];

    numerica=constantes(wn,wm,en,em);e2=numerica(2);e1=numerica(1);

    C=alfa*[M]+beta*[K];

    printf("\nCoeficientes de amortiguamiento modal equivalentes:\n\n",alfa)

    printf("Epsilon 1=%d \n",e1)

    printf("Epsilon 2=%d \n",e2)

  else
        
    C=zeros(size(K));

    alfa=0;beta=0;
    
    printf("\nSistema no amortiguado\n",alfa)
    
  endif

endfunction

  
