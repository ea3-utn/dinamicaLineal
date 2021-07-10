##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    ANALISIS LINEAL DINAMICO 
## /_____/_/  |_/___/___/___/     |
##                                |    integrador: Integracion numerica Duhamel
##                                |                     
##---------CICLO LECTIVO 2020----------------------------------------------------------------


function Y=integrador(duhamel,T)

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];

###########------------------------------

  
  PRESICION=1000;

  Y=zeros(size(T));

  Y(1)=0;
  
  for i=2:size(T,2)

    ## if (i==583)
    ##   keyboard
    ## endif
    
    dom=linspace(0,T(i),PRESICION);
    
    y=real(duhamel(T(i),dom));

    if (y==0)

      y=zeros(size(dom));

    endif
        
    Y(i)=trapz(dom,y);

  endfor  

  ## figure (1);clf;hold on;grid on;

  ## title ('Vista de integracion')

  ## plot(T,Y,["--" markStyle(2) color(2) ";Integracion;"]);

  ## hold off

  ## keyboard
endfunction

