##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    ANALISIS LINEAL DINAMICO
## /_____/_/  |_/___/___/___/     |
##                                |    Elimina componentes NaN de un vector simbolico arbitario 
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function [Vector]=eliminadoNaNvectorSimbolico(U);

  Vector=sym(zeros(sum(not(isnan(U))),1));

  j=1;
  
  for i=1:size(U,1)
    
    if (isnan(U(i))==0)

      Vector(j)=U(i);

      j++;
      
    endif
  endfor


endfunction


