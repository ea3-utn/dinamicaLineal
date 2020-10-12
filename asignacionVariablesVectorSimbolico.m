##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    ANALISIS LINEAL DINAMICO
## /_____/_/  |_/___/___/___/     |
##                                |    Asigna variables de forma masiva a un vector simbolico
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function [vectorModificado]=asignacionVariablesVectorSimbolico(vectorInicial,index,variables);
  
  for i=1:size(index,1)

    vectorInicial(index(i))=variables(i);
    
  endfor

  vectorModificado=vectorInicial;
  
endfunction


