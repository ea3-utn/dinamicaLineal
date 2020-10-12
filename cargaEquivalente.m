##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    ANALISIS LINEAL DINAMICO
## /_____/_/  |_/___/___/___/     |
##                                |    cargaEquivalente: Config. Desplazamientos, Cargas y matrices asoc. para analisis armonico 
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function [PxOK]=cargaEquivalente(Kxy,Mxy,CCx,CCy,CCw,GL,codigo);
  
  rigidezNula=all(~Kxy,2);
  
  indexGL=linspace(1,GL,GL)';

  U=sym(zeros(GL,1));
  
  switch (codigo)

	 case 1
	   
	   	   
	   CONDICIONES=double(sort([2*CCx(:,1)-1;2*CCy(:,1)]));

	   [U]=asignacionVariablesVectorSimbolico(U,2*CCx(:,1)-1,CCx(:,2));

	   [U]=asignacionVariablesVectorSimbolico(U,2*CCy(:,1),CCy(:,2));
	   
	 case 2

	   if (isempty(CCw)==0)
	     
	     CONDICIONES=double(sort([3*CCx(:,1)-2;3*CCy(:,1)-1;3*CCw(:,1)]));
	     
	     [U]=asignacionVariablesVectorSimbolico(U,3*CCx(:,1)-2,CCx(:,2));

	     [U]=asignacionVariablesVectorSimbolico(U,3*CCy(:,1)-1,CCy(:,2));

	     [U]=asignacionVariablesVectorSimbolico(U,3*CCw(:,1),CCw(:,2));

	   else

	     CONDICIONES=double(sort([3*CCx(:,1)-2;3*CCy(:,1)-1]));

	     [U]=asignacionVariablesVectorSimbolico(U,3*CCx(:,1)-2,CCx(:,2));

	     [U]=asignacionVariablesVectorSimbolico(U,3*CCy(:,1)-1,CCy(:,2));


	   endif
	   
	   
  endswitch


  ##-- Calculo matricial
  
  Mxy(CONDICIONES,:)=[];

  Kxy(CONDICIONES,:)=[];

  indexGL(CONDICIONES)=[];

  Mxy(:,indexGL)=[];

  Kxy(:,indexGL)=[];
  
  ##-- Adecuacion U
  
  U(indexGL)=NaN;

  [Uy]=eliminadoNaNvectorSimbolico(U);
  
  ###-- Calculo carga equivalente
  
  Px=-Mxy*diff(Uy,2)-Kxy*Uy;

  rigidezNula(CONDICIONES,:)=[];

  try
    
  Px(rigidezNula==1,:)=NaN;

  end_try_catch
  
  [PxOK]=eliminadoNaNvectorSimbolico(Px);
  
endfunction


