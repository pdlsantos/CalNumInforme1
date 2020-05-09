tol = 10^(-6);

tiemporala(1) = 0;      % Vector que guarda los tiempos de ejecución del algoritmo usando matrices ralas de distintos tamaños
tiempofull(1) = 0;      % Vector que guarda los tiempos de ejecución del algoritmo usando matrices ralas de distintos tamaños

funcionorala(1) = 1;    % Vector de 0 y 1, devuelve 1 si las primeras 6 cifras de la solución propuesta y de la solucion con LU coinciden
funcionofull(1) = 1;    % Vector de 0 y 1, devuelve 1 si las primeras 6 cifras de la solución propuesta y de la solucion con LU coinciden

%% Evalue el tiempo para distintos tamaños de matriz
%% Vea la diferencia entre usar matriz rala o normal

for n = 2:18
  N = 8 * (n^3);
  MAXIT = N;
  x0=zeros(N,1);
  
  % Construyo la matriz rala
  
  SD=sparse(1:N,1:N,6*ones(N,1),N,N);
  if1=setdiff(1:N-1, n:n:n*(8*n^2-1));
  SU1=sparse(if1,if1+1,-1*ones(1,(n-1)*8*(n^2)),N,N);
  [i, j]=find([ones(n*(2*n-1),4*n); zeros(n,4*n)]);
  if2=i+(j-1)*2*n^2;
  SU2=sparse(if2,if2+n,-1*ones(1,n*(2*n-1)*4*n),N,N);
  if3=1:2*n^2*(4*n-1);
  SU3=sparse(if3,if3+2*n^2,-1*ones(1,2*n^2*(4*n-1)),N,N);
  
  A=SD+SU1+SU1'+SU2+SU2'+SU3+SU3'; % Almacenamiento "ralo"
  
  b = ones(N,1)/(n + 1)^2;  % Vector del segundo miembro
  
  to = time();
  x = pcg(A, b, tol, MAXIT, [],  [], x0);
  t_calc = time() - to;
  
  tiemporala(n) = t_calc;   % Agrego el tiempo a mi vector de tiempos
  
  % Compruebo que funciono
  
  if floor((A*x-b)*(10^6)) == 0
    funcionorala(n) = 1;
  else
    funcionorala(n) = 0;
  endif
  
  clear('A');
  clear('x');
  
  A=full(SD+SU1+SU1'+SU2+SU2'+SU3+SU3'); % Por las dudas
  
  to = time();
  x = pcg(A, b, tol, MAXIT, [],[], x0);
  t_calc = time()-to;
  
  tiempofull(n) = t_calc;   % Agrego el tiempo a mi vector de tiempos
  
  % Compruebo que funciono
  
  if floor((A*x-b)*(10^6)) == 0
    funcionofull(n) = 1;
  else
    funcionofull(n) = 0;
  endif  
endfor

% https://octave.sourceforge.io/octave/function/pcg.html documentacion de pcg
% https://octave.sourceforge.io/octave/function/ilu.html documentacion de ilu
% https://octave.sourceforge.io/octave/function/lu.html  documentacion de lu