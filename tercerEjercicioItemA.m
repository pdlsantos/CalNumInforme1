tol = 10^(-6);

tiempoSP(1) = 0;  % Vector para guardar los tiempos de GCsp
tiempoCP(1) = 0;  % Vector para guardar los tiempos de GCp
tiempoLU(1) = 0;  % Vector para guardar los tiempos de LU

%% Evalue el tiempo para distintos tamaños de matriz
%% Vea la diferencia entre usar GCsp, GCp y LU

for n = 2:10 
  N = 8 * (n^3);
  maxit = N;
  x0 = zeros(N,1);
  
  % Genero la matriz rala
  
  SD  = sparse(1:N,1:N,6*ones(N,1),N,N);
  if1 = setdiff(1:N-1, n:n:n*(8*n^2-1));
  SU1 = sparse(if1,if1+1,-1*ones(1,(n-1)*8*(n^2)),N,N);
  [i, j] = find([ones(n*(2*n-1),4*n); zeros(n,4*n)]);
  if2 = i+(j-1)*2*n^2;
  SU2 = sparse(if2,if2+n,-1*ones(1,n*(2*n-1)*4*n),N,N);
  if3 = 1:2*n^2*(4*n-1);
  SU3 = sparse(if3,if3+2*n^2,-1*ones(1,2*n^2*(4*n-1)),N,N);
  
  A = SD+SU1+SU1'+SU2+SU2'+SU3+SU3'; % Almacenamiento "ralo"
  
  b = ones(N,1)/(n + 1)^2; % Vector del segundo miembro
  
  %ultimoresultado = [(A\b)(1), (A\b)(N)]; %Al final de cada bloque de procedimiento hay muestreos comentados
  
  %% Método pcg sin Precondicionamiento
  
  to = time();
  x  = pcg (A, b, tol, maxit, [], [], x0);
  t_calc = time() - to;
  
  tiempoSP(n) = t_calc;
  %ultimoSP = [x(1), x(N)];

  %% Método pcg con Precondicionamiento
  
  to = time();
  [L, U] = ilu(A);
  x  = pcg (A, b, tol, maxit, L, U, x0);
  t_calc = time()-to;
  
  tiempoCP(n) = t_calc;
  %ultimoCP = [x(1), x(N)];  

  %% Método LU

  to = time();
  [L, U, P, Q] = lu(A);
  y = L\P*b;
  x = Q*(U\y);
  t_calc = time()-to;
  
  tiempoLU(n) = t_calc;
  %ultimoLU = [x(1), x(N)];
endfor

% https://octave.sourceforge.io/octave/function/pcg.html documentacion de pcg
% https://octave.sourceforge.io/octave/function/ilu.html documentacion de ilu
% https://octave.sourceforge.io/octave/function/lu.html  documentacion de lu