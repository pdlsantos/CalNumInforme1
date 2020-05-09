%%% PRIMER PROBLEMA INCISO A

funciono(1,1) = 1; % Esta matriz va a dar 1 si se alcanzo la solucion con una tolerancia de 10^-10
TIEMPO(1,1) = 0;   % Cada columna de esta matriz representa un muestreo, en cada muestreo se prueba el método para distintos tamaños de matrices

cantMuestras = 1  % Cantidad de muestras a tomar
cantColumnas = 10  % Tamaño de la matriz mas grande a usar

for muestras = 1: cantMuestras
  for i = 2: cantColumnas
    n = 2**i;
    A = rand(n, n);               % Crea una matriz random
    x = rand(n, 1);               % Crea un vector random, este vector es la solución esperada

    b =A * x;                     % da el vector del segundo miembro de la igualdad Ax = b

    %%% Factorizacion LU  O(n3)
    
    U = A                         % Copiamos A en U, guardamos en U las matrices intermedias
    L = eye(n,n)                  % L comieza siendo una matriz identidad de n x n
    
    t0 = time();                  %Justo antes de empezar a resolver
    for k = 1:(n-1)               % Columanas a reducir
      for i = (k+1):n             % recorre los elementos por debajo de la diagonal
        L(i,k) = U(i,k) / U(k,k); % es necesario que U no se anule en la diagonal
        for j = (k+1):n           % opera sobre las i, recorre con j  toddos los coeficientes
          U(i,j) = U(i,j) - L(i,k)*U(k,j);
        endfor
        U(i,k) = 0;               % Prolijidad
      endfor
    endfor
    
    %%% Resolucion Ly=b     Algoritmo Forward - Sustitution     O(n2)

    for k = 1:n
      y(k) = b(k);                % Le asigna bi a yi
      for i = (k+1):n
        b(i) = b(i) - L(i,k)*y(k);% usa el y ya asignado para sacar el valor del nuevo y
      endfor
    endfor                        % Obtengo y
    
    %%% Resolucion Ux=y     Algoritmo Backward - Sustitution    O(n2)

    for k = n:-1:1
      X(k) = y(k) / U(k,k);       % lo mismo que el anterior pero en este caso la diagonal no es 1
      for i = 1:(k-1)
        y(i) = y(i) - U(i,k)*X(k);
      endfor
    endfor                        % X es la solución a mi sistema dada por el algoritmo
    
    t_calc = time() - t0;         % Justo al terminar de resolver
    
    TIEMPO(n, muestras) = t_calc;
    
    % Comprobar la tolerancia
    
    for w = 1:n
      if floor(x(w)*(10^10)) != floor(X(w)*(10^10)) % Sus primeros 10 decimales son iguales
        funciono(n, muestras) = 0;
        break
      else
        funciono(n, muestras) = 1;
      endif
    endfor
  endfor 
endfor
funciono(1,:) = 1;                %El tamaño más chico de matriz es 2x2





