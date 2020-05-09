%%% PRIMER PROBLEMA INCISO B

TIEMPO(1,1) = 0;   % Cada columna de esta matriz representa un muestreo, en cada muestreo se prueba el método para distintos tamaños de matrices
cantMuestras = 1  % Cantidad de muestras a tomar
cantColumnas = 20  % Tamaño de la matriz mas grande a usar

%funcionoSP = zeros(cantColumnas, cantMuestras); % Esta matriz va a dar la cantidad de decimales correctos que tiene la respuesta  sin preccondicionamiento
%funcionoCP = zeros(cantColumnas, cantMuestras); % Esta matriz va a dar la cantidad de decimales correctos que tiene la respuesta  con preccondicionamiento

errorCP = [];    %Este vector va a dar el error segun el tamaño de la ultima muestra precondicionada
errorSP = [];    %Este vector va a dar el error segun el tamaño de la ultima muestra sin precondicionar

for muestras = 1: cantMuestras
  for i = 2: cantColumnas
    n = i**2;
    A = rand(n, n);               % Crea una matriz random
    x = rand(n, 1);               % Crea un vector random, este vector es la solución esperada
    A(1,1) = (10^(-10));          % Le asigna al primer elemento de la primer columna 10^-10
    b = A * x;                    % Da el vector del segundo miembro de la igualdad Ax = b   
    b1 = b;                       % Guardo a b en la variable b1
    
    %%% Factorizacion LU  O(n3)
    
    U = A;                        % Copiamos A en U, guardamos en U las matrices intermedias
    L = eye(n,n);                 % L comieza siendo una matriz identidad de n x n
    
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
        b(i) = b(i) - L(i,k)*y(k);% usa el "y" ya asignado para sacar el valor del nuevo "y"
      endfor
    endfor                        % Obtengo y
    
    %%% Resolucion Ux=y     Algoritmo Backward - Sustitution    O(n2)

    for k = n:-1:1
      Z(k) = y(k) / U(k,k);       % lo mismo que el anterior pero en este caso la diagonal no es 1
      for i = 1:(k-1)
        y(i) = y(i) - U(i,k)*Z(k);
      endfor
    endfor                        % Z es la solución a mi sistema dada por el algoritmo
    
    clear('y');
    clear('L');                   % Por las dudas
    clear('U');
    
    %%% Factorizacion LU Con Pivoteo
    
    P = eye(n, n);                % Si no hay permutaciones P es una matriz identidad de n x n
    U = A;                        % Copiamos A en U, guardamos en U las matrices intermedias
    L = eye(n, n);                % L comieza siendo una matriz identidad de n x n 
    
    for k = 1:(n-1)               % Columanas a reducir
      
      % Antes de empezar el algoritmo de factorización, veo al elemento (k, k) de la matriz intermedia U
      
      if (abs(U(k,k)) <= (10^(-10))) %Valor absoluto porque los numeros negativos son menores que 10^-10
        max = U(k, k);
        wn = k;
        for w = (k+1):n
          if max < U(w, k)
            max = U(w, k);
            wn = w;
          endif
        endfor
        pivotP  = P(k,:);
        pivot = U(k,:);
        P(k, :) = P(wn, :);
        U(k, :) = U(wn, :);
        P(wn, :) = pivotP;
        U(wn, :) = pivot;
      endif
      
      % Factorización LU de la matriz pivoteada
      
      for i = (k+1):n             % recorre los elementos por debajo de la diagonal
        L(i,k) = U(i,k) / U(k,k); % es necesario que U no se anule en la diagonal
        for j = (k+1):n           % opera sobre las i, recorre con j  toddos los coeficientes
          U(i,j) = U(i,j) - L(i,k)*U(k,j);
        endfor
        U(i,k) = 0;               % Prolijidad
      endfor
    endfor
    
    b = P * b1;
    
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
    
    %%% Evaluación de la tolerancia
    
    for tol = 1:20
      conSP = 0;  % Quiero que le asigne la tolerancia solo si todos los elementos del vector tienen la misma tolerancia
      fac = 10^tol;
      for w = 1:n
        if floor(x(w)*fac) == floor(Z(w)*fac) % Sus primeros tol decimales son iguales
          conSP = conSP + 1;
        endif
      endfor
      if conSP == n
        funcionoSP(n, muestras) = tol;
      endif
    endfor
    
    clear('tol'); % Los puse separados para darles distintas cotas de tolerancia
    
    for tol = 1:30% Va a dejar el valor de la ultima tolerancia cumplida
      conCP = 0;  % Quiero que le asigne la tolerancia solo si todos los elementos del vector tienen la misma tolerancia
      fac = 10^tol;
      for w = 1:n
        if floor(x(w)*fac) == floor(X(w)*fac) % Sus primeros tol decimales son iguales
          conCP = conCP + 1;
        endif
      endfor
      if conCP == n
        funcionoCP(n, muestras) = tol;
      endif
    endfor
    
    % Agregar los errores a los vectores de errores
    sumCP = 0;
    sumSP = 0;
    for w = 1:n
      sumCP = (x(w)-X(w))**2 + sumCP;
      sumSP = (x(w)-Z(w))**2 + sumSP;
    endfor
    errorCP(i) = sqrt(sumCP);
    errorSP(i) = sqrt(sumSP);
  endfor 
endfor
funcionoSP(1,:) = 1;                %El tamaño más chico de matriz es 2x2
funcionoCP(1,:) = 1;                %El tamaño más chico de matriz es 2x2