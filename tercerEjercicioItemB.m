kappa1 = [];
kappa2 = [];

for n = 1:10   
  N = 8 * (n^3);
  
  % Genera matriz rala
  
  SD  = sparse(1:N,1:N,6*ones(N,1),N,N);
  if1 = setdiff(1:N-1, n:n:n*(8*n^2-1));
  SU1 = sparse(if1,if1+1,-1*ones(1,(n-1)*8*(n^2)),N,N);
  [i, j] = find([ones(n*(2*n-1),4*n); zeros(n,4*n)]);
  if2 = i+(j-1)*2*n^2;
  SU2 = sparse(if2,if2+n,-1*ones(1,n*(2*n-1)*4*n),N,N);
  if3 = 1:2*n^2*(4*n-1);
  SU3 = sparse(if3,if3+2*n^2,-1*ones(1,2*n^2*(4*n-1)),N,N);
  
  A = SD+SU1+SU1'+SU2+SU2'+SU3+SU3';
  
  [L, U] = ilu(A); % Factorización incompleta de A

  B = A'*A; % A transpuesta por A
  C = (inverse(M)*A*inverse(M))'*inverse(M)*A*inverse(M); %Debe haber alguna simplificación de esto
  
  SC = eig(B);
  CC = eig(C);
  
  kappa1(n) = max(SC)/min(SC);
  kappa2(n) = max(CC)/min(CC);
endfor
% https://octave.sourceforge.io/octave/function/pcg.html documentacion de pcg
% https://octave.sourceforge.io/octave/function/ilu.html documentacion de ilu
% https://octave.sourceforge.io/octave/function/lu.html  documentacion de lu
