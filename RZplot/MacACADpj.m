function P=MacACADpj(X)

global Acad

P = X;

if isfield(Acad,'ModP') 
   if Acad.ModP==1
   end
end

if isfield(Acad,'ModJ')
   if Acad.ModJ==1
      N  = round(length(X)/3);
      s  = X(1:N);
      Jo = X(2*N+1:end);
      Jn = Jo.*(Acad.ModJa1+Acad.ModJa2*s.^2).^Acad.ModJa3;
      P  = [X(1:2*N); Jn];      
   end
end
