function b = isintegervalue(x)
  b = isa(x,'integer') || (imag(x)==0 && mod(x,1)==0);
end