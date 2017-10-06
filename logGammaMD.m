function lgmd = logGammaMD(x,d)

xv=x-0.5*(0:1:d-1);

lgmd = 0.5*d*log(pi)+sum(gammaln(xv));