
// ESTA FUNCINO YA NO SE USA
// Dado uniforme con valores dimension n a m
// que va del valor entero vmin al vmaz     -- sin validacion de entradas
function ir=intrand(n,m,vmin,vmax)
    ir= floor(rand(n,m)*(vmax-vmin+1)+vmin);
endfunction// cambiar esta función por    grand(n,m,"uin",vmin,vmax)


// con un arreglo (arrayb) regresa un arreglo
// de tamaño n con una ditrifucion uniforme  --sv 
function ar=arrayrand(n,arrayb)
    indices=grand(1,n,"uin",1,length(arrayb) ) //intrand(1,n,1,length(arrayb));
    ar=arrayb( indices );
endfunction

// Calcula la energia dados los spines y los pesos 
function H=Energia(spines,J,pesos)
    n=length(spines);  ind=[n,1:n-1];
    H=-J*spines*spines(ind)' + spines*pesos';
endfunction


//Genera variante
function S1=variaS(spines)
    irand=grand(1,1,"uin",1, length(spines) ) //intrand(1,1,1,n);
    spines(irand)=-1*spines(irand);
    S1=spines
endfunction


//Cuenta los bloques de S
function nb=ndominios(S) 
    t=diff(S); indf=length(S);
    n=length(t(abs(t)==2))
    nb= max(n + bool2s( S(1)~=S(indf)),1)
//    if(nb==0);nb=1;end
endfunction


//Calcula el cambio de la energia variando un spin
function dE=CambioE(spines,J,pesos,ncambio)
    dE=2*J*(spiens(ncambio))
endfuction
