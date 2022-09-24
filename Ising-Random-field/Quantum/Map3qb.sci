

// ==================================================================
// Map3qb: caso general
// Kraus representation for the elemental 3 qubit operation
// ==================================================================
function [lst] = xMap3qb(bnk,T)   // modidificada 
  if (2*bnk<-4) then
  	P1  = getProjector(3, 0) + getProjector(3, 5 )+ ..
   			getProjector(3, 1) + getProjector(3, 4);
  	P0  = getProjector(3, 2) + getProjector(3, 3) + ..
        	getProjector(3, 6) + getProjector(3, 7);
//			disp(1)
  elseif (2*bnk>-4 & 2*bnk<0) then
  	P1  = getProjector(3, 2) + getProjector(3, 5)+ ..
  			getProjector(3, 1) + getProjector(3, 4);
  	P0  = getProjector(3, 0) + getProjector(3, 3) + ..
        	getProjector(3, 6) + getProjector(3, 7);
//			disp(2)
  elseif (2*bnk>0 & 2*bnk<4) then
  	P1  = getProjector(3, 2) + getProjector(3, 3)+ ..
  			getProjector(3,6) + getProjector(3, 5);
  	P0  = getProjector(3, 0) + getProjector(3, 1) + ..
        	getProjector(3, 4) + getProjector(3, 7);
//			disp(3)
  else
  	P1  = getProjector(3, 2) + getProjector(3, 7)+..
  			getProjector(3, 3) + getProjector(3, 6);
  	P0  = getProjector(3, 0) + getProjector(3, 1) + ..
        	getProjector(3, 4) + getProjector(3, 5);
//			disp(4)
  end

  
  id2 = eye(2, 2);
  //
  idSxid = kron( id2, kron(Sx(), id2) );  // esta coincide 
  idHid  = kron( id2, kron(H(),  id2) );  // esta no
  //idTid  = kron( id2, kron([al be ; -be al],  id2) );  // se cambio de halamar por  una matriz de transicion

  // caso general 
  t1=4; t2=2*bnk;
  be1=sqrt(1/(1+exp((t1+t2)/T)));  al1=sqrt(1-be1*be1);
  be2=sqrt(1/(1+exp(t2/T)));       al2=sqrt(1-be2*be2);
  be3=sqrt(1/(1+exp(-(t1+t2)/T))); al3=sqrt(1-be3*be3);
  be4=sqrt(1/(1+exp(-(t2)/T)));    al4=sqrt(1-be4*be4);
  be5=be2;                         al5=sqrt(1-be5*be5);
  be6=sqrt(1/(1+exp((-t1+t2)/T))); al6=sqrt(1-be6*be6);
  be7=be4;                         al7=sqrt(1-be7*be7);
  be8=sqrt(1/(1+exp((t1-t2)/T))) ; al8=sqrt(1-be8*be8);


 //P=exp()/(1+exp()) = 1/(exp(-)+exp())
  idTid  = [1  0   be3 0   0   0   0   0  ;
            0  al2 0   be4 0   0   0   0  ;
            0  0  -al3 0   0   0   0   0  ;
            0  be2 0  -al4 0   0   0   0  ;
            0  0   0   0   al5 0   be7 0  ;
            0  0   0   0   0   al6 0   0;
            0  0   0   0   be5 0  -al7 0  ;
            0  0   0   0   0   be6 0   1 ];
  
  //

   M1 = idSxid * P1;
   M2 = idHid * P0;
 // pause
  //
  lst = list(M1, M2);  
endfunction



// ==================================================================
// Map3qb: primer paso de cambio valores en la entrada
// Kraus representation for the elemental 3 qubit operation
// ==================================================================
function [lst] = iMap3qb(bnk,T)   // modidificada 
  
  // caso general 
  t1=4; t2=2*bnk;
  t3=1.0/sqrt(2);
  be1=sqrt(1/(1+exp((t1+t2)/T)));  al1=sqrt(1-be1*be1);
  be2=sqrt(1/(1+exp(t2/T)));       al2=sqrt(1-be2*be2);
  be3=sqrt(1/(1+exp(-(t1+t2)/T))); al3=sqrt(1-be3*be3);
  be4=sqrt(1/(1+exp(-(t2)/T)));    al4=sqrt(1-be4*be4);
  be5=be2;                         al5=sqrt(1-be5*be5);
  be6=sqrt(1/(1+exp((-t1+t2)/T))); al6=sqrt(1-be6*be6);
  be7=be4;                         al7=sqrt(1-be7*be7);
  be8=sqrt(1/(1+exp((t1-t2)/T))) ; al8=sqrt(1-be8*be8);

// P=exp()/(1+exp()) = 1/(exp(-)+exp())

  if  ( 2*bnk<0) then
      P1  = getProjector(3, 2) + getProjector(3, 5)+..
          getProjector(3, 1) + getProjector(3, 4);
      P0  = eye(8,8)-P1;
      // se cambio la columa 2 y 5 con temrminos t3  
      idTid  = [1  0   be3 0   0   0   0   0  ;
                0  t3  0   be4 0   0   0   0  ;
                0  0  -al3 0   0   0   0   0  ;
                0  t3  0  -al4 0   0   0   0  ;
                0  0   0   0   t3  0   be7 0  ;
                0  0   0   0   0   al6 0   0;
                0  0   0   0   t3  0  -al7 0  ;
                0  0   0   0   0   be6 0   1 ];
//    p1           x   x       x   x  


  else
      P1  = getProjector(3, 2) + getProjector(3, 5)+..
          getProjector(3,3) + getProjector(3, 6);
      P0  = eye(8,8)-P1;
      // se cambio la columna 4 y 7 con terminos t3
      idTid  = [1  0   be3 0   0   0   0   0  ;
                0 -al2 0   t3  0   0   0   0  ;
                0  0   al3 0   0   0   0   0  ;
                0  be2 0   t3  0   0   0   0  ;
                0  0   0   0  -al5 0   t3  0  ;
                0  0   0   0   0  -al6 0   0;
                0  0   0   0   be5 0   t3  0  ;
                0  0   0   0   0   be6 0   1 ];
//   p1                x   x       x   x 
  
  
  end

  id2 = eye(2, 2);
  //
  idSxid = kron( id2, kron(Sx(), id2) );  // esta coincide 
  //idHid  = kron( id2, kron(H(),  id2) );  // esta no


  //
   M1 =idSxid * P1;
   M2 = idTid * P0;
 // pause
  //
  lst = list(M1, M2);  
endfunction





// ==================================================================
// Map3qb: segundo paso de cambio singos 
// Kraus representation for the elemental 3 qubit operation
// ==================================================================
function [lst] = iiMap3qb(bnk,T)   // modidificada 
  
  // caso general 
  t1=4; t2=2*bnk;
  t3=1.0/sqrt(2);
  be1=sqrt(1/(1+exp((t1+t2)/T)));  al1=sqrt(1-be1*be1);
  be2=sqrt(1/(1+exp(t2/T)));       al2=sqrt(1-be2*be2);
  be3=sqrt(1/(1+exp(-(t1+t2)/T))); al3=sqrt(1-be3*be3);
  be4=sqrt(1/(1+exp(-(t2)/T)));    al4=sqrt(1-be4*be4);
  be5=be2;                         al5=sqrt(1-be5*be5);
  be6=sqrt(1/(1+exp((-t1+t2)/T))); al6=sqrt(1-be6*be6);
  be7=be4;                         al7=sqrt(1-be7*be7);
  be8=sqrt(1/(1+exp((t1-t2)/T))) ; al8=sqrt(1-be8*be8);

// P=exp()/(1+exp()) = 1/(exp(-)+exp())

  if  ( 2*bnk<0) then
      P1  = getProjector(3, 2) + getProjector(3, 5)+..
          getProjector(3, 1) + getProjector(3, 4);
      P0  = eye(8,8)-P1;
      // se cambio la columa 2 y 5 con temrminos t3  
      idTid  = [1  0   be3 0   0   0   0   0  ;
                0  t3  0   be4 0   0   0   0  ;
                0  0  -al3 0   0   0   0   0  ;
                0  t3  0  -al4 0   0   0   0  ;
                0  0   0   0   t3  0   be7 0  ;
                0  0   0   0   0   al6 0   0;
                0  0   0   0   t3  0  -al7 0  ;
                0  0   0   0   0   be6 0   1 ];
//    p1           x   x       x   x  


  else
      P1  = getProjector(3, 2) + getProjector(3, 5)+..
          getProjector(3,3) + getProjector(3, 6);
      P0  = eye(8,8)-P1;
      // se cambio la columna 4 y 7 con terminos t3
      idTid  = [1  0   be3 0   0   0   0   0  ;
                0  al2 0   t3  0   0   0   0  ;
                0  0  -al3 0   0   0   0   0  ;
                0  be2 0  -t3  0   0   0   0  ;
                0  0   0   0   al5 0   t3  0  ;
                0  0   0   0   0   al6 0   0;
                0  0   0   0   be5 0  -t3  0  ;
                0  0   0   0   0   be6 0   1 ];
//   p1                x   x       x   x 
  
  
  end

  id2 = eye(2, 2);
  //
  idSxid = kron( id2, kron(Sx(), id2) );  // esta coincide 
  //idHid  = kron( id2, kron(H(),  id2) );  // esta no


  //
   M1 =idSxid * P1;
   M2 = idTid * P0;
 // pause
  //
  lst = list(M1, M2);  
endfunction



// ==================================================================
// Map3qb: tercer paso cambio de proyectores a los de gustavo
// Kraus representation for the elemental 3 qubit operation
// ==================================================================
function [lst] = iiiMap3qb(bnk,T)   // modidificada 
  
  // caso general 
  t1=4; t2=2*bnk;
  t3=1.0/sqrt(2);
  be1=sqrt(1/(1+exp((t1+t2)/T)));  al1=sqrt(1-be1*be1);
  be2=sqrt(1/(1+exp(t2/T)));       al2=sqrt(1-be2*be2);
  be3=sqrt(1/(1+exp(-(t1+t2)/T))); al3=sqrt(1-be3*be3);
  be4=sqrt(1/(1+exp(-(t2)/T)));    al4=sqrt(1-be4*be4);
  be5=be2;                         al5=sqrt(1-be5*be5);
  be6=sqrt(1/(1+exp((-t1+t2)/T))); al6=sqrt(1-be6*be6);
  be7=be4;                         al7=sqrt(1-be7*be7);
  be8=sqrt(1/(1+exp((t1-t2)/T))) ; al8=sqrt(1-be8*be8);

// P=exp()/(1+exp()) = 1/(exp(-)+exp())

  if  ( 2*bnk<0) then
      P1  = getProjector(3, 2) + getProjector(3, 5);//+..
        //  getProjector(3, 1) + getProjector(3, 4);
      P0  = eye(8,8)-P1;
      // se cambio la columa 2 y 5 con temrminos t3  
      idTid  = [1  0   be3 0   0   0   0   0  ;
                0  t3  0   be4 0   0   0   0  ;
                0  0  -al3 0   0   0   0   0  ;
                0  t3  0  -al4 0   0   0   0  ;
                0  0   0   0   t3  0   be7 0  ;
                0  0   0   0   0   al6 0   0;
                0  0   0   0   t3  0  -al7 0  ;
                0  0   0   0   0   be6 0   1 ];
//    p1               x           x  


  else
      P1  = getProjector(3, 2) + getProjector(3, 5);//+..
        //  getProjector(3,3) + getProjector(3, 6);
      P0  = eye(8,8)-P1;
      // se cambio la columna 4 y 7 con terminos t3
      idTid  = [1  0   be3 0   0   0   0   0  ;
                0  al2 0   t3  0   0   0   0  ;
                0  0  -al3 0   0   0   0   0  ;
                0  be2 0  -t3  0   0   0   0  ;
                0  0   0   0   al5 0   t3  0  ;
                0  0   0   0   0   al6 0   0;
                0  0   0   0   be5 0  -t3  0  ;
                0  0   0   0   0   be6 0   1 ];
//   p1                x           x     
  
  
  end

  id2 = eye(2, 2);
  //
  idSxid = kron( id2, kron(Sx(), id2) );  // esta coincide 
  //idHid  = kron( id2, kron(H(),  id2) );  // esta no


  //
   M1 =idSxid * P1;
   M2 = idTid * P0;
 // pause
  //
  lst = list(M1, M2);  
endfunction





// ==================================================================
// Map3qb: cuarto se ajustan los valores en todas las columnas
// Kraus representation for the elemental 3 qubit operation
// ==================================================================
function [lst] = ivMap3qb(bnk,T)   // modidificada 
  
  // caso general 
  t1=4; t2=2*bnk;
  t3=1.0/sqrt(2);
  be1=sqrt(1/(1+exp((t1+t2)/T)));  al1=sqrt(1-be1*be1);
  be2=sqrt(1/(1+exp(t2/T)));       al2=sqrt(1-be2*be2);
  be3=sqrt(1/(1+exp(-(t1+t2)/T))); al3=sqrt(1-be3*be3);
  be4=sqrt(1/(1+exp(-(t2)/T)));    al4=sqrt(1-be4*be4);
  be5=be2;                         al5=sqrt(1-be5*be5);
  be6=sqrt(1/(1+exp((-t1+t2)/T))); al6=sqrt(1-be6*be6);
  be7=be4;                         al7=sqrt(1-be7*be7);
  be8=sqrt(1/(1+exp((t1-t2)/T))) ; al8=sqrt(1-be8*be8);

// P=exp()/(1+exp()) = 1/(exp(-)+exp())

  if  ( 2*bnk<0) then
      P1  = getProjector(3, 2) + getProjector(3, 5);//+..
        //  getProjector(3, 1) + getProjector(3, 4);
      P0  = eye(8,8)-P1;
      // se cambio la columa 2 y 5 con temrminos t3  
      idTid  = [1  0   be3 0   0   0   0   0  ;
                0  t3  0   t3  0   0   0   0  ;
                0  0  -al3 0   0   0   0   0  ;
                0  t3  0  -t3  0   0   0   0  ;
                0  0   0   0   t3  0   t3  0  ;
                0  0   0   0   0   al6 0   0;
                0  0   0   0   t3  0  -t3  0  ;
                0  0   0   0   0   be6 0   1 ];
//    p1               x           x  


  else
      P1  = getProjector(3, 2) + getProjector(3, 5);//+..
        //  getProjector(3,3) + getProjector(3, 6);
      P0  = eye(8,8)-P1;
      // se cambio la columna 4 y 7 con terminos t3
      idTid  = [1  0   be3 0   0   0   0   0  ;
                0  t3  0   t3  0   0   0   0  ;
                0  0  -al3 0   0   0   0   0  ;
                0  t3  0  -t3  0   0   0   0  ;
                0  0   0   0   t3  0   t3  0  ;
                0  0   0   0   0   al6 0   0;
                0  0   0   0   t3  0  -t3  0  ;
                0  0   0   0   0   be6 0   1 ];
//   p1                x           x     
  
  
  end

  id2 = eye(2, 2);
  //
  idSxid = kron( id2, kron(Sx(), id2) );  // esta coincide 
  //idHid  = kron( id2, kron(H(),  id2) );  // esta no


  //
   M1 =idSxid * P1;
   M2 = idTid * P0;
 // pause
  //
  lst = list(M1, M2);  
endfunction



// ==================================================================
// Map3qb: casos j 
// Kraus representation for the elemental 3 qubit operation
// ==================================================================
function [lst] = jMap3qb(bnk,T)   // modidificada 
  
  // caso general 
  t1=4; t2=2*bnk;
  t3=1.0/sqrt(2);
  be1=sqrt(1/(1+exp((t1+t2)/T)));  al1=sqrt(1-be1*be1);
  be2=sqrt(1/(1+exp(t2/T)));       al2=sqrt(1-be2*be2);
  be3=sqrt(1/(1+exp(-(t1+t2)/T))); al3=sqrt(1-be3*be3);
  be4=sqrt(1/(1+exp(-(t2)/T)));    al4=sqrt(1-be4*be4);
  be5=be2;                         al5=sqrt(1-be5*be5);
  be6=sqrt(1/(1+exp((-t1+t2)/T))); al6=sqrt(1-be6*be6);
  be7=be4;                         al7=sqrt(1-be7*be7);
  be8=sqrt(1/(1+exp((t1-t2)/T))) ; al8=sqrt(1-be8*be8);

// P=exp()/(1+exp()) = 1/(exp(-)+exp())

  if  ( 2*bnk<0) then
      P1  = getProjector(3, 2) + getProjector(3, 5)+getProjector(3, 4)+getProjector(3, 1);
      P0  = eye(8,8)-P1;
      idTid  = [1  0   t3  0   0   0   0   0  ;
                0  t3  0   t3  0   0   0   0  ;
                0  0  -t3  0   0   0   0   0  ;
                0  t3  0  -t3  0   0   0   0  ;
                0  0   0   0   t3  0   t3  0  ;
                0  0   0   0   0   t3  0   0;
                0  0   0   0   t3  0  -t3  0  ;
                0  0   0   0   0   t3  0   1 ];
//    p1  j            x       x   x  
//    p1  jj           x           x  


  else
      P1  = getProjector(3, 2) + getProjector(3, 5)+getProjector(3, 6)+getProjector(3, 3);
      P0  = eye(8,8)-P1;
      idTid  = [1  0   t3  0   0   0   0   0  ;
                0  t3  0   t3  0   0   0   0  ;
                0  0  -t3  0   0   0   0   0  ;
                0  t3  0  -t3  0   0   0   0  ;
                0  0   0   0   t3  0   t3  0  ;
                0  0   0   0   0   t3  0   0;
                0  0   0   0   t3  0  -t3  0  ;
                0  0   0   0   0   t3  0   1 ];
//   p1 j              x           x   x  
//   p1 jjj            x           x     
  
  
  end

  id2 = eye(2, 2);
  //
  idSxid = kron( id2, kron(Sx(), id2) );  // esta coincide 
  //idHid  = kron( id2, kron(H(),  id2) );  // esta no


  //
   M1 =idSxid * P1;
   M2 = idTid * P0;
 // pause
  //
  lst = list(M1, M2);  
endfunction




// ==================================================================
// Map3qb: caso2 G U, H P 
// Kraus representation for the elemental 3 qubit operation
// ==================================================================
function [lst] = xxxMap3qb(bnk,T)   // modidificada 
  if (2*bnk<0) then
  	P0  = getProjector(3, 2) + getProjector(3, 5)+getProjector(3,1);//+ ..
  		//	getProjector(3, 1) + getProjector(3, 4);
  	P1  = eye(8,8)-P0;// getProjector(3, 0) + getProjector(3, 3) + ..
         //	getProjector(3, 6) + getProjector(3, 7);
//			disp(2)
  else
  	P0  = getProjector(3, 2) + getProjector(3, 5);//+getProjector(3,3);//+ ..
  		//getProjector(3,6) + getProjector(3, 3);
  	P1  = eye(8,8)-P0;// getProjector(3, 0) + getProjector(3, 1) + ..
        	//getProjector(3, 4) + getProjector(3, 7);
//			disp(3)
  end

  
  id2 = eye(2, 2);
  //
  idSxid = kron( id2, kron(Sx(), id2) );  // esta coincide 
  //idHid  = kron( id2, kron(H(),  id2) );  // esta no
  //idTid  = kron( id2, kron([al be ; -be al],  id2) );  // se cambio de halamar por  una matriz de transicion
  
  // caso con B<<J
  t1=4; t2=2*bnk;
  be1=sqrt(1/(1+exp((t1)/T)));  al1=sqrt(1-be1*be1);
  be2=sqrt(1/(1+exp(t2/T)));       al2=sqrt(1-be2*be2);
  be3=sqrt(1/(1+exp(-(t1)/T))); al3=sqrt(1-be3*be3);
  be4=sqrt(1/(1+exp(-(t2)/T)));    al4=sqrt(1-be4*be4);
  be5=be2;                         al5=sqrt(1-be5*be5);
  be6=sqrt(1/(1+exp((-t1)/T))); al6=sqrt(1-be6*be6);
  be7=be4;                         al7=sqrt(1-be7*be7);
  be8=sqrt(1/(1+exp((t1)/T))) ; al8=sqrt(1-be8*be8);

//
// P=exp()/(1+exp()) = 1/(exp(-)+1)
  t3=1/sqrt(2.0);
  idTid  = [1   0   1   0   0   0   0   0  ;
  			0   t3  0   t3  0   0   0   0  ;
			0   0   0   0   0   0   0   0  ;
			0   t3  0  -t3  0   0   0   0  ;
			0   0   0   0   t3  0   t3  0  ;
			0   0   0   0   0   0   0   0  ;
			0   0   0   0   t3  0  -t3  0  ;
			0   0   0   0   0   1   0   1 ];
  
  //

   M1 = idSxid * P0;
   M2 = idTid * P1;
 // pause
  //
  lst = list(M1, M2);  
endfunction




// ==================================================================
// Map3qb: caso nuevo con nk=0
// Kraus representation for the elemental 3 qubit operation
// ==================================================================
function [lst] = Map3qb(bnk,T)   // modidificada 
 
  // caso general
  if (abs(bnk)>0) then
  t1=4; t2=2*bnk; t3=1/sqrt(2);
  be1=sqrt(1/(1+exp((t1+t2)/T)));  al1=sqrt(1-be1*be1);
  be2=sqrt(1/(1+exp(t2/T)));       al2=sqrt(1-be2*be2);
  be3=sqrt(1/(1+exp(-(t1+t2)/T))); al3=sqrt(1-be3*be3);
  be4=sqrt(1/(1+exp(-(t2)/T)));    al4=sqrt(1-be4*be4);
  be5=be2;                         al5=sqrt(1-be5*be5);
  be6=sqrt(1/(1+exp((-t1+t2)/T))); al6=sqrt(1-be6*be6);
  be7=be4;                         al7=sqrt(1-be7*be7);
  be8=sqrt(1/(1+exp((t1-t2)/T))) ; al8=sqrt(1-be8*be8);
  end
 
  id2 = eye(2, 2);

 P00=zeros(8,8);
 if (2*bnk<-4) then  // si entra aqui truena
  	P1  = getProjector(3, 0) + getProjector(3, 5 )+ ..
   			getProjector(3, 1) + getProjector(3, 4);
  	P0  = getProjector(3, 2) + getProjector(3, 3) + ..
        	getProjector(3, 6) + getProjector(3, 7);
 



  elseif (2*bnk>-4 & 2*bnk<0) then
  	P1  = getProjector(3, 2) + getProjector(3, 5);//+ ..
  	P0  = getProjector(3, 0) + getProjector(3, 7) ;
	P00 = getProjector(3, 6) + getProjector(3, 3);
 
 	P1 = P1+	getProjector(3, 1) + getProjector(3, 4);
//			disp(2)
    
	
	idTid=kron(id2,kron([al4,be4;-be4,al4],id2))
  
//   idTid  = [1  0   1   0   0   0   0   0  ;
//             0  al4 0   be4 0   0   0   0  ;
//             0  0   0   0   0   0   0   0  ;
//             0  be4 0   al4 0   0   0   0  ;
//             0  0   0   0   al7 0   be7 0  ;
//             0  0   0   0   0   0   0   0;
//             0  0   0   0   be7 0  al7 0  ;
//             0  0   0   0   0   1   0   1 ];
  
  M2 = idTid*P00+P0;



  elseif (2*bnk>0 & 2*bnk<4) then
  	P1  = getProjector(3, 2) + getProjector(3, 5);//+ ..
  	P0  = getProjector(3, 0) + getProjector(3, 7);
	P00 = getProjector(3, 4) + getProjector(3, 1);
	
  	P1  = P1 + getProjector(3,6) + getProjector(3, 3);
//			disp(3)

	idTid = kron(id2,kron([al2,-be2;be2,al2],id2)) ;
//  idTid  = [1  0   1   0   0   0   0   0  ;
//            0  al2 0   be2 0   0   0   0  ;
//            0  0   0   0   0   0   0   0  ;
//            0  be2 0   al2 0   0   0   0  ;
//            0  0   0   0   al5 0   be5 0  ;
//            0  0   0   0   0   0   0   0;
//            0  0   0   0   be5 0   al5 0  ;
//            0  0   0   0   0   1   0   1 ];
  
    M2 = idTid*P00+P0;

  elseif (2*bnk>4) // si entra aqui truena
  	P1  = getProjector(3, 2) + getProjector(3, 7)+..
  			getProjector(3, 3) + getProjector(3, 6);
  	P0  = getProjector(3, 0) + getProjector(3, 1) + ..
        	getProjector(3, 4) + getProjector(3, 5);
//			disp(4)\
  else(bnk==0)
  	P1  = getProjector(3, 2) + getProjector(3, 5);
  	P0  = getProjector(3, 0) + getProjector(3, 7)
	P00=  getProjector(3, 1) + getProjector(3, 3) +..
        	getProjector(3, 4) + getProjector(3, 6);
    
	
//   idTid= kron( id2, kron(H(),  id2) );
  idTid= kron( id2, kron([1 1 ; -1 1]*sqrt(1/2),  id2) );
    M2 = P0+idTid*P00;
   
  end

  
  //
  idSxid = kron( id2, kron([0, 1 ; 1 0], id2) );  // esta coincide 
  //idHid  = kron( id2, kron(H(),  id2) );  // esta no
  //idTid  = kron( id2, kron([al be ; -be al],  id2) );  // se cambio de halamar por  una matriz de transicion


 //P=exp()/(1+exp()) = 1/(exp(-)+exp())
  // idTid  = [1  0   be3 0   0   0   0   0  ;
  //          0  al2 0   be4 0   0   0   0  ;
  //          0  0  -al3 0   0   0   0   0  ;
  //          0  be2 0  -al4 0   0   0   0  ;
  //          0  0   0   0   al5 0   be7 0  ;
  //          0  0   0   0   0   al6 0   0;
  //          0  0   0   0   be5 0  -al7 0  ;
  //          0  0   0   0   0   be6 0   1 ];
  
  //

   M1 = idSxid * P1;
   //M2 = idHid * P0;
  //pause
  //
  lst = list(M1, M2);  
endfunction



