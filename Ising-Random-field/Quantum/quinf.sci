//
// Libreria para información cuántica                   iniciado: 28. Ene 2015
//                                                       revised: 21. Dic 2015
//                                                                23. Nov 2016
//
//   qb: single qubit (state, operation, ..)
//   sa: state array (register of qubits)

 function [qb]= b(i)
   if (i == 0) then
     qb= [1; 0]
   else
     qb= [0; 1]
   end
 endfunction


//
// Generate a random (real or complex) n-qubit(s) state
// Set the root of the random number generator in the main program with
// isee= 1234; grand('setsd',isee);

 function [qb]= Ranv(n,tstr)
   [lhs,rhs]= argn(0);                     // if tstr is missing, set output
   if rhs == 1 then                        // type to 'real'
     tstr = '';
   end

   nn= 2**n;
   if (tstr == 'cmplx') then               // complex random state
     rqb= grand(2*nn,1,'nor',0,2);
     s= sqrt(sum(rqb.*rqb));
     rqb= rqb/s
     qb= rqb(1:nn) + %i*rqb(nn+1:2*nn);
   else                                    // real random state
     rqb= grand(nn,1,'nor',0,1);
     s= sqrt(sum(rqb.*rqb));
     qb= rqb/s
   end
 endfunction


//
// Generate Pauli Matrices sigma_i for i=0,..,3    Luis Camacho
// Note, the Not gate is just sigma_x = Pauli(1)
 function [pm]= Pauli(i)
   select i
     case 0 then
       pm=[1,0;0,1];
     case 1 then
       pm=[0,1;1,0];
     case 2 then
       pm=[0,-%i;%i,0];
     case 3 then
       pm=[1,0;0,-1];
     else
       error("Wrong input");
   end
 endfunction


//
// Apply Pauli matrices to a qubit state           Luis Camacho

 function [pro]= aPauli(qb,i)
    pro=Pauli(i)*qb;
 endfunction


//
// Apply NOT gate to a qubit (we can use the Pauli matrix sigma_x)

 function [qbr]= aNot(qb)
   qbr= aPauli(qb,1);
 endfunction


//
// Generate single qubit Hadamard gate

 function [pm]= H()
   isq2= 1.0/sqrt(2.0);
   pm=[isq2, isq2; isq2, -isq2];
 endfunction


//
// Apply Hadamard gate

 function [qbr]= aH(qb)
   qbr= H()*qb;
 endfunction



//
// Here comes the first two-qubit gate -- a control not gate. In such a case
// the order of the qubits and the basis states becomes important. Let 
// (a_1,a_2) with a_k= {0,1} be a product basis. Then the standard ordering 
// will be as follows:
//
//               (0,0) (0,1) (1,0) (1,1)   -->   |00>  |01>  |10>  |11>
//    ordering :   1     2     3     4
//
// Then instead of refering to the first, second, .. qubit, we should refer to
// the high-order (leftmost) or low-order (rightmost) qubit. 

// Generate a control-not gate (the leftmost qubit is the control, the 
// low-order the target qubit.
//
  function [cn]= CNot()
    cn= [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
  endfunction

// Control-not gate, where the low-order qubit is the control qubit 
//
  function [cn]= NotC()
    cn= [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0];
  endfunction

// Swap gate
//
  function [sw]= Swap()
    sw= [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
  endfunction


// Tensor products:
// To obtain tensor products of states, or of matrices, we may simply use
// the scilab function "kron(A,B)" which generates a Matrix, with 
// block-entries B multiplied by scalar A(i,j).

//
// Generate a n-qubit register, with canonical basis state  1 <= bn <= 2^n
// This function is for educational purposes only. It really generates (in a
// complicated way) a canonical basis state with all components equal to zero,
// except for the bn'th component which is set equal to one.
//
 function [qa]= qrbas(n,bn)
   b= bn-1;             
   qa= 1;                    // n is the number of qubits of the register
   for i=n:-1:1              // bn is the number of the basis state starting
     if (b < 2**(i-1)) then  // with 1
       qb= [1; 0];
     else
       qb= [0; 1];
       b= b - 2**(i-1);
     end
     qa= kron(qa,qb);
   end
 endfunction

//
// Generate a unitary matrix which applies the single qubit gate umat, to the
// qubit 1 <= i <= n in a n-qubit register. The ordering is from right 
// (low-order) to left (high-order).

 function [qamat]= qrGate1(umat,i,n)
   L= eye(2**(n-i),2**(n-i));       // Identity in the n-i left-most qubits
   R= eye(2**(i-1),2**(i-1));       // Identity in the i-1 right-most qubits
   qamat= kron(L,kron(umat,R));
 endfunction

//
// Generate a unitary matrix which applies a two-qubit gate to the qubits
// i and i+1, in a n-qubit register

 function [qamat]= qrGate2(umat,i,n)
   L= eye(2**(n-i-1),2**(n-i-1));
   R= eye(2**(i-1),2**(i-1));
   qamat= kron(L,kron(umat,R));
 endfunction


// Apply arbitrary two-qubit gate between qubits it1 and it2 in a n-qubit
// register. In the standard situation, it2 is assumed to be the higher order 
// qubit, and it1 the lower order qubit, such that the qubit it2 appears to the
// left of the qubit it1, as follows: 
//
//    qb_n   qb_n-1   ..   qb_it2   ..   qb_it1   ..   qb_1
//
// If this is not the case, we apply a swap operation to the 2-qubit gate, to
// arrive at the standard situation
//
 function [qamat]= qrGate3(u2mat,it2,it1,n)

   if it1 > it2 then
     um= Swap()*u2mat*Swap();
   else
     um= u2mat;
   end

   ima= max(it1,it2);   // the high-order qubit to the left
   imi= min(it1,it2);   // the low-order qubit to the right

   if (ima > imi+1) then              // Generate a unitary transformation ums,
     ums= qrGate2(Swap(),imi,n);      // which applies Swap gates in order to
     for itj=imi+2:ima-1              // qubit ima to the position imi+1
       ums= qrGate2(Swap(),itj,n)*ums;
     end
   else
     ums= eye(2^n,2^n);
   end
   qamat= ums'*qrGate2(um,ima-1,n)*ums;
 endfunction


// Apply a partial trace over the last qubit
//
 function [rrho]= ptrace(rho)

   nn= size(rho,1);
   if size(rho,2) ~= nn then
     error("input matrix is not quadratic");
   elseif modulo(nn,2) == 1 then
     error("dimension of input matrix is not even");
   end 

   nnh= nn/2;
   rrho= rho(1:nnh,1:nnh) + rho(nnh+1:nn,nnh+1:nn);
 endfunction


// Perform a conditional measurement on the highest order qubit
// This function returns the normalized result of projecting the state qain 
// onto the subspace corresponding to the highest order qubit being in state 
// |b0>. It also gives the probability of obtaining this result.
//
 function [qa,p]= condmeas(qain,b0)

   qa= qain;

   nnh= length(qain)/2;
   if (b0 == 0) then
     p= sum(abs(qain(1:nnh)).^2);
     qa(nnh+1:2*nnh)= 0;
   else
     p= sum(abs(qain(nnh+1:2*nnh)).^2);
     qa(1:nnh)= 0;
   end

   qa= qa/sqrt(p);
 endfunction


// Similar to the above conditional measurement, perform a conditional 
// partial trace.
//
 function [qa,p]= cptrace(qain,b0)

   nnh= length(qain)/2;
   if (b0 == 0) then
     qa= qain(1:nnh);
   else
     qa= qain(nnh+1,2*nnh);
   end

   p= sum(abs(qa).^2);
   qa= qa/sqrt(p);
 endfunction
 
 
// ----------------------------------------------------------------------------
// Quantum optics section
// ----------------------------------------------------------------------------

// Generate matrix representation of the raising and lowering operators of the
// Harmonic oscillator. Note, n is the number of modes taken into account. 
// Thus, the highest mode taken into account is n-1.
//
 function [ama]= ho_lower(n)
   ama= zeros(n,n);
   for j=1:n-1
     ama(j,1+j)= sqrt(j);
   end
 endfunction

 
 
// ----------------------------------------------------------------------------
// RMT section
// ----------------------------------------------------------------------------

// Random SU(2) matrices, centered at the identity with normally distributed
// Euler angles

 function [U] = SU2_normal(sig_a, sig_b, sig_g)
   alh= 0.5*grand(1,1,'nor',0,sig_a);
   beh= 0.5*grand(1,1,'nor',0,sig_b); 
   gah= 0.5*grand(1,1,'nor',0,sig_g); 

   U= [ exp(-%i*(alh+gah))*cos(beh) , -exp(-%i*(alh-gah))*sin(beh); ..
        exp(%i*(alh-gah))*sin(beh) , exp(%i*(alh+gah))*cos(beh) ];
 endfunction


// Generate GOE matrix with semi-circle spectrum and different normalizations
// In general, the level density is
//                                   rho(x)= 2N/(pi a^2) sqrt(a^2 - x^2)
// Let <H_ij^2> = sd^2 , then
//
//   sd = 1           :  a= 2 sqrt(N) ,  rho(0)= sqrt(N)/pi
//        1/sqrt(N)   :  a= 2         ,
//        sqrt(N)/pi  :  a= 2N/pi     ,  rho(0)= 1

 function [H]= GOEmember(n,sd)             // <H_ij^2> = 1 for i<>j
   A= grand(n,n,'nor',0,sd/sqrt(2.0));     // level density with a= 2 sqrt(N)
   H= A + A';
 endfunction
 
 function [H]= GOE_unit_sd(n)
   H= GOEmember(n,1);
 endfunction
 
 function [H]= GOE_for_semicircle(n)  // Defines GOE matrices, such that the 
   H= GOEmember(n,1/sqrt(n));         // semi-circle spectrum lies in (-2,2)
 endfunction

 function [H]= GOE_unit_center(n)     // Defines GOE matrices, with unit level 
   H= GOEmember(n,sqrt(n)/%pi);       // level spacing in the center
 endfunction                   

  

 function [H]= GUEmember(n,sd)             // <H_ij^2> = 1 for i<>j
   A= grand(n,n,'nor',0,sd/sqrt(2.0));     // level density with a= 2 sqrt(N)
   for i=1:n
     for j=1:i-1
       H(i,j)= A(i,j) + %i*A(j,i);
       H(j,i)= A(i,j) - %i*A(j,i);
     end
     H(i,i)= sqrt(2.0)*A(i,i);
   end
 endfunction
 
 function [H]= GUE_unit_sd(n)
   H= GUEmember(n,1);
 endfunction
 
 function [H]= GUE_for_semicircle(n)  // Defines GUE matrices, such that the 
   H= GUEmember(n,1/sqrt(n));         // semi-circle spectrum lies in (-2,2)
 endfunction

 function [H]= GUE_unit_center(n)     // Defines GUE matrices, with unit level 
   H= GUEmember(n,sqrt(n)/%pi);       // level spacing in the center
 endfunction                   

  
