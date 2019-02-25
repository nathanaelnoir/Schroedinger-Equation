//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//=========================HEAVISIDE_FUNCTION================================
function u = heaviside(t)
   u(t<0) = 0;
   u(t>=0) = 1;
endfunction
//=========================HEAVISIDE_FUNCTION================================
function y=diracdelta(d)
    y=0;
    if d==0 then
      y=1
    end
endfunction
//=============================OSSZILATOR====================================
// INPUT: - x: vector for solving problem in the interval -L < x < L  .
//        - m: mass of the particle.
//        - omega : angular frequency
// OUTPUT: - Psi: Normalized Eigenfunctions stored in a Matrix.
//         - E: Eingenvalues stored in a diagonal Matrix.
//---------------------------------------------------------------------------
function [Psi, E] = osszilator(x,omega,m)
   N = length(x)
   dx = x(2) - x(1);
   // POTENTIAL
   U = 1/2*m*omega^2*x.^(2)
   // DISKRETISIERTER LAPLACE OPERATOR
   Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1)+ ...
                                               diag(ones((N-1),1),-1))/(dx^2);
   // Lap so that it is consistent with f(0) = f(L) = 0
   Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; //% So that f(0) = 0
   Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0; //% So that f(L) = 0
   //  Total Hamiltonian
   hbar = 1
   H = -1/2*(hbar^2/m) * Lap + diag(U,0);
   [Psi, E] = spec(H)

   for i = 1:3
      Norm(i) = intsplin(x,abs(Psi(:,i))^2)
      Psi(:,i) = Psi(:,i)./sqrt(Norm(i))
   end
endfunction
//===============ATTRACTIVE DELTA POTENTIAL (BOUNDSTATE)======================
// INPUT: - x: vector for solving problem in the interval -L < x < L  .
//        - m: mass of the particle.
//        - omega : angular frequency
// OUTPUT: - Psi: Normalized Eigenfunctions stored in a Matrix.
//         - E: Eingenvalues stored in a diagonal Matrix.
//---------------------------------------------------------------------------
function [Psi, E] = deltapot(x,omega,m)
   N = length(x)
   dx = x(2) - x(1);
   // POTENTIAL

   for i=1:length(x)
        f(i)=diracdelta(x(i));
   end

   U = -20*f
   Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1)+ ...
                                               diag(ones((N-1),1),-1))/(dx^2);
   // Lap so that it is consistent with f(0) = f(L) = 0
   Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; //% So that f(0) = 0
   Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0; //% So that f(L) = 0
   //  Total Hamiltonian
   hbar = 1
   H = -1/2*(hbar^2/m) * Lap + diag(U,0);
   [Psi, E] = spec(H)

   for i = 1:3
    Norm(i) = intsplin(x,abs(Psi(:,i))^2)
    Psi(:,i) = Psi(:,i)./sqrt(Norm(i))
   end
endfunction
//========================FINITE POTENTIAL WELL==============================
// INPUT: - x: vector for solving problem in the interval -L < x < L  .
//        - m: mass of the particle.
//        - omega : angular frequency
// OUTPUT: - Psi: Normalized Eigenfunctions stored in a Matrix.
//         - E: Eingenvalues stored in a diagonal Matrix.
//---------------------------------------------------------------------------
function [Psi, E] = well(x,omega,m)
   N = length(x)
   dx = x(2) - x(1);
   // POTENTIAL
   w = L/3;
   U = -500*(heaviside(x+w)-heaviside(x-w))+500;
   // DISKRETISIERTER LAPLACE OPERATOR
   Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1)+ ...
                                               diag(ones((N-1),1),-1))/(dx^2);
   // Lap so that it is consistent with f(0) = f(L) = 0
   Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; //% So that f(0) = 0
   Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0; //% So that f(L) = 0
   //  Total Hamiltonian
   hbar = 1
   H = -1/2*(hbar^2/m) * Lap + diag(U,0);
   [Psi, E] = spec(H)

   for i = 1:3
    Norm(i) = intsplin(x,abs(Psi(:,i))^2)
    Psi(:,i) = Psi(:,i)./sqrt(Norm(i))
   end
endfunction
//======================HYDROGEN ATOM BONUS==================================
// INPUT: - x: vector for solving problem in the interval 0 < x < R  .
//        - m: mass of the particle.
//        - n: # of state. (1 is the groundstate)
//             keep in mind e.g. n = 2; l = {0,1} .....
//        - omega : angular frequency
// OUTPUT: - Psi: Normalized Eigenfunctions stored in a Matrix.
//         - E: Eingenvalues stored in a diagonal Matrix.
//---------------------------------------------------------------------------
function [Psi, E] = hydrogen(x,omega,m,n)
   l = 0
   A = 0
   for i = 1:n

   N = length(x)
   dx = x(2) - x(1);
   // POTENTIAL
   U =  (1./x)-((1/2).* l * (l+1)./(2.*x.^2));
   U(1)= 0
   // DISKRETISIERTER LAPLACE OPERATOR
   Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1)+ ...
                                               diag(ones((N-1),1),-1))/(dx^2);
   // Lap so that it is consistent with f(0) = f(L) = 0
   Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; //% So that f(0) = 0
   Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0; //% So that f(L) = 0
   //  Total Hamiltonian
   hbar = 1
   H = -1/2*(hbar^2/m) * Lap - diag(U,0);
   [Psi, E] = spec(H)

    [A,ind] = gsort(diag(E),"g","i")
     A(:,1,i) = A



      Psi(:,i)=Psi(:,i)./x'
      Psi(1,i)=0
      Norm(i) = intsplin(x,x'.^2.*abs(Psi(:,i))^2)

      Psi(:,i) = Psi(:,i)./sqrt(Norm(i))

  //disp(Psi(:,2))
  //s(intsplin(x,x'.^2.*abs(Psi(:,2))^2))
  l = l + 1
  end
endfunction
//=============================Style-of-plot==================================
function set_my_line_styles(style, thickness)
    e = gce();
    e.children.line_style = style;
    e.children.thickness = thickness;
endfunction
//===============================OUTPUT=======================================
function S = repeat(s,g)
   S = repmat(s, [1,g]) ;
endfunction
//===========================================================================//
//                               MAIN                                        //
//===========================================================================//
function bsp18()
   // Parameters for solving problem in the interval -L < x < L
   clf
   clear
   clc
   hbar = 1 // red. planck constant
   omega = 1
   m = 1
   L = 10;
   x = -L:0.02:L
   //------------------------------------------------------------------------//
   //                          OSSZILATOR                                    //
   //------------------------------------------------------------------------//
   [Psi, E] = osszilator(x,omega,m)
   //------------------ FIRST PLOT NORMALIZED WAVEFUNCTION ------------------//
   clf
   subplot(331)
      xgrid
      for i = 1:3
        col = ["r","b","m"]
        c = col(i)
        plot(x',Psi(:,i),c)
        set_my_line_styles(1,3);
      end
      tx = 'Normalized wavefunction for 1d harmonic oscillator'
      title(tx,'fontsize',4);
      hl=legend(['n = 0';'n = 1';'n = 2']);

      [E_osz,ind] = gsort(diag(E),"g","i")
   //-------------SECOND PLOT NORMALIZED PROBABILITY DENSITY-----------------//
    subplot(334)
      xgrid
      for i = 1:3
         col = ["r","b","m"]
         c = col(i)
         plot(x',abs(Psi(:,i))^2,c)
         set_my_line_styles(1,3);
      end
      tx = 'Probability density for 1d harmonic oscillator'
      title(tx,'fontsize',4);
      hl=legend(['n = 0';'n = 1';'n = 2']);
   //----------------------- THIRD PLOT POTENTIAL ---------------------------//
    subplot(337)
      xgrid
      plot(x',1/2*m*omega^2*x.^(2),"k")
      set_my_line_styles(1,3);
      tx = 'Quadratic potential'
      title(tx,'fontsize',4);
   //------------------------------------------------------------------------//
   //                          FINITE POTENTIAL WELL                         //
   //------------------------------------------------------------------------//
    [Psi, E] = well(x,omega,m)
   //------------------ FIRST PLOT NORMALIZED WAVEFUNCTION ------------------//
    subplot(332)
      xgrid
      for i = 1:3
         col = ["r","b","m"]
         c = col(i)
         plot(x',Psi(:,i),c)
         set_my_line_styles(1,3);
      end
      tx = 'Normalized wavefunction in finite square well'
      title(tx,'fontsize',4);
      hl=legend(['n = 1';'n = 2';'n = 3']);
      [E,ind] = gsort(diag(E),"g","i")
   //-------------SECOND PLOT NORMALIZED PROBABILITY DENSITY-----------------//
    subplot(335)
      xgrid
      for i = 1:3
         col = ["r","b","m"]
         c = col(i)
         plot(x',abs(Psi(:,i))^2,c)
         set_my_line_styles(1,3);
      end
      tx = 'Probability density in finite square well'
      title(tx,'fontsize',4);
      hl=legend(['n = 1';'n = 2';'n = 3']);
   //----------------------- THIRD PLOT POTENTIAL ---------------------------//
    subplot(338)
      xgrid
      w = L/3;
      plot(x',-500 * (heaviside(x+w)-heaviside(x-w)) + 500,"k")
      set_my_line_styles(1,3)
      tx = 'Finite potential well'
      title(tx,'fontsize',4);
   //------------------------------------------------------------------------//
   //                          DELTA POTENTIAL                               //
   //------------------------------------------------------------------------//
    [Psi, E] = deltapot(x,omega,m)
   //------------------ FIRST PLOT NORMALIZED WAVEFUNCTION ------------------//
    subplot(333)
      xgrid
      for i = 1:1
         col = ["r","b","m"]
         c = col(i)
         plot(x',Psi(:,i),c)
         set_my_line_styles(1,3);
      end
      tx = 'Normalized Wavefunction in attractive delta-potential'
      title(tx,'fontsize',4);
      hl=legend(['n = 1']);

      [E,ind] = gsort(diag(E),"g","i")

   //-------------SECOND PLOT NORMALIZED PROBABILITY DENSITY-----------------//
    subplot(336)
      xgrid
      for i = 1:1
         col = ["r","b","m"]
         c = col(i)
         plot(x',abs(Psi(:,i))^2,c)
         set_my_line_styles(1,3);
      end
      tx = 'Probability density in attractive delta-potential'
      title(tx,'fontsize',4);
      hl=legend(['n = 1']);
   //----------------------- THIRD PLOT POTENTIAL ---------------------------//
    subplot(339)
      xgrid
      for i=1:length(x)
        f(i)=diracdelta(x(i));
      end

      plot(x',-20*f,"k")
      set_my_line_styles(1,3);
      tx = 'Attractive delta-potential (boundstate))'
      title(tx,'fontsize',4);
   //------------------------------------------------------------------------//
   //                          HYDROGEN BONUS                                //
   //------------------------------------------------------------------------//
    R = 20;
    x = 0:0.04:R
    n = 1 // GROUNDSTATE
   [Psi, E_hyd] = hydrogen(x,omega,m,n)
   //------------------ FIRST PLOT NORMALIZED WAVEFUNCTION ------------------//
   scf()
   subplot(131)
      xgrid
      for i = 1:n
        plot(x',abs(Psi(:,i)),"r")
        set_my_line_styles(1,3);
      end
      n = 2
      [Psi, E] = hydrogen(x,omega,m,n)
      col = ["b","g"]
      for i = 1:n
        plot(x',Psi(:,i),col(i))
        set_my_line_styles(1,3);
      end
      tx = '$\text{Radial part of the Hydrogen wavefunction}$'
      title(tx,'fontsize',4);
      hl=legend(['n = 1, l = 0';'n = 2, l = 1';'n = 2, l = 0']);

      [E_hyd,ind] = gsort(diag(E),"g","i")
      format_1 = "\n%5.2f\n";
   //-------------SECOND PLOT NORMALIZED PROBABILITY DENSITY-----------------//
    subplot(132)
    n = 1
   [Psi, E] = hydrogen(x,omega,m,n)
      xgrid
      for i = 1:n
         plot(x',x'.^2.*abs(Psi(:,i))^2,"r")
         set_my_line_styles(1,3);
      end
     n = 2
    [Psi, E] = hydrogen(x,omega,m,n)
    col = ["b","g"]
    for i = 1:n
         plot(x',x'.^2.*abs(Psi(:,i))^2,col(i))
         set_my_line_styles(1,3);
      end
      tx = '$\text{Probability density}: r^{2} |R_{nl}(r)|^{2}$'
      title(tx,'fontsize',4);
      hl=legend(['n = 1, l = 0';'n = 2, l = 1';'n = 2, l = 0']);

   //----------------------- THIRD PLOT POTENTIAL ---------------------------//
    subplot(133)
      xgrid
      l = 0
        U =  -(1./x)+((1/2).* l*(l+1)./(2.*x.^2))+10
        U(1)= 0
        plot(x',U,"k")
        set_my_line_styles(1,3);

      tx = '$\text{Potential U for n = 1}$'
      title(tx,'fontsize',4);
      hl=legend(['l = 0']);
   //=============================PRINT=======================================//
      S = repmat("=",58);
      M = repmat("-",58);
      format_1 = "%s\";
      format_2 = "\n%10s\%18.6f\t%10.0f\n";
      format_3 = "\n%36s";
      format_4 = "\n%32.2f\n";


      mprintf(format_1,S);
      mprintf("\n\tExample: %s\n", "EF and EV for the Schrödinger Equation");
      mprintf(format_1,S);
      mprintf(format_3,"[\hbar \omega]");
      mprintf(format_4, E_osz(1:10));
      mprintf(format_1,M);


endfunction
bsp18()
//================================END========================================//
