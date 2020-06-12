module vers
  character(len=50):: version='GIBBSTHUR - 30 April 2020'
  contains
    subroutine print_version(start)
    implicit none
    logical:: start
    
    
    print *,'-----------------'
    print *,trim(version)
    print *,'by L Varona, A Legarra'
    print *,'-----------------'
    if(start) print *,'started: '
    call printtime2
    end subroutine

      subroutine printtime2
      INTEGER  DATE_TIME (8)
      CHARACTER(LEN = 12) REAL_CLOCK (3)

      CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                      REAL_CLOCK (3), DATE_TIME)

      print *,'date: ',real_clock(1)(7:8),'/',real_clock(1)(5:6),&
      '/',real_clock(1)(1:4)
      print *,'time: ',&
      real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)

      end subroutine

  
end module vers

module kinds
! From Ignacy Misztal BLUPF90 distribution 
  integer, parameter :: single = SELECTED_REAL_KIND( 6, 37 )
  integer, parameter :: double = SELECTED_REAL_KIND( 15, 307 )
  integer, parameter :: extended = SELECTED_REAL_KIND( 18, 4931 )
  integer, parameter :: r4 = SELECTED_REAL_KIND( 6, 37 )
  integer, parameter :: r8 = SELECTED_REAL_KIND( 15, 307 )
  integer, parameter :: r16 = SELECTED_REAL_KIND( 18, 4931 )

  ! current precison for hash storage
  integer, parameter :: rh=r8
end module kinds

module modelo
  implicit none
      integer :: &
      maxrec=5000000,& !this is the number of non-zero elements in the MME
     iguard=1000000000 !never write continuation file
! total number of iterations : imue*icad
! burn-in: lap*icad
! thin interval or lag: icad
! number of samples used in estimating means and sd: (imue-lap)
! This is later read from par file

end module modelo

module listaligada
  use kinds
  USE modelo
  implicit none
  real(r8),allocatable:: zhz(:)
  integer,allocatable:: ifirst(:),ivcol(:),inext(:)
  integer :: nplace
 
  contains

SUBROUTINE Bubble_Sort(a,b,n)
  integer :: a(n,2)
  integer :: b(n)
  integer :: temp(3)
  INTEGER :: i, j, n
  LOGICAL :: swapped
  DO j = n-1, 1, -1
    swapped = .FALSE.
    DO i = 1, j
      IF (a(i,1).gt.a(i+1,1)) THEN 
	temp(1) = a(i,1)
	temp(2) = a(i,2)
	temp(3) = b(i)
        a(i,1) = a(i+1,1)
	a(i,2) = a(i+1,2)
	b(i) = b(i+1)
        a(i+1,1) = temp(1)
	a(i+1,2) = temp(2)
	b(i+1) = temp(3) 
        swapped = .TRUE.
      ENDIF
      IF (a(i,1).eq.a(i+1,1)) THEN
	IF (a(i,2).gt.a(i+1,2)) THEN   
	temp(1) = a(i,1)
	temp(2) = a(i,2)
	temp(3) = b(i)
        a(i,1) = a(i+1,1)
	a(i,2) = a(i+1,2)
	b(i) = b(i+1)
        a(i+1,1) = temp(1)
	a(i+1,2) = temp(2)
	b(i+1) = temp(3) 
        swapped = .TRUE. 
	END IF
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE 

!=======================================================================
      SUBROUTINE LINKS(IROW,ICOL,D)

!     CONSTRUCCION DE LISTAS LIGADAS POR FILAS
!     CON ORDEN CRECIENTE POR COLUMNAS
!     NO SE REQUIERE NINGUN ORDEN EN LOS DATOS


	  double precision:: D
	  integer:: ipre,iplace,irow,icol

      IPRE=0
      IPLACE=IFIRST(IROW)
4     IF(IPLACE.GT.0)THEN
        IF(IVCOL(IPLACE).GE.ICOL)THEN
          IF(IVCOL(IPLACE).EQ.ICOL)THEN
            ZHZ(IPLACE)=ZHZ(IPLACE)+D
            RETURN
          ELSE
            NPLACE=NPLACE+1
            IF(NPLACE.GT.MAXREC)THEN
              PRINT *,'There is no space in the linked list'
	      print *,'increasing linked list by 50%'
	      call increase_linked_list()
            END IF
            IF(IPRE.EQ.0)THEN
              INEXT(NPLACE)=IFIRST(IROW)
              IFIRST(IROW)=NPLACE
            ELSE
              INEXT(NPLACE)=INEXT(IPRE)
              INEXT(IPRE)=NPLACE
            END IF
            ZHZ(NPLACE)=D
            IVCOL(NPLACE)=ICOL
            RETURN
          END IF
        ELSE
          IPRE=IPLACE
          IPLACE=INEXT(IPLACE)
          GOTO 4
        END IF
      ELSE
        NPLACE=NPLACE+1
        IF(NPLACE.GT.MAXREC)THEN
              PRINT *,'There is no space in the linked list'
	      print *,'increasing linked list by 50%'
	      call increase_linked_list()
        END IF
        IF(IFIRST(IROW).GT.0)THEN
          INEXT(IPRE)=NPLACE
        ELSE
          IFIRST(IROW)=NPLACE
        END IF
        ZHZ(NPLACE)=D
        IVCOL(NPLACE)=ICOL
        INEXT(NPLACE)=0
      END IF
      END subroutine

!      -------------------------------------------------------------

 
  subroutine increase_linked_list()
  ! increase linked list structure by 50%
  ! simple because linked list is unordered
  ! AL 6/7/07
  implicit none
  integer,allocatable:: ivcol_temp(:),inext_temp(:)
  real(r8), allocatable:: zhz_temp(:)
  integer:: newmaxrec
    
  ! create temp storage
  allocate(ivcol_temp(size(ivcol)),&
           inext_temp(size(inext)),&
           zhz_temp(size(zhz))   )

  ivcol_temp=ivcol
  inext_temp=inext
  zhz_temp=zhz
         
  newmaxrec=floor(1.5*maxrec) !new size
  print *,'changing maxrec from: ', maxrec, 'to: ',newmaxrec
  deallocate(ivcol,&
             inext,&
        zhz)
  allocate(ivcol(0:newmaxrec), &
           inext(0:newmaxrec),&
           zhz(newmaxrec))	 
  ivcol=0
  inext=0
  zhz=0d0
  ivcol(0:maxrec)=ivcol_temp
  inext(0:maxrec)=inext_temp
  zhz(1:maxrec)=zhz_temp
  deallocate(ivcol_temp,&
             inext_temp,&
             zhz_temp)
  maxrec=newmaxrec	     
  
  end subroutine

end module listaligada

MODULE Ecuyer_random
! L'Ecuyer's 1996 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! N.B. This version is compatible with Lahey's ELF90
! http://www.ozemail.com.au/~milleraj
! Latest revision - 30 March 1999

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307) !=real(8) con DVF

! These are unsigned integers in the C version
INTEGER, SAVE :: s1 = 4234, s2 = -4267, s3 = 7850

CONTAINS

SUBROUTINE init_seeds(i1, i2, i3)
IMPLICIT NONE

INTEGER, INTENT(IN) :: i1, i2, i3

s1 = i1
s2 = i2
s3 = i3
IF (IAND(s1,-2) == 0) s1 = i1 - 1023
IF (IAND(s2,-8) == 0) s2 = i2 - 1023
IF (IAND(s3,-16) == 0) s3 = i3 - 1023

RETURN
END SUBROUTINE init_seeds


FUNCTION taus88() RESULT(random_numb)
! Generates a random number between 0 and 1.  Translated from C function in:
! Reference:
! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
! generators', Math. of Comput., 65, 203-213.

! The cycle length is claimed to be about 2^(88) or about 3E+26.
! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

IMPLICIT NONE
REAL (dp) :: random_numb

INTEGER   :: b

! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
!      to the left if j > 0, otherwise to the right.

b  = ISHFT( IEOR( ISHFT(s1,13), s1), -19)
s1 = IEOR( ISHFT( IAND(s1,-2), 12), b)
b  = ISHFT( IEOR( ISHFT(s2,2), s2), -25)
s2 = IEOR( ISHFT( IAND(s2,-8), 4), b)
b  = ISHFT( IEOR( ISHFT(s3,3), s3), -11)
s3 = IEOR( ISHFT( IAND(s3,-16), 17), b)
random_numb = IEOR( IEOR(s1,s2), s3) * 2.3283064365E-10_dp + 0.5_dp

RETURN
END FUNCTION taus88

END MODULE Ecuyer_random

MODULE mt19937
! A Fortran-program for MT19937: Real number version
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-11-26  Time: 17:09:23
! Latest revision - 5 February 2002
! A new seed initialization routine has been added based upon the new
! C version dated 26 January 2002.
! This version assumes that integer overflows do NOT cause crashes.
! This version is compatible with Lahey's ELF90 compiler,
! and should be compatible with most full Fortran 90 or 95 compilers.
! Notice the strange way in which umask is specified for ELF90.
 
!   genrand() generates one pseudorandom real number (double) which is
! uniformly distributed on [0,1]-interval, for each call.
! sgenrand(seed) set initial values to the working area of 624 words.
! Before genrand(), sgenrand(seed) must be called once.  (seed is any 32-bit
! integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.

! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Library General Public License as published by
! the Free Software Foundation; either version 2 of the License, or (at your
! option) any later version.   This library is distributed in the hope that
! it will be useful, but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General Public License
! along with this library; if not, write to the Free Foundation, Inc.,
! 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.

!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.

!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed

! This program uses the following standard intrinsics.
!   ishft(i,n): If n > 0, shifts bits in i by n positions to left.
!               If n < 0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.

!***********************************************************************

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

! Period parameters
INTEGER, PARAMETER :: n = 624, n1 = n+1, m = 397, mata = -1727483681
!                                    constant vector a
INTEGER, PARAMETER :: umask = -2147483647 - 1
!                                    most significant w-r bits
INTEGER, PARAMETER :: lmask =  2147483647
!                                    least significant r bits
! Tempering parameters
INTEGER, PARAMETER :: tmaskb= -1658038656, tmaskc= -272236544

!                     the array for the state vector
INTEGER, SAVE      :: mt(0:n-1), mti = n1
!                     mti==N+1 means mt[N] is not initialized

PRIVATE
PUBLIC :: dp, sgrnd, grnd, init_genrand

CONTAINS


SUBROUTINE sgrnd(seed)
! This is the original version of the seeding routine.
! It was replaced in the Japanese version in C on 26 January 2002
! It is recommended that routine init_genrand is used instead.

INTEGER, INTENT(IN)   :: seed

!    setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
!    [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]

mt(0)= IAND(seed, -1)
DO  mti=1,n-1
  mt(mti) = IAND(69069 * mt(mti-1), -1)
END DO

RETURN
END SUBROUTINE sgrnd
!***********************************************************************

SUBROUTINE init_genrand(seed)
! This initialization is based upon the multiplier given on p.106 of the
! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.

! This version assumes that integer overflow does NOT cause a crash.

INTEGER, INTENT(IN)  :: seed

INTEGER  :: latest

mt(0) = seed
latest = seed
DO mti = 1, n-1
  latest = IEOR( latest, ISHFT( latest, -30 ) )
  latest = latest * 1812433253 + mti
  mt(mti) = latest
END DO



RETURN
END SUBROUTINE init_genrand
!***********************************************************************

FUNCTION grnd() RESULT(fn_val)
REAL (dp) :: fn_val

INTEGER, SAVE :: mag01(0:1) = (/ 0, mata /)
!                        mag01(x) = x * MATA for x=0,1
INTEGER       :: kk, y

! These statement functions have been replaced with separate functions
! tshftu(y) = ISHFT(y,-11)
! tshfts(y) = ISHFT(y,7)
! tshftt(y) = ISHFT(y,15)
! tshftl(y) = ISHFT(y,-18)

IF(mti >= n) THEN
!                       generate N words at one time
  IF(mti == n+1) THEN
!                            if sgrnd() has not been called,
    CALL sgrnd(4357)
!                              a default initial seed is used
  END IF
  
  DO  kk = 0, n-m-1
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk+m), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  DO  kk = n-m, n-2
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk+(m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  y = IOR(IAND(mt(n-1),umask), IAND(mt(0),lmask))
  mt(n-1) = IEOR(IEOR(mt(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
  mti = 0
END IF

y = mt(mti)
mti = mti + 1
y = IEOR(y, tshftu(y))
y = IEOR(y, IAND(tshfts(y),tmaskb))
y = IEOR(y, IAND(tshftt(y),tmaskc))
y = IEOR(y, tshftl(y))

! old code AL
!IF(y < 0) THEN
!  fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0) 
!ELSE
!  fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
!END IF

! to make it (0-1) AL
IF(y < 0) THEN
  fn_val = (DBLE(y) + 2.0D0**32) 
ELSE
  fn_val = DBLE(y) 
END IF
fn_val=(fn_val+0.5d0)/4294967296.d0


RETURN
END FUNCTION grnd


FUNCTION tshftu(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-11)
RETURN
END FUNCTION tshftu


FUNCTION tshfts(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,7)
RETURN
END FUNCTION tshfts


FUNCTION tshftt(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,15)
RETURN
END FUNCTION tshftt


FUNCTION tshftl(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-18)
RETURN
END FUNCTION tshftl

END MODULE mt19937



module matrix_algebra
use kinds
double precision:: tol=1d-20

contains

! ----------------------------------------------------------------------


!Sacada de los programas de misztal EV. He quitado rank y he puesto n y he puesto entrada.
!    En principio, esta te devuelve la triangular inferior. Lo hemos cambiado para que 
!    te de la triangular superior y el programa sea coherente.

   SUBROUTINE CHOLS4(entrada,x,indiq)
    ! cholesky decomposition
    integer :: n,i,j,k,rank_o,indiq
    integer::rank
    double precision::diagsq
    double precision::x(:,:), entrada(:,:)
    indiq=0
    rank_o=0
    x=entrada
    n=size(x,1)
    do i=1,n
       diagsq=x(i,i)-dot_product(x(i,1:i-1),x(i,1:i-1))
       if (abs(diagsq).lt.tol) then
         x(i,:)=0;x(:,i)=0       !zero row and column
       elseif (diagsq.lt.0) then
         print*,' Matrix not semipositive-definite, row ',i
         indiq=1
	 exit  
       else
         rank_o=rank_o+1
         x(i,i)=sqrt(diagsq)
         do j=i+1,n     
           x(j,i)=(x(j,i)-dot_product(x(j,1:i-1),x(i,1:i-1)))/x(i,i)
           x(i,j)=x(j,i)
         enddo
       end if
    enddo

      ! zero upper-diagonals
      do i=1,n
         x(i,i+1:n)=0
      enddo   
    
!    Cambiamos a triangular superior 

    x=transpose(x)
      
    END SUBROUTINE CHOLS4

! ----------------------------------------------------------
      subroutine inver2(var,n)
        double precision VAR(:,:),VER(size(var,1)),VIR(size(var,1),size(var,1)),A,B
!        CALL TIME(LSEG)
!       ----------------------------------------------------------
!       LA MATRIZ CARGADA DE DECIMALES SE UTILIZA PARA DETECTAR 
!       ERRORES DE REDONDEO
!       ----------------------------------------------------------
       if(n>size(var,1)) then
         print *,'n>size of var in inver2',n,size(var,1)
         stop
       endif
       

       DO 4 I=1,N
!       INVERSION DE LA ESQUINA
         DO J=1,I-1
           VER(J)=0.0
           DO  K=1,I-1
             VER(J)=VER(J)+VAR(I,K)*VAR(K,J)
           enddo
           VAR(I,I)=VAR(I,I)-VER(J)*VAR(I,J)
         enddo
         VAR(I,I)=1.0/VAR(I,I)

!        INVERSION DEL MARGEN
         DO  J=1,I-1
           VAR(I,J)=-VAR(I,I)*VER(J)
         enddo

!        INVERSION DE LA MATRIZ BLOQUE
         DO  J=1,I-1
           DO  K=1,I-1
             VAR(J,K)=VAR(J,K)-VER(J)*VAR(I,K)
           enddo
	 enddo

!        INVERSION DEL OTRO MARGEN
         DO  J=1,N-1
           VAR(J,I)=VAR(I,J)
         enddo
4      CONTINUE

       END subroutine
    

!Sacada de los programas de misztal EV
     SUBROUTINE INVCS4(l)
! calculates inverse of LL', where L is cholesky decomposition 
      integer  ::n,i,j,rank_o,indiq
      double precision :: l(:,:)
      double precision ::w(size(l,1))
    double precision ::kkita(size(l,1),size(l,1))
    
    n=size(l,1)
    call chols4(l,kkita,indiq)     ! cholesky factorization
!      Modificado por EV
    l=kkita

    do i=1,n
      w(i:n)=0
      if (abs(l(i,i)).gt.tol) w(i)=1/l(i,i)
  ! forward substitution
      do j=i+1,n
      if (abs(l(j,j)).gt.tol) &
      w(j)=-dot_product(l(j,i:j-1),w(i:j-1))/l(j,j)
      enddo
    !backward substitution
      do j=n,i,-1
      if (abs(l(j,j)).gt.tol) &
     w(j)=(w(j)-dot_product(l(j+1:n,j),w(j+1:n)))/l(j,j)
      enddo
      l(i:n,i)=w(i:n)
      l(i,i:n)=w(i:n)
      enddo   


    END SUBROUTINE INVCS4


    

! --------------------------------------
      subroutine ginv1(a,n,m,tol,irank)
! returns generalized inverse of matrix x of size n x n declared
! as m x m. tol is working zero and irank returns the rank of
! the matrix. w is a work vector of size m,
! by rohan fernando, slightly structured by i. misztal 05/05/87

      double precision:: a(m,m),w(m),re,sum,tol
      irank=n
      do 10 i=1,n
         do 20 j=1,i-1
              re=a(i,j)
              do 20 ii=i,n
20                 a(ii,i)=a(ii,i)-re*a(ii,j)
         if (a(i,i).lt.tol) then
              a(i,i)=0.0
              do 45 ii=i+1,n
45                 a(ii,i)=0.
           irank=irank-1
           else
              a(i,i)=sqrt(a(i,i))
              do 40 ii=i+1,n
40                a(ii,i)=a(ii,i)/a(i,i)
         endif
10    continue
 
      do 100 i=1,n
         if (a(i,i).eq.0.) then
              do 150 ii=i+1,n
150                a(ii,i)=0.
           else
              a(i,i)=1.0/ a(i,i)
              do 200 ii=i+1,n
200               w(ii)=0.0
              do 300 ii=i+1,n
                  iim1=ii-1
                  re=a(iim1,i)
                  do 400 iii=ii,n
400                   w(iii)=w(iii)-a(iii,iim1)*re
                  if (a(ii,ii).eq.0.) then
                      a(ii,i)=0.
                    else
                      a(ii,i)=w(ii)/a(ii,ii)
                  endif
300           continue
          endif
100     continue
 
      do 110 j=1,n
         do 110 i=j,n
              sum=0
              do 130 ii=i,n
130                sum=sum+a(ii,j)*a(ii,i)
110           a(i,j)=sum
      do 600 i=1,n
          do 600 j=i,n
600           a(i,j)=a(j,i)

      end subroutine


end module matrix_algebra


module random_generator
  use kinds
  use matrix_algebra
  use Ecuyer_random
  use mt19937
  
  contains


!     ---------------------
      subroutine unif(ix,u)
!      use kinds

!----------------------------------------------
! SUBRUTINA PARA MUESTREO DE UNIFORME- NUEVA 
! La semilla debe ser un entero comprendido entre:
! ix=2140000001      ! 0<ix<2147483647
!----------------------------------------------
      implicit none
      integer ix,k1
      double precision u
!      k1=ix/127773
!      ix=16807*(ix-k1*127773)-k1*2836
!      if(ix.lt.0)ix=ix+2147483647
!      u=ix*4.656612875E-10
      u=taus88()
!      u=grnd() !Mersenne-Twister
      end subroutine


      function fnormal() result(z)
!      generacion de un numero normal z -> n(0,1)
!      x1 es la semilla
      implicit double precision(a-h,o-z)
      double precision:: z,u1,u2
      call unif(0,u1)
      call unif(0,u2)
      z=((-2.*log(u1))**0.5)*cos(2.*3.1416*u2)
      end function

      subroutine normal(x1,z)
!      generacion de un numero normal z -> n(0,1)
!      x1 es la semilla
      implicit double precision(a-h,o-z)
    integer x1
      double precision:: z,u1,u2
      call unif(x1,u1)
      call unif(x1,u2)
      z=((-2.*log(u1))**0.5)*cos(2.*3.1416*u2)
      end subroutine

!     ================Nuevas==============

! -----------------------------------------------------
      subroutine gen_trunc_normal_old(cota1,cota2,z,x1)
!      generacion de un numero normal z -> normal truncada entre cota1 y cota2
!     (cotas estandarizadas a varianza=1)
!      x1 es la semilla
!     Metodo brutal a base de prueba y error 
!     Existe un metodo mas fino en Sorensen, p.146
!     AL 15/4/04
      implicit double precision(a-h,o-z)
      double precision:: z,xnor
      double precision:: cota1,cota2
      integer i,x1

      i=0
      do 
        i=i+1
!        if (i.gt.10000000) then
!          print *,i
!          stop
!        endif
!        z=xnor(x1)
        call normal(x1,z)
        if ((z.gt.cota1).and.(z.lt.cota2)) goto 101
      enddo
 101  continue
      end  subroutine

! ----------------------------------------------------------------------
      subroutine chin(se,ne,x,u)
!      se = varianza a priori.
!      ne = grados de credibilidad.
      implicit double precision (a-h,o-z)
      double precision ne
      integer x
1     call unif(x,u1)
      rg=se*14/(ne**.5)
      t=2*rg*u1+(se-rg)
      if (t.lt.0) goto 1
      o=(ne/(ne+2))*se
      chiot=-((ne+2)/2)*log(o)+(-(ne*se/2)/o)
      chitt=-((ne+2)/2)*log(t)+(-(ne*se/2)/t)
      ref=dexp(chitt-chiot)
      call unif(x,u1)
      if (u1.lt.ref) then
        u=t
      else
        goto 1
      endif

      end subroutine


      subroutine wish(nrank,se,ve,ne,x1)
       implicit double precision(a-h,o-z)
       double precision se(:,:),ve(:,:),ne
       double precision t(size(se,1),size(se,1)),a(size(se,1),size(se,1)),l(size(se,1),size(se,1))
       double precision b(size(se,1),size(se,1))
       integer x1,indiq
!      double precision se(nrank,nrank),ve(nrank,nrank)
!      double precision t(nrank,nrank),a(nrank,nrank),l(nrank,nrank)
!       double precision b(nrank,nrank)
         ia=ne/2

       do 1 i=1,nrank
           au=real(ne-i+1)/2.
           fio=gamdev(au,x1)
         t(i,i)=sqrt(2*fio)
         do 2 j=1,i-1
!           u=xnor(x1)
           call normal(x1,u)
            t(i,j)=0.
            t(j,i)=u
2         continue
1       continue
       do 3 i=1,nrank
          do 4 j=1,nrank
             a(i,j)=0.
             do 5 k=1,nrank
             a(i,j)=a(i,j)+t(i,k)*t(j,k)
5             continue
4          continue
3       continue
       do 6 i=1,nrank
         do 6 j=1,nrank
           l(i,j)=0.
6      continue
       call chols4(se,l,indiq)
       do 30 i=1,nrank
          do 40 j=1,nrank
             b(i,j)=0.
             do 50 k=1,nrank
               b(i,j)=b(i,j)+l(k,i)*a(k,j)
50           continue
40        continue
30     continue
       do 31 i=1,nrank
          do 41 j=1,nrank
             ve(i,j)=0.
             do 51 k=1,nrank
               ve(i,j)=ve(i,j)+b(i,k)*l(k,j)
51             continue
41        continue
31     continue
       do 101 i=1,nrank
         do 101 j=1,nrank
           ve(i,j)=ve(i,j)/ne
101    continue

       end subroutine


    


! ---------------------------------------------------      
      subroutine inv_con_wish(ncar,se,ve,ne,x1)
!      muestrea una wishart invertida para una suma de cuadrados (NO para su inversa)
!      condicionada a que un elemento 
!      (el del car\E1cter umbral) sea igual a 1.
!      Teor\EDa de Korsgaard, de 
!       http://www.csc.fi/ttn/ccb99/articles/IKorsgaard.pdf
!      Solo trabaja con el ultimo caracter
!      AL, 15/4/2004
!      10/6/04 -> corregido error en calculo de meant2
!      (no hacia multiplicacion matricial)

      double precision:: se(:,:),ve(:,:)
      integer x1
      integer i,j,ncar
!     auxiliares
      double precision:: se_old(size(se,1),size(se,1)),ne
      double precision:: V11(size(se,1),size(se,1)),V11_old(size(se,1),size(se,1))
      double precision:: t2(size(se,1)),meant2(size(se,1)),vart2(size(se,1),size(se,1))

! se como en el paper (sumas de cuadrados y no sumas de cuadrados entre el n de datos)
      do i=1,ncar
        do j=1,ncar
          se(i,j)=se(i,j)*real(ne)
        enddo
      enddo
! lo invierto para q sea como en el paper
      call inver2(se,ncar)
! y lo guardo
      do i=1,ncar
        do j=1,ncar
          se_old(i,j)=se(i,j)
        enddo
      enddo

! 1- V11
      call wish(ncar-1,se,V11,ne-1,x1)
      do i=1,ncar-1
        do j=1,ncar-1
! por la parametrizacion de la subrutina Wishart del paper donde E(W|se,ne)=se*ne
! wish da E(wish)=se
! por tanto lo reescalo
          V11(i,j)=V11(i,j)*real(ne-1)
          V11_old(i,j)=V11(i,j)
        enddo
      enddo
! inv(V11)
      call inver2(V11,ncar-1)

! 2- t2
! inv(se)
      call inver2(se,ncar)
      do i=1,ncar-1
        do j=1,ncar-1
! para sacar se_sup(22)
          vart2(i,j)=V11(i,j)*(1.d0/se(ncar,ncar))
        enddo
      enddo

! se de nuevo bien
      call inver2(se,ncar)
! para sacar inv(se_11)
      call inver2(se,ncar-1)
      do i=1,ncar-1
        meant2(i)=0.d0
        do j=1,ncar-1
          meant2(i)=meant2(i)+se(i,j)*se_old(j,ncar)
        enddo
      enddo

      call mvn(meant2,vart2,t2,ncar-1,x1)
      
      do i=1,ncar-1
        do j=1,ncar-1
          ve(i,j)=V11(i,j)+t2(i)*t2(j)
        enddo
        ve(i,ncar)=-t2(i)
        ve(ncar,i)=ve(i,ncar)
      enddo
      ve(ncar,ncar)=1.d0

      end subroutine
      
!     -------------------------------------------------------
      subroutine inv_con_wish_multiple(ncar,nres,se,ve,ne,x1)
!  muestrea una wishart invertida para una suma de cuadrados (NO para su inversa)
!  condicionada a que una submatriz 
! (la de los caracteres umbrales) sea igual a I.
! Teor\EDa de Korsgaard, de http://www.csc.fi/ttn/ccb99/articles/IKorsgaard.pdf
! AL, 10/6/04
! Ultima correccion 20/6/04
! 

! Problema2: puede producir matrices no positivas-definidas... no se si esto 
! choca con la teoria o no.
! Modified 6/7/07 to tale out "maxcar"
      ! be aware that most of this sizes anre not correct and make use of inver2(a,n) inverting n just in 
      ! its nxn upper left matrix

      double precision:: se(:,:),ve(:,:),ne 
! se -> sigma
      integer x1
      integer i,j,ncar,k1,k2
!     auxiliares
      double precision:: se_old(size(se,1),size(se,1)),se_221(size(se,1),size(se,1))
      double precision:: V11(size(se,1),size(se,1)),V11_old(size(se,1),size(se,1))
      double precision:: t2((ncar-nres)*nres),meant2((ncar-nres)*nres)
      double precision:: t2new(size(se,1),size(se,1))
      double precision:: vart2((ncar-nres)*nres,(ncar-nres)*nres)
      integer pos1,pos2
      integer ncont,nres
! ncont -> n caracteres continuos + umbrales con var. residual estimada
!          (p.ej. con varias categorias)
! nres numero de caracteres umbral con var res restringida a 1
      ncont=ncar-nres

! se como en el paper (sumas de cuadrados y no sumas de cuadrados 
! entre el n de datos)
      do i=1,ncar
        do j=1,ncar
          se(i,j)=se(i,j)*real(ne)
        enddo
      enddo
! lo invierto para q sea como en el paper
      call inver2(se,ncar)
! y lo guardo
      do i=1,ncar
        do j=1,ncar
          se_old(i,j)=se(i,j)
        enddo
      enddo

! 1- V11
      call wish(ncont,se,V11,ne-nres,x1)
      do i=1,ncont
        do j=1,ncont
! por la parametrizacion de la subrutina Wishart del paper donde E(W|se,ne)=se*ne
! mientras que la subrutina wish da E(wish)=se
! por tanto lo reescalo
          V11(i,j)=V11(i,j)*real(ne-nres)
          V11_old(i,j)=V11(i,j)
        enddo
      enddo
! inv(V11) -> v_11exp-1
      call inver2(V11,ncont)

! 2- t2
! inv(se)-> sigma-1
      call inver2(se,ncar)
! extraigo de sigma-1 el bloque correspondiente a se_22.1    
      do i=1,nres
        do j=1,nres
          se_221(i,j)=se(ncont+i,ncont+j)
        enddo
      enddo
! y lo invierto -> se_22.1
      call inver2(se_221,nres)

      do i=1,ncont
        do j=1,ncont
          do k1=1,nres
            do k2=1,nres
              pos1=(i-1)*nres+k1
              pos2=(j-1)*nres+k2
              vart2(pos1,pos2)=V11(i,j)*se_221(k1,k2)
            enddo
          enddo
        enddo
      enddo

! se de nuevo bien ->sigma
      call inver2(se,ncar)
! para sacar inv(se_11)
      call inver2(se,ncont)
      do i=1,ncont
        do k1=1,nres
          pos1=(i-1)*nres+k1
          meant2(pos1)=0.d0
          do j=1,ncont
            meant2(pos1)=meant2(pos1)+se(i,j)*se_old(j,ncont+k1)
          enddo
        enddo
      enddo

      call mvn(meant2,vart2,t2,nres*ncont,x1)
      
!      recoloco t2
      do i=1,ncont
        do k1=1,nres
          pos1=(i-1)*nres+k1
          t2new(i,k1)=t2(pos1)
        enddo
      enddo

! v11-1+t2t2'
      do i=1,ncont
        do j=1,ncont
          ve(i,j)=V11(i,j)
          do k1=1,nres
            do k2=1,nres
              ve(i,j)=ve(i,j)+t2new(i,k1)*t2new(j,k2)
            enddo
          enddo
        enddo
      enddo

! -t2
      do i=1,ncont
        do k1=1,nres
          ve(i,ncont+k1)=-t2new(i,k1)
          ve(ncont+k1,i)=ve(i,ncont+k1)
        enddo
      enddo

! Im4
      do i=1,nres
        do j=1,nres
          if (i.eq.j) then
            ve(ncont+i,ncont+j)=1.d0
          else 
            ve(ncont+i,ncont+j)=0.d0
          endif
        enddo
      enddo
      
! Chequeo si la inversa no tiene negativos
      do i=1,ncar
        do j=1,ncar
          se_old(i,j)=ve(i,j)
        enddo
      enddo
      call inver2(se_old,ncar)
      do i=1,ncar
          if (se_old(i,i).lt.0.d0) then
            print*,'negativo en inv(ve)',i,i,se_old(i,i)
            do j=1,ncar
              print *,(ve(j,k),k=1,ncar)
            enddo
          endif
      enddo

      end subroutine


! -------------------------------------------------
      subroutine inv_wish_special(rank,se,ve,ne,x1)
! Devuelve wishart invertidas para matrices con elementos operacionalmente 0
! Para no tener problemas numericos
! AL
      double precision:: se(:,:),ve(:,:),zero,ne
      integer rank,kk,x1
! Cero operacional
      zero=1e-15     

      call ginv1(se,rank,size(se,1),zero,kk)
      call wish(rank,se,ve,ne,x1)
      call ginv1(ve,rank,size(se,1),zero,kk)

      end subroutine


! ------------------------------
      subroutine mvn(a,v,b,n,x1)
! Genera normales multivariadas b~MVN(a,v) de rango n con semilla x1
! A partir de una subrutina de Luis Varona
      integer i,n,j,x1,indiq
      double precision a(:),v(:,:),cv(size(v,1),size(v,1)),b(n)
      double precision u(n),xnor
      call chols4(v,cv,indiq) 
      do i=1,n
        call normal(x1,u(i))
!        u(i)=xnor(x1)
      enddo
      do i=1,n
        b(i)=a(i)
        do j=1,i
          b(i)=b(i)+u(j)*cv(j,i)
        enddo
      enddo
      end subroutine


!---- Subrutinas para obtener normales truncadas; principalmente de programas
!       de Misztal

! ---------------------------------
      function normal_invcdf(p) 
      ! return inverse of CDF, i.e., such x: p=cdf(x)
      !   resticted to |invcdf(x)|<10

!      use kinds
      double precision:: p,x,eps,normal_invcdf
      integer i

      if (p.lt.0. .or. p .gt. 1.) then
         print*,'normal_invcdf: arguments outside of 0..1'
         stop
      endif

      eps=10

      x=0
      do i=1,50
         if (normalcdf(x) .lt. p) then     
           x=x+eps
         else
           x=x-eps
         endif    
         if (abs(x) .gt. 10.1) exit
         eps=eps/2.
      enddo      
      normal_invcdf=x
      end function

! ---------------------------
      function normalcdf(x)
      ! returns cumulative normal density function
      implicit none
      double precision:: normalcdf,x,alnorm,q,pdf
      !
      !normalcdf=alnorm(x,.false.)
       call  NORMP(x, normalcdf, Q, PDF)

      end function


!-----------------------------------
      SUBROUTINE NORMP(Z, P, Q, PDF)
!
!    Normal distribution probabilities accurate to 1.e-15.
!    Z = no. of standard deviations from the mean.
!    P, Q = probabilities to the left & right of Z.   P + Q = 1.
!       PDF = the probability density.
!
!       Based upon algorithm 5666 for the error function, from:
!       Hart, J.F. et al, 'Computer Approximations', Wiley 1968
!
!       Programmer: Alan Miller
!
!    Latest revision - 30 March 1986
!
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA P0, P1, P2, P3, P4, P5, P6/220.2068679123761D0, &
           221.2135961699311D0, 112.0792914978709D0, &
           33.91286607838300D0, 6.373962203531650D0, &
          .7003830644436881D0, .3526249659989109D-01/, &
           Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7/440.4137358247522D0, &
           793.8265125199484D0, 637.3336333788311D0, &
           296.5642487796737D0, 86.78073220294608D0, &
           16.06417757920695D0, 1.755667163182642D0, &
           .8838834764831844D-1/, &
           CUTOFF/7.071D0/, ROOT2PI/2.506628274631001D0/
!
      ZABS = ABS(Z)
!
!      |Z| > 37.
!
      IF (ZABS .GT. 37.D0) THEN
        PDF = 0.D0
        IF (Z .GT. 0.D0) THEN
          P = 1.D0
          Q = 0.D0
        ELSE
          P = 0.D0
          Q = 1.D0
        END IF
        RETURN
      END IF
!
!      |Z| <= 37.
!
      EXPNTL = EXP(-0.5D0*ZABS**2)
      PDF = EXPNTL/ROOT2PI
!
!      |Z| < CUTOFF = 10/sqrt(2).
!
      IF (ZABS .LT. CUTOFF) THEN
        P = EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS + &
             P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS + &
             Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS + &
             Q0)
!
!      |Z| >= CUTOFF.
!
      ELSE
        P = PDF/(ZABS + 1.D0/(ZABS + 2.D0/(ZABS + 3.D0/(ZABS + 4.D0/ &
         (ZABS + 0.65D0)))))
      END IF
!
      IF (Z .LT. 0.D0) THEN
        Q = 1.D0 - P
      ELSE
        Q = P
        P = 1.D0 - Q
      END IF
      END subroutine

! ---------------------------
      function dnormal(x)
      !
      ! Returns probability density of normal density function
      !
      implicit none
      double precision:: dnormal,x
      double precision:: pi
      
      pi=3.14159265358979
      dnormal=1.d0/sqrt(2.d0*pi)*dexp(-x*x/2.d0)
      end function

      subroutine gen_trunc_normal(a,b,x,x1)
! Generates normal distribution N(0,1) truncated to <a,b>

      double precision:: a,b,cdfa,cdfb
      double precision:: un,x
      integer x1


      if (b .lt. a) then
         Print*,'GEN_TRUNC_NORMAL: bound b =',b,' < boound a =',a
         stop
      endif   

      cdfa=normalcdf(a)
      cdfb=normalcdf(b)
1     continue        
      call unif(x1,un)  !uniform random number generator UN(0,1)
      x=normal_invcdf(cdfa+(cdfb-cdfa)*un)       
      if ((x.lt.a).or.(x.gt.b)) then
        print *,a,x,b
	print *,'I can''t sample this truncated normal, subroutine gen_trunc_normal'
!        goto 1
      endif
      end subroutine
      
! ----------------------------------
       FUNCTION XNORMD(PROB)
!       SOLO DEFINIDA PARA PROB EN EL ESPACIO 0-.5
       IMPLICIT double precision (A-H,O-Z)
       double precision,save:: limit=tiny(1d0)
       if(prob<limit) prob=limit !AL to avoid log(0)
       T=SQRT(LOG(1./(PROB*PROB)))
       XNORMD=T-(2.515517+T*(.802853+T*.010328))/(1.+T* &
       (1.432788+T*(.189269+T*.001308)))
       END function

! -------------------------------------

      function xnor_i(x1) result(u)
      integer :: x1
      double precision:: u
      
      u=fnormal()
      end function
      
      
      FUNCTION XNOR(x1)

       IMPLICIT double precision (A-H,O-Z)
       integer x1
       double precision:: xnor

!       CALCULA LA ABCISA (X) PARA CUALQUIER INCIDENCIA (P) 
! Genera numeros de distribucion normal mediante transformacion a partir
! de una uniforme (x1 es el valor de esa uniforme)
! Segun Luis no es tan fina como (normal) pero es mejor en el modelo umbral
! para muestrear normales truncadas "a pelo"
! Modificada AL para usar unif() ya que la otra subrutina dio problemas
! de ciclicidad  

! Generacion de muestra de la uniforme

!       divis=2.**31.-1.
!       trans=7.**5.
!       divid=trans*x1
!       lsol=int(divid/divis)
!       x1=divid-lsol*divis
!       prob=x1/divis

       call unif(x1,prob)

!      Transformacion
       IF(PROB.LT..5) THEN
          XNOR=-XNORMD(PROB)
       ELSE IF(PROB.GT..5)THEN
          XNOR=XNORMD(1.-PROB)
       ELSE
          XNOR=0.
       ENDIF
       END function

!     ----------------------------------------------------------------
      function gamdev(a,x1)
       implicit double precision(a-h,o-z)
       integer x1
       if(a.lt.1)then
         print *,'illegal alpha parameter in gamdev'
       endif
       if(a.lt.6)then
         x=1.
           ia=int(a)
         do 11 j=1,ia
            call unif(x1,u)
            x=x*u
11        continue
         x=-log(x)
       else
1         call unif(x1,u)
         v1=2.*u-1.
         call unif(x1,u)
         v2=2.*u-1.
         if(v1**2+v2**2.gt.1.)goto 1
         y=v2/v1
           am=a-1
         s=sqrt(2.*am+1.)
         x=s*y+am
         if (x.le.0.)goto 1
            e=(1.+y**2)*exp(am*log(x/am)-s*y)
         call unif(x1,u)
         if (u.gt.e)goto 1
       endif
       gamdev=x
       end function


end module random_generator


module aux1

contains

      subroutine printtime
      INTEGER  DATE_TIME (8)
      CHARACTER(LEN = 12) REAL_CLOCK (3)

      CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                      REAL_CLOCK (3), DATE_TIME)

      print *,real_clock(1)(7:8),'/',real_clock(1)(5:6),&
      '/',real_clock(1)(1:4),' ', &
      real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)

      end subroutine

      function sqrt1(a)
!     returns square root of a if a>0, otherwise 0
      double precision:: a,sqrt1
      if(a.le.0.d0) sqrt1=0.d0
      if(a.gt.0.d0) sqrt1=sqrt(a)
      end function

end module aux1

module aux_options
! from IM's blupps1.f90
implicit none
character(len=20):: xc(3)
integer:: nitem
real:: xcnum(3)

contains

  subroutine nums2(a,n,x,xc)
  ! separates array a into items delimited by blanks. character elements are
  ! put into optional character vector xc, decoded numeric values
  ! into optional real vector x, and n contains the number of items. The
  ! dimension of x and xc can be lower than n.
  ! A modification of nums() from f77 to f90
  ! Now accepts real numbers
  ! 2/23/2000

  character (*)::a
  character (*),optional::xc(:)
  real,optional::x(:)
  integer::n,curr,first,last,lena,stat,i

  curr=1;lena=len(a);n=0

  do
    ! search for first nonspace
    first=0
    do i=curr,lena
       if (a(i:i) /= ' ') then
          first=i
          exit
       endif
    enddo
    if (first == 0) exit


    ! search for first space
    curr=first+1
    last=0
    do i=curr,lena
        if (a(i:i) == ' ') then
          last=i
          exit
        endif
    enddo

    if (last == 0) last=lena

    n=n+1
    if (present(xc)) then
       if (size(xc) >= n) then
             xc(n)=a(first:last)
           else
             stop "size of xc in nums2 too small"
       endif
    endif
    if (present(x)) then
       if (size(x) >= n) then
              read(a(first:last),'(f12.0)',iostat=stat)x(n)
              if (stat /=0) x(n)=0
           else
              stop "size of x in nums2 too small"
       endif
    endif

    curr=last+1
  enddo
  end subroutine
  subroutine getoption_unit(io_unit,name,n,x,xc)
  ! In unit io_unit, locates line:
  ! OPTION name str1 str2
  ! where str1, str2 are strings separated by spaces.
  ! Then, it assigns: xc(1)=str1,  xc(2)=str2,...
  ! and attempts to decode strings into real values: x1=value(str1),....
  !
  ! n contains the number of strings.  x and xc are optional and their
  ! dimensions may be smaller than n in which case some strings/values are
  ! not stored.
  !
  ! Upon exit, unit io_p=40 points to line next to the one located.
  !
  ! If the line cannot be located, n=-1
  !

  character (*)::name
  integer::n,io_unit
  real,optional::x(:)
  character (*),optional::xc(:)

  real::x1(100)
  integer::stat,m
  character (50)::xc1(100)
  character (200)::a

  rewind io_unit
  n=-1

  do
     read (io_unit,'(a)',iostat=stat)a
     if (stat /= 0) exit
        call nums2(a,m,x1,xc1)
        if ( m>1 .and. xc1(1) == 'OPTION' .and. xc1(2) == name) then
           n=m-2
           if (present(xc)) xc=xc1(3:size(xc)+2)
           if (present(x)) x=x1(3:size(x)+2)
           exit
        endif
  enddo
  end subroutine
  
end module




!        ----------
         Program GIBBSTHUR     
!        ----------

! use denseop temporary disabled AL
use kinds
use modelo
use listaligada
use matrix_algebra
use random_generator
use aux1
use aux_options
use vers

implicit none

!
real(r8):: uu, ee, landa
!
real(r8),allocatable::  b(:,:)
real(r8):: u,var
real(r8),allocatable:: vara(:,:),vare(:,:),varp(:,:,:),vark(:)
real(r8),allocatable:: vara_a(:,:),vare_a(:,:),varp_a(:,:,:)
real(r8),allocatable:: se(:,:),ve(:,:),va(:,:),sa(:,:),ve2(:,:),ve3(:,:)
real(r8):: ne
real(r8),allocatable:: y(:,:)
real(r8),allocatable:: reg(:,:)
real(r8),allocatable:: thr(:,:)
real(r8):: cota1,cota2,xcor
real(r8):: total,coef
real(r8),allocatable:: mint(:,:),maxt(:,:)
real(r8):: xm
real(r8):: xmin,xmi,xmax,xd,a,b1,xa,xm1
real(r8),allocatable:: xmm(:)
!ev_mad
real(r8),allocatable::valor(:,:)
!
integer:: npdve,npdvp,npdva ! numero de veces q Ve Vp y Va son no semipositivo definido
logical::  stat
integer::  nfac
integer:: ncarthur
integer:: iold,ifen,icar,ic
!ev_mad
integer:: ncov
!
integer::  ngrup
integer::  nefani
integer::  i,j,k,ijk,kji,jk,ji,ij
logical::  kk
integer::  nrow
integer::  mm
integer::  nanim
integer::  nangru
integer::  mami
integer::  ia,ip,im
integer:: imue,lap,icad
integer::  is,iss,imgs
integer::  ndat1
integer::  iplace
integer::  ilev,rank_o
integer::  irango
integer::  l
integer :: na,n
integer::  itercor, itermas 
integer::  ncar,indiq
integer::  x1,ifentest,ifentest2
!     AL 
integer::  presenthr,nper,iper
integer::  ncarumb,ncont,poscar,nres,maxthr
integer::  pos1,pos2,iefani,iefani2
integer::  posani
integer::  ia1,ip1,im1,IS1,ISS1,IMGS1

logical,allocatable :: kk2(:,:)
integer, allocatable :: ipru(:,:,:),iloc(:),indth(:,:),iaux(:,:),indaux(:)
integer,allocatable :: ifac(:),nfix(:,:),mfac(:),ind(:,:),nthr(:), &
                       posvara1(:,:),posvara2(:,:)
real(r8),allocatable :: nzz(:),dia(:),yy(:,:)

real(r8),allocatable:: sumvara(:,:),ssvara(:,:),sumvarp(:,:,:),ssvarp(:,:,:), &
                     sumvare(:,:),ssvare(:,:)
real(r8),allocatable:: cora(:,:),corp(:,:,:),core(:,:),h2(:)
real(r8),allocatable:: sumcora(:,:),sumcore(:,:),sumcorp(:,:,:),sscora(:,:), &
                     sscore(:,:),sscorp(:,:,:)
real(r8),allocatable:: sumb(:,:),ssb(:,:)
real(r8):: x,x2,sum,diagsq
real:: seeds(3)=0d0
character(len=40) fichero,datos,genea
character(len=100) rotulo(100)
character(len=40) tipomodelo,recursividad
character(10)::kkk,task
logical :: VCE=.true.
!EV
real(r8),allocatable:: d(:),v(:,:) !Eigenvalues y eigenvectors
real(r8) :: liabilitybound=4.d0 !other option is 4.d0 to constrain it

! to get contrasts
integer:: efftest,sizecontrast,io
integer,allocatable:: levefftest(:)
character(5):: contrast


call print_version(start=.true.)


npdve=0.d0; npdvp=0.d0; npdva=0.d0
itercor=0; itermas=0;
total=0.d0
u=0.d0
presenthr=0
iper=0
poscar=0
pos1=0
pos2=0
iefani=0
iefani2=0
posani=0
ia1=0
ip1=0
im1=0
is1=0
iss1=0
imgs1=0
nfac=0;ngrup=0;nefani=0;l=0;i=0;j=0;ijk=0;kji=0;jk=0;ji=0;mm=0;nanim=0;nangru=0
ia=0;ip=0;im=0
is=0;iss=0;imgs=0
iplace=0
irango=0
ilev=0
na=0
! Pone a 0 
NPLACE=0
!EV



! --> Ficheros de salida
open(15,file='samples.txt',position='append')


!-----------------------------------------
! --> Parameter file reading <--
!-----------------------------------------

   
OPEN(11,FILE='parameters.txt',status='old')
read(11,*)datos
ndat1=0
open(13,file=datos,status='old')
do 
  read(13,*,iostat=io)
  if(io/=0) exit
  ndat1=ndat1+1
enddo  
print *,'number of records',ndat1
close(13)

read(11,*)genea
tipomodelo='animal'
! Numero de efectos, numero de grupos geneticos (al menos uno), numero de caracteres, 
! Se asume que el ultimo caracter es el umbral
!ev-mad: number of efects (nfac) =number of covariates (ncov) + number of crossclassified effects (ncrossclassified)


READ(11,*)NFAC
READ(11,*)NCOV
READ(11,*)NGRUP
READ(11,*)NCAR
READ(11,*)ncarumb

IF(NCARumb.GT.nCAR)THEN
  PRINT *,'more threshold traits than total traits',NCARumb,ncar
  STOP
END IF
ncont=ncar-ncarumb
nres=0

allocate( nthr(ncarumb) )
nthr=0

!  numero de categorias (orden:primero 3 o mas, segundo 2)
maxthr=0
READ(11,*)(nthr(i),i=1,ncarumb)
do i=1,ncarumb
  if (nthr(i).eq.2) nres=nres+1
enddo 
  maxthr=maxval( nthr )

allocate( thr(ncar,maxthr),mint(ncar,maxthr) )
allocate( maxt(ncar,maxthr))


! numero de caracteres thrustonianos
read(11,*)ncarthur
nres=nres+ncarthur
print *,'number of traits with var(e) constrained to 1 -->',nres
IF(ncarthur.gt.ncar)THEN
  PRINT *,'more trhustonian traits than total traits',NCARumb,ncar
  STOP
END IF
IF((ncarthur+ncarumb).GT.nCAR)THEN
  PRINT *,'more trhustonian and threshold traits than total traits',ncarthur,NCARumb,ncar
  STOP
END IF
ncont=ncont-ncarthur
allocate(iloc(ncarthur))
read(11,*)(iloc(i),i=1,ncarthur)

! numero de efectos permanentes
read(11,*) nper

! numero de efectos animales (correlacionados)
read(11,*) nefani
allocate( ifac(nfac+1)); ifac=0

! Numero de niveles de cada efecto (en los animales NO incluye grupo genetico)
READ(11,*)(IFAC(I),I=1,NFAC)
! Animales son los \FAltimos efectos 
! Permanente ("fijo aleatorio") los penultimos
NANIM=IFAC(NFAC)
do i=1,nefani
  ifac(nfac-nefani+i)=ifac(nfac-nefani+i)+ngrup
enddo

! Prepare address array (ifac)
NROW=0

! Ahora ifac(i) tiene la posicion del efecto (i-1)-esimo en su ultimo nivel
DO  I=1,NFAC


!ev-mad Modified to account for the levels of each covariate (1)


  if (i .le. ncov) then
    MM=1+nrow
    ifac(i)=nrow
    nrow=mm
  else
!
    MM=NROW+IFAC(I)
    IFAC(I)=NROW
    NROW=MM
  endif
enddo

! Numero total de ecuaciones
! USAR NROW PARA ALLOCATAR
IFAC(NFAC+1)=NROW

NANGRU=NANIM+NGRUP

PRINT *,'Animals, unknown parent groups =',NANIM,NGRUP
MAMI=0


allocate( b(nrow,ncar), &
     nzz(nrow),nfix(nfac,ndat1),y(ndat1,nCAR),MFAC(nfac), &
     vara(ncar*nefani,ncar*nefani), &
     varp(ncar,ncar,nper), &
     vare(nCAR,nCAR),vark(ncar) &
     ,se(nCAR,nCAR),ve(nCAR,nCAR)&
     ,sa(ncar*nefani,ncar*nefani),va(ncar*nefani,ncar*nefani)&
     ,reg(ndat1,ncar),IND(nCAR,nfac),KK2(nrow,nCAR) )
allocate(vara_a(ncar*nefani,ncar*nefani), &
     varp_a(ncar,ncar,nper), &
     vare_a(nCAR,nCAR))
allocate(ipru(ndat1,2,ncarthur))
allocate(iaux(ndat1,2))
allocate(indth(ndat1,ncarthur))
allocate(indaux(ndat1))
! Outputs
allocate( sumvara(ncar*nefani,ncar*nefani), &
     ssvara(ncar*nefani,ncar*nefani) )
allocate( sumvarp(ncar,ncar,nper), &
       ssvarp(ncar,ncar,nper) )
allocate( sumvare(ncar,ncar),ssvare(ncar,ncar) )

allocate( cora(ncar*nefani,ncar*nefani), &
     corp(ncar,ncar,nper) )
allocate( core(ncar,ncar) )
allocate( sumcora(ncar*nefani,ncar*nefani),sscora(ncar*nefani,ncar*nefani) )
allocate( sumcorp(ncar,ncar,nper), &
        sscorp(ncar,ncar,nper) )
allocate( sumcore(ncar,ncar),sscore(ncar,ncar) )
allocate( sumb(nrow,ncar),ssb(nrow,ncar) )
allocate(h2(ncar*nefani))
allocate(posvara1(nefani,ncar),posvara2(nefani,ncar))

! allocate linked list
allocate(zhz(maxrec),ifirst(nrow),ivcol(0:maxrec), &
     inext(0:maxrec) )

allocate( dia(nrow),yy(nrow,ncar) )
allocate(xmm(ndat1))
!EV
allocate(d(ncar), v(ncar,ncar))
!ev_mad para la cov
allocate(valor(ndat1,nfac))
valor=0d0


!     Set to zero (f90 style)
sumvara=0.d0
ssvara=0.d0
sumvarp=0.d0
ssvarp=0.d0
sumvare=0.d0
ssvare=0.d0
sumcora=0.d0
sscora=0.d0
sumcorp=0.d0
sscorp=0.d0
sumcore=0.d0
sscore=0.d0
sumb=0.d0
ssb=0.d0
d=0.d0;v=0.d0
var=0.d0
se=0.d0
ve=0.d0
sa=0.d0
va=0.d0
cora=0.d0
corp=0.d0
core=0.d0
h2=0.d0
b=0.d0
posvara1=0
posvara2=0
DIA=0.E+00
YY=0.E+00
reg=0.d0
y=0.d0
nfix=0
nzz=0.d0
IFIRST=0


!EV

! different models for each trait

DO I=1,NCAR
  READ(11,*)(IND(I,J),J=1,NFAC)
ENDDO

! Reading features of the Gibbs sampling
read(11,*) imue
read(11,*) lap
read(11,*) icad
imue=imue/icad
lap=lap/icad
print *,'Total n of iterations:',imue*icad
print *,'Burn-in: discarded for results.txt and solutions.txt',lap*icad
print *,'Thin interval',icad
print *,'Total n of samples in the Gibbs sampler:',imue
! semillas de varianzas

! lectura de prior variances
! Reading initial variances
! Genetic
vara=0
do i=1,ncar*nefani
   vara(i,i)=1.
enddo
varp=0
do iper=1,nper
  do i=1,ncar
   varp(i,i,iper)=1.
  enddo
enddo
! Residual
vare=0
do i=1,ncar
  vare(i,i)=1.
enddo

call getoption_unit(11,'RandomSeeds',pos1,x=seeds)
if(pos1/=-1) then
	if(pos1==3)then
		print *,'putting random seeds',int(seeds)
		call init_seeds(int(seeds(1)),int(seeds(2)),int(seeds(3)))
	endif
endif


close(11)


!--------------
!-> Preliminares
!---------------


! funcion para la posicion de las varianzas
do iefani=1,nefani
  do j=1,ncar
    posvara1(iefani,j)=(iefani-1)*ncar+j
    posvara2(iefani,j)=(iefani-1)*ncar+j
  enddo
enddo

!--------------------------------
! --> Loop Lectura genealogia <--
!--------------------------------

! El fichero debe estar ordenado por animal para chequeos

!     SE LEE FICHERO GEN Y SE CONSTRUYE LA INVERSA DE A
!       SIN CONSANGUINIDAD Y GPD


  print *,'animal model with genetic groups '
  OPEN(12,FILE=GENEA,status='old')
  3 READ(12,*,END=4)IA,IP,IM

! VERIFICACION DE LA ESTRUCTURA DE LOS DATOS
  MAMI=MAMI+1
  if(mod(mami,10000)==0) print *,'pedigree file, record',mami


  IF((IA.GT.NANIM).OR.(IP.GT.NANGRU).OR.(IM.GT.NANGRU).OR. &
    (IA.NE.MAMI).OR.(IP.LT.1).OR.(IM.LT.1))THEN
    PRINT *,'non-sorted or incorrect genealogy file?',IA,IP,IM
    STOP
  ENDIF

! Almacenamos A-1 repetidamente
  do i=1,nefani
    ! Si hay genealogia desconocida 
    XM=0.d0
    IF(IP.GT.NANIM) XM=XM+1.
    IF(IM.GT.NANIM) XM=XM+1.
    X=4./(XM+2.)
    ! Posiciones
    IA1=IA+IFAC(NFAC-nefani+i)
    IP1=IP+IFAC(NFAC-nefani+i)
    IM1=IM+IFAC(NFAC-nefani+i)

    ! Almacenamos el valor del animal en la diagonal (es la diagonal de A-1)
    DIA(IA1)=DIA(IA1)+X

    ! Animal - Padre y madre
    X2=X/2.
      CALL LINKS(IP1,IA1,-X2)
      CALL LINKS(IA1,IP1,-X2)
      CALL LINKS(IM1,IA1,-X2)
      CALL LINKS(IA1,IM1,-X2)

    ! Padre-madre
    X2=X/4.
    DIA(IP1)=DIA(IP1)+X2
    DIA(IM1)=DIA(IM1)+X2
    IF(IM1/=IP1)THEN
      CALL LINKS(IM1,IP1,X2)
      CALL LINKS(IP1,IM1,X2)
    ELSE
      DIA(IM1)=DIA(IM1)+2*X2
    ENDIF
  enddo
  GOTO 3
4 CLOSE(12)
!----------------------------------
! -- Fin Loop Lectura genealogia --
!----------------------------------

!-------------------------------------
! --> Lectura del fichero de datos <--
!-------------------------------------

!     LECTURA DEL FICHERO DE DATOS
OPEN(13,FILE=datos,status='old')
NDAT1=0

	iold=0
	ifen=0
5 NDAT1=NDAT1+1
  if(mod(ndat1,10000)==0) print *,'Data file, record',ndat1


  !  Lectura de datos: covariables(mod_ev_mad), efectos fijos, fijo aleatorio,animal, caracteres
  !  leo todo en mfac (a diferencia de antes q eran dos vectores, ia y mfac)
  !      READ(13,*,END=6)(MFAC(I),I=nfac-nefani+1,NFAC),
  !       animales
  !     +(MFAC(I),I=1,NFAC-nefani)
  !       fijos y permanentes
  !     1,(REG(NDAT1,I),i=1,ncar)
!ev_mad the vector valor takes the value of the covariates in their positions
!mfac vector takes the value of each crossclassified effect


  READ(13,*,END=6) (valor(ndat1,i),i=1,ncov),(MFAC(I),I=ncov+1,NFAC), & !efectos
            (REG(NDAT1,I),i=1,ncar)  !registros

!	if (ncarthur.eq.1) then
!	if (mfac(nfac-nefani-nper).lt.iold) then
!	print *,'non ordered fixed effect for thrustonian trait'
!	stop
!	else
!	if (mfac(nfac-nefani-nper).eq.iold) then
!	    ifen=ifen+1
!	    ifentest=int(reg(ndat1,ncar))
!	    if ((ifentest.ne.ifen).and.(ifentest.ne.99)) then
!	       print *,'non ordered data for thrustonian trait 1'
!	       stop
!	    endif	
!	else
!	    iold=mfac(nfac-nefani-nper)
!	    ifen=1
!	    ifentest=int(reg(ndat1,ncar))
!	     if ((ifentest.ne.ifen).and.(ifentest.ne.99)) then
!	      print *,'non ordered data for thrustonian trait 2'
!	      stop
!	    endif
!	endif
!	endif
!	endif	

	do i=1,ncarthur
	      ipru(ndat1,1,i)=mfac(iloc(i))
	      ipru(ndat1,2,i)=int(reg(ndat1,ncont+ncarumb+i))
!	      print *,ipru(ndat1,1,i),ipru(ndat1,2,i),ndat1
	enddo
  do i=nfac-nefani+1,nfac
    IF((mfac(i).LT.1).OR.(mfac(i).GT.NANIM))THEN
     PRINT *,'animal=0 or higher number than expected',mfac(i),NANIM
     STOP                       
    END IF
  enddo
!ev_mad: I modified it in order to skip the covariates
  DO  I=1+NCOV,NFAC-nefani
    IF((MFAC(I).LT.1).OR.(MFAC(I).GT.IFAC(I+1)-IFAC(I)))THEN
      PRINT *,'level of effect',I,' beyond bounds in record',  &
       ndat1,' : ',IFAC(I+1)-IFAC(I),' --->',MFAC(I)
      STOP
    END IF
  enddo
!ev_mad
!valor vector takes the value of the covariate in their corresponding positions and 1 in the corresponding positions of
!the crossclassified effects
!mfac vector takes the value of 1 in the positions corresponding to the covariates and the level of each crossclassified
!effect in ther corresponding positions for the crossclassified effects
!Example: data 1: 278 250 5 (cov1 cov2 ef1), valor will contain (270 250 1), and mfac (1 1 5)


        do j=1,nfac
                if (j .le. ncov) then
                mfac(j)=1
                else
                valor(ndat1,j)=1.d0
                endif
        enddo


  ! Nfix almacena toda la estructura de datos para todos los registros. O sea [X;Z]
  ! En realidad es [Z;X] modificarlo a [Z1:Zn;X]
  ! loop a lo largo de todos los ef animales
  do i=1,nefani
    NFIX(i,NDAT1)=mfac(nfac-nefani+i)
  enddo
  DO  I=1,NFAC-nefani
    NFIX(I+nefani,NDAT1)=MFAC(I)
  enddo
  
  ! SE CONSTRUYE RHS: XX, ZX 

  do i=1,nfac
    pos1=ifac(i)+mfac(i)
    ! ef fijo/permanente -> X'X
!ev_mad *valor(ndat1,i) was added. for covariates, we multiply by the value of the covariate**2,that is, valor*valor)
    if (i.le.(nfac-nefani)) dia(pos1)=dia(pos1)+valor(ndat1,i)*valor(ndat1,i)
        ! ef aditivo -> Z'Z
    if (i.gt.(nfac-nefani)) nzz(pos1)=nzz(pos1)+valor(ndat1,i)*valor(ndat1,i)
    do j=1,nfac
      pos2=ifac(j)+mfac(j)
      ! -> X'Z, Z1'Z2,etc
          !ev_mad valor(ndat1,i)*valor(ndat1,j) was added. If pos1 corresponds to a covariate and pos2 corresponds to another effect 
          !(not covariate), the value of the covariate*1 is added. If pos1 and pos2 are not covariates, 1*1 is added
          if (pos1/=pos2) CALL LINKS(pos1,pos2,(valor(ndat1,i)*valor(ndat1,j)))
        enddo
  enddo

GOTO 5
6 CLOSE(13)
NDAT1=NDAT1-1

!---------------------------------------
! -- Fin Lectura del fichero de datos --
!---------------------------------------
!	pause

       PRINT*,'number of records, non-null elements = ',NDAT1,NPLACE*2

! preparaciÃ³n de los ficheros para el modelo thurstoniano.......

	do i=1,ncarthur
	  do j=1,ndat1
	    iaux(j,1)=ipru(j,1,i)
	    iaux(j,2)=ipru(j,2,i)
	    indaux(j)=j
!	    print *,iaux(j,1),iaux(j,2),indaux(j)
	  enddo
	  call Bubble_Sort(iaux,indaux,ndat1)
	  do j=1,ndat1
!	     print *,iaux(j,1),iaux(j,2),indaux(j)
	     ipru(j,1,i)=iaux(j,1)
	     ipru(j,2,i)=iaux(j,2)
	     indth(j,i)=indaux(j)
	  enddo
!	  pause
	enddo
	
!===========================================================================
!
!   B.-GIBBS


! \C1 couple of preliminaries...

! kk2 es un indicador de si el muestreo en 
!  la posicion (efecto, caracter) es v\E1lida (.true.) o no (.false.)

      KK2=.false.

      DO I=1,NFAC
        DO J=1,NCAR
          IF (IND(J,I).EQ.1) THEN
            DO K=IFAC(I)+1,IFAC(I+1)
               KK2(k,j)=.true.
            ENDDO
          ENDIF
        ENDDO
      ENDDO

       nanim=nanim+ngrup
!      Posicion de partida de los animales en las MME
      ilev=ifac(nfac-nefani+1)


! Varianzas iniciales


!- Inicializo los umbrales
! Primer umbral es 0, segundo (si >=2 categorias) 1, y ultimo cuasi-infinito 
! loop para todos los car. umbrales
      do j=1,ncarumb
        do i=1,nthr(j)-1
          thr(j,i)=1.d0*(i-1)
        enddo
        thr(j,nthr(j))=999.d0
	!AL 6/04/06 to avoid jumpins in binary traits
	if (nthr(j)==2) thr(j,nthr(j))=liabilitybound !for binary traits only!!!

        mint(j,1)=-999.d0
        maxt(j,1)=0.d0
        if(nthr(j).gt.2) then
          mint(j,2)=0.d0
!          maxt(j,2)=10.d0
          maxt(j,2)=1.d0
        endif
        maxt(j,nthr(j))=thr(j,nthr(j))

!  si >2 categorias escapo los 2 primeros y el ultimo;
! si no, el 1o y el ultimo
        if (nthr(j).gt.2) then
          do i=3,nthr(j)-1
            mint(j,i)=thr(j,i-1)+1.d-2
            maxt(j,i)=thr(j,i)-1.d-2
          enddo
        else
          do i=2,nthr(j)-1
            mint(j,i)=thr(j,i-1)+1.d-2
            maxt(j,i)=thr(j,i)-1.d-2
          enddo
        endif
        mint(j,nthr(j))=thr(j,nthr(j)-1)+1.d-2
      enddo


! inicializo el thrustoniano
	do ij=1,ncarthur

	   poscar=ncont+ncarumb+ij
	   do i=1,ndat1
	   xmm(i)=-1000.
	   enddo

	   icar=0
	   xd=0
	   do i=1,ndat1
	   if (ipru(i,1,ij).ne.icar) then
	       y(indth(i,ij),poscar)=0.
               icar=ipru(i,1,ij)
	       xd=0.
	   else
	       xd=xd-0.3
	       y(indth(i,ij),poscar)=xd
	   endif
!	   print *,i,ipru(i,1,ij),y(indth(i,ij),poscar),indth(i,ij)
	   enddo
!	   pause
	enddo
!................................................................................................
!     Muestreo de Gibbs
      do 9999 ijk=1,imue
      do 9998 kji=1,icad
!      pause


! ----- Invierto matrices de covarianzas -----------------

        call ginv1(vare,ncar,ncar,tol,irango)
        call ginv1(vara,ncar*nefani,ncar*nefani,tol,irango)

        do i=1,nper
          se=varp(:,:,i)
          call ginv1(se,ncar,ncar,tol,irango)
          varp(:,:,i)=se
        enddo

! ------- aumento de datos para datos faltantes o censurados ---------------
!         y para caracteres discontinuos


!  Caracteres lineales
! Es un poco embarullado y me pierdo un poco
        do j=1,ncont
        do i=1,ndat1

!EV inicializo vark
!         vark=0.d0
         if (reg(i,j).lt.0.01) then
              xm=0.
              do k=1,ncar
! vark es una variable de trabajo
               vark(k)=y(i,k)
              enddo

!              do jk=1,nfac-1
              do jk=1,nfac-nefani
! b es el vector de soluciones
! por tanto esto es Xb
!ev_mad he anhadido *valor(i,jk)
                  xm=xm+b(nfix(jk+nefani,i)+ifac(jk),j)*valor(i,jk)
                  do k=1,ncar
                    if (k.ne.j) then
                      vark(k)=vark(k)-b(nfix(jk+nefani,i)+ifac(jk),k)*valor(i,jk)
                    endif
                  enddo
              enddo
! y esto es Zb (son los animales)
! loop para efectos maternos 
!ev_mad he anhadido *valor(i,nfac-nefani+jk)
              do jk=1,nefani
                xm=xm+b(nfix(jk,i)+ifac(nfac-nefani+jk),j)*valor(i,nfac-nefani+jk)
!              xm=xm+b(nfix(1,i)+ifac(nfac),j)
                do k=1,ncar
                  if (k.ne.j) then
!ev_mad he anhadido *valor(i,nfac-nefani+jk)
                    vark(k)=vark(k)-b(nfix(jk,i)+ifac(nfac-nefani+jk),k)*valor(i,nfac-nefani+jk)
!                    vark(k)=vark(k)-b(nfix(1,i)+ifac(nfac),k)
                  endif
                enddo
            enddo
! e~MVN(e,R)
              do k=1,ncar
                if (k.ne.j) then
                  xm=xm-vark(k)*vare(k,j)/vare(j,j)
                endif
              enddo
              var=1.d0/vare(j,j)

!este para missing value
			  if (reg(i,j).gt. -0.01) then
				  call normal(x1,u)
				  y(i,j)=xm+u*sqrt(var)
			  elseif (reg(i,j) .lt. -0.01) then
!EV  para datos censurados, se codifica la cota inferior como -cota
				  cota1=-reg(i,j) -xm
				  cota2=999.d0 -xm
				  cota1=cota1/sqrt(var)
				  cota2=cota2/sqrt(var)
!				  if ((abs(cota1).gt.8).and.(abs(cota2).gt.8)) goto 891
				  call gen_trunc_normal(cota1,cota2,u,x1)
				  y(i,j)=xm+u*sqrt(var)
!891		                  continue
			  endif
! Var usa una formulacion equivalente q se basa en las 
! propiedades de la inversa de MVN.
         else
              y(i,j)=reg(i,j)
         endif
        enddo
        enddo

! -->   aumento de los caracteres umbrales (segundos)

        do j=1,ncarumb 
          poscar=ncont+j
          do i=1,ndat1

              xm=0.
              do k=1,ncar
               vark(k)=y(i,k)
              enddo

! Xb
!              do jk=1,nfac-1
              do jk=1,nfac-nefani
!ev_mad He anhadido *valor(i,jk)
                  xm=xm+b(nfix(jk+nefani,i)+ifac(jk),poscar)*valor(i,jk)
                                  
                  do k=1,ncar
                    if (k.ne.poscar) then
!ev_mad He anhadido *valor(i,jk)
                      vark(k)=vark(k)-b(nfix(jk+nefani,i)+ifac(jk),k)*valor(i,jk)
                    endif
                  enddo
              enddo
! Zu
              do jk=1,nefani
!ev_mad He anhadido *valor(i,nfac-nefani+jk)
                xm=xm+b(nfix(jk,i)+ifac(nfac-nefani+jk),poscar)*valor(i,nfac-nefani+jk)
                do k=1,ncar
                  if (k.ne.poscar) then
!ev_mad He anhadido *valor(i,nfac-nefani+jk)
                    vark(k)=vark(k)-b(nfix(jk,i)+ifac(nfac-nefani+jk),k)*valor(i,nfac-nefani+jk)
                                  endif
                enddo
            enddo
! e~MVN(e,R)
              do k=1,ncar
                if (k.ne.poscar) then
                  xm=xm-vark(k)*vare(k,poscar)/vare(poscar,poscar)
                endif
              enddo

              var=1.d0/vare(poscar,poscar)
          
!       Bajo que umbral esta?

              presenthr=reg(i,poscar)
! --> generacion de liability si 'missing value' : no hay limites, es igual que un
!  dato lineal
! ahora lo limito entre -999 y +999 (umbrales max y min)
              if (dble(presenthr).lt.0.01d0) then !missing or censored
	        if (dble(presenthr).gt.-0.01d0) then
		  !real missing value, no lower bound
                  cota1=(-999.d0-xm)
                  cota2=(999.d0-xm)
		  !AL 6/4/06 For binary traits to avoid going to infinite the 
  		  !upper and lower bound are 4 (4 s.d. from the threshold)
	          if (nthr(j)==2) then
		    cota1=-liabilitybound-xm
		    cota2=liabilitybound-xm	
		  endif
                  cota1=cota1/sqrt(var)
                  cota2=cota2/sqrt(var)
!		  if ((abs(cota1).gt.8).and.(abs(cota2).gt.8)) goto 993
                  call gen_trunc_normal(cota1,cota2,u,x1)
!                 call normal(x1,u)
                  y(i,poscar)=xm+u*sqrt(var)
!993		  continue
	        else !censored categorical trait - has to be more than 3 categories, otherwise censoring is meaningless
                  if (presenthr.eq.-1) then
                    cota1=-999.d0-xm
                    !AL 6/4/06 For binary traits to avoid going to infinite
		    if (nthr(j)==2) then
		    cota1=-liabilitybound-xm
		  endif
                  cota2=999.d0-xm
                  else
                    cota1=(thr(j,abs(presenthr)-1)-xm)
                    cota2=(999.d0-xm)
		    !AL 6/4/06 For binary traits to avoid going to infinite
		    if (nthr(j)==2) then
		      cota2=+liabilitybound-xm
		    endif
                  endif
                  cota1=cota1/sqrt(var)
                  cota2=cota2/sqrt(var)
!                  if ((abs(cota1).gt.8).and.(abs(cota2).gt.8)) goto 992
		  call gen_trunc_normal(cota1,cota2,u,x1)
                  y(i,poscar)=xm+u*sqrt(var)
!992		  continue
	        endif
		
		  
! --> generacion de liability cuando si se tiene el fenotipo
              else 
                if (presenthr.eq.1) then
                  cota1=-999.d0-xm
		  !AL 6/4/06 For binary traits to avoid going to infinite
		  if (nthr(j)==2) then
		    cota1=-liabilitybound-xm
		  endif
                  cota2=thr(j,presenthr)-xm
                else
                  cota1=(thr(j,presenthr-1)-xm)
                  cota2=(thr(j,presenthr)-xm)
		  !AL 6/4/06 For binary traits to avoid going to infinite
		  if (nthr(j)==2) then
		    cota2=+liabilitybound-xm
		  endif
                endif
                cota1=cota1/sqrt(var)
                cota2=cota2/sqrt(var)
!		if ((abs(cota1).gt.8).and.(abs(cota2).gt.8)) goto 991
                call gen_trunc_normal(cota1,cota2,u,x1)
                y(i,poscar)=xm+u*sqrt(var)
!	        print *,i,y(i,poscar),xm,var
!991		continue
!-- Almaceno minimos y maximos para el muestreo de los umbrales
!    solo cuando son datos "verdaderos"
                if (y(i,poscar).lt.mint(j,presenthr)) then
                  mint(j,presenthr)=y(i,poscar)
                endif
                if (y(i,poscar).gt.maxt(j,presenthr)) then
                  maxt(j,presenthr)=y(i,poscar)
                endif
                if((y(i,poscar).lt.0).and.(presenthr.eq.2)) then
 !                 print *,y(i,poscar),reg(i,poscar),cota1,cota2,u,xm,var
 !                 print *,'kk1'
 !                 print *,vare(poscar,poscar)
                  stop
                endif
              endif
          enddo
        enddo


!	aumento de carÃ¡cteres thustonianos
	do j=1,ncarthur
              poscar=ncont+ncarumb+j     
	      icar=0
!	      print *,j,poscar,ncont,ncarumb
	      do ij=1,ndat1
	        i=indth(ij,j)	      
		xm=0.
                do k=1,ncar
                   vark(k)=y(i,k)
                enddo
                do jk=1,nfac-nefani
                  xm=xm+b(nfix(jk+nefani,i)+ifac(jk),poscar)*valor(i,jk)
                  do k=1,ncar
                    if (k.ne.poscar) then
                      vark(k)=vark(k)-b(nfix(jk+nefani,i)+ifac(jk),k)*valor(i,jk)
                    endif
                  enddo
                enddo
                do jk=1,nefani
                  xm=xm+b(nfix(jk,i)+ifac(nfac-nefani+jk),poscar)*valor(i,nfac-nefani+jk)
                  do k=1,ncar
                    if (k.ne.poscar) then
                      vark(k)=vark(k)-b(nfix(jk,i)+ifac(nfac-nefani+jk),k)*valor(i,nfac-nefani+jk)
                    endif
                  enddo
                enddo
                do k=1,ncar
                  if (k.ne.poscar) then
                    xm=xm-vark(k)*vare(k,poscar)/vare(poscar,poscar)
                   endif
                enddo
                var=1.d0/vare(poscar,poscar)
 
!		print *,i,xm,var,j,ij,ipru(ij,1,j),ipru(ij,2,j)
		ifentest=ipru(ij,2,j)
!		print *,ifentest,ipru(ij,1,j),j
								
!	        algoritmo de los caracters thustonianos
         if (ifentest.eq.0) then
		       call normal(x1,u)
		       y(i,poscar)=xm+u*sqrt(var)
		       goto 8131
         endif

		     if (icar.ne.ipru(ij,1,j)) then
             y(i,poscar)=0
	           icar=ipru(ij,1,j)
	           xmax=y(i,poscar)
	           xmin=xmm(icar)
	           xmm(icar)=-1000.
	           goto 8131
          endif
		      if (icar.eq.ipru(ij,1,j)) then


!	            caso 1		    
 	            if ((ipru(ij+1,1,j).eq.icar).and.(ifentest.ne.99)) then
!	                 print *,'caso1'
			             ifentest2=ipru(ij+1,2,j)
		               if (ifentest2.ne.99) then
	                     xmi=y(indth(ij+1,j),poscar)
	                 else
	                     xmi=xmin
	                 endif
!			 print *,xmi,xmax,1
	                 a=(xmi-xm)/sqrt(var)
	                 b1=(xmax-xm)/sqrt(var)
			             if (abs(a).gt.7.) then
			             if (abs(b1).gt.7.) then
				               goto 777
			             endif
	                 endif
!  			print *,a,b1,1,xm,var
	                call gen_trunc_normal(a,b1,xa,x1)
	                y(i,poscar)=xa*sqrt(var)+xm
777	                xmax=y(i,poscar)	       
	              endif

!	       caso 2
 	            if ((ipru(ij+1,1,j).ne.icar).and.(ifentest.ne.99)) then
!	               print *,'caso2'
                 a=-1000.
	               b1=(xmax-xm)/sqrt(var)
!  		       print *,a,b1,2,xm,var
		             if ((abs(a).gt.7.).and.(abs(b1).gt.7.)) goto 778
	                call gen_trunc_normal(a,b1,xa,x1)
	                 y(i,poscar)=xa*sqrt(var)+xm
778	               xmax=y(i,poscar)	       
	               endif

!		caso 3
	             if (ifentest.eq.99) then
!			print *,'caso3'
!	       PRINT *,XMAX
	                a=-1000.
	                b1=(xmax-xm)/sqrt(var)
!		        PRINT *,a,b1,3,xm,var
			            if ((abs(a).gt.7.).and.(abs(b1).gt.7.)) goto 779
   	                call gen_trunc_normal(a,b1,xa,x1)
	                  y(i,poscar)=xa*sqrt(var)+xm	       
779	                if (y(i,poscar).gt.xmm(icar)) xmm(icar)=y(i,poscar)
	                endif
		            endif	 
8131	          continue
!	       print *,i,y(i,poscar),reg(i,poscar),ipru(i)
	       
	     enddo
	   enddo

!        Fin del loop de aumento de datos

! --------- Muestreo de umbrales ------------------

! Escapo el primer umbral porque es 0 y el ultimo pq es +infinito
! Escapo los dos primeros umbrales (0,1) si hay > 2 categorias

        do j=1,ncarumb 
! Para el umbral 1 (=0)
          maxt(j,1)=0.d0              
! Para el umbral 2 (=10)
          if(nthr(j).gt.2) then
            mint(j,2)=0.d0
            maxt(j,2)=1.d0      
          endif
          if (nthr(j).gt.2) then
            do i=3,nthr(j)-1
              cota1=maxt(j,i)
              cota2=mint(j,i+1)
              if (cota2.lt.cota1) then
                print *,cota1,cota2,'cotas'
                stop
              endif
              call unif(x1,u)
              thr(j,i)=cota1+u*(cota2-cota1)
            enddo

! Borro minimos y maximos para el muestreo de los umbrales
            do i=1,nthr(j)
              maxt(j,i)=-1000.d0
              mint(j,i)=+1000.d0
            enddo
            mint(j,1)=-999.d0
            maxt(j,nthr(j))=999.d0
          endif
        enddo

! ___________________________________________________________________


! ------------- contruccion de rhs --------------------------------

        yy=0.d0

        do i=1,ndat1
           do k=1,ncar
               do j=1,nfac-nefani
! Construye X'y en una columna por caracter
!ev_mad He anhadido valor(i,j)
                 yy(nfix(j+nefani,i)+ifac(J),k)=  &
                 yy(nfix(j+nefani,i)+ifac(J),k)+y(i,k)*valor(i,j)
               enddo               
! Contruye Z'y id.
!               yy(nfix(1,i)+ilev,k)=yy(nfix(1,i)+ilev,k)+y(i,k)
               do j=1,nefani
                 pos1=ifac(nfac-nefani+j)
!ev_mad He anhadido valor(i,j+nfac-nefani)
                 yy(nfix(j,i)+pos1,k)=yy(nfix(j,i)+pos1,k)+y(i,k)*valor(i,j+nfac-nefani)
               enddo
           enddo
        enddo


! ---------- muestreo de efectos fijos -----------------------------------------

         do 1005 i=1,ifac(nfac-nper-nefani+1)
           do 5001 j=1,ncar
               b(i,j)=0.d0
! Primero asigna a b el valor del RHS (xR-1y)
               do 2001 k=1,ncar
                 b(i,j)=b(i,j)+yy(i,k)*vare(j,k)
                 if (j.ne.k) then
! Y quita el efecto de correlacion de los otros caracteres
                   b(i,j)=b(i,j)-b(i,k)*vare(j,k)*dia(i)
                 endif
2001           continue

! Le resta las soluciones de los demas efectos q estan acumuladas en ZHZ
               iplace=ifirst(i)
1006           if(iplace.gt.0)then
                 ji=ivcol(iplace)               
                   do 2002 k=1,ncar
                     b(i,j)=b(i,j)-b(ji,k)*zhz(iplace)*vare(j,k)             
2002               continue
                 iplace=inext(iplace)
                 goto 1006
               endif
! Sistema para que muestree los efectos segun el modelo
               KK=KK2(I,J)
               IF (KK) THEN
                 b(i,j)=b(i,j)/(dia(i)*vare(j,j))
                 var=1./(dia(i)*vare(j,j))
                u=xnor(x1)
!               call normal(x1,u)
                b(i,j)=b(i,j)+u*sqrt(var)
               ELSE
                 B(I,J)=0.d0
               ENDIF
5001       continue
1005     continue

!-------muestreo de efectos fijos aleatorios---------------------------
! loop a lo largo de todos los efectos permanentes
        do iper=1,nper
!        do i=ifac(nfac-nper+iper-1)+1,ifac(nfac-nper+iper)
!        do i=ifac(nfac-nper+iper-nefani)+1,ifac(nfac-nper+iper-nefani+1)
        pos1=ifac(nfac-nper+iper-nefani)+1
      pos2=ifac(nfac-nper+iper-nefani+1)
        do i=pos1,pos2
!           PRINT *,I,DIA(I)
          do j=1,ncar
             b(i,j)=0.d0
             do k=1,ncar
               b(i,j)=b(i,j)+yy(i,k)*vare(j,k)
               if (j.ne.k) then
                 b(i,j)=b(i,j)-b(i,k)*vare(j,k)*dia(i)
                 b(i,j)=b(i,j)-b(i,k)*varp(j,k,iper)
               endif
             enddo    
!            PRINT *,IFIRST(I)
             iplace=ifirst(i)
1426         if(iplace.gt.0)then
               ji=ivcol(iplace)
         
               do k=1,ncar
                 b(i,j)=b(i,j)-b(ji,k)*zhz(iplace)*vare(j,k)
               ENDDO
               iplace=inext(iplace)
               goto 1426
             endif
! Sistema para que muestree los efectos segun el modelo

             KK=KK2(I,J)
             IF(KK) then
               b(i,j)=b(i,j)/(dia(i)*vare(j,j)+varp(j,j,iper))
               var=1./(dia(i)*vare(j,j)+varp(j,j,iper))
               u=xnor(x1)
!              call normal(x1,u)
               b(i,j)=b(i,j)+u*sqrt(var)
             else
               b(i,j)=0.d0
             endif
        
          ENDDO
        ENDDO
        enddo

!------Muestreo de aleatorios -------------------------------------------

       do iefani=1,nefani
       pos1=ifac(nfac-nefani+iefani)+1
       pos2=ifac(nfac-nefani+iefani+1)
       do 1007 i=pos1,pos2
         do 5002 j=1,ncar
!           posvara1=(iefani-1)*ncar+j
           b(i,j)=0.
           do 2003 k=1,ncar
             b(i,j)=b(i,j)+yy(i,k)*vare(j,k)
             if (j.ne.k) then
!              covarianzas residuales
               b(i,j)=b(i,j)-b(i,k)*vare(j,k)*nzz(i)
             endif
2003       continue
!         esto es la correlacion de los otros efectos aditivos
!         del mismo individuo
           do l=1,nefani
             do k=1,ncar
               if (posvara1(iefani,j).ne.posvara2(l,k)) then
!                 posicion de la solucion del animal
!                 b(i,j)=b(i,j)-b(nanim*(l-1)+i,k)*
!                 posani=i-ifac(nfac-nefani+iefani)+ifac(nfac-nefani+l)
                 b(i,j)=b(i,j)- &
                b(i-ifac(nfac-nefani+iefani)+ifac(nfac-nefani+l),k)* &
                vara(posvara1(iefani,j),posvara2(l,k))*dia(i) 
               endif
             enddo
           enddo

!     barrido del resto 
           iplace=ifirst(i)
1008       if(iplace.gt.0)then
             ji=ivcol(iplace)
             do k=1,ncar
! fijos / permanentes/Z1'Z2
!               if((ji.le.pos1).or.(ji.gt.pos2))then
               if((ji.lt.pos1).or.(ji.gt.pos2))then
                 b(i,j)=b(i,j)-b(ji,k)*zhz(iplace)*vare(j,k)
               else 
! A-1*G-1
                 do l=1,nefani
!     posicion de la solucion del animal
!                   b(i,j)=b(i,j)-b(nanim*(l-1)+ji,k)*zhz(iplace)*
!                 posani=ji-ifac(nfac-nefani+iefani)+ifac(nfac-nefani+l)
                   b(i,j)=b(i,j)- &
      b(ji-ifac(nfac-nefani+iefani)+ifac(nfac-nefani+l),k)*zhz(iplace)* &
                        vara(posvara1(iefani,j),posvara2(l,k))
                 enddo
               endif
             enddo
             iplace=inext(iplace)
             goto 1008
           endif
! muestreo segun modelo...
           KK=KK2(I,J)
           IF(KK) then
             var=1./(nzz(i)*vare(j,j)+dia(i)* &
             vara(posvara1(iefani,j),posvara1(iefani,j))) 
             b(i,j)=b(i,j)/(nzz(i)*vare(j,j)+ &
             dia(i)*vara(posvara1(iefani,j),posvara1(iefani,j)))
             u=xnor(x1)
!            call normal(x1,u)
             b(i,j)=b(i,j)+u*sqrt(var)
           else
             b(i,j)=0.d0
           endif

5002     continue
1007   continue
       

! Pone a 0 el ultimo grupo genetico
        do i=1,ncar 
          b(ifac(nfac-nefani+iefani+1),i)=0.d0
        enddo
       enddo
!        fin loop nefani

!-- SAMPLING VARIANCE COMPONENTS ONLY IF DESIRED --


!------Muestreo de varianzas residuales----------------------------------
        call ginv1(vare,ncar,ncar,tol,irango)
!	print *,'vare=',vare
! Utiliza flat priors para todas las varianzas
         se=0.d0
! Calculo de residuos
        do 1011 i=1,ndat1
           do 2009 k=1,ncar
              vark(k)=y(i,k)
!  Fijos y permanentes
              do 1019 j=1,nfac-nefani
!ev_mad he anhadido valor(i, j)
                 vark(k)=vark(k)-b(nfix(j+nefani,i)+ifac(j),k)*valor(i,j)
1019          continue
!  Aleatorios
              do j=1,nefani
!ev_mad he anhadido valor(i, j+nfac-nefani)
                vark(k)=vark(k)-b(nfix(j,i)+ifac(nfac-nefani+j),k)*valor(i,j+nfac-nefani)
              enddo
!                vark(k)=vark(k)-b(nfix(1,i)+ifac(nfac),k)
2009      continue
!   Suma de cuadrados
          do 2010 k=1,ncar
            do 2010 l=1,ncar
             se(k,l)=se(k,l)+vark(k)*vark(l)
2010      continue   
1011    continue

! Se ~ IW

        se=se/real(ndat1-NCAR-1)
       
         call ginv1(se,ncar,ncar,tol,irango)
         call wish(ncar,se,ve,dble(ndat1-NCAR-1),x1)
         call ginv1(ve,ncar,ncar,tol,irango)       

	vare=ve
!	pause
!------Muestreo de varianzas permanentes----------------------------------
! loop para todos los permanentes
        do iper=1,nper
          se=0.d0

! Suma de cuadrados
!          do i=ifac(nfac-nper+iper-1)+1,ifac(nfac-nper+iper)
          pos1=ifac(nfac-nper+iper-nefani)+1
          pos2=ifac(nfac-nper+iper-nefani+1)
          do i=pos1,pos2
            do k=1,ncar
              do l=1,ncar
                 se(k,l)=se(k,l)+b(i,k)*b(i,l)
              enddo
            enddo
          enddo

! Muestreo

              se=se/real(pos2-pos1+1-NCAR-1)

!          do
            call ginv1(se,ncar,ncar,tol,irango)
            call wish(ncar,se,ve,dble(pos2-pos1+1-NCAR-1),x1)
            call ginv1(ve,ncar,ncar,tol,irango)

          varp(:,:,iper)=ve
        enddo

!----------muestreo de varianzas aditivos-------------------------------
! 
        sa=0.d0

! Suma a'A-1a
! A-1 esta almacenado en ZHZ

!        trabajamos con A-1 abajo del todo
        do i=ifac(nfac)+1,ifac(nfac+1)
          do iefani=1,nefani
            pos1=i-ifac(nfac)+ifac(nfac-nefani+iefani)
! Elementos diagonales
            do iefani2=1,nefani
            pos2=i-ifac(nfac)+ifac(nfac-nefani+iefani2)
             do  k=1,ncar
!             posvara1=(iefani-1)*ncar+k
              do  l=1,ncar
!              posvara2=(iefani2-1)*ncar+l
                sa(posvara1(iefani,k),posvara2(iefani2,l))= &
                 sa(posvara1(iefani,k),posvara2(iefani2,l))+ &
                 b(pos1,k)*b(pos2,l)*dia(i) 
              enddo
             enddo
            enddo
          enddo
! elementos no diagonales
          iplace=ifirst(i)
! Salta elementos que no son A-1
1014      if((ivcol(iplace).le.ifac(nfac)).and.(iplace.gt.0))then
            iplace=inext(iplace)
            goto 1014
          endif
1015      if(iplace.gt.0)then
            do iefani=1,nefani
              pos1=i-ifac(nfac)+ifac(nfac-nefani+iefani)
              do iefani2=1,nefani
                pos2=ivcol(iplace)-ifac(nfac)+ifac(nfac-nefani+iefani2)

                do  k=1,ncar
                  do  l=1,ncar
                    sa(posvara1(iefani,k),posvara2(iefani2,l))= &
                     sa(posvara1(iefani,k),posvara2(iefani2,l))+ &
                     b(pos1,k)*b(pos2,l)*zhz(iplace)
                  enddo
                enddo
              enddo
            enddo
            iplace=inext(iplace)
            goto 1015
          endif
        enddo

        na=nanim-ngrup


        sa=sa/real(na-NCAR*nefani-1)       

!      do
        call ginv1(sa,ncar*nefani,ncar*nefani,tol,irango)
        call wish(ncar*nefani,sa,va,dble(na-NCAR*nefani-1),x1)
        call ginv1(va,ncar*nefani,ncar*nefani,tol,irango)



		vara=va


9998  continue 

! reparametrizo las variables......
!  varianza residual
      do i=1,ncar-ncarthur-ncarumb
         vare_a(i,i)=vare(i,i)
      enddo
      do i=ncar-ncarthur-ncarumb+1,ncar
        vare_a(i,i)=1
      enddo

!     modificación de las covarianzas residuales
      do i=1,ncar
         do j=1,ncar
         if (i.ne.j) then
         coef=vare(i,j)/sqrt(vare(i,i)*vare(j,j))
         vare_a(i,j)=sqrt(vare_a(i,i)*vare_a(j,j))*coef
         endif
         enddo
      enddo

!     modificación de las varianzas permanentes
      do k=1,nper
         do i=1,ncar
              coef=vare_a(i,i)/vare(i,i)
              varp_a(i,i,k)=varp(i,i,k)*coef
         enddo
      enddo
      
!     modificación de las covarianzas permanentes
      do k=1,nper
         do i=1,ncar
           do j=1,ncar
              if (i.ne.j) then
                coef=varp(i,j,k)/sqrt(varp(i,i,k)*varp(j,j,k))
                varp_a(i,j,k)=sqrt(varp_a(i,i,k)*varp_a(j,j,k))*coef
               endif
           enddo
         enddo
      enddo

!     modificación de las varianzas aditivas
      ic=0
      do i=1,ncar
      do j=1,nefani
           ic=ic+1
           coef=vare_a(i,i)/vare(i,i)
           vara_a(ic,ic)=vara(ic,ic)*coef
      enddo
      enddo

!    modificación de la covarianzas aditivas

      do i=1,ncar*nefani
           do j=1,ncar*nefani
              if (i.ne.j) then
                coef=vara(i,j)/sqrt(vara(i,i)*vara(j,j))
                vara_a(i,j)=sqrt(vara_a(i,i)*vara_a(j,j))*coef
               endif
         enddo
      enddo

      if ((mod(ijk,50)).eq.0) then
!      do i=1,ncar*nefani
!        print *,(varA(i,j),j=1,ncar*nefani)
!      enddo
      do i=1,ncar*nefani
        print *,(vara_a(i,j),j=1,ncar*nefani)
      enddo

      do iper=1,nper
      do i=1,ncar
        print *,(varp_a(i,j,iper),j=1,ncar)
      enddo
      enddo

!      do i=1,ncar
!        print *,(vare(i,j),j=1,ncar)
!      enddo
      
      do i=1,ncar
        print *,(vare_a(i,j),j=1,ncar)
      enddo

      print *,' imue ',ijk
      call printtime
      endif
      do i=1,ncar
        total=vara(i,i)+vare(i,i)
        do j=1,nper
          total=total+varp(i,i,j)
        enddo
        h2(i)=vara(i,i)/total
      enddo

! Write h2 if animal model only	
	 write(15,'(40f15.8)')  &
	    ((vara_a(i,j),j=i,ncar*nefani),i=1,ncar*nefani), &
	    (((varp_a(i,j,k),j=i,ncar),i=1,ncar),k=1,nper),  &
	    ((vare_a(i,j),j=i,ncar),i=1,ncar)
	
! ---------	
! Contrasts
! ---------
!    uncomment next line if you want contrasts
!    write(20,'(20f15.8)') b(1:3,1),b(1:3,2)
! -------------
! end contrasts
! -------------

!      write(15,'(20f15.8)')((vara(i,j),j=i,ncar),i=1,ncar),
!     +(h2(i),i=1,ncar)
!      write(16,'(20f15.8)')((thr(j,i),i=1,nthr(j)),j=1,ncarumb)

!-- Preparacion de la salida        
      AfterBurnin: if (ijk.gt.lap) then
      
        VarianceComponentEstimation: if(VCE) then
	  do i=1,ncar
            do j=1,ncar
              sumvare(i,j)=sumvare(i,j)+vare_a(i,j)
              do k=1,nper
        	sumvarp(i,j,k)=sumvarp(i,j,k)+varp_a(i,j,k)
              enddo
              ssvare(i,j)=ssvare(i,j)+vare_a(i,j)**2
              do k=1,nper
        	ssvarp(i,j,k)=ssvarp(i,j,k)+varp_a(i,j,k)**2
              enddo
              if(i.ne.j)then
        	do k=1,nper
   !             para evitar problemas en efectos permanentes que no se usan
        	  if ((varp(i,i,k)*varp(j,j,k)).le.0.d0) then
                    corp(i,j,k)=0.d0
        	  else
                    corp(i,j,k)=varp(i,j,k)/sqrt(varp(i,i,k)*varp(j,j,k))
        	  endif
        	enddo
        	core(i,j)=vare(i,j)/sqrt(vare(i,i)*vare(j,j))
              endif
              if(i.eq.j)then
        	total=vare(i,i)
        	do k=1,nefani
                   total=total+vara(posvara1(k,i),posvara1(k,i))
        	enddo
        	do k=1,nper
        	  total=total+varp(i,i,k)
        	enddo
        	do k=1,nefani
        	  cora(posvara1(k,i),posvara1(k,i))= &
        	  vara(posvara1(k,i),posvara1(k,i))/total
        	enddo
        	do k=1,nper
        	  corp(i,j,k)=varp(i,i,k)/total
        	enddo
        	core(i,j)=vare(i,i)/total
              endif
              do k=1,nper
        	sumcorp(i,j,k)=sumcorp(i,j,k)+corp(i,j,k)
              enddo
              sumcore(i,j)=sumcore(i,j)+core(i,j)
              do k=1,nper
        	sscorp(i,j,k)=sscorp(i,j,k)+corp(i,j,k)**2
              enddo
              sscore(i,j)=sscore(i,j)+core(i,j)**2
            enddo
	  enddo
	  do i=1,ncar*nefani
            do j=1,ncar*nefani
              sumvara(i,j)=sumvara(i,j)+vara_a(i,j)
              ssvara(i,j)=ssvara(i,j)+vara_a(i,j)**2
              if(i.ne.j)then
        	if ( (vara(i,i)*vara(j,j)).le.0.d0 )then
        	  cora=0.d0
        	else  
        	  cora(i,j)=vara(i,j)/sqrt(vara(i,i)*vara(j,j))
        	endif
              endif
              sumcora(i,j)=sumcora(i,j)+cora(i,j)
              sscora(i,j)=sscora(i,j)+cora(i,j)**2
            enddo
	  enddo
        endif VarianceComponentEstimation
	

	
	!mean and sd of solutions
	do i=1,nrow
	  do j=1,ncar
	    sumb(i,j)=sumb(i,j)+b(i,j)
	    ssb(i,j)=ssb(i,j)+b(i,j)**2
	  enddo
	enddo
	

        SavingResults:  if (mod(ijk-lap,100).eq.0) then
          
            open(unit=14,file='results.txt',status='replace')
            write(14,'(a,i10)')' Iteration number:',ijk*icad
            write(14,'(a,i10)')' Burn-in:',lap*icad
            write(14,*)
            write(14,*)'Average additive variance'
            do i=1,ncar*nefani
              write(14,'(20f15.8)')((sumvara(i,j)/(ijk-lap)), &
                     j=1,ncar*nefani)
            enddo
            write(14,*)' Sd Additive variance'
            do i=1,ncar*nefani           
              write(14,'(20f15.8)') (sqrt1((ssvara(i,j)-(sumvara(i,j)**2)/ &
               (ijk-lap)) / (ijk-lap-1)),j=1,ncar*nefani)
            enddo
            do k=1,nper
              write(14,*)
              write(14,'(a,i3,a)')' Average environmental variance  ',k,'-th'
              do i=1,ncar
        	write(14,'(20f15.8)')((sumvarp(i,j,k)/(ijk-lap)),j=1,ncar)
              enddo
              write(14,*)' Sd environmental variance'
              do i=1,ncar
        	write(14,'(20f15.8)') (sqrt1((ssvarp(i,j,k)- &
        	 (sumvarp(i,j,k)**2)/ &
               (ijk-lap)) / (ijk-lap-1)),j=1,ncar)
              enddo
            enddo
            write(14,*)
            write(14,*)' Average residual variance'
            do i=1,ncar
              write(14,'(20f15.8)')((sumvare(i,j)/(ijk-lap)),j=1,ncar)
            enddo
            write(14,*)' Sd residual variance '
            do i=1,ncar
              write(14,'(20f15.8)') (sqrt1((ssvare(i,j)-(sumvare(i,j)**2)/ &
             (ijk-lap)) / (ijk-lap-1)),j=1,ncar)
            enddo
	    
	         if(nefani==1) then
	    
            write(14,*)
            write(14,*)' Average h2 and additive correlation '
            do i=1,ncar*nefani
              write(14,'(20f15.8)')((sumcora(i,j)/(ijk-lap)), &
                  j=1,ncar*nefani)
            enddo
            write(14,*)' Sd h2 and additive correlation'
            do i=1,ncar*nefani
              write(14,'(20f15.8)') (sqrt1((sscora(i,j)-(sumcora(i,j)**2)/ &
              (ijk-lap)) / (ijk-lap-1)),j=1,ncar*nefani)
            enddo
            do k=1,nper
              write(14,*)
              write(14,'(a,i3,a)')' Average c2 and environmental cor ',k,'-th'
              do i=1,ncar
        	write(14,'(20f15.8)')((sumcorp(i,j,k)/(ijk-lap)),j=1,ncar)
              enddo
              write(14,*)' Sd c2 and environmental cor'
              do i=1,ncar
        	write(14,'(20f15.8)') (sqrt1((sscorp(i,j,k)- &
        	 (sumcorp(i,j,k)**2)/ &
        	 (ijk-lap)) / (ijk-lap-1)),j=1,ncar)
              enddo
            enddo
            write(14,*)
            write(14,*)' Average he2 and residual cor'
            do i=1,ncar
              write(14,'(20f15.8)')((sumcore(i,j)/(ijk-lap)),j=1,ncar)
            enddo
            write(14,*)' Sd he2 and residual cor'
            do i=1,ncar
              write(14,'(20f15.8)') (sqrt1((sscore(i,j)-(sumcore(i,j)**2)/ &
               (ijk-lap)) / (ijk-lap-1)),j=1,ncar)
            enddo
	    ELSE
	    
	      WRITE(14,*) '--------------------- NOTE -----------------------------------'
	      WRITE(14,*) 'Heritabilities and r_g have to be calculated from ''samples.txt'' '
	      WRITE(14,*) '--------------------------------------------------------------'
	      
	    ENDIF 

            close(14)
            
          open(17,file='solutions.txt',status='replace')
          do i=1,ifac(nfac)
              write(17,'(20f15.8)')(sumb(i,j)/(ijk-lap), &
     sqrt1( (ssb(i,j)-(sumb(i,j)**2)/ &
            (ijk-lap)) / (ijk-lap-1) ),j=1,ncar)
          enddo
          close(17)
          open(18,file='breeding.txt',status='replace')
          do i=ifac(nfac)+1,nrow
              write(18,'(i15,20f15.8)')i-ifac(nfac),(sumb(i,j)/(ijk-lap), &
     sqrt1( (ssb(i,j)-(sumb(i,j)**2)/ &
            (ijk-lap)) / (ijk-lap-1) ),j=1,ncar)
          enddo
          close(18)


        endif SavingResults


      endif Afterburnin

9999  continue 
!-- Fin del gibbs sampler
      stop 
      end 

