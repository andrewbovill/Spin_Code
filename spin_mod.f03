      module spin_mod
!
!     This module supports the main program "spin.f03".
!     Contains subroutines called in main program 
!
!     -A. J. Bovill, 2022.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use mqc_algebra
      use iso_fortran_env

!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
      logical::DEBUG=.false.

      CONTAINS

!
!PROCEDURE Sz**2 value
!
      subroutine SpinsquaredSz(NElectrons,NElAlpha,NElBeta,SZsquared)
!
!     Calculates Spinsquared Sz value. Equation from 
!     Szabo and Ostlund page 107
!
      implicit none

      integer(kind=int64),intent(in)::NElectrons,NElAlpha,NElBeta
      type(MQC_Variable),intent(out)::SZsquared
      
!     Need float values from all integers

      SZsquared = (float(NElAlpha - NElBeta)/float(2))*((float(NElAlpha - NElBeta)/float(2))+float(1))

      end subroutine SpinsquaredSz
!
!PROCEDURE S**2 value
!
      subroutine Spinsquaredtot(NElectrons,NElAlpha,NElBeta,SZsquared,Calpha,Cbeta,Smatrix,Stotsquared)

      implicit none

      integer(kind=int64),intent(in)::NElectrons,NElAlpha,NElBeta
      integer::i,j
      real(kind=real64)::Overlaptmp
      type(MQC_Variable),intent(in)::SZsquared,Calpha,Cbeta,Smatrix
      type(MQC_Variable),intent(out)::Stotsquared

!     Temp matrices for Matmul operations
!     MQC does not incorporate slicing with object types
!     So you need intrinsic fortran matrices.

      type(MQC_Variable)::tmpMQCvar,tmpMQCvar1,tmpMQCvar2
      real(kind=real64),allocatable,dimension(:,:)::tmpMatrixS,  &
        tmpMatrixCa,tmpMatrixCb,Overlaptot,Sz

      tmpMatrixS = Smatrix
      tmpMatrixCa = Calpha
      tmpMatrixCb = Cbeta

      Overlaptot = MatMul(Transpose(tmpMatrixCa(:,1:NElAlpha)),  &
        MatMul(tmpMatrixS,tmpMatrixCb(:,1:NElBeta)))

      Overlaptmp = 0
      do i = 1,NElAlpha
        do j = 1,NElBeta
          Overlaptmp = Overlaptmp + abs(Overlaptot(i,j)*Overlaptot(i,j))
        end do
      end do

      tmpMQCvar = NElBeta
      tmpMQCvar1 = tmp
!     Sz = Szsquared

      Stotsquared = Szsquared + tmpMQCvar - tmpMQCvar1
!     write(*,*) 'Sz**2', Sz

      end subroutine Spinsquaredtot

      end module spin_mod
