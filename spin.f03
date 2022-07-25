INCLUDE 'spin_mod.f03'
      program spin
!
!     This program calculates Spin Squared values from matrix files.
!     This code's structure is heavily derived from H.P.Hratchian
!     NIO code https://github.com/hphratchian/nio.git
!
!     The purpose of this code is just a 'proof of concept' nothing more
!
!     -A. J. Bovill, 2022.
!
!
!     USE Connections
!
      use spin_mod
!
!     Variable Declarations
!
      implicit none

      integer(kind=int64)::nCommands,i,nAtoms,nAtoms2,  &
        nBasis,nBasis2,nBasisUse,nBasisUse2,nEl1,nEl2,nElAlpha1,  &
        nElBeta1,nElAlpha2,NElBeta2
      character(len=512)::matrixFilename1
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile1

      type(MQC_Variable)::PMatrixAlpha,PMatrixBeta
      type(MQC_Variable)::CAlpha,CBeta
      type(MQC_Variable)::SMatrixAO,SMatrixEVecs,SMatrixEVals
      type(MQC_Variable)::SSquaredtot, SSquaredz

!
!     Format Statements !
 1000 Format(1x,'Enter Program NIO.')
 1010 Format(1x,'Matrix File 1: ',A,/)
 1020 Format(1x,'nAtoms=',I4,3x,'nBasis   =',I4,3x,'nBasisUse=',I4,/,  &
             1x,'nEl1  =',I4,3x,'nElAlpha1=',I4,3x,'nElBeta  =',I4,/)  
 1030 Format(1x,'Hello World!!')
 8999 Format(/,1x,'END OF SPIN PROGRAM')

      write(IOut,1030)
      call mqc_version_print(iOut)

!     Open the Gaussian matrix file and load the number of atomic centers.
      nCommands = command_argument_count()
      if(nCommands.ne.1)  &
        call mqc_error('Just one Matrix file in command line.')
      call get_command_argument(1,matrixFilename1)
      call GMatrixFile1%load(matrixFilename1)
      write(IOut,1010) TRIM(matrixFilename1)
!
!     Do some consistency checks and load the number of atoms, basis functions,
!     and linearly independent basis functions.
!
      nAtoms  = GMatrixFile1%getVal('nAtoms')
      nBasis  = GMatrixFile1%getVal('nBasis')
      nBasisUse  = GMatrixFile1%getVal('nBasisUse')
      nEl1      = GMatrixFile1%getVal('nElectrons')
      nElAlpha1 = GMatrixFile1%getVal('nAlpha')
      nElBeta1  = GMatrixFile1%getVal('nBeta')

      write(IOut,1020) nAtoms,nBasis,nBasisUse,nEl1,nElAlpha1,nElBeta1  

!
!     Obtain Overlap, Calpha, Cbeta matrices for calculations
!

      call GMatrixFile1%getArray('OVERLAP',mqcVarOut=SMatrixAO)
      call GMatrixFile1%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha)
      call GMatrixFile1%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=Calpha)
      call GMatrixFile1%getArray('BETA MO COEFFICIENTS',mqcVarOut=Cbeta)

!
!     Calculate Sz**2 
!
      call SpinsquaredSz(nEl1,nElAlpha1,nElBeta1,SSquaredZ)
      call SSquaredZ%print(header='Sz**2')
!
!     Calculate S**2 
!

      call Spinsquaredtot(nEl1,nElAlpha1,nElBeta1,SSquaredZ,CAlpha,CBeta,SMatrixAO,SSquaredtot)
      call SSquaredtot%print(header='S**2')

  999 Continue
      write(IOut,8999)

      End program spin
