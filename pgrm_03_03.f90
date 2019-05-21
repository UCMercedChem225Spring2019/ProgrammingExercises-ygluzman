      program pgrm_03_03
!
!     This program computes components of the Hartree-Fock energy and the number
!     of electrons using contraction of various matrices loaded from
!     user-provided files (that presumably come from Gaussian calculations).
!
!     At run time, the program expects 5 command line arguments:
!       (1) the number of electrons;
!       (2) the number of atomic orbital basis functions;
!       (3) an input file containing core-Hamiltonian matrix elements (in
!           symmetric upper/column storage form);
!       (4) an input file containing the overlap matrix elements (in
!           symmetric upper/column storage form); and
!       (5) an input file containing Fock matrix elements (in symmetric
!           upper/column storage form).
!
!     At run time, the program outputs 6 items:
!       (1) the MO energies;
!       (2) the MO coefficients;
!       (3) the density matrix;
!       (4) the one-electron contribution to the total electronic energy;
!       (5) the two-electron contribution to the total electronic energy; and
!       (6) the tr(PS), which should equal the number of electrons.
!
!
!
      implicit none
      integer,parameter::unitIn=10
      integer::i,iError,nElectrons,nOcc,nBasis,lenSym,j,k
      real::oneElectronEnergy,twoElectronEnergy,tracePS
      real,dimension(:),allocatable::symFock,symCoreHamiltonian,& 
        symOverlap,moEnergies,tempSymMatrix
      real,dimension(:,:),allocatable::sqFock,fockTilde,& 
        sqCoreHamiltonian,InvSqrtOverlap,moCoefficients,densityMatrix,& 
        tempSqMatrix,SSPEV_Scratch
      character(len=256)::cmdlineArg
!
!
!     Begin by reading the number of basis functions, allocating array memory,
!     and loading the symmetric Fock and overlap matrices from input files
!     provided on the command line.
!
      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nElectrons
      nOcc = nElectrons/2
!
      call Get_Command_Argument(2,cmdlineArg)
      read(cmdlineArg,'(I)') nBasis
      lenSym = (nBasis*(nBasis+1))/2
      allocate(symFock(lenSym),symCoreHamiltonian(lenSym),  &
        symOverlap(lenSym),tempSymMatrix(lenSym))
      allocate(moEnergies(nBasis))
      allocate(sqFock(nBasis,nBasis),fockTilde(nBasis,nBasis),  &
        sqCoreHamiltonian(nBasis,nBasis),invSqrtOverlap(nBasis,nBasis),&
        moCoefficients(nBasis,nBasis),densityMatrix(nBasis,nBasis),  &
        tempSqMatrix(nBasis,nBasis),SSPEV_Scratch(nBasis,3))

!
!
!
      Call Get_Command_Argument(3,cmdlineArg)
      Open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=IError)
      If(IError.ne.0) then
        Write(*,*)' Error opening input file.1'
        STOP
      endIf
      DO i = 1, lenSym !Read triangle elements
         read(UnitIn,*) symCoreHamiltonian(i)
      END DO
      Close(Unit=UnitIn)
        

      Call Get_Command_Argument(4,cmdlineArg)
      Open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=IError)
      If(IError.ne.0) then
        Write(*,*)' Error opening input file.2'
        STOP
      endIf
      DO i = 1, lenSym !Read triangle elements
         read(unitIn,*) symOverlap(i)
      END DO
      Close(Unit=unitIn)

      Call Get_Command_Argument(5,cmdlineArg)
      Open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=IError)
      If(IError.ne.0) then
        Write(*,*)' Error opening input file.3'
        STOP
      endIf
      DO i = 1, lenSym !Read triangle elements
         read(unitIn,*) symFock(i)
      END DO
      Close(Unit=unitIn)

      call SymmetricPacked2Matrix_UpperPacked(nBasis,symOverlap,& 
        invSqrtOverlap)
!
!     Form the square-root of the overlap matrix.
      tempSymMatrix = symOverlap
      call InvSQRT_SymMatrix(nBasis,tempsymMatrix,invSqrtOverlap)
!     Form fTilde and solve for the MO energies and coefficients.
!
      call SymmetricPacked2Matrix_UpperPacked(nBasis,symFock,sqFock)
      tempSqMatrix = MatMul(invSqrtOverlap,sqFock)
      fockTilde = MatMul(tempSqMatrix,invSqrtOverlap)
      call Sq2SymMatrix(nBasis,fockTilde,tempSymMatrix)
      call SSPEV('V','U',nBasis,tempSymMatrix,moEnergies,  &
        moCoefficients,nBasis,SSPEV_Scratch,iError)
      If(iError.ne.0) then
        write(*,*)' Failure in SSPEV.'
        STOP
      endIf
      write(*,*)' MO Energies:'
      call Print_Matrix_Full_Real(Reshape(moEnergies,[nBasis,1]),nBasis,1)
      tempSqMatrix = moCoefficients
      moCoefficients = MatMul(invSqrtOverlap,tempSqMatrix)
      write(*,*)' MO Coefficients:'
      call Print_Matrix_Full_Real(moCoefficients,nBasis,nBasis)
      
      
      densityMatrix=2*matmul(mocoefficients,transpose(moCoefficients))
      write(*,*) "Density Matrix: "
      call Print_Matrix_Full_Real(densitymatrix,nbasis,nbasis)
      call SymmetricPacked2Matrix_UpperPacked(nBasis, &
           symCoreHamiltonian,sqCoreHamiltonian) 
      tempSqMatrix = MatMul(densityMatrix,sqCoreHamiltonian)
      
      do i = 1,nbasis
        oneElectronEnergy = oneElectronEnergy + tempSqMatrix(i,i)
      enddo

      write(*,*)' One electron energy contribution: ',oneElectronEnergy
      tempSqMatrix = sqFock - sqCoreHamiltonian  
      tempSqMatrix = (MatMul(densityMatrix,tempSqMatrix))/2
      Do i=1, nbasis
        twoElectronEnergy = twoElectronEnergy + tempSqMatrix(i,i)
      enddo

      write(*,*)'Two electron energy contribution: ', &
        twoelectronenergy
      call SymmetricPacked2Matrix_UpperPacked(nBasis, &
           symOverlap, tempSqMatrix)
      tempSqMatrix = MatMul(densityMatrix,tempSqMatrix)
      
      
      
      do i=1, nBasis 
        tracePS = tracePS + tempSqMatrix(i,i)
      enddo
      write(*,*)' The number of electrons is: ',tracePS
      end program pgrm_03_03

      Subroutine InvSQRT_SymMatrix(N,inputSymMatrix,Mat)

        Implicit None
        Integer,Intent(In)::N
        Real,Dimension((N*(N+1))/2),Intent(In)::inputSymMatrix
        Real,Dimension(N,N),Intent(Out)::Mat

        Integer::i,j,k,IError,s
        Real,Dimension(:,:),Allocatable::EVecs,Temp_Matrix,Evalmat
        Real,Dimension(:),Allocatable::EVals,Temp_Vector,Temp2
        Allocate(EVals(N),EVecs(N,N),Temp_Vector(3*N),Evalmat(N,N))
        Allocate(Temp_Matrix(N,N),Temp2((N*(N+1))/2))
        Do i = 1,(N*(N+1))/2
          Temp2(i)= inputSymMatrix(i)
        EndDo
        Call SSPEV("V","U",N,Temp2,EVals,EVecs,N,  &
        Temp_Vector,IError)

        Do i= 1,N
         Evalmat(i,i) =1/(sqrt(Evals(i)))
        Enddo

        Mat = Matmul(Matmul(EVecs,Evalmat),transpose(EVecs))

        Return
      End Subroutine
      Subroutine SymmetricPacked2Matrix_UpperPacked(N,ArrayIn,AMatOut)
!
!     This subroutine accepts an array, ArrayIn, that is (N*(N+1))/2
!     long.
!     It then converts that form to the N-by-N matrix AMatOut taking
!     ArrayIn to be in upper-packed storage form. Note: The storage mode
!     also assumes the upper-packed storage is packed by columns.
!
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension((N*(N+1))/2),Intent(In)::ArrayIn
      Real,Dimension(N,N),Intent(Out)::AMatOut
!
      Integer::i,j,k,s
!
!     Loop through the elements of AMatOut and fill them appropriately
!     from
!     Array_Input.
!
!
! ********************************************************
         Do j=1,N
          Do i = 1,N
            AMatOut(i,j) = ArrayIn(j*(j-1)/2+i)
          EndDo
        EndDo
         Do  i = 1,N
          Do  j = 1,i
             AMatOut(i,j) = AMatOut(j,i)
          EndDo
        EndDo

! *************************************************************************
!
!
      Return
      End Subroutine SymmetricPacked2Matrix_UpperPacked
      subroutine Sq2SymMatrix(N,Sq,Sym)

        Implicit None
        Integer,Intent(In)::N
        Real,Dimension(N,N),Intent(in)::Sq
        Real,Dimension((N*(N+1))/2),Intent(out)::Sym


        Integer::i,j

        Do j=1,N
          Do i=1,j
            Sym(j*(j-1)/2+i) = Sq(i,j)
        enddo
          enddo
      end subroutine Sq2SymMatrix
      Subroutine Print_Matrix_Full_Real(AMat,M,N)
!
!     This subroutine prints a real matrix that is fully dimension -
!     i.e.,
!     not stored in packed form. AMat is the matrix, which is
!     dimensioned
!     (M,N).
!
!     The output of this routine is sent to unit number 6 (set by the
!     local
!     parameter integer IOut).
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::M,N
      real,dimension(M,N),intent(in)::AMat
!
!     Local variables
      integer,parameter::IOut=6,NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
!
      Return
      End Subroutine Print_Matrix_Full_Real


