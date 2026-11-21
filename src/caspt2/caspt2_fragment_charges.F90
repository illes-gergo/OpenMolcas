!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
!***********************************************************************

subroutine CASPT2_Fragment_Charges(nSym,nBas,bName,Dtri,nDens)

  use Index_Functions, only: iTri, nTri_Elem
  use OneDat, only: sNoNuc, sNoOri
  use Molcas, only: LenIn
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp, u6

  implicit none
  integer(kind=iwp), intent(in) :: nSym, nDens
  integer(kind=iwp), intent(in) :: nBas(nSym)
  character(len=LenIn+8), intent(in) :: bName(*)
  real(kind=wp), intent(in) :: Dtri(nDens)

  integer(kind=iwp) :: i, j, iFrag, iAtom, iBas, iComp, iSyLbl, iOpt, iRc
  integer(kind=iwp) :: nFrag, nCount, nAtom, nPair, iOff
  integer(kind=iwp) :: lWork, info
  integer(kind=iwp) :: nFound
  logical :: Found
  real(kind=wp) :: electron
  real(kind=wp), allocatable :: S(:,:), D(:,:), Shalf(:,:), Tmp(:,:), Eig(:), Work(:), Buf(:)
  real(kind=wp), allocatable :: AtomCharge(:), AtomElectron(:), Zeff(:)
  integer(kind=iwp), allocatable :: AtomCounts(:), Center(:)
  character(len=LenIn), allocatable :: AtomNames(:)
  character(len=8) :: Label

  if (nSym /= 1) then
    write(u6,'(/,6X,A)') 'CASPT2 fragment Loewdin charges require C1 (NoSym) calculations.'
    return
  end

  call qpg_iScalar('nCoordFiles',Found)
  if (.not. Found) return
  call Get_iScalar('nCoordFiles',nFrag)
  if (nFrag < 1) return

  call qpg_iArray('Fragment Atom Counts',Found,nCount)
  if (.not. Found) return
  if (nCount /= nFrag) return
  call mma_allocate(AtomCounts,nFrag,Label='CASPT2 Fragment Atom Counts')
  call Get_iArray('Fragment Atom Counts',AtomCounts,nFrag)

  call Get_iScalar('Unique atoms',nAtom)
  if (nAtom /= sum(AtomCounts)) then
    write(u6,'(/,6X,A)') 'CASPT2 fragment Loewdin charges: fragment atom counts do not match the molecular system.'
    call mma_deallocate(AtomCounts)
    return
  end

  iBas = nBas(1)
  nPair = nTri_Elem(iBas)
  if (nDens /= nPair) then
    call mma_deallocate(AtomCounts)
    return
  end

  call mma_allocate(S,iBas,iBas,Label='CASPT2 Loewdin overlap')
  call mma_allocate(D,iBas,iBas,Label='CASPT2 Loewdin density')
  call mma_allocate(Shalf,iBas,iBas,Label='CASPT2 Loewdin square root')
  call mma_allocate(Tmp,iBas,iBas,Label='CASPT2 Loewdin work')
  call mma_allocate(Eig,iBas,Label='CASPT2 Loewdin eigenvalues')
  call mma_allocate(Buf,nPair+4,Label='CASPT2 overlap integral buffer')
  call mma_allocate(Zeff,nAtom,Label='CASPT2 nuclear charges')
  call mma_allocate(AtomElectron,nAtom,Label='CASPT2 atom populations')
  call mma_allocate(AtomCharge,nAtom,Label='CASPT2 atom charges')

  processing: do
  iRc = -1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iComp = 1
  iSyLbl = 1
  Label = 'Mltpl  0'
  call RdOne(iRc,iOpt,Label,iComp,Buf,nPair,iSyLbl)
  if (iRc /= 0) then
    write(u6,'(/,6X,A)') 'CASPT2 fragment Loewdin charges: failed to read the overlap matrix.'
    exit processing
  end
  call Square(Buf,S,iBas,1,iBas)

  D(:,:) = Zero
  iOff = 0
  do i=1,iBas
    do j=1,i
      D(i,j) = Dtri(iOff+iTri(i,j))
      if (i /= j) D(i,j) = D(i,j)*0.5_wp
      D(j,i) = D(i,j)
    end do
  end do

  call mma_allocate(Work,1,Label='CASPT2 Loewdin query work')
  call DSYEV_('V','U',iBas,S,iBas,Eig,Work,-1_iwp,info)
  if (info /= 0) then
    write(u6,'(/,6X,A)') 'CASPT2 fragment Loewdin charges: failed to size the overlap diagonalization.'
    exit processing
  end
  lWork = max(1_iwp,int(Work(1),kind=iwp))
  call mma_deallocate(Work)
  call mma_allocate(Work,lWork,Label='CASPT2 Loewdin diagonalization work')
  call DSYEV_('V','U',iBas,S,iBas,Eig,Work,lWork,info)
  if (info /= 0) then
    write(u6,'(/,6X,A)') 'CASPT2 fragment Loewdin charges: failed to diagonalize the overlap matrix.'
    exit processing
  end

  Tmp(:,:) = Zero
  do i=1,iBas
    if (Eig(i) <= 1.0e-12_wp) then
      write(u6,'(/,6X,A)') 'CASPT2 fragment Loewdin charges: overlap matrix is linearly dependent.'
      exit processing
    end
    do j=1,iBas
      Tmp(j,i) = S(j,i)*sqrt(Eig(i))
    end do
  end do
  call DGEMM_('N','T',iBas,iBas,iBas,One,Tmp,iBas,S,iBas,Zero,Shalf,iBas)

  call DGEMM_('N','N',iBas,iBas,iBas,One,Shalf,iBas,D,iBas,Zero,Tmp,iBas)
  call DGEMM_('N','N',iBas,iBas,iBas,One,Tmp,iBas,Shalf,iBas,Zero,D,iBas)

  call mma_allocate(Center,iBas,Label='CASPT2 basis center map')
  call qpg_iArray('Center Index',Found,nFound)
  if (Found .and. (nFound >= iBas)) then
    call Get_iArray('Center Index',Center,iBas)
  else
    call mma_allocate(AtomNames,nAtom,Label='CASPT2 atom names')
    call Get_cArray('Unique Atom Names',AtomNames,LenIn*nAtom)
    do i=1,iBas
      Center(i) = 0
      do iAtom=1,nAtom
        if (trim(AtomNames(iAtom)) == trim(adjustl(bName(i)(1:LenIn)))) then
          Center(i) = iAtom
          exit
        end
      end do
    end do
    call mma_deallocate(AtomNames)
  end

  call Get_dArray('Effective nuclear charge',Zeff,nAtom)
  AtomElectron(:) = Zero
  do i=1,iBas
    iAtom = Center(i)
    if ((iAtom >= 1) .and. (iAtom <= nAtom)) AtomElectron(iAtom) = AtomElectron(iAtom)+D(i,i)
  end do
  AtomCharge(:) = Zeff(:)-AtomElectron(:)

  write(u6,'(/,6X,A)') 'CASPT2 fragment Loewdin charge analysis'
  write(u6,'(6X,A)')   '---------------------------------------'
  iAtom = 0
  do iFrag=1,nFrag
    electron = Zero
    do j=1,AtomCounts(iFrag)
      iAtom = iAtom+1
      electron = electron+AtomCharge(iAtom)
    end do
    write(u6,'(8X,A,I0,A,F12.6,A)') 'Fragment ',iFrag,' charge: ',electron,' e'
  end do
  end do processing

  if (allocated(AtomCounts)) call mma_deallocate(AtomCounts)
  if (allocated(Center)) call mma_deallocate(Center)
  if (allocated(Zeff)) call mma_deallocate(Zeff)
  if (allocated(AtomElectron)) call mma_deallocate(AtomElectron)
  if (allocated(AtomCharge)) call mma_deallocate(AtomCharge)
  if (allocated(Work)) call mma_deallocate(Work)
  if (allocated(Buf)) call mma_deallocate(Buf)
  call mma_deallocate(S)
  call mma_deallocate(D)
  call mma_deallocate(Shalf)
  call mma_deallocate(Tmp)
  call mma_deallocate(Eig)

end subroutine CASPT2_Fragment_Charges
