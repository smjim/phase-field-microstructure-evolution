module initialize_conditions
  implicit none
  private

  public :: print_hi
  public :: initialize_polycrystal 

contains

  ! Subroutine for print hi
  subroutine print_hi(Nx, Ny, Nz)
    integer, intent(in) :: Nx, Ny, Nz
    write(*, '(A, 3I8)') "Nx, Ny, Nz: ", Nx, Ny, Nz
  end subroutine print_hi

  subroutine initialize_polycrystal(Nx, Ny, Nz, phi, con, ist, ien, con_0_ppt, var, myid)
    integer, intent(in) :: Nx, Ny, Nz, ist(3), ien(3), var, myid
    real *8, intent(in) :: con_0_ppt    
    real *8, intent(out) :: phi(:,:,:,:), con(:,:,:)
    integer :: i, j, k, ii, jj, kk, ivar

    write(*,'(A, 7I8)') "myid, ist, ien", myid, ist, ien
    write(*,'(4I8)') size(phi, dim=1), size(phi, dim=2), size(phi, dim=3), size(phi, dim=4)

    open(3,file='/scratch/jroger87/phase-field-microstructure-evolution/inputs/polycrystal_configs/mc_ivar_fin',status='old')
 
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx
 
      read(3,*) ii, jj, kk, ivar 
      if(ii.ge.ist(1).and.ii.le.ien(1).AND. &
         jj.ge.ist(2).and.jj.le.ien(2).AND. &
         kk.ge.ist(3).and.kk.le.ien(3)) then 
!         phi(ivar,ii,jj,kk) = 1.d0 
         phi(ivar,ii,jj,kk)=0.0
         if(ivar.eq.var) con(ii,jj,kk) = con_0_ppt
      end if
    end do
    end do
    end do
 
    close(3)
    print *, "finished polycrystal"

  end subroutine initialize_polycrystal


end module initialize_conditions

