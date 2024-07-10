module initialize_conditions
  implicit none
  private

  public :: initialize_polycrystal 
  public :: initialize_plane
!  public :: initialize_circle
!  public :: assign_additional_precipitates
!  public :: assign_concentration

contains

  ! Subroutine for initializing polycrystal from file
  subroutine initialize_polycrystal(Nx, Ny, Nz, phi, con, ist, ien, var)
    integer, intent(in) :: Nx, Ny, Nz, ist(3), ien(3), var
    logical, intent(in) :: bonus_present
    real *8, intent(in) :: con_0_ppt
    real *8, intent(out) :: phi(1:var,ist(1):ien(1),ist(2):ien(2),ist(3):ien(3)), con(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
    integer :: i, j, k, ii, jj, kk, ivar

    open(3,file='/scratch/jroger87/phase-field-microstructure-evolution/inputs/polycrystal_configs/mc_ivar_fin',status='old')
 
    do k = 1, Nz
    do j = 1, Ny
    do i = 1, Nx
 
      read(3,*) ii, jj, kk, ivar 
      if(ii.ge.ist(1).and.ii.le.ien(1).AND. &
         jj.ge.ist(2).and.jj.le.ien(2).AND. &
         kk.ge.ist(3).and.kk.le.ien(3)) then 
         phi(ivar,ii,jj,kk) = 1.d0 
      end if
    end do
    end do
    end do
 
    close(3)

  end subroutine initialize_polycrystal

  subroutine initialize_plane(Nx, Ny, Nz, phi, ist, ien, var, i_ppt, j_ppt, k_ppt, kount, ppt_rad)
    integer, intent(in) :: Nx, Ny, Nz, var
    integer, intent(inout) :: kount, i_ppt, j_ppt, k_ppt
    integer, intent(in) :: ist(3), ien(3), ppt_rad(3)
    real *8, intent(out) :: phi(1:var,ist(1):ien(1),ist(2):ien(2),ist(3):ien(3))

    integer :: i, j, k, ii, jj, kk, ivar, j1
    real *8 :: term_1, term_2, term_3, sum_ijk

    do kk = ist(3), ien(3)
    do jj = ist(2), ien(2)
    do ii = ist(1), ien(1)

       phi(1,ii,jj,kk) = 1.0

    end do
    end do
    end do

    do kk = ist(3), ien(3)
    do jj = ist(2), ien(2)
    do ii = ist(1), ien(1)

    if(jj.gt.j_ppt) then
       phi(:,ii,jj,kk) = 0.0
       phi(2,ii,jj,kk) = 1.d0
    end if

    end do
    end do
    end do


    do kk = ist(3), ien(3)
    do jj = ist(2), ien(2)
    do ii = ist(1), ien(1)

    j1 = Ny - j_ppt + 1
    if(jj.le.j1) then
       phi(:,ii,jj,kk) = 0.0
       phi(2,ii,jj,kk) = 1.d0
    end if

    end do
    end do
    end do

!   Introduce a precipitate phase at the center 
!   and assign concentrations

    do k = -ppt_rad(3), ppt_rad(3)
    do j = -ppt_rad(2), ppt_rad(2)
    do i = -ppt_rad(1), ppt_rad(1)

!      j_ppt = Ny/2
       term_1 = (float(i)/float(ppt_rad(1)))**2
       term_2 = (float(j)/float(ppt_rad(2)))**2
       term_3 = (float(k)/float(ppt_rad(3)))**2
       sum_ijk = term_1 + term_2 + term_3

       if(sum_ijk.le.1.0) then
         kount = kount + 1
         ii = i + i_ppt
         jj = j + j_ppt
         kk = k + k_ppt
         if(ii.ge.ist(1).AND.ii.le.ien(1).AND. &
         jj.ge.ist(2).AND.jj.le.ien(2).AND. &
         kk.ge.ist(3).AND.kk.le.ien(3)) then
         phi(3,ii,jj,kk) = 1.d0
         phi(1,ii,jj,kk) = 0.d0
         phi(2,ii,jj,kk) = 0.d0
!         con(ii,jj,kk) = con_0_ppt
         end if
       end if

    end do
    end do
    end do

    j_ppt = Ny - j_ppt + 1

    do k = -ppt_rad(3), ppt_rad(3)
    do j = -ppt_rad(2), ppt_rad(2)
    do i = -ppt_rad(1), ppt_rad(1)

       term_1 = (float(i)/float(ppt_rad(1)))**2
       term_2 = (float(j)/float(ppt_rad(2)))**2
       term_3 = (float(k)/float(ppt_rad(3)))**2
       sum_ijk = term_1 + term_2 + term_3

       if(sum_ijk.le.1.0) then
         kount = kount + 1
         ii = i + i_ppt
         jj = j + j_ppt
         kk = k + k_ppt
         if(ii.ge.ist(1).AND.ii.le.ien(1).AND. &
         jj.ge.ist(2).AND.jj.le.ien(2).AND. &
         kk.ge.ist(3).AND.kk.le.ien(3)) then
         phi(3,ii,jj,kk) = 1.d0
         phi(1,ii,jj,kk) = 0.d0
         phi(2,ii,jj,kk) = 0.d0
!         con(ii,jj,kk) = con_0_ppt
         end if
       end if

    end do
    end do
    end do

  end subroutine initialize_plane

!  subroutine initialize_circle(Nx, Ny, Nz, phi, ist, ien, var, kount, ppt_rad, grn_rad, num_ppt, ierr, myid, iseed)
!    integer, intent(in) :: Nx, Ny, Nz, var, grn_rad, num_ppt, myid, iseed
!    integer, intent(in) :: ist(3), ien(3), ppt_rad(3)
!    integer, intent(inout) :: kount, ierr
!    real *8, intent(out) :: phi(1:var,ist(1):ien(1),ist(2):ien(2),ist(3):ien(3))
!
!    integer :: i, j, k, ii, jj, kk, ix, jy, kz, ivar, j1, nppt, phi_count, phi_tot, rad_ppt, radius
!    real *8 :: term_1, term_2, term_3, sum_ijk, ran_2
!
!    phi(1,:,:,:) = 1.0
!
!!  Introduce a circular grain  at the center 
!
!     do k = -grn_rad, grn_rad
!     do j = -grn_rad, grn_rad
!     do i = -grn_rad, grn_rad
!
!        term_1 = (float(i)/float(grn_rad))**2
!        term_2 = (float(j)/float(grn_rad))**2
!        term_3 = (float(k)/float(grn_rad))**2
!        sum_ijk = term_1 + term_2 + term_3
!
!        if(sum_ijk.le.1.0) then
!          kount = kount + 1
!          ii = i + Nx/2
!          jj = j + Ny/2
!          kk = k + Nz/2
!          if(ii.ge.ist(1).AND.ii.le.ien(1).AND. &
!          jj.ge.ist(2).AND.jj.le.ien(2).AND. &
!          kk.ge.ist(3).AND.kk.le.ien(3)) then
!            phi(2,ii,jj,kk) = 1.d0
!            phi(1,ii,jj,kk) = 0.d0
!            phi(3,ii,jj,kk) = 0.d0
!          end if
!        end if
!
!     end do
!     end do
!     end do
!
!!  Introduce multiple precipitates within the circular grain
!
!     nppt = 1
!     rad_ppt = (ppt_rad(1)**2 + ppt_rad(2)**2 + ppt_rad(3)**2)*1.0
!
!     do while(nppt.le.num_ppt)
!
!13     if(myid.eq.0) then
!
!16      ix = ran_2(iseed)*Nx + 1.0
!        jy = ran_2(iseed)*Ny + 1.0
!        kz = ran_2(iseed)*Nz + 1.0
!
!       end if
!
!       call MPI_Bcast(ix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!       call MPI_Bcast(jy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!       call MPI_Bcast(kz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!
!       phi_count = 0
!       do k  = ist(3), ien(3)
!       do j  = ist(2), ien(2)
!       do i  = ist(1), ien(1)
!
!         radius = (i - ix)**2 + (j - jy)**2 + (k-kz)**2
!         if(radius.le.rad_ppt.AND.phi(var,i,j,k).gt.0.0) &
!         phi_count = phi_count +1
!
!       end do
!       end do
!       end do
!
!       call MPI_Allreduce(phi_count, phi_tot, 1, MPI_INTEGER, &
!           MPI_SUM, MPI_COMM_WORLD, ierr)
!
!       if(phi_tot.gt.0) go to 13
!
!       if(phi_tot.eq.0) then
!
!         do k  = ist(3), ien(3)
!         do j  = ist(2), ien(2)
!         do i  = ist(1), ien(1)
!
!             radius = (i - ix)**2 + (j - jy)**2 + (k-kz)**2
!             if(radius.le.rad_ppt) then
!               phi(3,i,j,k) = 1.0
!               phi(2,i,j,k) = 0.0
!               phi(1,i,j,k) = 0.0
!!               con(i,j,k) = con_0_ppt
!             end if
!
!         end do
!         end do
!         end do
!         nppt = nppt + 1
!       end if
!     end do
!     if(myid.eq.0) write(*,*) 'nppt=', nppt
!
!  end subroutine initialize_circle           

!  ! Assign nppt additional precipitates on specified boundary types for wetting analysis
!  subroutine assign_additional_precipitates(Nx, Ny, Nz, phi, con, ist, ien, wetting_param, num_ppt, ppt_rad, myid, var, ierr)
!        implicit none
!
!        integer, intent(in) :: Nx, Ny, Nz, ist(3), ien(3), myid, ierr, var
!        integer, intent(in) :: wetting_param, num_ppt, ppt_rad(3)
!        real *8, intent(inout) :: phi(1:var,ist(1):ien(1),ist(2):ien(2),ist(3):ien(3)), con(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
!
!        integer :: i, j, k, ii, jj, kk, ll, nppt
!
!        integer, allocatable :: boundary(:,:,:), boundary_locations(:,:), precipitate_locations(:,:)
!        integer :: phi_above_threshold, num_boundary_locations, loc_count
!        integer :: ix2, jy2, kz2, rad_ppt
!        real *8 :: threshold, ran_2
!
!        if(myid.eq.0) print *, "assigning additional precipitates"
!        
!        threshold = 0.0 ! 0.0 is great
!        allocate(boundary(ist(1):ien(1),ist(2):ien(2),ist(3):ien(3)))
!
!  !     Introduce multiple precipitates on boundaries 
!  !++++++++++++++++++++++++++
!          ! Determine possible locations for additional precipitates
!          ! Necessary to determine boundary locations using gradient
!          ! method, since phi initialized sharply without transition
!          ! ===============================
!          num_boundary_locations = 0
!  
!          !$omp parallel do collapse(3) &
!          !$omp& private(ii, i, j, k, phi_above_threshold) & 
!          !$omp& shared(ist, ien, var, threshold, phi, boundary, wetting_param) &
!          !$omp& reduction(+:num_boundary_locations)
!!          if (myid.eq.3 .or. myid.eq.14 .or. myid.eq.15) write(*,"(A, 7I8)") "myid, ist, ien", myid, ist, ien
!        
!          do k = ist(3), ien(3) 
!            do j = ist(2), ien(2)
!              do i = ist(1), ien(1)  
!          
!                phi_above_threshold = 0 
!          
!                do ii = 1, var ! loop through all precipitates 
!                  if (phi(ii,i,j,k) > threshold) phi_above_threshold = phi_above_threshold + 1 
!                end do
!
!                boundary(i,j,k) = phi_above_threshold ! if more than one phi is nonzero, that signifies a boundary
!                if (phi_above_threshold == 2) num_boundary_locations = num_boundary_locations + 1
!          
!              end do
!            end do
!          end do
!        
!          !$omp end parallel do
!          ! ===============================
!  
!          if (num_boundary_locations > 0) then
!            
!            allocate(boundary_locations(num_boundary_locations, 3))
!            loc_count = 0
!             
!            ! Second pass to store boundary locations
!            do k = 1, Nz
!              do j = 1, Ny
!                do i = 1, Nx
!                  if (boundary(i,j,k) == wetting_param) then
!                    loc_count = loc_count + 1
!                    boundary_locations(loc_count, :) = (/ i, j, k /)
!                  end if
!                end do
!              end do
!            end do
!          else
!            print *, "No boundary locations found."
!            stop
!          end if 
!          call MPI_Barrier(MPI_COMM_WORLD, ierr)
!          print *, "finished boundary location thing: ", myid
!      
!        ! ===============================
!          nppt = 1
!          rad_ppt = (ppt_rad(1)**2 + ppt_rad(2)**2 + ppt_rad(3)**2)*1.0
!  
!          call MPI_Barrier(MPI_COMM_WORLD, ierr)
!          allocate(precipitate_locations(num_ppt, 3))
!  
!          do while (nppt.le.num_ppt)
!            ! Choose additional precipitates
!            if (myid == 0) then
!              ll = int(ran_2(iseed)*num_boundary_locations) + 1 ! integer for the ll_th boundary voxel (1, num_boundary_locations)
!              ix = boundary_locations(ll, 1)
!              jy = boundary_locations(ll, 2)
!              kz = boundary_locations(ll, 3)
!            end if      
!  
!            call MPI_Bcast(ix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!            call MPI_Bcast(jy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!            call MPI_Bcast(kz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!    
!            ! Ensure additional precipitates do not overlap
!            phi_count = 0
!            if (nppt>1) then
!              do ll = 1, nppt-1
!                write(*,'(A, 3I8, 6I8)') "= myid, nppt, ll &
!                    (x2, x1, y2, y1, z2, z1)", myid, nppt, ll, ix2, ix, jy2, jy, kz2, kz
!    
!                ix2 = precipitate_locations(ll, 1)
!                jy2 = precipitate_locations(ll, 2)
!                kz2 = precipitate_locations(ll, 3)
!    
!                radius = (ix - ix2)**2 + (jy - jy2)**2 + (kz - kz2)**2
!                !if(radius.le.rad_ppt.AND.phi(var,i,j,k).gt.0.0) &
!    
!                if(radius.le.rad_ppt) &
!                phi_count = phi_count +1
!              end do
!            end if
!            if(myid.eq.0) print *, "here 2"
!  
!            call MPI_Allreduce(phi_count, phi_tot, 1, MPI_INTEGER, &
!                   MPI_SUM, MPI_COMM_WORLD, ierr)
!  
!            if (phi_tot .gt. 0) cycle
!            if (phi_tot .eq. 0) then  
!              ! Assign additional precipitates
!              rad_ppt = (ppt_rad(1)**2 + ppt_rad(2)**2 + ppt_rad(3)**2)*1.0
!              do k  = ist(3), ien(3)
!                do j  = ist(2), ien(2)
!                  do i  = ist(1), ien(1)
!      
!                    radius = (i - ix)**2 + (j - jy)**2 + (k-kz)**2
!                    if(radius.le.rad_ppt) then
!                      phi(3,i,j,k) = 1.0
!                      phi(2,i,j,k) = 0.0
!                      phi(1,i,j,k) = 0.0
!                      con(i,j,k) = con_0_ppt
!  
!                      ! Add precipitate location to list for overlap checking
!                      ll = nppt
!                      precipitate_locations(ll, :) = (/ i, j, k /)
!                    end if
!      
!                  end do
!                end do
!              end do
!              nppt = nppt + 1
!            end if
!  
!          end do
!        ! ===============================
!          if(myid.eq.0) print *, "here 3"
!      
!          deallocate(boundary)
!          deallocate(boundary_locations)
!          deallocate(precipitate_locations)
!
!          if (myid.eq.0) write(*,*) 'nppt=', nppt
!  !++++++++++++++++++++++++++
!
!  end subroutine assign_additional_precipitates

end module initialize_conditions

