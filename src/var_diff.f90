      program var_diff 
      !$acc routine (eldis) gang

      use p3dfft
      use mpi
      
!      use initialize_conditions
      implicit none

      integer Nx, Ny, Nz, ndim, dims(2), kount , N_step, step
      integer ist(3), ien(3), isize(3), fst(3), fen(3), fsize(3)
      integer i, j, k, l, ii, jj, kk, var, ivar, jvar, iseed, ix, jy, kz
      integer ppt_rad(3), grn_rad, ifreq, mat_var, phi_count, phi_tot
      integer rank, nprocs, myid, ierr , i_ppt, j_ppt, k_ppt, nppt, num_ppt
      integer j1, j2, nrun, rad_ppt, radius, phi_num, ppt_int
      real *8 dx, dy, dz, t_step , c_av, time(6), run_time
      real *8 term_1, term_2, term_3, sum_ijk, term
      real *8 grad_coeff_phi, grad_coeff_c, am, ap, mob_phi 
      real *8 con_0_mat, con_0_ppt, gb_force, phi_mat, phi_ppt
      real *8 fi_i, fi_j, gb_en, fch_mat, fch_ppt, beta_en , beta_en_tot
      real *8 num_ijk, den_ijk, phi_1, phi_2, phi_3, con_ijk, den, phi_3_max
      real *8 f_int_tot, f_eta_tot, f_ch_tot, e_tot
      real *8 f_int, f_eta, f_ch, c_min, c_min_glob, sigma_1, sigma_2
      real *8 D_gb, D_mean, phi_sum, fi, D_vol, D_S, D_L, glob_sum 
      real *8 rad_test, c_mat, c_ppt
      real *8 c_max, c_max_glob, phi_max, phi_max_glob
      real *8  ran_2, phi_min, phi_min_glob, phi_sq
      real *8  mob_min, mob_max, min_glob, max_glob, ave_mob, glob_mob
      character start*10, file_num*4, dir_name*256 
      character input_file_path*256, output_dir*256

      real *8, dimension(:,:), allocatable :: gamma
      real *8, dimension(:,:,:),  allocatable :: con
      real *8, dimension(:,:,:),  allocatable :: df_dc
      real *8, dimension(:,:,:),  allocatable :: mob_c
      real *8, dimension(:,:,:),  allocatable :: dummy 
      real *8, dimension(:,:,:),  allocatable :: dummy2 
      real *8, dimension(:,:,:),  allocatable :: dummy_x, dummy_y, dummy_z
      real *8, dimension(:,:,:),  allocatable :: numr
      real *8, dimension(:,:,:),  allocatable :: denr 
      real *8, dimension(:,:,:),  allocatable :: grad_x_c 
      real *8, dimension(:,:,:),  allocatable :: grad_y_c 
      real *8, dimension(:,:,:),  allocatable :: grad_z_c 
      real *8, dimension(:,:,:,:),  allocatable :: phi 
      real *8, dimension(:,:,:,:),  allocatable :: df_dphi 
      real *8, dimension(:,:,:,:),  allocatable :: grad_x 
      real *8, dimension(:,:,:,:),  allocatable :: grad_y 
      real *8, dimension(:,:,:,:),  allocatable :: grad_z 
      real *8, dimension(:,:,:,:),  allocatable :: term1 
      real *8, dimension(:,:,:,:),  allocatable :: term2 
      real *8, dimension(:,:,:,:),  allocatable :: term3 
      real *8, dimension(:,:,:,:),  allocatable :: term4 
      real *8, dimension(:,:,:,:),  allocatable :: term5 
      real *8, dimension(:),  allocatable :: D_bulk

      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_con
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_df_dc
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_dummy
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_dummy_x
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_dummy_y
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_dummy_z
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_dummy2
      complex(p3dfft_type), dimension(:,:,:,:),  allocatable :: dft_grad_x 
      complex(p3dfft_type), dimension(:,:,:,:),  allocatable :: dft_grad_y 
      complex(p3dfft_type), dimension(:,:,:,:),  allocatable :: dft_grad_z 
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_grad_x_c 
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_grad_y_c 
      complex(p3dfft_type), dimension(:,:,:),  allocatable :: dft_grad_z_c 

      double complex, dimension(:,:,:,:), allocatable :: kf, kf_sq, kf_4
      double complex  k_sq, k_4, kf_sum, k_cube
      complex(p3dfft_type) term11, term22 , term33, term44

      real *8 :: loop_time(14), loop_duration(14), avg_loop_time(14)

      ! Input filename specification
      call get_command_argument(1, input_file_path, status=ierr)
    
      ! Output directory specification
      call get_command_argument(2, output_dir, status=ierr)

      open(93,file=trim(output_dir)//'time_step.dat')
      write(93, '(A)') "# step number, iteration duration, &
                              run time elapsed" 

      open(95,file=trim(output_dir)//'timing_data.dat')
      write(95, '(A)') "# step number, loop 1 avg, loop 2 avg, loop 3 avg, loop 4 avg, loop 5 avg, loop 6 avg, loop 7 avg, &
                                       loop 8 avg, loop 9 avg, loop 10 avg, loop 11 avg, loop 12 avg, loop 13 avg, loop 14 avg, &
                                       loop duration, run time elapsed"

      open(97,file=trim(output_dir)//'gamma_matrix.dat')

      open(98,file=trim(output_dir)//'summary.dat')
      write(98, '(A)') "# step number, f_int_tot, f_ch_tot, &
        f_eta_tot, e_tot, c_max_glob, c_min_glob, phi_min_glob, phi_3_max" 


!     MPI Initializations
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

!     Read data from input file
      open(2, file=input_file_path, status='old')

      read(2,*) nrun, ppt_int
      read(2,*) ppt_rad(1), ppt_rad(2), ppt_rad(3), N_step, ifreq
      read(2,*) i_ppt, j_ppt, k_ppt, iseed
      read(2,*) grn_rad, num_ppt, c_mat, c_ppt
      read(2,*) Nx, Ny, Nz, var, dx, dy, dz, t_step
      read(2,*) grad_coeff_phi, grad_coeff_c, D_S, D_L, D_gb, mob_phi
      read(2,*) sigma_1, sigma_2
      read(2,*) am, ap, con_0_mat, con_0_ppt, gb_force
      read(2,*) ndim
      if ( ndim .eq. 1 ) then
         dims(1) = 1
         dims(2) = nprocs
      else if ( ndim .eq. 2 ) then
         read(2,*) dims(1), dims(2)
      end if
      read(2,*) start

!     Default value of free energy from thermocalc is J/mole
!     Convert to J/m^3 using molar value of Ti-6Al-4V

!     Set up work structures for P3DFFT
!     -- note use of TRANSPOSED arrays
      call p3dfft_setup (dims, Nx, Ny, Nz, MPI_COMM_WORLD,Nx,Ny,Nz,.true.)
      if(myid.eq.0) write(*,*) dims(1), dims(2)

!     Get dimensions for the original array of real numbers, X-pencils
      call p3dfft_get_dims(ist, ien, isize, 1)
      if(myid.eq.0) write(*,*) ist, ien 

!     Get dimensions for the R2C-forward-transformed array of complex numbers
!     Z-pencils
      call p3dfft_get_dims(fst, fen, fsize, 2)
      if(myid.eq.0) write(*,*) fst, fen 

!     Determine local array sizes and allocate arrays
      allocate (gamma(var,var))
      allocate (D_bulk(var))
      allocate( con(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( mob_c(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( df_dc(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( dummy(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( dummy2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( dummy_x(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( dummy_y(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( dummy_z(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( numr(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( denr(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )

      allocate( phi(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( grad_x(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( grad_y(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( grad_z(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( grad_x_c(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( grad_y_c(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( grad_z_c(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( term1(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( term2(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( term3(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( term4(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( term5(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate( df_dphi(var,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )

      allocate( dft_con(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_df_dc(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_dummy(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_dummy_x(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_dummy_y(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_dummy_z(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_dummy2(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_grad_x(var, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_grad_y(var, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_grad_z(var, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_grad_x_c(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_grad_y_c(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( dft_grad_z_c(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )

      allocate( kf(3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( kf_sq(3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate( kf_4(3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )

      write(file_num,'(i4.4)') myid+1

      do ivar = 1, var
        if(ivar.lt.var) then
          D_Bulk(ivar) = D_S
        else
          D_Bulk(ivar) = D_L
        end if
      end do

      phi = 0.d0
      con = c_mat
      gamma = 0.0
      mat_var = var-1
      iseed = iseed - myid
!     maximum mc number for the matrix grains is mat_var
!     mc number equal to var is used for the precipitate phase
!     in this example we use 1 precipitate phase


      do j = 1, var
      do i = 1, var
        if (i.ne.j) then
           if(i.le.mat_var.and.j.le.mat_var) gamma(i,j) = sigma_1 
        end if
        if(j.eq.var) then
           if(i.lt.var) gamma(i,j) = sigma_2
        end if
        if(i.eq.var) then
           if(j.lt.var) gamma(i,j) = sigma_2
        end if
      end do 
      end do 

      if(myid.eq.0) then
        do i = 1,var
          write(97,*) i, (gamma(i,j),j=1,var)
        end do
      end if
      call flush(97)

!     set initial order parameters

      IF(nrun.eq.1) THEN

      ! Subroutine for initializing polycrystal with input file
      if(start.eq.'polycryst') then
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

      end if

      if(start.eq.'plane') then

      j1 = Ny - j_ppt + 1
      if(myid.eq.0) write(*,*) 'j1=', j1
      
      do kk = ist(3), ien(3)
      do jj = ist(2), ien(2)
      do ii = ist(1), ien(1)
       
         phi(2,ii,jj,kk) = 1.0

      end do  
      end do  
      end do  

      do kk = ist(3), ien(3)
      do jj = ist(2), ien(2)
      do ii = ist(1), ien(1)

      if(jj.lt.j1) then
         phi(:,ii,jj,kk) = 0.0
         phi(1,ii,jj,kk) = 1.d0
      end if

      end do  
      end do  
      end do  

      
      do kk = ist(3), ien(3)
      do jj = ist(2), ien(2)
      do ii = ist(1), ien(1)

      if(jj.ge.j_ppt) then
         phi(:,ii,jj,kk) = 0.0
         phi(1,ii,jj,kk) = 1.d0
      end if

      end do  
      end do  
      end do  

!     Introduce a precipitate phase at the center 
!     and assign concentrations

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
           phi(3,ii,jj,kk) = 1.0
           phi(1,ii,jj,kk) = 0.d0
           phi(2,ii,jj,kk) = 0.d0
           end if
         end if

      end do
      end do
      end do

      j_ppt = Ny - j_ppt +1

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
           phi(3,ii,jj,kk) = 1.0
           phi(1,ii,jj,kk) = 0.d0
           phi(2,ii,jj,kk) = 0.d0
           end if
         end if

      end do
      end do
      end do

      end if

      if(start.eq.'circle') then

         phi(1,:,:,:) = 1.0

!     Introduce a circular grain  at the center 

        do k = -grn_rad, grn_rad
        do j = -grn_rad, grn_rad
        do i = -grn_rad, grn_rad

           term_1 = (float(i)/float(grn_rad))**2
           term_2 = (float(j)/float(grn_rad))**2
           term_3 = (float(k)/float(grn_rad))**2
           sum_ijk = term_1 + term_2 + term_3

           if(sum_ijk.le.1.0) then
             kount = kount + 1
             ii = i + Nx/2
             jj = j + Ny/2
             kk = k + Nz/2
             if(ii.ge.ist(1).AND.ii.le.ien(1).AND. &
             jj.ge.ist(2).AND.jj.le.ien(2).AND. &
             kk.ge.ist(3).AND.kk.le.ien(3)) then
               phi(2,ii,jj,kk) = 1.d0
               phi(1,ii,jj,kk) = 0.d0
               phi(3,ii,jj,kk) = 0.d0
             end if
           end if

        end do
        end do
        end do

!     Introduce multiple precipitates within the circular grain

        nppt = 1
        rad_ppt = (ppt_rad(1)**2 + ppt_rad(2)**2 + ppt_rad(3)**2)*1.0

        do while(nppt.le.num_ppt)  

13        if(myid.eq.0) then

18         ix = ran_2(iseed)*Nx + 1.0
           jy = ran_2(iseed)*Ny + 1.0
           kz = ran_2(iseed)*Nz + 1.0

          end if

          call MPI_Bcast(ix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          call MPI_Bcast(jy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          call MPI_Bcast(kz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

          phi_count = 0
          do k  = ist(3), ien(3)
          do j  = ist(2), ien(2)
          do i  = ist(1), ien(1)

            radius = (i - ix)**2 + (j - jy)**2 + (k-kz)**2
            if(radius.le.rad_ppt.AND.phi(var,i,j,k).gt.0.0) &
            phi_count = phi_count +1

          end do
          end do
          end do

          call MPI_Allreduce(phi_count, phi_tot, 1, MPI_INTEGER, &
          MPI_SUM, MPI_COMM_WORLD, ierr)

          if(phi_tot.gt.0) go to 13

          if(phi_tot.eq.0) then

            do k  = ist(3), ien(3)
            do j  = ist(2), ien(2)
            do i  = ist(1), ien(1)

                radius = (i - ix)**2 + (j - jy)**2 + (k-kz)**2
                if(radius.le.rad_ppt) then
                  phi(3,i,j,k) = 1.0
                  phi(2,i,j,k) = 0.0
                  phi(1,i,j,k) = 0.0
                  con(i,j,k) = con_0_ppt
                end if

            end do
            end do
            end do
            nppt = nppt + 1
          end if
          if(myid.eq.0) write(*,*) 'nppt=', nppt
        end do 

      end if 

      ELSE

      open(5,file='data_rerun.'//file_num, status='old', &
      form='unformatted')

      read(5) ist, ien
      read(5) phi, con

      close(5) 

      END IF

      if(myid.eq.0) write(*,*) "precipitates assigned"

!     Set up k vectors in Fourier space
      call k_space(Nx, Ny, Nz, fst, fen, kf, kf_sq, kf_4)
      phi_min_glob = 1.0

      run_time = 0.0

! ========================================
! Start time loop
! ========================================
      do step = 1, N_step
      time(1) = MPI_Wtime()

!     Introduce multiple precipitates within the polycrystal 

      if(step.eq.ppt_int) then

        nppt = 1
        rad_ppt = (ppt_rad(1)**2 + ppt_rad(2)**2 + ppt_rad(3)**2)

        do while(nppt.le.num_ppt)

19        if(myid.eq.0) then

16         ix = ran_2(iseed)*Nx + 1.0
           jy = ran_2(iseed)*Ny + 1.0
           kz = ran_2(iseed)*Nz + 1.0

          end if

           call MPI_Bcast(ix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
           call MPI_Bcast(jy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
           call MPI_Bcast(kz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

          phi_sum = 0.0
          if(ix.ge.ist(1).and.ix.le.ien(1).AND. &
             jy.ge.ist(2).and.jy.le.ien(2).AND. &
             kz.ge.ist(3).and.kz.le.ien(3)) then
              do ivar = 1, mat_var
                fi = phi(ivar,ix,jy,kz)
                phi_sum =phi_sum + fi**2
              end do
          end if

          call MPI_Allreduce(phi_sum, glob_sum, 1, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ierr)

          !if(myid.eq.0) write(94,*) ix, jy, kz, glob_sum

          if(glob_sum.gt.0.8) go to 19

          phi_count = 0
          do k  = ist(3), ien(3)
          do j  = ist(2), ien(2)
          do i  = ist(1), ien(1)

            radius = (i - ix)**2 + (j - jy)**2 + (k-kz)**2
            if(radius.le.rad_ppt.AND.phi(var,i,j,k).gt.0.0) &
            phi_count = phi_count +1

          end do
          end do
          end do

          call MPI_Allreduce(phi_count, phi_tot, 1, MPI_INTEGER, &
          MPI_SUM, MPI_COMM_WORLD, ierr)

          if(phi_tot.gt.0) go to 19

          if(phi_tot.eq.0) then

            do k  = ist(3), ien(3)
            do j  = ist(2), ien(2)
            do i  = ist(1), ien(1)

                radius = (i - ix)**2 + (j - jy)**2 + (k-kz)**2
                if(radius.le.rad_ppt) then
                  phi(var,i,j,k) = 1.0
                  do ivar = 1, mat_var
                  phi(ivar,i,j,k) = 0.0
                  end do
                end if

            end do
            end do
            end do
          end if
          nppt = nppt + 1
          write(*,*) 'nppt=', nppt
        end do
      end if

      c_min = 10.0

!     Calculate the total system energy
!     First gradient energy due to phi gradients

      loop_time(1) = MPI_Wtime()
! ===============================
      !$omp parallel do collapse(4) &
      !$omp& private(i, j, k, ivar, dummy, dft_dummy) &
      !$omp& shared(kf, grad_x, grad_y, grad_z, dft_grad_x, dft_grad_y, dft_grad_z)
      DO ivar = 1, var

      dummy(:,:,:) = phi(ivar,:,:,:)
      call f_trans(dummy(ist(1),ist(2),ist(3)),  &
                  dft_dummy(fst(1),fst(2),fst(3)),  &
                  Nx, Ny, Nz, ist, ien, fst, fen)

      do k = fst(3), fen(3)
      do j = fst(2), fen(2)
      do i = fst(1), fen(1)

         dft_grad_x(ivar,i,j,k) = dft_dummy(i,j,k)*kf(1,i,j,k)
         dft_grad_y(ivar,i,j,k) = dft_dummy(i,j,k)*kf(2,i,j,k)
         dft_grad_z(ivar,i,j,k) = dft_dummy(i,j,k)*kf(3,i,j,k)

      end do
      end do
      end do

      dft_dummy(:,:,:) = dft_grad_x(ivar,:,:,:)
      call inv_trans(dft_dummy(fst(1),fst(2),fst(3)), &
                    dummy(ist(1),ist(2),ist(3)), ist, ien, fst, fen)
      grad_x(ivar,:,:,:) = dummy(:,:,:)


      dft_dummy(:,:,:) = dft_grad_y(ivar,:,:,:)
      call inv_trans(dft_dummy(fst(1),fst(2),fst(3)), &
                    dummy(ist(1),ist(2),ist(3)), ist, ien, fst, fen)
      grad_y(ivar,:,:,:) = dummy(:,:,:)


      dft_dummy(:,:,:) = dft_grad_z(ivar,:,:,:)
      call inv_trans(dft_dummy(fst(1),fst(2),fst(3)), &
                    dummy(ist(1),ist(2),ist(3)), ist, ien, fst, fen)
      grad_z(ivar,:,:,:) = dummy(:,:,:)

      END DO
      !$omp end parallel do
! ===============================
      loop_time(1) = MPI_Wtime() - loop_time(1)

      loop_time(2) = MPI_Wtime()
      f_int = 0.d0
      DO  ivar=1,var

! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k) &
      !$omp& shared(ivar, grad_coeff_phi, grad_x, grad_y, grad_z) &
      !$omp& reduction(+:f_int)
      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

      f_int = f_int + 0.5* grad_coeff_phi *  &
                     (grad_x(ivar,i,j,k)**2 + & 
                     grad_y(ivar,i,j,k)**2 + grad_z(ivar,i,j,k)**2)
      end do
      end do
      end do
      !$omp end parallel do
! ===============================

      END DO
      loop_time(2) = MPI_Wtime() - loop_time(2)

!     Next, gradient energy due to concentration gradients

      loop_time(3) = MPI_Wtime()
      call f_trans(con, dft_con, Nx, Ny, Nz, ist, ien, fst, fen)

! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k) &
      !$omp& shared(dft_grad_x_c, dft_grad_y_c, dft_grad_z_c, dft_con, kf) 
      do k = fst(3), fen(3)
      do j = fst(2), fen(2)
      do i = fst(1), fen(1)

         dft_grad_x_c(i,j,k) =  dft_con(i,j,k) * kf(1,i,j,k)
         dft_grad_y_c(i,j,k) =  dft_con(i,j,k) * kf(2,i,j,k)
         dft_grad_z_c(i,j,k) =  dft_con(i,j,k) * kf(3,i,j,k)

      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      loop_time(3) = MPI_Wtime() - loop_time(3)

      loop_time(4) = MPI_Wtime()
      call inv_trans(dft_grad_x_c, grad_x_c, ist, ien, fst, fen)
      call inv_trans(dft_grad_y_c, grad_y_c, ist, ien, fst, fen)
      call inv_trans(dft_grad_z_c, grad_z_c, ist, ien, fst, fen)

! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k) &
      !$omp& shared(grad_coeff_c, grad_x_c, grad_y_c, grad_z_c) &
      !$omp& reduction(+:f_int)
      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

      f_int = f_int + 0.5 * grad_coeff_c*(grad_x_c(i,j,k)**2 + &
                      grad_y_c(i,j,k)**2 + grad_z_c(i,j,k)**2)

      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      loop_time(4) = MPI_Wtime() - loop_time(4)

!     Total system chemical energy consisting of matrix grains and
!     precipitate

!     First matrix grains interfacial energy

      loop_time(5) = MPI_Wtime()
      f_eta = 0.d0

! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k, ivar, jvar, fi_i, fi_, gb_enj) &
      !$omp& shared(grad_coeff_c, grad_x_c, grad_y_c, grad_z_c) &
      !$omp& reduction(+:f_eta)
      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

      do ivar = 1, var

        fi_i = phi(ivar,i,j,k)
        f_eta = f_eta -0.5* fi_i**2 + 0.25* fi_i**4 

      do jvar = 1, var

        fi_j = phi(jvar,i,j,k)
        if(jvar.gt.ivar) then
!        if(jvar.ne.ivar) then

          gb_en = gamma(ivar,jvar)
          f_eta = f_eta + gb_en*fi_i**2*fi_j**2

        end if

      end do
      end do

      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      loop_time(5) = MPI_Wtime() - loop_time(5)

      loop_time(6) = MPI_Wtime()
      f_ch = 0.0

! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k, num_ijk, den_ijk, con_ijk, phi_mat, phi_ppt, fch_mat, fch_ppt) &
      !$omp& shared(df_dc) &
      !$omp& reduction(+:f_ch)
      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

        num_ijk = 0.0
        den_ijk = 0.0
        do ivar = 1, var
          if (ivar.lt.var) num_ijk = num_ijk + phi(ivar,i,j,k)**4 
          den_ijk = den_ijk + phi(ivar,i,j,k)**4
        end do
        con_ijk = con(i,j,k)

        phi_mat = num_ijk / den_ijk
        phi_ppt = 1.0 - phi_mat 

        fch_mat = phi_mat*am*(con_ijk-con_0_mat)**2
        fch_ppt = phi_ppt*ap*(con_ijk-con_0_ppt)**2

        f_ch = f_ch + fch_mat + fch_ppt 

        df_dc(i,j,k) = 2.d0*phi_mat*am*(con_ijk-con_0_mat) + &
                       2.d0*phi_ppt*ap*(con_ijk-con_0_ppt)

      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      loop_time(6) = MPI_Wtime() - loop_time(6)

!    Calculate the gradients of the total energy with respect to phi and
!    c required for solving the evolution equations

!    gradient wrt phi consists of four terms (Eqn 4 of paper by Clang
!    and Moelans)

      loop_time(7) = MPI_Wtime()
      term2 = 0.d0
      term1 = 0.0

! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i,j,k,ivar,jvar,fi_i,fi_j,gb_en) &
      !$omp& shared (term1, term2, phi)
      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

      do ivar = 1, var

        fi_i = phi(ivar,i,j,k)
        term1(ivar,i,j,k)  = fi_i**3 - fi_i
     
      do jvar = 1, var

       if(jvar.ne.ivar) then
!       if(jvar.gt.ivar) then

          fi_j = phi(jvar,i,j,k)

          gb_en = gamma(ivar,jvar)
          term2(ivar,i,j,k) = term2(ivar,i,j,k) + 2.d0*gb_en*fi_i*fi_j**2

       end if

      end do
      end do


      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      loop_time(7) = MPI_Wtime() - loop_time(7)


!     calculating term4 and term5

      loop_time(8) = MPI_Wtime()
! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i,j,k,ivar,con_ijk,den_ijk,fch_mat,fch_ppt) &
      !$omp& shared (term1, term2, phi)

      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

       den_ijk = 0.0
       con_ijk = con(i,j,k)

       do ivar = 1, var
         den_ijk = den_ijk + phi(ivar,i,j,k)**4
       end do

       fch_mat = am*(con_ijk-con_0_mat)**2
       fch_ppt = ap*(con_ijk-con_0_ppt)**2
       term_2 = 4.0*phi(var,i,j,k)**3*(den_ijk - phi(var,i,j,k)**4)

       do ivar = 1, var

       term_1 = 4.0*phi(ivar,i,j,k)**3*phi(var,i,j,k)**4 
       if(ivar.le.mat_var) term4(ivar,i,j,k) = term_1*(fch_mat - fch_ppt)/den_ijk**2
       if(ivar.gt.mat_var) term4(ivar,i,j,k) = term_2*(fch_ppt - fch_mat)/den_ijk**2

       end do

      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      loop_time(8) = MPI_Wtime() - loop_time(8)

      loop_time(9) = MPI_Wtime()
! ===============================
      !$omp parallel do collapse(4) &
      !$omp& private(i, j, k, ivar) &
      !$omp& shared (term1, term2, term4, df_dphi)
      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

      do ivar = 1,var

      df_dphi(ivar,i,j,k)= term1(ivar,i,j,k) + term2(ivar,i,j,k) + &
                           term4(ivar,i,j,k) 
      end do

      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      loop_time(9) = MPI_Wtime() - loop_time(9)

      loop_time(10) = MPI_Wtime()
      call f_trans(con, dft_con, Nx, Ny, Nz, ist, ien, fst, fen)
      call f_trans(df_dc, dft_df_dc, Nx, Ny, Nz, ist, ien, fst, fen)

!     Set up position-dependent diffusion coefficients

! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k, ivar, phi_sum, fi) &
      !$omp& shared(phi, mob_c, D_bulk, D_gb)
      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

!      phi_num = 0
!      do ivar = 1, var
!       fi = phi(ivar,i,j,k)
!       if(fi.gt.0.0) phi_num = phi_num + 1 
!      end do
!      if(phi_num.gt.1) then
!       mob_c(i,j,k) = D_gb
!      else
!       mob_c(i,j,k) = D_bulk
!      end if

       phi_sum = 0.0
       do ivar = 1, var
         fi = phi(ivar,i,j,k)
         phi_sum = phi_sum + fi*fi
       end do
       mob_c(i,j,k) = D_bulk(var) + 2.0*(1.0-phi_sum)*D_gb

      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      loop_time(10) = MPI_Wtime() - loop_time(10)

      !D_mean = (D_bulk + D_gb) / 2.0
      !D_mean = D_gb
      !D_mean = D_bulk

      !call MPI_Allreduce(ave_mob, glob_mob, 1, MPI_DOUBLE_PRECISION, &
      ! MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(sum(mob_c), glob_mob, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, ierr) 

      D_mean = D_L

      !D_mean = glob_mob/float(Nx*Ny*Nz) 
      !D_mean = D_L*4.0*D_gb

!     Set up additional terms in Fourier space due to variable mobility

      loop_time(11) = MPI_Wtime()
! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k, kf_sum, k_sq) &
      !$omp& shared(dft_dummy, dft_df_dc, grad_coeff_c, dft_con, kf, kf_sq)
      do k = fst(3), fen(3)
      do j = fst(2), fen(2)
      do i = fst(1), fen(1)

       kf_sum = kf(1,i,j,k) + kf(2,i,j,k) + kf(3,i,j,k)
       k_sq = kf_sq(1,i,j,k)+kf_sq(2,i,j,k)+kf_sq(3,i,j,k)
       k_cube = kf(1,i,j,k)*kf_sq(1,i,j,k) + kf(2,i,j,k)*kf_sq(2,i,j,k) &
        + kf(3,i,j,k)*kf_sq(3,i,j,k)
       dft_dummy(i,j,k) = kf_sum*dft_df_dc(i,j,k) - grad_coeff_c*dft_con(i,j,k)*k_cube

      end do
      end do
      end do
      !$omp end parallel do
! ===============================

      call inv_trans(dft_dummy, dummy, ist, ien, fst, fen)
      loop_time(11) = MPI_Wtime() - loop_time(11)

      loop_time(12) = MPI_Wtime()
! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k) &
      !$omp& shared(dummy, mob_c, D_mean)
      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)

      dummy(i,j,k) = dummy(i,j,k)*(mob_c(i,j,k)-D_mean)

      end do
      end do
      end do
      !$omp end parallel do
! ===============================

      call f_trans(dummy, dft_dummy, Nx, Ny, Nz, ist, ien, fst, fen)
      loop_time(12) = MPI_Wtime() - loop_time(12)

!     Solve C-H Equation in Fourier Space

      loop_time(13) = MPI_Wtime()
! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private (i, j, k, kf_sum, k_4, k_sq, term11, term22, term33) &
      !$omp& shared(kf, kf_4, kf_sq, dft_dummy, dft_con, D_mean, t_step, grad_coeff_c, dft_df_dc)
      do k = fst(3), fen(3)
      do j = fst(2), fen(2)
      do i = fst(1), fen(1)

        kf_sum = kf(1,i,j,k) + kf(2,i,j,k) + kf(3,i,j,k)
        !kf_sum = sqrt(kf(1,i,j,k)**2 + kf(2,i,j,k)**2 + kf(3,i,j,k)**2)
        k_4  = kf_4(1,i,j,k) + kf_4(2,i,j,k) + kf_4(3,i,j,k)
        k_sq = kf_sq(1,i,j,k)+kf_sq(2,i,j,k)+kf_sq(3,i,j,k)

        term11 = dft_dummy(i,j,k)*t_step*kf_sum
        term22 = 1.0 + D_mean*t_step*grad_coeff_c*k_4
        term33 = D_mean*k_sq*t_step*dft_df_dc(i,j,k)

        dft_con(i,j,k) = (dft_con(i,j,k) + term11 + term33) / term22

      end do
      end do
      end do
      !$omp end parallel do
! ===============================
      call inv_trans(dft_con, con, ist, ien, fst, fen)
      loop_time(13) = MPI_Wtime() - loop_time(13)

!    Solve time-dependent G-L equation

      loop_time(14) = MPI_Wtime()
! ===============================
      !$omp parallel do collapse(3) &
      !$omp& private(i, j, k, k_sq, kf_sum, term11, term22, dummy, dummy2, dft_dummy, dft_dummy2) &
      !$omp& shared(phi, df_dphi, kf_sq, kf, mob_phi, t_step, grad_coeff_phi)
      do ivar = 1, var

      dummy(:,:,:) = phi(ivar,:,:,:)
      call f_trans(dummy(ist(1),ist(2),ist(3)),  &
                  dft_dummy(fst(1),fst(2),fst(3)),  &
                  Nx, Ny, Nz, ist, ien, fst, fen)
      dummy2(:,:,:) = df_dphi(ivar,:,:,:)
      call f_trans(dummy2(ist(1),ist(2),ist(3)),  &
                  dft_dummy2(fst(1),fst(2),fst(3)),  &
                  Nx, Ny, Nz, ist, ien, fst, fen)

      do k = fst(3), fen(3)
      do j = fst(2), fen(2)
      do i = fst(1), fen(1)

        k_sq = kf_sq(1,i,j,k)+kf_sq(2,i,j,k)+kf_sq(3,i,j,k)
        !kf_sum = kf(1,i,j,k) + kf(2,i,j,k) + kf(3,i,j,k)
        !kf_sum = sqrt(kf(1,i,j,k)**2 + kf(2,i,j,k)**2 + kf(3,i,j,k)**2)

!       Implicit

        term11 = mob_phi*t_step*dft_dummy2(i,j,k)
        term22 = 1.0 - mob_phi*grad_coeff_phi*t_step*k_sq
        dft_dummy(i,j,k) = ( dft_dummy(i,j,k) - term11 ) / term22

!       Explicit

!       term11 = mob_phi*t_step*dft_dummy2(i,j,k)
!       term22 = mob_phi*grad_coeff_phi*t_step*k_sq
!       dft_dummy(i,j,k) = dft_dummy(i,j,k) -term11 - term22*dft_dummy(i,j,k)

      end do
      end do
      end do

!     perfrom inverse transformation to get the new phi values
!     This new phi which will be used in the Cahn-Hilliard equation
!     in the next step

      call inv_trans(dft_dummy(fst(1),fst(2),fst(3)), &
                    dummy(ist(1),ist(2),ist(3)), ist, ien, fst, fen)
      phi(ivar,:,:,:) = dummy(:,:,:)

      END DO
      !$omp end parallel do
! ===============================
      loop_time(14) = MPI_Wtime() - loop_time(14)

! get the run time for each step and the cumulative time

      time(1) = MPI_Wtime() - time(1)
      run_time = run_time + time(1)
      if(myid.eq.0) write(93,*) step, time(1), run_time

      if ( mod(step,2) .eq. 0 ) then
         call MPI_Reduce(f_int, f_int_tot, 1, MPI_DOUBLE_PRECISION, &
              MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_Reduce(f_ch, f_ch_tot, 1, MPI_DOUBLE_PRECISION, &
              MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_Reduce(f_eta, f_eta_tot, 1, MPI_DOUBLE_PRECISION, &
              MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_Reduce(beta_en, beta_en_tot, 1, MPI_DOUBLE_PRECISION, &
              MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         if ( myid .eq. 0 )  then
            !e_tot = f_int_tot + f_ch_tot + f_eta_tot + beta_en_tot
            e_tot = f_int_tot + f_ch_tot + f_eta_tot
            write(98,10) step, f_int_tot, f_ch_tot, f_eta_tot, e_tot, & 
                        c_max_glob, c_min_glob, phi_min_glob, phi_3_max
         end if
      end if
      call flush(98)

      phi_min = 1.0
      phi_3 = 0.0
      c_max = 0.0
      c_min = 1.0

! ===============================
      do k = ist(3), ien(3) 
      do j = ist(2), ien(2) 
      do i = ist(1), ien(1) 
         if(phi(var,i,j,k).gt.phi_3) phi_3 = phi(var,i,j,k)
         if(con(i,j,k).le.c_min) c_min = con(i,j,k)
         if(con(i,j,k).gt.c_max) c_max = con(i,j,k)
         phi_sum = 0.0
         phi_sq = 0.0
         do ivar = 1, var
          phi_sum = phi_sum + phi(ivar,i,j,k)
          phi_sq = phi_sq + phi(ivar,i,j,k)**2
         end do
         if(phi_sq.lt.phi_min) phi_min=phi_sq
!         if(phi_sum.gt.1.0) then
!           do ivar=1,var
!             phi(ivar,i,j,k) = phi(ivar,i,j,k)/phi_sum
!           end do
!           if(con(i,j,k).lt.1.e-7) con(i,j,k)=1.e-7
!           if(con(i,j,k).gt.1.0) con(i,j,k)=1.0
!         end if
           
      end do
      end do
      end do
! ===============================

      call MPI_Reduce(c_min, c_min_glob, 1, MPI_DOUBLE_PRECISION, &
              MPI_MIN, 0, MPI_COMM_WORLD, ierr)

      call MPI_Reduce(phi_min, phi_min_glob, 1, MPI_DOUBLE_PRECISION, &
              MPI_MIN, 0, MPI_COMM_WORLD, ierr)

      call MPI_Reduce(c_max, c_max_glob, 1, MPI_DOUBLE_PRECISION, &
              MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      call MPI_Reduce(phi_3, phi_3_max, 1, MPI_DOUBLE_PRECISION, &
              MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      if(mod(step,ifreq).eq.0) then

        write(dir_name,'(A, a5,i7.7)') trim(output_dir),'step_',step
        !write(dir_name,'(a5,i7.7)') 'step_',step

        if(myid.eq.0) then
        call system('mkdir -p '//trim(dir_name))
          open(9,file=trim(dir_name)//'/data_extract_mult.in')
          write(9,*) Nx, Ny, Nz, nprocs, var
          write(9,*) dims
          close(9)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        open (4, file=trim(dir_name)//'/data_save.'//file_num, &
              form='unformatted')
        rewind(4)
        write(4) ist, ien
        write(4) phi, con
        call flush(4)
        close (4)

      end if

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      ! Print average loop durations across processes 
      call MPI_Allreduce(loop_time, avg_loop_time, 14, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      avg_loop_time = avg_loop_time / real(nprocs)
      if (myid.eq.0) &
        write(95, '(I8, 16(A, F8.6))') step, &
              ' ', avg_loop_time(1), ' ', avg_loop_time(2), ' ', avg_loop_time(3), ' ', avg_loop_time(4), ' ', avg_loop_time(5), ' ', avg_loop_time(6), ' ', avg_loop_time(7), &
              ' ', avg_loop_time(8), ' ', avg_loop_time(9), ' ', avg_loop_time(10), ' ', avg_loop_time(11), ' ', avg_loop_time(12), ' ', avg_loop_time(13), ' ', avg_loop_time(14), &
              ' ', time(1), ' ', run_time


      END DO

10    format(i6, 8e15.7)
9     format(3i6, 4e15.7)


      call p3dfft_clean

      call MPI_Finalize(ierr)

      end

