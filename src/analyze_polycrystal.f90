      program paraview

      implicit none
      
      include 'mpif.h'
      ! --
      ! Index integers and MPI Variables 
      integer Nx, Ny, Nz, Nxyz, nprocs, dims(2), myid, ist(3), ien(3), Nxy
      integer i, j, k, kk, var, ivar, step, var_max, ii, ll
      integer ix, jy, kz
      integer ierr
      ! --
      ! GL arrays for saving output 
      integer, dimension(:), allocatable :: var_count
      integer, dimension(:), allocatable :: num_var_gl
      real *8, dimension(:,:,:,:), allocatable :: phi
      real *8, dimension(:,:,:), allocatable :: con
      real *8, dimension(:), allocatable ::  con_gl
      real *8, dimension(:), allocatable ::  phi_sq_gl
      real *8, dimension(:), allocatable ::  phi_var_gl, grad_phi_var_gl
      real *8, dimension(:), allocatable :: boundary_gl 
      character(len=256), dimension(:), allocatable :: precipitate_gl 
      ! --
      ! Calculation of precipitate coverage of boundary 
      real *8 threshold, phi_above_threshold, precipitate_fraction
      integer boundary_voxels, precipitate_boundary_voxels
      integer :: my_start_step, my_end_step
      ! --
      ! Reading and calculating from input file
      real *8 term, phi_sq, phi_sum, phi_tot, phi_max, phi_temp
      character(len=256) :: input_dir, output_dir, step_dir
      character(len=256) :: infile, outfile_phi_sq, outfile_con, outfile_boundary
      character(len=256) :: outfile_var_num, outfile_phi_var, outfile_grad_phi_var
      character(len=256) :: command
      character(len=64), dimension(:), allocatable :: step_strings
      character(len=64) :: step_line, threshold_str
      character(len=32) :: file_num_str
      character(len=7) :: step_string
      integer :: num_steps, unit, ios, step_i, file_num
      ! --
      ! Calculation of mean thickness estimate
      integer :: intersections, intersections_glob, iseed, phi_var_sum, phi_var_sum_glob, N, total_length
      real *8 :: N_L, ran_2, thickness_est, phi_var_vol_frac

      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
      
!! DO THIS WITH MYID==0

      ! Process 0 reads the command arguments and broadcasts them
      if (myid == 0) then
        call get_command_argument(1, input_dir, status=ierr)
        call get_command_argument(2, output_dir, status=ierr)
        call get_command_argument(3, threshold_str, status=ierr)
        read(threshold_str, *) threshold
        print *, "threshold = ", threshold
      endif

      ! Broadcast threshold and directory names to all processes
      call MPI_Bcast(threshold, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(input_dir, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(output_dir, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      ! Process 0 determines steps and broadcasts to step_strings
      if (myid == 0) then
        ! Determine step directories
        command = 'ls ' // trim(input_dir) // &
                    ' | grep "^step_" | cut -d"_" -f2 > steps.txt'
        call execute_command_line(command, wait=.true., exitstat=ierr)
  
        if (ierr /= 0) then
              print *, "Error executing command to list steps"
              stop
        end if
  
        ! Open the file with step numbers
        open(unit=99, file='steps.txt', status='old', action='read')
  
        num_steps = 0
        do
              read(99, '(A)', iostat=ios) step_line
              if (ios /= 0) exit
              num_steps = num_steps + 1
        end do

        allocate(step_strings(num_steps))
        allocate(precipitate_gl(num_steps))
  
        rewind(99)
        do i = 1, num_steps
              read(99, '(A)', iostat=ios) step_line
              if (ios /=0) then
                    print *, "Error reading steps.txt at line ", i
                    stop
              end if
  
              ! Extract the numerical part of the step directory anme
              step_line = adjustl(trim(step_line))   ! Remove leading and trailing whitespace
              
              step_strings(i) = step_line
        end do
  
        close(99)
        !call execute_command_line('rm steps.txt')
  
        open(93, file=trim(output_dir) // '/thickness_calculation.dat', status='unknown')
        write(93,'(A)') "# step, avg thickness (estimated)"
        write(*,'(A)') "# step, precipitate_voxels, &
                     boundary_voxels, precipitate_fraction"
      end if

      ! Broadcast number of steps and step_strings to all processes
      call MPI_Bcast(num_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (myid /= 0) then
        allocate(step_strings(num_steps))
        allocate(precipitate_gl(num_steps))
        !open(96, file=trim(output_dir) // '/precipitate_fraction.dat', status='unknown', action='write', position='append')
      end if
      call MPI_Bcast(step_strings, num_steps*64, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      ! Determine range of steps each process will handle
      my_start_step = (num_steps / nprocs) * myid + 1
      if (myid == nprocs - 1) then
        my_end_step = num_steps
      else
        my_end_step = (num_steps / nprocs) * (myid + 1)
      end if

!! DO THIS LOOP WITH MPI DISTRIBUTED PROCESSES
! There will be approximately 20-50 num_steps
! There will be approximately 16 data_save.xxxx files

      ! Loop through step directories assigned to this process
      do step_i = my_start_step, my_end_step 
            step_string = trim(step_strings(step_i))
            print *, trim(step_string)
            step_dir = trim(input_dir) // '/step_' & 
                  //trim(adjustl(step_string)) // '/'
            infile = trim(step_dir) // 'data_extract_mult.in'

            !print *, "Analyzing directory: ", step_dir

            !open(10, file=trim(step_dir)//'var_dist', status='unknown')
            open(11, file=trim(infile), status='old')
            read(11,*) Nx, Ny, Nz, nprocs, var
            read(11,*) dims(1), dims(2)
            close(11)
      
            Nxyz = Nx*Ny*Nz
            Nxy = Nx*Ny

            allocate (phi_sq_gl(1:Nxyz))
            allocate (phi_var_gl(1:Nxyz))
            allocate (grad_phi_var_gl(1:Nxyz))
            allocate (con_gl(1:Nxyz))
            allocate (num_var_gl(1:Nxyz))
            allocate (var_count(var))
            allocate (boundary_gl(1:Nxyz))
      
            phi_tot = 0.0
            var_count = 0

            boundary_voxels = 0
            precipitate_boundary_voxels = 0

            phi_var_sum = 0
            file_num = 0

            ! Loop through data_save.xxxxx files according to dims
            do kk = 1, dims(2)
              do j = 1, dims(1)
                file_num = file_num + 1
                write(file_num_str, '(i4.4)') file_num
                !write(*,'(A, I4, 4(A))') "myid, step_dir, file_num: ", myid, " ", step_dir, " ", file_num_str ! Debug MPI
                open(12,file=trim(step_dir)//'data_save.'//file_num_str, & 
                    status='old', form='unformatted')
                read(12) ist, ien
                if (j == 1) allocate( phi(var, Nx, Ny, ist(3):ien(3)) )
                if (j == 1) allocate( con(Nx, Ny, ist(3):ien(3)) )
                read(12) phi(1:var,1:Nx,ist(2):ien(2),ist(3):ien(3)),  &
                         con(1:Nx,ist(2):ien(2),ist(3):ien(3))
                close(12)

              end do
              
              ! ===============================
              !$omp parallel do collapse(3) &
              !$omp& private(i, j, j, ll, phi_sq, phi_max, var_max, phi_above_threshold) & 
              !$omp& shared(con_gl, con, phi_var_gl, phi, boundary_gl) &
              !$omp& reduction(+:boundary_voxels, +:phi_var_sum, +:precipitate_boundary_voxels)
  
              do k = ist(3), ien(3)
              do j = 1, Ny
              do i = 1, Nx 
  
                ll = (k-1) * Nxy + (j-1) * Nx + i
                con_gl(ll) = con(i,j,k)
                phi_sq = 0.0
                phi_max = 0.8
                var_max = 0
   
                phi_above_threshold = 0 
   
                do ii = 1, var ! loop through all precipitates 
                  if (ii == var) phi_var_gl(ll) = phi(var,i,j,k) ! var is the "precipitate phase"?
                  !write(*, '(A, 4I8, 4e15.7)') "ivar, i, j, k, phi", ii, i, j, k, phi(ivar,i,j,k)
  
                  phi_sq = phi_sq + phi(ii,i,j,k)*phi(ii,i,j,k)
                  if (phi(ii,i,j,k) > phi_max) var_max = ii
                  if (phi(ii,i,j,k) > threshold) then
                    phi_above_threshold = phi_above_threshold + 1
                    !print *, "threshold, phi above threshold: ", phi_above_threshold

                    ! It is confusing here why the last precipitate is
                    ! chosen to be the important one, I am unsure if
                    ! it is correct to say that this means it is "the"
                    ! precipitate
                  end if
                end do

                ! To calculate volume fraction
                if (phi(var,i,j,k) > threshold) &
                  phi_var_sum = phi_var_sum + 1
   
                !write(95,'(3(I6), 4(F6.2))') i, j, k, phi(1,i,j,k), phi(var,i,j,k), phi(3,i,j,k), phi_sq
   
                phi_sq_gl(ll) = phi_sq
                num_var_gl(ll) = var_max
   
                ! if more than one phi is nonzero, that signifies a boundary
                if (phi_above_threshold > 1) then 
                  boundary_gl(ll) = 1.0 !phi_above_threshold
                  boundary_voxels = boundary_voxels + 1

                  ! precipitate voxels within boundary region
                  ! phi(var) is above threshold means precipitate present
                  if (phi(var,i,j,k) > threshold) precipitate_boundary_voxels = precipitate_boundary_voxels + 1 
                else
                  boundary_gl(ll) = 0.0
                end if
  
              end do
              end do
              end do
  
              !$omp end parallel do
              ! ===============================

              deallocate (phi)
              deallocate (con)
  
              !write(*,'(A, I4, 3(A), 2(I8), F6.2)') "myid file_num prec, bound, frac ", myid, " ", file_num_str, " ", precipitate_boundary_voxels, boundary_voxels, real(precipitate_boundary_voxels, kind=8)/ real(boundary_voxels, kind=8) 
            end do
            if (myid == 0) print *, "=========== FINISHED GL ANALYSIS ==========="

            ! ====================
            ! Statistical estimate of mean thickness of second phase
            ! ====================
            ! Calculate volume fraction of second phase
            iseed = -1
            N = 16 

            file_num = myid * dims(1) ! TODO requires nproc==dims(2)
            do j = 1, dims(1)
              file_num = file_num + 1
              write(file_num_str, '(i4.4)') file_num
              open(12,file=trim(step_dir)//'data_save.'//file_num_str, & 
                  status='old', form='unformatted')
              read(12) ist, ien
              if (j == 1) allocate( phi(var, Nx, Ny, ist(3):ien(3)) )
              if (j == 1) allocate( con(Nx, Ny, ist(3):ien(3)) )
              read(12) phi(1:var,1:Nx,ist(2):ien(2),ist(3):ien(3)),  &
                       con(1:Nx,ist(2):ien(2),ist(3):ien(3))
              close(12)

            end do

            intersections = 0

            ! N random rays in x,z 
            do ii = 1, N
              ix = ran_2(iseed)*Nx + 1.0
              kz = ran_2(iseed)*(ien(3) - ist(3)) + ist(3) 
              ! Calculate intersections of rays
              do jy = 1, Ny
                !write(*,'(A,6I8,2F8.4)') "myid, ii, ix, jy, kz, intersections, phi_var, threshold: ", myid, ii, ix, jy, kz, intersections, phi(var,ix,jy,kz), threshold
                if (phi(var,ix,jy,kz) >= threshold) &
                  intersections = intersections + 1
              end do
            end do
            ! N random rays in y,z 
            do ii = 1, N
              jy = ran_2(iseed)*Nx + 1.0
              kz = ran_2(iseed)*(ien(3) - ist(3)) + ist(3) 
              ! Calculate intersections of rays
              do ix = 1, Nx
                !write(*,'(A,6I8,2F8.4)') "myid, ii, ix, jy, kz, intersections, phi_var, threshold: ", myid, ii, ix, jy, kz, intersections, phi(var,ix,jy,kz), threshold
                if (phi(var,ix,jy,kz) >= threshold) &
                  intersections = intersections + 1
              end do
            end do
            ! N random rays in x,y 
            do ii = 1, N*nprocs
              ix = ran_2(iseed)*Nx + 1.0
              jy = ran_2(iseed)*Ny + 1.0
              ! Calculate intersections of rays
              do kz = ist(3), ien(3) 
                !write(*,'(A,6I8,2F8.4)') "myid, ii, ix, jy, kz, intersections, phi_var, threshold: ", myid, ii, ix, jy, kz, intersections, phi(var,ix,jy,kz), threshold
                if (phi(var,ix,jy,kz) >= threshold) &
                  intersections = intersections + 1
              end do
            end do

            deallocate (phi)
            deallocate (con)

            call MPI_Allreduce(phi_var_sum, phi_var_sum_glob, 1, MPI_INTEGER, &
                MPI_SUM, MPI_COMM_WORLD, ierr)
            phi_var_vol_frac = real(phi_var_sum_glob, kind=8)/ real(Nx * Ny * Nz, kind=8)

            ! Average intersections per unit length
            call MPI_Allreduce(intersections, intersections_glob, 1, MPI_INTEGER, &
                MPI_SUM, MPI_COMM_WORLD, ierr)
            total_length = nprocs * ( (N * (Nx + Ny)) + (nprocs * N * (ien(3) - ist(3))) )
            N_L = real(intersections_glob, kind=8)/ real(total_length, kind=8)

            ! Calculate average mean thickness prediction
            thickness_est = 4*phi_var_vol_frac/ 2*N_L
            !write(93,'(2A,I8,3F8.4,I8)') "step, myid, thickness_est, phi_var_vol_frac, N_L, intersections: ", trim(step_string), myid, thickness_est, phi_var_vol_frac, N_L, intersections

            ! ====================

            ! -----------
            ! Write outputs
            ! -----------

            ! Write mean thickness estimate
            write(93,*) trim(step_string), thickness_est
            !write(93,'(I8, F8.4)') step, mean_thickness_est
 
            ! Write VTK files from each MPI Process
            outfile_phi_sq = trim(output_dir) // '/phi_sq_' // trim(adjustl(step_string)) // '.vtk'
            outfile_con = trim(output_dir) // '/con_' // trim(adjustl(step_string)) // '.vtk'
            outfile_var_num = trim(output_dir) // '/var_num_' // trim(adjustl(step_string)) // '.vtk'
            outfile_phi_var = trim(output_dir) // '/phi_var_' // trim(adjustl(step_string)) // '.vtk'
            outfile_grad_phi_var = trim(output_dir) // '/grad_phi_var_' // trim(adjustl(step_string)) // '.vtk'
            outfile_boundary = trim(output_dir) // '/boundary_' // trim(adjustl(step_string)) // '.vtk'

            open(21, file=trim(outfile_phi_sq), status='unknown')
            open(22, file=trim(outfile_con), status='unknown')
            open(23, file=trim(outfile_var_num), status='unknown')
            open(24, file=trim(outfile_phi_var), status='unknown')
            open(25, file=trim(outfile_grad_phi_var), status='unknown')
            open(26, file=trim(outfile_boundary), status='unknown')
      
            ! Write header
            do ii = 21,26
              write(ii,'(a26)') "# vtk DataFile Version 3.0"
              write(ii,'(a11)') "Sample Data"
              write(ii,'(a5)') "ASCII"
              write(ii,'(a25)') "DATASET STRUCTURED_POINTS"
              write(ii,'(a10,2x,3i6)') "DIMENSIONS",nx,ny,nz
              write(ii,'(a24)') "ASPECT_RATIO 1.0 1.0 1.0"
              write(ii,'(a18)') "ORIGIN 0.0 0.0 0.0"
              write(ii,'(a10,2x,i12)') "POINT_DATA",Nxyz
              write(ii,'(a22)') "SCALARS Field double 1"
              write(ii,'(a20)') "LOOKUP_TABLE default"
            end do

            ! Calculate gradient of phi field
            call calculate_phi_grad(Nx, Ny, Nz, phi_var_gl, grad_phi_var_gl)
      
            ! Write outputs
            write(21,*) phi_sq_gl
            write(22,*) con_gl
            write(23,*) num_var_gl
            write(24,*) phi_var_gl
            write(25,*) grad_phi_var_gl
            write(26,*) boundary_gl
      
            close(21)
            close(22)
            close(23)
            close(24)
            close(25)
            close(26)
  
            ! Each process needs to write the following:
            if (boundary_voxels > 0) then
              precipitate_fraction = real(precipitate_boundary_voxels, kind=8)/ real(boundary_voxels, kind=8)
            else
              precipitate_fraction = 0
            end if
            !write(precipitate_gl, '(A, 2(I8), F6.2)') trim(step_string), precipitate_boundary_voxels, boundary_voxels, precipitate_fraction
            !write(96, '(A, 2(I8), F6.2)') trim(step_string), precipitate_boundary_voxels, boundary_voxels, precipitate_fraction
            write(*, '(A, 2(I8), F6.2)') trim(step_string), precipitate_boundary_voxels, boundary_voxels, precipitate_fraction

            deallocate(phi_sq_gl)
            deallocate(phi_var_gl)
            deallocate(grad_phi_var_gl)
            deallocate(con_gl)
            deallocate(num_var_gl)
            deallocate(var_count)
            deallocate(boundary_gl)

      end do
      deallocate(step_strings)

      call MPI_Finalize(ierr)
     
      contains

      ! Subroutine for calculation of gradients of phi
      subroutine calculate_phi_grad(Nx, Ny, Nz, phi, grad_phi_mag)
            implicit none

            integer, intent(in) :: Nx, Ny, Nz
            real(8), intent(in) :: phi(:)
            real(8), dimension(:,:,:), allocatable :: phi_var, grad_phi_var_mag
            real(8), dimension(:,:,:), allocatable :: grad_x, grad_y, grad_z
            real(8), intent(out) :: grad_phi_mag(:)
            integer :: i, j, k

            allocate(grad_x(Nx, Ny, Nz))
            allocate(grad_y(Nx, Ny, Nz))
            allocate(grad_z(Nx, Ny, Nz))

            allocate(phi_var(Nx, Ny, Nz))
            allocate(grad_phi_var_mag(Nx, Ny, Nz))

            ! Reshape phi array
            phi_var = reshape(phi, [Nx, Ny, Nz])

            ! Initialize gradients to zero
            grad_x = 0.0
            grad_y = 0.0
            grad_z = 0.0

            ! Compute gradients using central differences
            do k = 2, Nz-1
            do j = 2, Ny-1
            do i = 2, Nx-1
              grad_x(i,j,k) = (phi_var(i+1,j,k) - phi_var(i-1,j,k)) / 2.0
              grad_y(i,j,k) = (phi_var(i,j+1,k) - phi_var(i,j-1,k)) / 2.0
              grad_z(i,j,k) = (phi_var(i,j,k+1) - phi_var(i,j,k-1)) / 2.0 

              ! Calculate grad magnitude array
              grad_phi_var_mag(i,j,k) = sqrt(grad_x(i,j,k)**2 + grad_y(i,j,k)**2 + grad_z(i,j,k)**2)
            end do
            end do
            end do

            ! Reshape grad magnitude array
            grad_phi_mag = reshape(grad_phi_var_mag, [Nx*Ny*Nz])
            
            deallocate(phi_var)
            deallocate(grad_phi_var_mag)

            deallocate(grad_x)
            deallocate(grad_y)
            deallocate(grad_z)

      end subroutine calculate_phi_grad

      end program paraview

