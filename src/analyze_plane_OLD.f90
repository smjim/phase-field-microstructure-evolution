      ! Determine max diameter of ppt introduced on boundaries of plane ic
      program analyze_plane

      use mpi

      implicit none

      ! --
      ! Index integers and MPI Variables 
      integer Nx, Ny, Nz, Nxyz, nprocs, dims(2), myid, ist(3), ien(3), Nxy
      integer i, j, k, kk, var, ivar, step, ii
      integer ierr
      ! --
      ! Arrays for input data 
      real *8, dimension(:,:,:,:), allocatable :: phi
      real *8, dimension(:,:,:), allocatable :: con
      ! --
      ! Calculation of precipitate coverage of boundary 
      real *8 threshold, phi_above_threshold, precipitate_fraction
      integer boundary_voxels, precipitate_boundary_voxels
      integer :: my_start_step, my_end_step
      ! --
      ! Reading and calculating from input file
      character(len=256) :: input_dir, output_dir, step_dir
      character(len=256) :: infile
      character(len=256) :: command
      character(len=64), dimension(:), allocatable :: step_strings
      character(len=64) :: step_line, threshold_str
      character(len=32) :: file_num_str
      character(len=7) :: step_string
      integer :: num_steps, unit, ios, step_i, file_num
      ! --
      ! Calculating precipitate coverage 
      integer, dimension(:), allocatable :: displs, recv_counts
      integer, dimension(:,:), allocatable :: boundary_voxels_list, all_boundary_voxels
      integer :: p, dx, dy, dz, max_boundary_voxels, phi_contained, j_ppt, j1
      real *8 :: distance, max_distance

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
        call execute_command_line('rm steps.txt')
      end if

      ! Broadcast number of steps and step_strings to all processes
      call MPI_Bcast(num_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (myid /= 0) allocate(step_strings(num_steps))
      call MPI_Bcast(step_strings, num_steps*64, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      ! Determine range of steps each process will handle
      my_start_step = (num_steps / nprocs) * myid + 1
      if (myid == nprocs - 1) then
        my_end_step = num_steps
      else
        my_end_step = (num_steps / nprocs) * (myid + 1)
      end if

!! DO THIS LOOP WITH MPI DISTRIBUTED PROCESSES

      ! Loop through step directories assigned to this process
      do step_i = my_start_step, my_end_step 
            step_string = trim(step_strings(step_i))
            step_dir = trim(input_dir) // '/step_' & 
                  //trim(adjustl(step_string)) // '/'
            infile = trim(step_dir) // 'data_extract_mult.in'

            !print *, "Analyzing directory: ", step_dir

            open(10, file=trim(step_dir)//'var_dist', status='unknown')
            open(11, file=trim(infile), status='old')
            read(11,*) Nx, Ny, Nz, nprocs, var
            read(11,*) dims(1), dims(2)
            close(11)
      
            Nxyz = Nx*Ny*Nz
            Nxy = Nx*Ny

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
              j_ppt = 70
              j1 = Ny - j_ppt + 1

              ! ====================
              ! Analyze the first plane
              j = j1
              ii = 0
              
              do k = ist(3), ien(3)
              do i = 1, Nx 

                ! Determine precipitate boundary
                phi_contained = 0
                do ivar = 1, var
                  if (0.1 <= phi(ivar,i,j,k) .and. phi(ivar,i,j,k) <= 0.11) &
                    phi_contained = phi_contained + 1
                    end do

                if (phi_contained.eq.3) then
                  ! Determine number of triple boundary voxels  
                  ii = ii + 1
                end if

              end do
              end do

              allocate(boundary_voxels_list(ii,3))

              ii = 0
              do k = ist(3), ien(3)
              do i = 1, Nx 

                ! Determine precipitate boundary
                phi_contained = 0
                do ivar = 1, var
                  if (0.1 <= phi(ivar,i,j,k) .and. phi(ivar,i,j,k) <= 0.11) &
                    phi_contained = phi_contained + 1
                end do

                if (phi_contained.eq.3) then
                  ! Save ppt boundary voxel indices to shared data
                  boundary_voxels_list(ii,:) = (/ i, j, k /)
                  ii = ii + 1
                end if

              end do
              end do

              ! Share memory across mpi processes to process 0
              allocate(recv_counts(nprocs))
              call MPI_Gather(ii, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

              ! Calculate displacements for data gathering
              allocate(displs(nprocs))
              displs(1) = 0
              do p = 2, nprocs
                displs(p) = displs(p-1) + recv_counts(p-1)
              end do

              ! Allocate memory to gather all boundary voxels to process 0
              if (myid == 0) then
                max_boundary_voxels = maxval(recv_counts)
                allocate(all_boundary_voxels(max_boundary_voxels,3))
              end if

              ! Gather all boundary voxels to process 0
              call MPI_Gatherv(boundary_voxels, ii, MPI_INTEGER, &
                     all_boundary_voxels, recv_counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

              if (myid.eq.0) then
                ! Find max diameter of precipitate (farthest apart)
                max_distance = 0.d0

                ! Loop through all precipitate boundary voxels
                do i = 1, max_boundary_voxels
                  do j = 1, max_boundary_voxels
                    ! Calculate distance between boundary voxels i and j
                    dx = all_boundary_voxels(j, 1) - all_boundary_voxels(i, 1)
                    dy = all_boundary_voxels(j, 2) - all_boundary_voxels(i, 2)
                    dz = all_boundary_voxels(j, 3) - all_boundary_voxels(i, 3)
      
                    distance = (dx**2 + dy**2 + dz**2)**0.5

                    if (distance > max_distance) &
                      max_distance = distance 

                  end do
                end do

                print *, "step, max_distance plane 1: ", step_i, max_distance

                deallocate(all_boundary_voxels)
              end if

              deallocate(recv_counts)
              deallocate(displs)
              deallocate(boundary_voxels_list)
              ! ====================

              ! ====================
              ! Analyze the second plane
              j = j_ppt
              ii = 0
              
              do k = ist(3), ien(3)
              do i = 1, Nx 

                ! Determine precipitate boundary
                phi_contained = 0
                do ivar = 1, var
                  if (0.1 <= phi(ivar,i,j,k) .and. phi(ivar,i,j,k) <= 0.11) &
                    phi_contained = phi_contained + 1
                end do

                if (phi_contained.eq.3) then
                  ! Determine number of triple boundary voxels  
                  ii = ii + 1
                end if

              end do
              end do

              allocate(boundary_voxels_list(ii,3))

              ii = 0
              do k = ist(3), ien(3)
              do i = 1, Nx 

                ! Determine precipitate boundary
                phi_contained = 0
                do ivar = 1, var
                  if (0.1 <= phi(ivar,i,j,k) .and. phi(ivar,i,j,k) <= 0.11) &
                    phi_contained = phi_contained + 1
                end do

                if (phi_contained.eq.3) then
                  ! Save ppt boundary voxel indices to shared data
                  boundary_voxels_list(ii,:) = (/ i, j, k /)
                  ii = ii + 1
                end if

              end do
              end do

              ! Share memory across mpi processes to process 0
              allocate(recv_counts(nprocs))
              call MPI_Gather(ii, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

              ! Calculate displacements for data gathering
              allocate(displs(nprocs))
              displs(1) = 0
              do p = 2, nprocs
                displs(p) = displs(p-1) + recv_counts(p-1)
              end do

              ! Allocate memory to gather all boundary voxels to process 0
              if (myid == 0) then
                max_boundary_voxels = maxval(recv_counts)
                print *, "step, total boundary voxels: ", step_i, max_boundary_voxels
                allocate(all_boundary_voxels(max_boundary_voxels,3))
              end if

              ! Gather all boundary voxels to process 0
              call MPI_Gatherv(boundary_voxels, ii, MPI_INTEGER, &
                     all_boundary_voxels, recv_counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

              if (myid.eq.0) then
                ! Find max diameter of precipitate (farthest apart)
                max_distance = 0.d0
                print *, "myid, step, max_boundary_voxels", myid, step_i, max_boundary_voxels

                ! Loop through all precipitate boundary voxels
                do i = 1, max_boundary_voxels
                  do j = 1, max_boundary_voxels
                    ! Calculate distance between boundary voxels i and j
                    dx = all_boundary_voxels(j, 1) - all_boundary_voxels(i, 1)
                    dy = all_boundary_voxels(j, 2) - all_boundary_voxels(i, 2)
                    dz = all_boundary_voxels(j, 3) - all_boundary_voxels(i, 3)
      
                    distance = (dx**2 + dy**2 + dz**2)**0.5

                    if (distance > max_distance) &
                      max_distance = distance

                  end do
                end do

                print *, "step, max_distance plane 2: ", step_i, max_distance

                deallocate(all_boundary_voxels)
              end if

              deallocate(recv_counts)
              deallocate(displs)
              deallocate(boundary_voxels_list)
              ! ====================

              deallocate(phi)
              deallocate(con)

            end do
      end do

      deallocate(step_strings)

      end program analyze_plane

