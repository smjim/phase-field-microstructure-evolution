      ! Determine max diameter of ppt introduced on boundaries of plane ic
      program analyze_plane

      implicit none
      integer Nx, Ny, Nz, nxyz, nprocs, dims(2), myid, ist(3), ien(3), nxy
      integer i, j, k, kk, var, ivar, step, var_max, ii, ll
      integer ix, jy, kz
      integer, dimension(:), allocatable :: var_count
      integer, dimension(:), allocatable :: num_var_gl
      real *8, dimension(:,:,:,:), allocatable :: phi
      real *8, dimension(:,:,:), allocatable :: con
      real *8, dimension(:), allocatable ::  phi_sq_gl
      real *8, dimension(:), allocatable ::  phi_3_gl
      real *8, dimension(:), allocatable ::  con_gl
      real *8 term, phi_sq, phi_sum, phi_tot, phi_max, phi_temp
      ! --
      ! Reading and calculating from input file
      character(len=256) :: input_dir, output_dir, step_dir
      character(len=256) :: infile
      character(len=256) :: command
      character(len=64), dimension(:), allocatable :: step_strings
      character(len=64) :: step_line, threshold_str
      real *8 :: threshold
      character(len=32) :: file_num_str
      character(len=7) :: step_string
      integer :: num_steps, unit, ios, step_i, file_num
      ! --
      ! Calculating precipitate coverage
      integer, dimension(:,:,:), allocatable :: boundary_locations
      integer :: dx, dy, dz, phi_contained, j_ppt, j1, ierr, p
      real *8 :: distance, max_distance

      ! Read the command arguments
      call get_command_argument(1, input_dir, status=ierr)
      call get_command_argument(2, output_dir, status=ierr)
      call get_command_argument(3, threshold_str, status=ierr)
      read(threshold_str, *) threshold
      print *, "threshold = ", threshold

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

      ! Write header to output file
      open(91, file=trim(output_dir) // "ppt_wetting_p1.dat", status='unknown')
      write(91, *) "# step, ppt max diameter, num_boundary_locations"
      open(92, file=trim(output_dir) // "ppt_wetting_p2.dat", status='unknown')
      write(92, *) "# step, ppt max diameter, num_boundary_locations"

      ! Loop through all step directories
      do step_i = 1, num_steps 
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

            open(3, file='shape_phi', status='unknown')
            open(4, file='phi_var', status='unknown')
      
            nxyz = Nx*Ny*Nz
            nxy = Nx*Ny
      
            phi_tot = 0.0
            var_count = 0
      
            !write(1,*) 'zone', ' ', 'i=', Nx, ' ', 'j=', Ny, ' ', 'k=', Nz, 'f=', 'point'
            !write(3,*) 'zone', ' ', 'i=', Nx, ' ', 'j=', Ny, ' ', 'k=', Nz, 'f=', 'point'

            ! Determine boundary locations
            ii = 0
            allocate(boundary_locations(nxyz,2,3))

            file_num = 0

            ! Loop through data_save.xxxxx files according to dims
            do kk = 1, dims(2)
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

               ! Determine boundary locations on planes
               j_ppt = 180 
               !j_ppt = 70
               j1 = Ny - j_ppt + 1

               ! ====================
               ! Analyze the first plane
               j = j1
               p = 1
 
               do k = ist(3), ien(3)
               do i = 1, Nx

                  phi_sq = 0.d0
                  do ivar = 1, var
                    phi_sq = phi_sq + phi(ivar,i,j,k)*phi(ivar,i,j,k)
                  end do
 
!                 ! Determine precipitate boundary as 3 or more overlapping phases containing ppt
!                 do ivar = 1, var
!                   if (phi(ivar,i,j,k) > threshold) &
!                     phi_contained = phi_contained + 1
!                 end do

                 !if (phi(var,i,j,k) > threshold .and. phi_contained >= 3) then

                 ! Determine precipitate boundary as phi(var)>threshold and phi_sq<threshold(?)
                 if (phi(var,i,j,k) > threshold .and. phi_sq < threshold) then
                   ii = ii + 1
                   boundary_locations(ii,p,:) = (/ i, j, k /)
                 end if

                 ! Determine precipitate boundary as exactly 3 overlapping phases
                 !phi_contained = 0
                 !do ivar = 1, var
                 !  if (0.1 <= phi(ivar,i,j,k) .and. phi(ivar,i,j,k) <= 0.11) &
                 !    phi_contained = phi_contained + 1
                 !end do
 
                 !if (phi_contained.eq.3) then
                 !  ii = ii + 1
                 !  boundary_locations(ii,:) = (/ i, j, k /)
                 !end if
 
               end do
               end do

               ! ====================
               ! Analyze the second plane
               j = j_ppt
               p = 2 
 
               do k = ist(3), ien(3)
               do i = 1, Nx

                  phi_sq = 0.d0
                  do ivar = 1, var
                    phi_sq = phi_sq + phi(ivar,i,j,k)*phi(ivar,i,j,k)
                  end do
 
!                 ! Determine precipitate boundary as 3 or more overlapping phases containing ppt
!                 do ivar = 1, var
!                   if (phi(ivar,i,j,k) > threshold) &
!                     phi_contained = phi_contained + 1
!                 end do

                 !if (phi(var,i,j,k) > threshold .and. phi_contained >= 3) then

                 ! Determine precipitate boundary as phi(var)>threshold and phi_sq<threshold(?)
                 if (phi(var,i,j,k) > threshold .and. phi_sq < threshold) then
                   ii = ii + 1
                   boundary_locations(ii,p,:) = (/ i, j, k /)
                 end if

                 ! Determine precipitate boundary as exactly 3 overlapping phases
                 !phi_contained = 0
                 !do ivar = 1, var
                 !  if (0.1 <= phi(ivar,i,j,k) .and. phi(ivar,i,j,k) <= 0.11) &
                 !    phi_contained = phi_contained + 1
                 !end do
 
                 !if (phi_contained.eq.3) then
                 !  ii = ii + 1
                 !  boundary_locations(ii,:) = (/ i, j, k /)
                 !end if
 
               end do
               end do
      
               deallocate (phi)
               deallocate (con)
      
            end do

            ! ====================
            ! Find max extent for each ppt
            ! ====================
            do p = 1, 2
              ! Find max diameter of precipitate (farthest apart)
              max_distance = 0.d0
  
              ! Loop through all precipitate boundary voxels on each plane separately
              do i = 1, ii 
                if (boundary_locations(i, p, 1) /= 0 .and. boundary_locations(i, p, 2) /= 0 .and. boundary_locations(i, p, 3) /= 0) then
                  do j = 1, ii 
                    if (boundary_locations(j, p, 1) /= 0 .and. boundary_locations(j, p, 2) /= 0 .and. boundary_locations(j, p, 3) /= 0) then
                      ! Calculate distance between boundary voxels i and j
                      dx = boundary_locations(j, p, 1) - boundary_locations(i, p, 1)
                      dy = boundary_locations(j, p, 2) - boundary_locations(i, p, 2)
                      dz = boundary_locations(j, p, 3) - boundary_locations(i, p, 3)
      
                      distance = (dx**2 + dy**2 + dz**2)**0.5
      
                      if (distance > max_distance) &
                        max_distance = distance
                    end if
    
                  end do
                end if
              end do
              
              if (p.eq.1) &
                write(91, '(I8, F6.2, 2I8)') step_i, max_distance, ii
              if (p.eq.2) &
                write(92, '(I8, F6.2, 2I8)') step_i, max_distance, ii

            end do
            ! ====================

            ! ====================
            ! Statistical estimate of mean thickness of second phase
            ! ====================
            ! TODO: Calculate volume fraction of second phase
            ! TODO: N random rays in xz plane intersecting both planes 
              ! TODO: Calculate intersections of rays
              ! TODO: Calculate mean thickness prediction
            ! TODO: Calculate average mean thickness prediction

            ! ====================

            deallocate(boundary_locations)

      end do

      deallocate(step_strings)

      end program analyze_plane
