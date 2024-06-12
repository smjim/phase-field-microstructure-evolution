      program paraview
      implicit none
      integer Nx, Ny, Nz, Nxyz, nprocs, dims(2), myid, ist(3), ien(3), Nxy
      integer i, j, k, kk, var, ivar, step, var_max, ii, ll
      integer ix, jy, kz
      integer ierr
      integer, dimension(:), allocatable :: var_count
      integer, dimension(:), allocatable :: num_var_gl
      real *8, dimension(:,:,:,:), allocatable :: phi
      real *8, dimension(:,:,:), allocatable :: con
      real *8, dimension(:), allocatable ::  con_gl
      real *8, dimension(:), allocatable ::  phi_sq_gl
      real *8, dimension(:), allocatable ::  phi_3_gl, grad_phi_3_gl
      real *8, dimension(:), allocatable :: boundary_gl 
      real *8 threshold, phi_above_threshold

      real *8 term, phi_sq, phi_sum, phi_tot, phi_max, phi_temp
      character(len=256) :: input_dir, output_dir, step_dir
      character(len=32) :: file_num
      character(len=256) :: infile, outfile_phi_sq, outfile_con 
      character(len=256) :: outfile_var_num, outfile_phi_3, outfile_grad_phi_3
      character(len=256) :: outfile_boundary
      character(len=256) :: command
      character(len=64), dimension(:), allocatable :: step_strings
      character(len=64) :: step_line
      character(len=7) :: step_string
      integer :: num_steps, unit, ios, step_i

      call get_command_argument(1, input_dir, status=ierr)
      call get_command_argument(2, output_dir, status=ierr)

      threshold = 0.2 ! Threshold for boundary calculation
      
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

      ! Allocate array for step_strings
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

      ! Loop through step directories
      do step_i = 1, num_steps 
            step_string = trim(step_strings(step_i))
            step_dir = trim(input_dir) // '/step_' & 
                  //trim(adjustl(step_string)) // '/'

            print *, "Analyzing directory: ", step_dir

            infile = trim(step_dir) // 'data_extract_mult.in'

            open(10, file=trim(step_dir)//'var_dist', status='unknown')
            open(11, file=trim(infile), status='old')
      
            read(11,*) Nx, Ny, Nz, nprocs, var
            read(11,*) dims(1), dims(2)
            close(11)
            !open(3, file=trim(output_dir)//'/shape_phi', status='unknown')
            !open(4, file=trim(output_dir)//'/phi_var', status='unknown')
      
            Nxyz = Nx*Ny*Nz
            Nxy = Nx*Ny

      
            allocate (phi_sq_gl(1:Nxyz))
            allocate (phi_3_gl(1:Nxyz))
            allocate (grad_phi_3_gl(1:Nxyz))
            allocate (con_gl(1:Nxyz))
            allocate (num_var_gl(1:Nxyz))
            allocate (var_count(var))
            allocate (boundary_gl(1:Nxyz))
      
            ! Write VTK files
            outfile_phi_sq = trim(output_dir) // '/phi_sq_' // trim(adjustl(step_string)) // '.vtk'
            outfile_con = trim(output_dir) // '/con_' // trim(adjustl(step_string)) // '.vtk'
            outfile_var_num = trim(output_dir) // '/var_num_' // trim(adjustl(step_string)) // '.vtk'
            outfile_phi_3 = trim(output_dir) // '/phi_3_' // trim(adjustl(step_string)) // '.vtk'
            outfile_grad_phi_3 = trim(output_dir) // '/grad_phi_3_' // trim(adjustl(step_string)) // '.vtk'
            outfile_boundary = trim(output_dir) // '/boundary_' // trim(adjustl(step_string)) // '.vtk'

            open(21, file=trim(outfile_phi_sq), status='unknown')
            open(22, file=trim(outfile_con), status='unknown')
            open(23, file=trim(outfile_var_num), status='unknown')
            open(24, file=trim(outfile_phi_3), status='unknown')
            open(25, file=trim(outfile_grad_phi_3), status='unknown')
            open(26, file=trim(outfile_boundary), status='unknown')
      
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

            print *, "writing file for step ", step_string
            open(95, file=trim(output_dir) // '/phi_values' // trim(adjustl(step_string)) // '.dat', status='unknown')
            write(95,'(A)') "# x, y, z, phi 1, phi 2, phi 3, phi_sq"

            phi_tot = 0.0
            var_count = 0
      
            !write(1,*) 'zone', ' ', 'i=', Nx, ' ', 'j=', Ny, ' ', 'k=', Nz, 'f=', 'point'
            !write(3,*) 'zone', ' ', 'i=', Nx, ' ', 'j=', Ny, ' ', 'k=', Nz, 'f=', 'point'
      
            myid = 0
            do kk=1,dims(2)
               do j=1,dims(1)
                  myid = myid + 1
                  write(file_num,'(i4.4)') myid
                  open(12,file=trim(step_dir)//'/data_save.'//file_num, status='old', &
                      form='unformatted')
                  read(12) ist, ien
                  if ( j.eq. 1 ) allocate( phi(var, Nx, Ny, ist(3):ien(3)) )
                  if ( j.eq. 1 ) allocate( con(Nx, Ny, ist(3):ien(3)) )
                  read(12) phi(1:var,1:Nx,ist(2):ien(2),ist(3):ien(3)),  &
                           con(1:Nx,ist(2):ien(2),ist(3):ien(3))
                  close(12)
               end do
                  !write(95,*) Nx, Ny, ist, ien
                  !call flush(95)
      
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
                    if (ii.eq.var) phi_3_gl(ll) = phi(3,i,j,k) ! var is the "precipitate phase"?
                    phi_sq = phi_sq + phi(ii,i,j,k)*phi(ii,i,j,k)
                    if(phi(ii,i,j,k).gt.phi_max) var_max = ii
                    if(phi(ii,i,j,k).gt.threshold) phi_above_threshold = phi_above_threshold + 1
                  end do

                  write(95,'(3(I6), 4(F6.2))') i, j, k, phi(1,i,j,k), phi(2,i,j,k), phi(3,i,j,k), phi_sq

                  phi_sq_gl(ll) = phi_sq
                  num_var_gl(ll) = var_max

                  ! if more than one phi is nonzero, that signifies a boundary
                  !if (phi_above_threshold .gt. 1) then 
                  !  boundary_gl(ll) = 1.0 !phi_above_threshold
                  !else
                  !  boundary_gl(ll) = 0.0
                  !end if
                  boundary_gl(ll) = phi_above_threshold
               end do
               end do
               end do
      
               deallocate (phi)
               deallocate (con)
            end do

            ! Calculate gradient of phi field
            call calculate_phi_grad(Nx, Ny, Nz, phi_3_gl, grad_phi_3_gl)
      
            ! Write outputs

               write(21,*) phi_sq_gl
               write(22,*) con_gl
               write(23,*) num_var_gl
               write(24,*) phi_3_gl
               write(25,*) grad_phi_3_gl
               write(26,*) boundary_gl
      
            close(21)
            close(22)
            close(23)
            close(24)
            close(25)
            close(26)

            close(95)

            deallocate(phi_sq_gl)
            deallocate(phi_3_gl)
            deallocate(grad_phi_3_gl)
            deallocate(con_gl)
            deallocate(num_var_gl)
            deallocate(var_count)
            deallocate(boundary_gl)

      end do
     
      contains

      ! Subroutine for calculation of gradients of phi
      subroutine calculate_phi_grad(Nx, Ny, Nz, phi, grad_phi_mag)
            implicit none

            integer, intent(in) :: Nx, Ny, Nz
            real(8), intent(in) :: phi(:)
            real(8), dimension(:,:,:), allocatable :: phi_3, grad_phi_3_mag
            real(8), dimension(:,:,:), allocatable :: grad_x, grad_y, grad_z
            real(8), intent(out) :: grad_phi_mag(:)
            integer :: i, j, k

            allocate(grad_x(Nx, Ny, Nz))
            allocate(grad_y(Nx, Ny, Nz))
            allocate(grad_z(Nx, Ny, Nz))

            allocate(phi_3(Nx, Ny, Nz))
            allocate(grad_phi_3_mag(Nx, Ny, Nz))

            ! Reshape phi array
            phi_3 = reshape(phi, [Nx, Ny, Nz])

            ! Initialize gradients to zero
            grad_x = 0.0
            grad_y = 0.0
            grad_z = 0.0

            ! Compute gradients using central differences
            do k = 2, Nz-1
            do j = 2, Ny-1
            do i = 2, Nx-1
              grad_x(i,j,k) = (phi_3(i+1,j,k) - phi_3(i-1,j,k)) / 2.0
              grad_y(i,j,k) = (phi_3(i,j+1,k) - phi_3(i,j-1,k)) / 2.0
              grad_z(i,j,k) = (phi_3(i,j,k+1) - phi_3(i,j,k-1)) / 2.0 

              ! Calculate grad magnitude array
              grad_phi_3_mag(i,j,k) = sqrt(grad_x(i,j,k)**2 + grad_y(i,j,k)**2 + grad_z(i,j,k)**2)
            end do
            end do
            end do

            ! Reshape grad magnitude array
            grad_phi_mag = reshape(grad_phi_3_mag, [Nx*Ny*Nz])
            
            deallocate(phi_3)
            deallocate(grad_phi_3_mag)

            deallocate(grad_x)
            deallocate(grad_y)
            deallocate(grad_z)

      end subroutine calculate_phi_grad

      end program paraview

