      program analyze_phase 

      implicit none
      integer Nx, Ny, Nz, nxyz, nprocs, dims(2), myid, ist(3), ien(3), nxy
      integer i, j, k, kk, var, ivar, step, var_max, ii, ll
      integer ix, jy, kz
      integer ierr
      integer, dimension(:), allocatable :: var_count
      integer, dimension(:), allocatable :: num_var_gl
      real *8, dimension(:,:,:,:), allocatable :: phi
      real *8, dimension(:,:,:), allocatable :: con
      real *8, dimension(:), allocatable ::  phi_sq_gl
      real *8, dimension(:), allocatable ::  phi_3_gl
      real *8, dimension(:), allocatable ::  con_gl
      real *8 :: term, phi_sq, phi_sum, phi_tot, phi_max, phi_temp, boundary_val
      character(len=256) :: input_dir, output_dir, step_dir
      character(len=32) :: file_num
      character(len=256) :: infile, data_filename
      character(len=256) :: command
      character(len=64), dimension(:), allocatable :: step_strings
      character(len=64) :: step_line, boundary_val_str
      character(len=7) :: step_string
      integer :: num_steps, unit, ios, step_i

      call get_command_argument(1, input_dir, status=ierr)
      call get_command_argument(2, output_dir, status=ierr)
      call get_command_argument(3, boundary_val_str, status=ierr)

      read(boundary_val_str, *) boundary_val
      print *, boundary_val
      
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

      data_filename = trim(output_dir) // '/boundary_data.dat'
      print *, data_filename 
      open(20, file=data_filename, status='unknown')
      write(20, '(A)') "# time_step, boundary_voxels, &
                  boundary_phase_fraction" 
      close(20)

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
            open(3, file=trim(output_dir)//'/shape_phi', status='unknown')
            open(4,   file=trim(output_dir)//'/phi_var', status='unknown')
      
            nxyz = Nx*Ny*Nz
            nxy = Nx*Ny

      
            allocate (phi_sq_gl(1:nxyz))
            allocate (phi_3_gl(1:nxyz))
            allocate (con_gl(1:nxyz))
            allocate (num_var_gl(1:nxyz))
            allocate (var_count(var))
      
            phi_tot = 0.0
            var_count = 0
      
            myid = 0
            do kk=1,dims(2)
               do j=1,dims(1)
                  myid = myid + 1
                  write(file_num,'(i4.4)') myid
                  open(12,file=trim(step_dir)//'data_save.'//file_num, status='old', &
                      form='unformatted')
                  read(12) ist, ien
                  if ( j.eq. 1 ) allocate( phi(var, Nx, Ny, ist(3):ien(3)) )
                  if ( j.eq. 1 ) allocate( con(Nx, Ny, ist(3):ien(3)) )
                  read(12) phi(1:var,1:Nx,ist(2):ien(2),ist(3):ien(3)),  &
                           con(1:Nx,ist(2):ien(2),ist(3):ien(3))
                  close(12)
               end do
                  !write(95,*) Nx, Ny, ist, ien
                  call flush(95)
      
               do k = ist(3), ien(3)
               do j = 1, Ny
               do i = 1, Nx 
      
                  ll = (k-1) * nxy + (j-1) * Nx + i
                  con_gl(ll) = con(i,j,k)
                  phi_sq = 0.0
                  phi_max = 0.8
                  var_max = 0
                  do ii = 1, var
                  if (ii.eq.var) phi_3_gl(ll) = phi(3,i,j,k)
                  phi_sq = phi_sq + phi(ii,i,j,k)*phi(ii,i,j,k)
                  if(phi(ii,i,j,k).gt.phi_max) var_max = ii
                  end do
                  phi_sq_gl(ll) = phi_sq
                  num_var_gl(ll) = var_max
               end do
               end do
               end do
      
               deallocate (phi)
               deallocate (con)
            end do
      
            ! Analyze data here for timestep trim(adjustl(step_string))
            ! The following data are ready to be analyzed
            ! phi_sq_gl, con_gl, num_var_gl, phi_3_gl
            ! save boundary data to data_filename
            call analyze_data(Nx, Ny, Nz, phi_3_gl, boundary_val, data_filename, output_dir, step_string)

            deallocate(phi_sq_gl)
            deallocate(phi_3_gl)
            deallocate(con_gl)
            deallocate(num_var_gl)
            deallocate(var_count)

      end do

      contains 

      ! Subroutine for calculation of gradients of phi
      subroutine calculate_phi_grad(Nx, Ny, Nz, phi, grad_x, grad_y, grad_z)
            implicit none

            integer, intent(in) :: Nx, Ny, Nz
            real(8), intent(in) :: phi(:,:,:)
            real(8), intent(out) :: grad_x(:,:,:), grad_y(:,:,:), grad_z(:,:,:)
            integer :: i, j, k

            ! Initialize gradients to zero
            grad_x = 0.0
            grad_y = 0.0
            grad_z = 0.0
        
            ! Compute gradients using central differences
            do k = 2, Nz-1
            do j = 2, Ny-1
            do i = 2, Nx-1
              grad_x(i,j,k) = (phi(i+1,j,k) - phi(i-1,j,k)) / 2.0
              grad_y(i,j,k) = (phi(i,j+1,k) - phi(i,j-1,k)) / 2.0
              grad_z(i,j,k) = (phi(i,j,k+1) - phi(i,j,k-1)) / 2.0
            end do
            end do
            end do
      end subroutine calculate_phi_grad

      ! Subroutine for analysis of .vtk data
      subroutine analyze_data(Nx, Ny, Nz, phi_3_gl, boundary_val, output_file, output_dir, step_string)
            implicit none

            integer, intent(in) :: Nx, Ny, Nz
            real *8, dimension(:), intent(in) :: phi_3_gl
            real *8, intent(in) :: boundary_val

            integer :: nxy
            character(len=*), intent(in) :: output_file, output_dir, step_string
            integer :: i, j, k, boundary_voxels, phase_1_boundary_voxels
            real(8) :: fraction_boundary_voxels
            real(8), dimension(:,:,:), allocatable :: grad_phi_3_x, grad_phi_3_y, grad_phi_3_z
            real(8), dimension(:,:,:), allocatable :: phi_3(:,:,:)

            character(len=256) :: vtk_filename
            integer :: unit

            ! Reshape 1D array to 3D
            allocate(phi_3(Nx, Ny, Nz))
            phi_3 = reshape(phi_3_gl, [Nx, Ny, Nz])
    
            allocate(grad_phi_3_x(Nx, Ny, Nz))
            allocate(grad_phi_3_y(Nx, Ny, Nz))
            allocate(grad_phi_3_z(Nx, Ny, Nz))

            ! Calculate gradients
            call calculate_phi_grad(Nx, Ny, Nz, phi_3, grad_phi_3_x, grad_phi_3_y, grad_phi_3_z) 
    
            boundary_voxels = 0
            phase_1_boundary_voxels = 0
    
            ! Determine boundary voxels
            nxy = Nx*Ny
            !open(14, file=data_filename, status='unknown')
            do k = 2, Nz-1
            do j = 2, Ny-1
            do i = 2, Nx-1
              !write(14, '(F16.8)') sqrt(grad_phi_3_x(i,j,k)**2 + grad_phi_3_y(i,j,k)**2 + grad_phi_3_z(i,j,k)**2)
              if (sqrt(grad_phi_3_x(i,j,k)**2 + grad_phi_3_y(i,j,k)**2 + grad_phi_3_z(i,j,k)**2) > boundary_val) then
                boundary_voxels = boundary_voxels + 1
                if (num_var_gl((k-1)*nxy + (j-1)*Nx + i) == 1) then
                  phase_1_boundary_voxels = phase_1_boundary_voxels + 1
                end if
              end if
            end do
            end do
            end do
            !close(14)
    
            fraction_boundary_voxels = real(phase_1_boundary_voxels) / real(boundary_voxels)
            if (real(boundary_voxels)==0) fraction_boundary_voxels = 0
    
            ! Write VTK file for boundary voxels visualization
            vtk_filename = trim(output_dir) // '/phi_boundary_' // trim(adjustl(step_string)) // '.vtk'
            print *, "vtk filename", vtk_filename
            open(10, file=vtk_filename, status='unknown')
            write(10, '(A)') "# vtk DataFile Version 3.0"
            write(10, '(A)') "Points Example"
            write(10, '(A)') "ASCII"
            write(10, '(A)') "DATASET POLYDATA"
            write(10, '(A, I8, A)') "POINTS ", boundary_voxels, " float"
    
            ! Write points data
            do k = 1, Nz
            do j = 1, Ny
            do i = 1, Nx
              if (sqrt(grad_phi_3_x(i,j,k)**2 + grad_phi_3_y(i,j,k)**2 + grad_phi_3_z(i,j,k)**2) > boundary_val) then
                write(10, '(3I8)') i, j, k
              end if
            end do
            end do
            end do

            ! Write lines data, this is necessary for paraview
            ! visualization unfortunately
            write(10, '(A, I8, I8)') "LINES ", INT(boundary_voxels/2)+1, boundary_voxels+2
            do i = 0, INT(boundary_voxels/2)
              write(10, '(A, I8, A, I8)') "2 ", 2*i, " ", 2*i+1
            end do
    
            close(10)
    
            ! Write data file with boundary voxel information
            print *, "output file: ", output_file
            open(unit=20, file=output_file, status="OLD", &
                        action="write", position="append")
            write(20,'(A, I8, F16.8)') trim(adjustl(step_string)), boundary_voxels, fraction_boundary_voxels
            close(20)
    
            deallocate(grad_phi_3_x)
            deallocate(grad_phi_3_y)
            deallocate(grad_phi_3_z)
      end subroutine analyze_data
     
      end program analyze_phase
