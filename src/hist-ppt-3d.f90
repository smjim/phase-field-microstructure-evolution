program histogram_ppt_sizes

    use mpi

    implicit none

    ! I/O files
    character(len=256) :: input_dir, infile, outfile

    ! Step tracking and mpi organization
    integer :: ierr, myid, ios, nprocs, num_steps, i
    integer :: my_start_step, my_end_step, step_i
    character(len=256) :: command, step_line
    character(len=7) :: step_string
    character(len=64), dimension(:), allocatable :: step_strings


    ! Initialize mpi
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  
    ! Process 0 reads the command arguments and broadcasts them
    if (myid == 0) then
      call get_command_argument(1, input_dir, status=ierr)
    endif

    call MPI_Bcast(input_dir, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

          ! Process 0 determines steps and broadcasts to step_strings
      if (myid == 0) then
        ! Determine step directories
        command = 'ls ' // trim(input_dir) // &
                    ' | grep "^ppt_loc_" | cut -d"_" -f3 | cut -d"." -f1 > steps.txt'
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
        !call execute_command_line('rm steps.txt')
      end if

      ! Broadcast number of steps and step_strings to all processes
      call MPI_Bcast(num_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (myid /= 0) then
        allocate(step_strings(num_steps))
      end if
      call MPI_Bcast(step_strings, num_steps*64, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      ! Determine range of steps each process will handle
      my_start_step = (num_steps / nprocs) * myid + 1
      if (myid == nprocs - 1) then
        my_end_step = num_steps
      else
        my_end_step = (num_steps / nprocs) * (myid + 1)
      end if
      !print *, "myid, num_steps, my_start_step, my_end_step: ", myid, num_steps, my_start_step, my_end_step

      ! Each process gets its own section of steps to analyze
      do step_i = my_start_step, my_end_step
        step_string = trim(step_strings(step_i))
        print *, myid, " analyzing step ", step_string

        infile = trim(input_dir) // 'ppt_loc_'  // trim(step_string) // '.dat'
        outfile = trim(input_dir) // 'ppt_hist_'  // trim(step_string) // '.dat'
        call analyze_step(infile, outfile)
      end do
    
contains

    subroutine analyze_step(infile, outfile)

        implicit none

        character(len=256), intent(in) :: infile, outfile
        
        integer, parameter :: max_locations = 100000    ! Adjust as needed
        integer, parameter :: max_region_sizes = 100000 ! Adjust as needed
        
        integer :: ppt_locations(max_locations, 3)      ! Array to store ppt locations (x, y, z)
        logical :: checked_locations(max_locations)     ! Boolean array to track checked locations
        integer :: region_sizes(max_region_sizes)       ! Array to store sizes of each ppt region
        integer, allocatable :: stack(:,:)              ! Stack for flood fill (stores x, y, z)
        integer :: num_locations, num_checked_locations
        integer :: region_size, num_regions, stack_size
        integer :: i, j, k, m, n, Nx, Ny, Nz
        integer :: random_ppt_loc(3), current_ppt_loc(3), neighbor(3)
    
        logical :: in_stack
        
        ! Read ppt locations from file
        call read_ppt_locations(ppt_locations, num_locations, Nx, Ny, Nz, infile)
        
        if (num_locations > 0) then
          ! Initialize checked_locations array
          checked_locations = .false.
          
          ! Initialize variables
          num_checked_locations = 0
          num_regions = 0
      
          allocate(stack(num_locations,3))
          
          ! Loop until all locations are checked
          do while (num_checked_locations < num_locations)
          
              ! Find unchecked ppt location
              do i = 1, num_locations
                  if (.not. checked_locations(i)) then
                      random_ppt_loc = ppt_locations(i, :)
                      exit
                  end if
              end do
              !write(*,'(A,5I8)') "num_regions, num_checked_locations, ppt_loc(3):", num_regions, num_checked_locations, random_ppt_loc
              
              ! Perform flood fill starting from random_ppt_loc
              num_regions = num_regions + 1
              stack = 0
              stack_size = 1
              stack(stack_size, :) = random_ppt_loc
              
              region_size = 0
      
              do while (stack_size > 0)
                  !write(93,'(I8)') stack_size
      
                  ! Pop stack element
                  current_ppt_loc = stack(stack_size, :)
                  stack(stack_size, :) = 0
      
                  ! Determine the index of `current_ppt_loc` in `ppt_locations`
                  do i = 1, num_locations
                      if (all(ppt_locations(i, :) == current_ppt_loc) .and. &
                          .not. checked_locations(i)) then
                          ! Mark the ppt location as checked
                          checked_locations(i) = .true.
                          num_checked_locations = num_checked_locations + 1
                          exit ! Exit loop once found
                      end if
                  end do
      
                  ! Decrement stack_size to move to pop stack element
                  stack_size = stack_size - 1
                  
                  ! Increment region size
                  region_size = region_size + 1
      
                  ! Find neighboring ppt locations
                  do k = -1, 1
                    do j = -1, 1
                      do i = -1, 1
                        if (i == 0 .and. j == 0 .and. k == 0) then
                          cycle
                        end if
      
                        neighbor(1) = current_ppt_loc(1) + i
                        neighbor(2) = current_ppt_loc(2) + j
                        neighbor(3) = current_ppt_loc(3) + k
      
                        ! Check if neighbor is within bounds
                        if (neighbor(1) >= 0 .and. neighbor(1) <= Nx .and. &
                            neighbor(2) >= 0 .and. neighbor(2) <= Ny .and. &
                            neighbor(3) >= 0 .and. neighbor(3) <= Nz) then
                          ! Check if neighbor is a ppt location
                          do m = 1, num_locations
                            if (ppt_locations(m, 1) == neighbor(1) .and. &
                                ppt_locations(m, 2) == neighbor(2) .and. &
                                ppt_locations(m, 3) == neighbor(3) .and. &
                                .not. checked_locations(m)) then
                              in_stack = .false.
                              do n = 1, stack_size
                                if (stack(n,1) == neighbor(1) .and. &
                                    stack(n,2) == neighbor(2) .and. &
                                    stack(n,3) == neighbor(3)) then 
                                  in_stack = .true. 
                                  exit
                                end if
                              end do
      
                              ! Push neighbor onto stack
                              if (.not. in_stack) then
                                stack_size = stack_size + 1
                                !print *, "neighbor is added to stack, stacksize: ", neighbor(:), stack_size
                                stack(stack_size, :) = neighbor(:)
                                ! HERE IS WHERE THE ERROR IS!!
                              end if
      
                            end if
                          end do
      
                        end if
      
                      end do
                    end do
                  end do
              end do 
              
              ! Store region size
              !print *, "filename, num_regions, region_size: ", infile, num_regions, region_size
              region_sizes(num_regions) = region_size
          end do 
      
          deallocate(stack)
  
          !write(*, '(2A)') "finished analysis of ", trim(infile)
          
          ! Write region sizes to a .dat file
          !write(*, '(2A,I8)') "num_regions: ", trim(infile), num_regions
          if (num_regions > 0) &
            call write_region_sizes(region_sizes, num_regions, outfile)
        end if

    end subroutine

    subroutine read_ppt_locations(ppt_locations, num_locations, Nx, Ny, Nz, filename)
        implicit none
        integer, intent(out) :: ppt_locations(:, :)
        integer, intent(out) :: num_locations
        integer, intent(out) :: Nx, Ny, Nz 
        character(len=*), intent(in) :: filename
        integer :: iunit, i, j
        
        ! Open the file
        open(unit=10, file=filename, status='old', action='read', iostat=iunit)
        if (iunit /= 0) then
            write(*, '(2A)') "Error opening file: ", filename
            stop
        end if
        
        ! Read number of locations
        read(10, *) num_locations, Nx, Ny, Nz
        
        ! Read ppt locations
        do i = 1, num_locations
            read(10, *) ppt_locations(i, :)
        end do
        
        ! Close the file
        close(10)
        
        ! Display ppt locations
        !write(*,'(2A, I8)') "Filename, Number of ppt locations read: ", filename, num_locations
        !write(*,'(I8)') num_locations
        !do i = 1, num_locations
        !    print *, "PPT Location", i, ":", ppt_locations(i, :)
        !end do
    end subroutine read_ppt_locations
    
    subroutine write_region_sizes(region_sizes, num_regions, filename)
        implicit none
        integer, intent(in) :: region_sizes(:)
        integer, intent(in) :: num_regions
        character(len=256) :: filename
        integer :: i
        
        ! Open the file
        open(unit=20, file=filename, status='replace', action='write')
        
        ! Write region sizes
        do i = 1, num_regions
            write(20, '(I8)') region_sizes(i)
        end do
        
        ! Close the file
        close(20)
        
        ! Display message
        !write(*,'(2A)') "Region sizes saved to ", trim(filename)
    end subroutine write_region_sizes
    
end program histogram_ppt_sizes

