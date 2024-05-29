!     Parallelized with OpenMP 5-28-2024
!     subroutine for calculating the k vectors in Fourier space

      subroutine k_space(Nx, Ny, Nz, fst, fen, kf, kf_sq, kf_4)
      implicit none
      integer Nx, Ny, Nz, fst(3), fen(3), i, j, k, ii, jj, kk
      real *8 pi
      double complex lamda1, lamda2, lamda3, unit1, unit2, unit3
      double complex unit, kf(3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
      double complex kf_sq(3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
      double complex kf_4(3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

      pi = datan(1.d0)*4.0
      unit = (0.0,1.0)
      unit = 2.d0*pi*unit
      unit1 = unit/float(Nx)
      unit2 = unit/float(Ny)
      unit3 = unit/float(Nz)

      !$omp parallel do private(i, j, k, ii, jj, kk, lambda1, lambda2, lambda3)
      do i = fst(3), fen(3)
      ii = i-1

      if(ii.le.Nx/2) then
        lamda1 = float(ii) * unit1
      else
        lamda1 = float(ii-Nx) * unit1
      end if

      do j = fst(2), fen(2)
      jj = j-1

      if(jj.le.Ny/2) then
        lamda2 = float(jj) * unit2
      else
        lamda2 = float(jj-Ny) * unit2
      end if

      do k = fst(1), fen(1)
      kk = k-1

      if(kk.le.Nz/2) then
        lamda3 = float(kk) * unit3
      else
        lamda3 = float(kk-Nz) * unit3
      end if

      kf(1,k,j,i) = lamda1
      kf(2,k,j,i) = lamda2
      kf(3,k,j,i) = lamda3
      if (2*ii.eq.Nx) kf(1,k,j,i) = 0.d0
      if (2*jj.eq.Ny) kf(2,k,j,i) = 0.d0
      if (2*kk.eq.Nz) kf(3,k,j,i) = 0.d0
      kf_sq(1,k,j,i) = lamda1*lamda1
      kf_sq(2,k,j,i) = lamda2*lamda2
      kf_sq(3,k,j,i) = lamda3*lamda3
      kf_4(1,k,j,i) = lamda1*lamda1*lamda1*lamda1
      kf_4(2,k,j,i) = lamda2*lamda2*lamda2*lamda2
      kf_4(3,k,j,i) = lamda3*lamda3*lamda3*lamda3

      end do
      end do
      end do
      !$omp end parallel do

      return 
      end
