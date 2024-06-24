      subroutine f_trans(in, out, Nx, Ny, Nz, ist, ien, fst, fen)

      use p3dfft
      implicit none

      integer Nx, Ny, Nz, ist(3), ien(3), fst(3), fen(3)
      real *8 in( ist(1):ien(1), ist(2):ien(2), ist(3):ien(3) )
      double complex out( fst(1):fen(1), fst(2):fen(2), fst(3):fen(3) )

!     write(*,*) 'entered f_trans'
      call p3dfft_ftran_r2c(in, out, 'fft')
      out = out/float(Nx*Ny*Nz)

      return
      end
