      subroutine inv_trans(in, out, ist, ien, fst, fen)

      use p3dfft
      implicit none

      integer ist(3), ien(3), fst(3), fen(3)
      double complex in( fst(1):fen(1), fst(2):fen(2), fst(3):fen(3) )
      real *8 out( ist(1):ien(1), ist(2):ien(2), ist(3):ien(3) )

      call p3dfft_btran_c2r(in, out, 'tff')

      return
      end
