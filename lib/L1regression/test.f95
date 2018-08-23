
!! This test-program produces a number of random test-problems and calls
!! the L1-regression-routine for each and checks the optimality of the 
!! calculated regression.
!!
!! Compile with gfortran (for example) with
!!
!! gfortran L1_regression_module.f95 test.f95
!!
!! that produces an executable a.out under Linux. a.exe under Windows.
!! gfortran is the GNU Fortran compiler. Downloadable on the WWW for all sorts of
!! platforms.

PROGRAM test
   use L1_regression_module

   real(realPara), dimension(:,:), allocatable :: x
   real(realPara), dimension(:), allocatable :: y 
   real(realPara), dimension(:), allocatable :: a, a2
   real(realPara) :: b, f
   integer :: n, m, i 
   logical :: ok
            
   do i = 1, 10
      
      n = randomInteger( 20000 )
      m = randomInteger( 100 )
         
      if(n <= m) n = n + m
         
      allocate( x(m,n), y(n), a(m), a2(m) )
         
      call random_number(a2)
      a2 = 20*a2 - 10
         
      call random_number(x)
         
      call random_number(y)
         
      y = 2*y + matmul(a2, x)
      
      deallocate( a2 )
         
      call L1_regression( i_x = x, i_y = y, &
                          o_a = a, &
                          o_b = b, &
                          o_f = f )
                          
      print*,"Testproblem ", i, " with n, m= ", n, m

      call checkOptimality( i_x = x, i_y = y, i_a = a, i_b = b, &
                               i_f = f, i_delta = 0.01_realPara,  &
                               o_opt_ok = ok )
                         
      deallocate( x, y, a )
                         
      print*,"ok status= ", ok
                         
      if(.not.ok) stop
         
   end do

   
 contains


   SUBROUTINE checkOptimality( i_x, i_y, i_a, i_b, i_f, i_delta, o_opt_ok )
   
      real(realPara), intent(in), dimension(:,:) :: i_x
      real(realPara), intent(in), dimension(:) :: i_y, i_a 
      real(realPara), intent(in) :: i_b, i_f, i_delta
      logical, intent(out), optional :: o_opt_ok
      
      integer :: i
      real(realPara) :: f, b
      real(realPara), dimension(size(i_a)) :: a
      real(realPara), dimension(size(i_y)) :: ax
      logical :: ok
      
      ok = .true.
      
      do i = 1, size(i_a)
      
         a = i_a 
         a(i) = a(i) - i_delta
         
         ax = matmul(a, i_x)
         
         b = calcMedian2( i_y - ax )
         
         f = sum(abs( i_y - ax - b ))
         
         if( f < i_f ) then
         
            print*,"i= ", i, "a(i) minus delta with ", f
            ok = .false.
            
         end if
         
         a(i) = a(i) + 2*i_delta
         
         ax = matmul(a, i_x)
         
         b = calcMedian2( i_y - ax )
         
         f = sum(abs( i_y - ax - b ))
         
         if( f < i_f ) then
         
            print*,"i= ", i, "a(i) plus delta with ", f
            ok = .false.
            
         end if
           
      end do
            
      if(present(o_opt_ok)) o_opt_ok = ok
      
   END SUBROUTINE checkOptimality
   
   !! Produces a random integer in the interval [1,i_n].
   
   FUNCTION randomInteger( i_n ) result(res)
      implicit none
      integer, intent(in) :: i_n
      integer :: res

      real :: harvest
      intrinsic :: random_number, floor

      call random_number(harvest)

      res = floor( harvest * i_n ) + 1

   END FUNCTION randomInteger
 
END PROGRAM test