
!! This module provides routine L1_regression, which
!! calculates the L1-norm-regression coefficients to data points x(:,i),y(i).
!!
!! L1-norm regression is similar to the well known minimization of square-error
!! regression, i.e. L2-norm regression.
!! There are the following main differences between them:
!! - 1) L2-regression allows for a solution by solving a set of linear
!!      equations. i.e. It is easy and therefor the standard.
!!      L1-regression does not allow for such a simple solution.
!!      Instead a convex unconstraint minimization problem must be solved.
!!      The convex objective function is built up of hyperplane pieces.
!! - 2) L2-regression is unstable with respect to outliers, L1-regression is
!!      stable. That is why L1-regression is also called robust regression.
!!      The reason is that the solution of L2-regression is in some sense the
!!      mean of the slopes, whereas for L1-regression it is the median.
!!      When some points move outwards to become outliers it affects the mean
!!      but not the median.
!!
!! This software was written by Erich W. Steiner.
!! It comes with NO warranty / guarantee. Everybody using it is responsible
!! for the results himself. The author does not accept any liability for
!! the results of this software or errors that maybe contained in the software.
!! This software is provided as is under the 
!! GNU-Lesser-Public-Licence. This note must not be removed.

MODULE L1_regression_module

   !! realPara is the kind-type for reals. 8 is usually double precision.
   integer, parameter :: realPara = 8  
   
   !! g_tolerance defines a small positive number for considering
   !! values to be zero. Round-off tolerance. Set it to be roughly to
   !! square root of the precision of realPara real.
   
   real(realPara), parameter :: g_tolerance = 1.0E-8 
   
   !! V*U=S  *=matmul.  Hyperplane i: (x-Xs(:,i)*V(:,i) = 0. *=dot_product.
   !! The components of GrSchmidtType are treated read-only in the code
   !! apart from those routine having GrSchmidt in their name which are
   !! intended to modify such an object.
   
   type GrSchmidtType

      integer :: n !! The present number of hyperplanes stored.
   
      real(realPara), dimension(:,:), allocatable :: Xs  !! Stuetzvektoren
      real(realPara), dimension(:,:), allocatable :: V   !! gradients of hyperplanes
      real(realPara), dimension(:,:), allocatable :: U   !! Upper-triangular matrix.
      real(realPara), dimension(:,:), allocatable :: S !!orthogonal Gr-Schmidt vectors.
   
   end type GrSchmidtType
   
   private
   public :: L1_regression, realPara, g_tolerance, calcMedian2
   
 CONTAINS
 
      
   !! Routine calculates the linear regression for the points (i_x(:,i), i_y(i))
   !! minimizing the L1-norm. I.e the i-th column of i_x gives the x-coordinates
   !! of the i-th data point having y-value i_y(i). 
   !! The linear regression equation is:
   !!
   !! a*x + b = y     
   !!
   !! * meaning scalar dot_product. o_f is the minimal absolute deviation 
   !! between the approximated y and actual y values of the solution o_a, o_b.
   
   SUBROUTINE L1_regression( i_x, i_y, o_a, o_b, o_f )
   
      real(realPara), intent(in), dimension(:,:) :: i_x
      real(realPara), intent(in), dimension(size(i_x, dim=2)) :: i_y
      real(realPara), intent(out), dimension(size(i_x,dim=1)) :: o_a
      real(realPara), intent(out) :: o_b, o_f
      
      real(realPara), dimension(size(o_a)) :: grad 
      real(realPara), dimension(size(o_a),2) :: nextA
      real(realPara), dimension(size(o_a)+1) :: g1, a1, objDir, reducedObjDir
      real(realPara), dimension(2) :: nextB, nextF, last_f
      type(GrSchmidtType) :: GrS
      integer :: i, k
            
      o_a = 0
      o_b = calcMedian2(i_y)
      
      GrS = GrSchmidtConstructorFn( i_n = size(o_a) + 1 )
      
      grad = calcGradient( i_x = i_x, &
                        i_left = (/ (o_b,i=1,size(i_y)) /), &
                           i_y = i_y )
                           
      o_f = sum(abs(i_y - o_b))
      last_f = o_f * 2
      
      g1(:size(grad)) = grad
      g1(size(g1))    = -1
      
      a1(:size(o_a)) = o_a
      a1(size(a1))   = o_f
                           
      call addHyperplaneGrSchmidt( c_x = GrS, &
                                i_vect = g1, &
                                   i_x = a1 )
                 
      objDir = 0
      objDir(size(objDir)) = -1

      loop: do
               
         if(all(abs(grad) < g_tolerance)) exit loop
      
         call lineSearch( i_x = i_x, &
                          i_y = i_y, &
                          i_a = o_a, &
                          i_b = o_b, &
                 i_searchDirA = -grad, &
                      o_nextA = nextA, &
                      o_nextB = nextB, &
                      o_nextF = nextF )
                      
         if(nextF(1) < nextF(2)) then
            k = 1
         else
            k = 2
         end if
         
         o_a = nextA(:,k)
         o_b = nextB(k)
         o_f = nextF(k)
         
         a1(:size(o_a)) = nextA(:,k)
         a1(size(a1))   = nextF(k)
         
         call releaseHyperplIfNotTightGrSchmidt( c_x = GrS, &
                                                 i_x = a1 )
         
         do i = 1, 2
         
            g1(:size(grad)) = calcGradient( i_x = i_x, &
                              i_left = matmul(nextA(:,i), i_x) + nextB(i) , &
                                 i_y = i_y )
                                                                  
            call addHyperplaneGrSchmidt( c_x = GrS, &
                                      i_vect = g1, &
                                         i_x = a1 )
            if( i == 2 ) then
            
               call removeHyperplIfKTcoeffNegGrSchmidt( c_x = GrS, &
                                                 i_objVect = objDir, &
                                    o_opt_reducedObjVect = reducedObjDir )
                                    
            else if( GrS%n == size(GrS%V,dim=1) ) then
            
               call removeHyperplIfKTcoeffNegGrSchmidt( c_x = GrS, &
                                                 i_objVect = objDir )
            end if

         end do
                  
         if( GrS%n == size(GrS%V, dim=1) ) exit loop
         
         if(all(abs(reducedObjDir) < g_tolerance )) exit loop
         
         if( last_f(1)*(1.0-g_tolerance) <= o_f ) exit loop  ! no progress stop
                                                             ! condition
         last_f(1) = last_f(2)
         last_f(2) = o_f
         
         grad = reducedObjDir(:size(grad))
         
      end do loop
      
      call freeSpaceGrSchmidt( GrS )
      
   END SUBROUTINE L1_regression
          
 
   !! Solves 2-dimensional L1-norm regression problem with bisection approach
   !! of minimizing the sum of absolute deviations.
   !! o_a(1), o_b(1) are the left solutions of the bisection with
   !! negative or = 0 gradient o_g(1) and sum of absolute deviations o_f(1),
   !! o_a(2), o_b(2), o_g(2), o_f(2) are the corresponding ones of the right
   !! bisection point solution with strictly positive gradient o_g(2).
   !! Stop when o_a(2)-o_a(1) <= g_tolerance or when o_g(1) < g_tolerance.
   !! The regression equation is:
   !!
   !! a*x + b = y.
   !!
   !! The used bisection incorporates binary bisection steps and 
   !! Newton-Raphson like subgradients intersection steps.
 
   SUBROUTINE regression_2D( i_x, i_y, o_a, o_b, o_g, o_f, i_opt_step )
   
      real(realPara), intent(in), dimension(:) :: i_x
      real(realPara), intent(in), dimension(size(i_x)) :: i_y 
      real(realPara), intent(out), dimension(2) :: o_a, o_b, o_g, o_f
      real(realPara), intent(in), optional :: i_opt_step
      
      integer :: i, outL
      real(realPara) :: step, a, b, g, f, proj_f, intvLen_old, intvLen
      real(realPara), dimension(size(i_x)) :: v, yMinusV, yPlusTol, yMinusTol
      logical :: bisect, checkOpti, isOptimal
      
      yPlusTol  = i_y + g_tolerance
      yMinusTol = i_y - g_tolerance
      
      o_a(1) = 0
      
      o_b(1) = calcMedian2( i_y )
            
      o_g(1) = gradient( i_left = (/ (o_b(1), i=1,size(i_x)) /) )
      
      !! Find o_a(2), o_b(2), o_g(2), o_f(2).
      
      step = 1.0
      if(present(i_opt_step)) step = abs(i_opt_step)
      
      outerLoop: do outL = 1, 5      ! nur zur Sicherheit
                                     !um unendlich-loop zu vermeiden.
      
         if( o_g(1) > 0.0 ) step = -step     
      
         loop: do
      
            o_a(2) = o_a(1) + step
         
            v = o_a(2) * i_x
         
            yMinusV = i_y - v
         
            o_b(2) = calcMedian2( yMinusV )
         
            o_g(2) = gradient( i_left = v + o_b(2) )
         
            if( o_g(2) * o_g(1) <= 0 ) exit loop
         
            o_a(1) = o_a(2)
            o_b(1) = o_b(2)
            o_g(1) = o_g(2)
         
            step = 2*step
         
         end do loop
      
         o_f(1) = sum( abs( i_y - o_a(1)*i_x - o_b(1) ) ) 
         o_f(2) = sum( abs( i_y - o_a(2)*i_x - o_b(2) ) )
      
         if( step < 0 ) then
            o_a = o_a(2:1:-1)
            o_b = o_b(2:1:-1)
            o_g = o_g(2:1:-1)
            o_f = o_f(2:1:-1)
         end if
      
      ! do the bisection: bisect or subgradient-bisect
      
         bisect = .true.
      
         loop2: do
            
            if( bisect ) then
         
               a = sum(o_a) / 2
                        
            else
         
            ! calc projected a, proj_f  as intersection of the left and right
            ! subgradient lines:
            
               a = ( o_a(2)*o_g(2) - o_a(1)*o_g(1) + o_f(1) - o_f(2) ) / &
                   ( o_g(2) - o_g(1) )
                                
               if( a < o_a(1) .or. a > o_a(2) ) then
                        
                  bisect = .true.
               
                  cycle loop2
               
               end if
                
               proj_f = o_g(1) * (a - o_a(1)) + o_f(1)
            
            end if
         
            v = a * i_x
         
            yMinusV = i_y - v
         
            b = calcMedian2( yMinusV )
         
            g = gradient( i_left = v + b )
         
            f = sum( abs( yMinusV - b ) )
                     
            intvLen_old = o_a(2) - o_a(1)
         
            if( g <= 0 ) then
               o_a(1) = a
               o_b(1) = b
               o_g(1) = g
               o_f(1) = f
            else
               o_a(2) = a
               o_b(2) = b
               o_g(2) = g
               o_f(2) = f
            end if
          
         ! Iteration stop conditions:
          
            intvLen = o_a(2) - o_a(1)
          
            checkOpti = .false.
          
            if( intvLen <= g_tolerance ) then
         
               checkOpti = .true.
            
            else if( o_g(2) - o_g(1) <= g_tolerance ) then
         
               checkOpti = .true.
         
            else if(.not.bisect) then
                     
               if( abs(proj_f - f) <= g_tolerance ) checkOpti = .true.
            
            end if
         
            if(checkOpti) then
         
               call checkForOptimality( o_isOptimal = isOptimal )
            
               if( isOptimal ) exit outerLoop
            
               cycle outerLoop
            
            end if
         
         
            if(bisect) then
         
               bisect = .false.
            
            else
         
               bisect = intvLen > 0.5 * intvLen_old
            
            end if
         
         end do loop2
         
      end do outerLoop
      
   
      contains
      
         function gradient( i_left ) result(res)
            real(realPara), intent(in), dimension(:) :: i_left
            real(realPara) :: res 
            
            res = sum( i_x, mask = i_left > yPlusTol )
            
            res = res - sum( i_x, mask = i_left < yMinusTol )
            
         end function gradient
         
         
         subroutine checkForOptimality( o_isOptimal )
            logical, intent(out) :: o_isOptimal
            
            if( o_f(1) < o_f(2) ) then
               a = o_a(1) - 100*g_tolerance
            else if( o_f(1) > o_f(2) ) then
               a = o_a(2) + 100*g_tolerance
            else
               o_isOptimal = .true.
               return
            end if
            
            v = a * i_x
            yMinusV = i_y - v
            b = calcMedian2( yMinusV )
            f = sum( abs( yMinusV - b ) )
                        
            if( a < o_a(1) ) then
               
               if( f >= o_f(1) ) then
                  o_isOptimal = .true.
                  return
               end if
               
               o_isOptimal = .false.
               o_g(1) = 1
               step   = 1000*g_tolerance
               
            else
            
               if( f >= o_f(2) ) then
                  o_isOptimal = .true.
                  return
               end if
               
               o_isOptimal = .false.
               o_a(1) = o_a(2)
               o_b(1) = o_b(2)
               o_g(1) = -1
               o_f(1) = o_f(2)
               step   = 1000*g_tolerance
               
            end if
            
         end subroutine checkForOptimality
               
   END SUBROUTINE regression_2D
      
      
   !! Allocates the components of res to n*n matrices and assigns 0 as values to them.
   
   FUNCTION GrSchmidtConstructorFn( i_n ) result(res)
      integer, intent(in) :: i_n
      type(GrSchmidtType) :: res 
      
      allocate( res%Xs(i_n,i_n) )
      allocate( res%V(i_n,i_n) )
      allocate( res%U(i_n,i_n) )
      allocate( res%S(i_n,i_n) )
      
      res%n = 0
      res%Xs= 0
      res%V = 0
      res%U = 0
      res%S = 0
      
   END FUNCTION GrSchmidtConstructorFn
   
   
   !! Deallocates the matrix components of c_x.
   
   SUBROUTINE freeSpaceGrSchmidt( c_x )
      type(GrSchmidtType), intent(inout) :: c_x
      
      if(allocated(c_x%Xs )) deallocate(c_x%Xs)
      if(allocated(c_x%V )) deallocate(c_x%V)
      if(allocated(c_x%U )) deallocate(c_x%U)
      if(allocated(c_x%S )) deallocate(c_x%S)
      
   END SUBROUTINE freeSpaceGrSchmidt
   
   
   !! Adds the hyperplane  (x-i_x)*i_vect=0 to c_x. *=dot_product.
   !! If i_vect is linearly dependent to presently stored hyperplane normal vectors
   !! then i_vect, i_x is not added.
   
   SUBROUTINE addHyperplaneGrSchmidt( c_x, i_vect, i_x )
      type(GrSchmidtType), intent(inout) :: c_x
      real(realPara), intent(in), dimension(size(c_x%V,dim=1)) :: i_vect, i_x
      
      real(realPara), dimension(size(i_x)) :: s
      real(realPara), dimension(c_x%n) :: z
      real(realPara) :: norm
      integer :: i
      intrinsic :: dot_product, sqrt, matmul
      
      s = i_vect
      
      do i = 1, c_x%n
      
         z(i) = -dot_product(i_vect, c_x%S(:,i))
      
         s = s + z(i)*c_x%S(:,i)
         
      end do
      
      if(all(abs(s) < g_tolerance)) return   ! i_vect is linearly dependent on
                                            !presently stored %S, %V vectors.

      c_x%n = c_x%n + 1
      
      c_x%V(:,c_x%n) = i_vect
      
      c_x%Xs(:,c_x%n) = i_x 

      norm = sqrt(dot_product(s,s))
      
      s = s / norm
      
      z = z / norm
      
      c_x%S(:,c_x%n) = s
      
      i = c_x%n - 1
      
      c_x%U(c_x%n,c_x%n) = 1.0 / norm
      
      c_x%U(:i,c_x%n) = matmul( c_x%U(:i,:i), z )
      
   END SUBROUTINE addHyperplaneGrSchmidt
      
      
   SUBROUTINE updateHyperplaneGrSchmidt( c_x, i_vect, i_x, i_objVect, &
                                          o_isOptimal, o_reducedObjVect )
      type(GrSchmidtType), intent(inout) :: c_x
      real(realPara), intent(in), dimension(size(c_x%V,dim=1)) :: i_vect, &
                                                            i_x, i_objVect
      logical, intent(out) :: o_isOptimal
      real(realPara), intent(out), dimension(size(c_x%V,dim=1)) :: o_reducedObjVect
      
      call releaseHyperplIfNotTightGrSchmidt( c_x, i_x )
      
      call addHyperplaneGrSchmidt( c_x, i_vect, i_x )
      
      call removeHyperplIfKTcoeffNegGrSchmidt( c_x, i_objVect, o_reducedObjVect )
      
      if( c_x%n == size(c_x%V, dim=1) ) then
      
         o_isOptimal = .true.
      
      else 
      
         o_isOptimal = all(abs(o_reducedObjVect) < g_tolerance )
         
      end if
      
   END SUBROUTINE updateHyperplaneGrSchmidt
   
   
   SUBROUTINE releaseHyperplIfNotTightGrSchmidt( c_x, i_x )
   
      type(GrSchmidtType), intent(inout) :: c_x
      real(realPara), intent(in), dimension(size(c_x%V,dim=1)) :: i_x
      
      logical, dimension(c_x%n) :: maske
      integer :: i

      intrinsic :: abs, dot_product
      
      forall( i = 1:size(maske) ) maske(i) = &
         abs(dot_product( i_x - c_x%Xs(:,i), c_x%V(:,i) )) <= g_tolerance
         
      call keepOnlyMaskedHyperplGrSchmidt( c_x = c_x, i_maske = maske )
                     
   END SUBROUTINE releaseHyperplIfNotTightGrSchmidt
      
      
   SUBROUTINE keepOnlyMaskedHyperplGrSchmidt( c_x, i_maske )
   
      type(GrSchmidtType), intent(inout) :: c_x
      logical, intent(in), dimension(:) :: i_maske
      
      integer :: i
      integer, dimension(:), allocatable :: indxs
      real(realPara), dimension(:,:), allocatable :: tmpXs, tmpV
      intrinsic :: all
      
      if(all(i_maske)) return
      
      if(all(.not.i_maske)) then
      
         c_x%n = 0
         
         return
         
      end if
      
      i = firstIndxTrue( .not.i_maske )
      
      c_x%n = i - 1
      
      call getIndxSet3( indxSet = indxs, maske = i_maske(i:) )
      
      if(size(indxs) == 0) goto 100 
      
      indxs = indxs + i - 1
         
      allocate( tmpXs( size(c_x%V,dim=1), size(indxs) ) )
      allocate(  tmpV( size(c_x%V,dim=1), size(indxs) ) )
         
      tmpXs = c_x%Xs(:,indxs)
      tmpV  =  c_x%V(:,indxs)
         
      do i = 1, size(indxs)
      
         call addHyperplaneGrSchmidt( c_x = c_x,  &
                                   i_vect = tmpV(:,i), &
                                      i_x = tmpXs(:,i) )
      end do
   
      deallocate( tmpXs, tmpV )
  100 if(allocated(indxs)) deallocate(indxs)
                     
   END SUBROUTINE  keepOnlyMaskedHyperplGrSchmidt   
      
      
   !! Removes those hyperplanes with negative KT-coefficient.
   !! If i_opt_iterativeRemove given and = true then does such remove
   !! steps iteratively until all KT-coefficients are positive.
      
   SUBROUTINE removeHyperplIfKTcoeffNegGrSchmidt( c_x, i_objVect, &
                                                  o_opt_reducedObjVect, &
                                                  i_opt_iterativeRemove )
      type(GrSchmidtType), intent(inout) :: c_x
      real(realPara), intent(in), dimension(size(c_x%V,dim=1)) :: i_objVect
      real(realPara), intent(out), dimension(size(c_x%V,dim=1)), optional :: &
                                       o_opt_reducedObjVect
      logical, intent(in), optional :: i_opt_iterativeRemove
                                       
      real(realPara), dimension(:), allocatable :: k
      integer :: n_old
      logical :: iterativeRemove 
      intrinsic :: matmul
      
      if( c_x%n == 0 ) then    !nothing to be removed.
      
         if(present(o_opt_reducedObjVect)) &
            o_opt_reducedObjVect = i_objVect
         
         return
         
      end if
      
      iterativeRemove = .false.
      if(present(i_opt_iterativeRemove)) &
         iterativeRemove = i_opt_iterativeRemove
      
      loop: do
      
         if(allocated(k)) deallocate(k)
         allocate( k(c_x%n) )
      
         k = matmul( i_objVect, c_x%S(:,:c_x%n) )
      
         k = matmul( c_x%U(:c_x%n,:c_x%n), k )
      
         n_old = c_x%n
               
         call keepOnlyMaskedHyperplGrSchmidt( c_x = c_x, &
                                          i_maske = k >= 0.0 )
                                          
         if(.not.iterativeRemove) exit loop
         if(n_old == c_x%n) exit loop
         if(c_x%n == 0) exit loop
         
      end do loop
                      
      if(.not.present(o_opt_reducedObjVect)) goto 100
                      
      if( c_x%n == n_old ) then
      
         o_opt_reducedObjVect = i_objVect - matmul( c_x%V(:,:c_x%n), k )
         
      else if( c_x%n > 0 ) then
      
         if(allocated(k)) deallocate(k)
         allocate( k(c_x%n) )
         
         k = matmul( i_objVect, c_x%S(:,:c_x%n) )
         
         o_opt_reducedObjVect = i_objVect - matmul( c_x%S(:,:c_x%n), k )
               
      else
      
         o_opt_reducedObjVect = i_objVect
         
      end if
      
 100  if(allocated(k)) deallocate(k)
         
   END SUBROUTINE removeHyperplIfKTcoeffNegGrSchmidt
      
             
         
   FUNCTION calcGradient( i_x, i_left, i_y ) result(g)
      real(realPara), intent(in), dimension(:,:) :: i_x 
      real(realPara), intent(in), dimension(size(i_x,dim=2)) :: i_left, i_y
      real(realPara), dimension(size(i_x,dim=1)) :: g
      
      logical, dimension(size(i_y)) :: maske1, maske2
      integer :: i
      
      maske1 = i_left > i_y
      maske2 = i_left < i_y
      
      forall( i=1:size(g) ) g(i) = sum(i_x(i,:), mask=maske1) - &
                                   sum(i_x(i,:), mask=maske2)
      
   END FUNCTION calcGradient

   
   SUBROUTINE lineSearch( i_x, i_y, i_a, i_b, i_searchDirA, & 
                          o_nextA, o_nextB, o_nextF )
                          
      real(realPara), intent(in), dimension(:,:) :: i_x
      real(realPara), intent(in), dimension(size(i_x,dim=2)) :: i_y
      real(realPara), intent(in), dimension(size(i_x,dim=1)) :: i_a
      real(realPara), intent(in) :: i_b
      real(realPara), intent(in), dimension(size(i_a)) :: i_searchDirA
      real(realPara), intent(out), dimension(size(i_a),2) :: o_nextA
      real(realPara), intent(out), dimension(2) :: o_nextB, o_nextF
      
      real(realPara), dimension(size(i_y)) :: a_x, sa_x
      real(realPara), dimension(2) :: lt, lb, lg
      integer :: i
      real(realPara), dimension(size(i_searchDirA)) :: searchDirA 
      intrinsic :: matmul
            
      searchDirA = i_searchDirA / sqrt(dot_product(i_searchDirA, i_searchDirA))
      
      a_x = matmul(i_a, i_x)
      sa_x = matmul(searchDirA, i_x)
      
      call regression_2D( i_x = sa_x, &
                          i_y = i_y - i_b - a_x, &
                          o_a = lt, &
                          o_b = lb, &
                          o_g = lg, &
                          o_f = o_nextF )
                                                    
      forall(i = 1:2)
      
         o_nextA(:,i) = i_a + lt(i)*searchDirA
         
         o_nextB(i)   = i_b + lb(i)
                  
      end forall
      
   END SUBROUTINE lineSearch
         
   
      !! Let n=size(i_x). The function sorts i_x then returns sortedX(t) where
      !! t = (n+1)/2 if n is odd, or returns
      !! (sortedX(t) + sortedX(t+1)) / 2 where t=n/2 if n is even.

    PURE FUNCTION calcMedian2( i_x )  result(res)
         implicit none
         real(realPara), intent(in), dimension(:) :: i_x
         real(realPara) :: res

         integer, dimension(size(i_x)) :: perm
         integer :: n, t
         intrinsic :: modulo
      
         call real_sorting( i_arr = i_x, o_perm = perm )
         n = size(i_x)
         
         if(modulo(n,2) == 0) then
         
            t = n / 2
            res = (i_x(perm(t)) + i_x(perm(t+1))) / 2
            
         else
         
            t = (n+1) / 2
            res = i_x(perm(t))
            
         end if
         
    END FUNCTION calcMedian2
 
 
      !! Routine sorts i_arr such that i_arr(o_perm) is an increasing
      !! Array. Assumption: size(i_arr) = size(o_perm).
      !! Implements heap-sort. complexity n*log(n) with n=size(i_arr).

    PURE SUBROUTINE real_sorting(i_arr, o_perm)
      implicit none
      real(realPara), intent(in),  dimension(:) :: i_arr
      integer,        intent(out), dimension(:) :: o_perm

      integer, dimension(size(i_arr)) :: heap
      integer :: i, j, parent, child, child2, n

      !make heap: (parent is smaller than children.)

      forall(i = 1 : size(heap)) heap(i) = i

      loop1: do n = 2, size(heap)

         child = n

         loop2: do

            parent = child / 2    !swap child and parent, if parent bigger.

            if( i_arr(heap(parent)) > i_arr(heap(child)) ) then

               call swapping( heap(parent), heap(child) )

               if(parent == 1) exit loop2

               child = parent

            else

               exit loop2

            end if

         end do loop2

      end do loop1

      !heap order is established.

      n = size(o_perm)    !n = number of heap elements

      loop3: do i = 1, size(o_perm)

         o_perm(i) = heap(1)   !Take top element of heap to define o_perm(i)

         heap(1) = heap(n)     !Put last heap element on top

         n = n - 1

         !Restore heap order. Swap with the smallest child recursively, if
         !                    it is smaller than parent.

         parent = 1

         loop4: do

            child = 2 * parent

            if( child < n ) then

               child2 = child + 1    !Define j = smaller of the children

               if( i_arr(heap(child2)) < i_arr(heap(child)) ) then

                  j = child2

               else

                  j = child

               end if          !swap parent and j if child j is smaller.

               if( i_arr(heap(j)) < i_arr(heap(parent)) ) then

                  call swapping( heap(j), heap(parent) )

                  parent = j

               else

                  exit loop4

               end if

            else if( child > n ) then

               exit loop4

            else    !child = n

               if( i_arr(heap(child)) < i_arr(heap(parent)) ) &
                    call swapping( heap(child), heap(parent) )

               exit loop4

            end if

         end do loop4

      end do loop3

    END SUBROUTINE real_sorting

 
    ELEMENTAL SUBROUTINE swapping(a, b)
        implicit none
        integer, intent(inout) :: a, b

        integer :: temp

        temp = a

        a = b

        b = temp

    END SUBROUTINE swapping


        !! Routine gibt res zur�ck, so dass maske(res) = .true.
        !! und maske(i) = .false. f�r alle i < res.
        !! Wenn alle maske=.false. dann res = 0.

    PURE FUNCTION firstIndxTrue( maske )  result(res)
        implicit none
        logical, intent(in), dimension(:) :: maske
        integer :: res

        integer :: i

        res = 0

        loop: do i = 1, size(maske)

           if(maske(i)) then

              res = i

              exit loop

           end if

        end do loop

    END FUNCTION firstIndxTrue

        !! Routine alloziert indxSet und definiert indxSet, so dass
        !! maske(indxSet)==.true., indxSet maximal ist, indxSet der Groesse
        !! nach geordnet ist mit dem kleinsten Element vorne, kein Element
        !! in indxSet mehr als einmal vorkommt.

    PURE SUBROUTINE getIndxSet3(indxSet, maske)
        implicit none
        integer, intent(out), allocatable, dimension(:) :: indxSet
        logical, intent(in), dimension(:) :: maske

        integer :: i, j
        intrinsic :: count

        if(allocated(indxSet)) deallocate(indxSet)

        allocate( indxSet(count(maske)) )

        if(size(indxSet)==0) return


        j = 0

        loop: do i=1, size(maske)

           if(maske(i)) then

              j = j + 1

              indxSet(j) = i

              if(j==size(indxSet)) exit loop

           end if

        end do loop

    END SUBROUTINE getIndxSet3
      
END MODULE L1_regression_module
