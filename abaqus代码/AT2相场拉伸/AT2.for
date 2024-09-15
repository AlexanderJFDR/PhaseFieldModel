!**********************************************************************************************************
!
      module NumKind
!
!**********************************************************************************************************
        implicit none
        integer (kind(1)), parameter :: ikind = kind(1), 
     &                                  rkind = kind(0.D0), 
     &                                  lkind = kind(.true.)
!       
      end module Numkind
      

!**********************************************************************************************************
!
      module ModelParam
!
!**********************************************************************************************************
        use NumKind
        implicit none

        ! Constants

        ! Flag of initilization
        logical (lkind) :: bInitialized = .false.

        ! Tolerance
        real    (rkind), parameter :: TOL = 1.0d-12 
        ! number of guass points
        integer (ikind), parameter :: ngp = 2
        
        ! geometric function parameter
        real(rkind), parameter :: c0 = 2.d0
        
        ! 
        real(rkind) :: thk, A11, A12, A22, A66, Gf, lb
        real(rkind) :: De(3, 3)        
        real(rkind) :: gp(ngp), gw(ngp)
        real(rkind) :: QQ(12,12) 

        !
        contains

        !===================================
          subroutine Initialize(props, nprops)
        !===================================

            integer (ikind), intent (in) :: nprops
            real    (rkind), intent (in) :: props(nprops)

            !********************************************
            integer(ikind) :: indexq(12), i

            ! material properties
            A11  =  props(1) ! props(1) -- A11
            A12  =  props(2) ! props(2) -- A12
            A22  =  props(3) ! props(3) -- A22
            A66  =  props(4) ! props(4) -- A66
            Gf   =  props(5) ! props(5) -- fracture energy
            lb   =  props(6) ! props(6) -- length scale            
            thk = 1.0
            write (*, *) 'props = ', props
            ! elastic stiffness matrix
            De(:,1) = (/ A11,  A12, 0.D0/)
            De(:,2) = (/ A12,  A22, 0.D0/)
            De(:,3) = (/0.D0, 0.D0,  A66/)
            
            ! integration points
            gp = (/ -1.d0, 1.d0 /) / dsqrt(3.d0)
            gw = (/  1.d0, 1.d0 /)
            
            ! dof interchange
            indexq = (/ 1,2,9, 3,4,10, 5,6,11, 7,8,12 /)
            ! interchange the locations of dofs
            QQ = 0.d0
            do i = 1, 12
              QQ(indexq(i),i) = 1.d0
            end do             
            
            bInitialized = .true.
            
            return
          end subroutine Initialize
      !========================================================================= 
      end module ModelParam

!**********************************************************************************************************
!
      module FEM
!
!**********************************************************************************************************
        use NumKind
        implicit none

        contains      
          !==================shape function and its derivative with xi and eta======================    
          subroutine shapefuc(n, dn_xieta, xi, eta)
          
            implicit none      
            real(rkind) :: n(4), dn_xieta(2, 4), xi, eta

            n(1) = 0.25d0*(1.d0 - xi)*(1.d0 - eta)
            n(2) = 0.25d0*(1.d0 + xi)*(1.d0 - eta)
            n(3) = 0.25d0*(1.d0 + xi)*(1.d0 + eta)
            n(4) = 0.25d0*(1.d0 - xi)*(1.d0 + eta)
            
            dn_xieta(1, 1) = -0.25d0*(1.d0 - eta)
            dn_xieta(1, 2) =  0.25d0*(1.d0 - eta)
            dn_xieta(1, 3) =  0.25d0*(1.d0 + eta)
            dn_xieta(1, 4) = -0.25d0*(1.d0 + eta)
            
            dn_xieta(2, 1) = -0.25d0*(1.d0 - xi)
            dn_xieta(2, 2) = -0.25d0*(1.d0 + xi)
            dn_xieta(2, 3) =  0.25d0*(1.d0 + xi)
            dn_xieta(2, 4) =  0.25d0*(1.d0 - xi)
            
            return 
          end subroutine shapefuc

          !===============traditional b matrix==============================================      
          subroutine b_matrix(nd,bd,b,det_jacb, coords,xi,eta)
          
            implicit none
            real(rkind) :: nd(4), bd(2,4), b(3,8)
            real(rkind) :: jacb(2,2), inv_jacb(2,2), coords(2, 4)
            real(rkind) :: det_jacb, xi, eta
            
            !local varibles
            real(rkind) :: n(4), dn_xieta(2,4), dn_x(4), dn_y(4)
            integer(ikind) :: i, j
             
            ! shape functions 
            call shapefuc(n,dn_xieta,xi,eta)
            nd = n
            
            ! jacob matrix
            jacb = matmul(dn_xieta, transpose(coords))            
            det_jacb = jacb(1,1)*jacb(2,2) - jacb(1,2)*jacb(2,1)
            inv_jacb(1, 1) = jacb(2, 2)
            inv_jacb(1, 2) =-jacb(1, 2)
            inv_jacb(2, 1) =-jacb(2, 1)
            inv_jacb(2, 2) = jacb(1, 1)
            inv_jacb = 1.d0/det_jacb*inv_jacb            
            
            !initialize varibles
            do i = 1,4
              dn_x(i) = inv_jacb(1,1)*dn_xieta(1,i)
     &                + inv_jacb(1,2)*dn_xieta(2,i)
              dn_y(i) = inv_jacb(2,1)*dn_xieta(1,i)
     &                + inv_jacb(2,2)*dn_xieta(2,i)
            end do
            
            ! B matrix for displacement
            b = 0.d0
            do j = 1, 4
              b(1, 2*(j-1) + 1) = dn_x(j)
              b(2, 2*(j-1) + 2) = dn_y(j)
              b(3, 2*(j-1) + 1) = dn_y(j)
              b(3, 2*(j-1) + 2) = dn_x(j)
            end do
            
            ! B matrix for damage
            do j = 1,4
              bd(1,j) = dn_x(j)
              bd(2,j) = dn_y(j)
            end do
          
            return
          end subroutine b_matrix
      
        !********************************************************************
        ! define the dyadic function
          function dyadic(vector1,vector2, vlen)
        !********************************************************************
            integer (ikind) :: vlen, i, j
            real    (rkind) :: vector1(vlen),vector2(vlen)
            real    (rkind) :: dyadic(vlen,vlen)
          
            do i = 1, vlen
              do j = 1, vlen
                dyadic(i,j) = vector1(i) * vector2(j)
              end do
            end do

            return
          end function dyadic

      end module FEM
      

!**********************************************************************************************************
!
      subroutine pfczm(rhs,amatrix,coords,u,svars)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        use FEM
        implicit none

        real(rkind):: rhs(12), amatrix(12,12), coords(2,4)
        real(rkind):: svars(4), u(12)
       
        ! local varibles
        real(rkind):: b(3,8), nd(4), bd(2,4)
        real(rkind):: uu(8), dd(4), rd(4), ru(8), kdd(4,4), kuu(8,8)
        real(rkind):: rr(12), kk(12,12)
        real(rkind):: strain(3), stressEff(3)
        real(rkind):: det_jacb, energy_crk, phi, DH, DDH
        real(rkind):: omega, domega, ddomega
        real(rkind):: dalpha, ddalpha, phi_source, dphi_source, dvol
        integer(ikind):: i, j, k
        real(rkind):: savg, sdif, sdev, smax, smin
        
        ! extrat nodal displacement and damage dofs
        do i = 1, 4
          uu(2*i - 1) = u(3*i - 2)
          uu(2*i    ) = u(3*i - 1)
          dd(i)       = u(3*i)
        end do
        
        ! initialize varibles
        rd  = 0.d0
        kdd = 0.d0
        kuu = 0.d0
        do i = 1, ngp
          do j = 1, ngp      
            call b_matrix(nd,bd,b,det_jacb, coords,gp(i),gp(j))
              
            strain = matmul(b, uu)  ! strain field
            stressEff = matmul(De, strain) ! effective stress
            phi  = dot_product(nd,dd) ! crack phase-field
                        
            call geometricFunc(dalpha,ddalpha,phi) ! geometric function
            call historyFunc(DH,DDH,strain,phi)
     
            phi_source  = DH + Gf/(c0*lb)*dalpha
            dphi_source = DDH + Gf/(c0*lb)*ddalpha

            ! residual for damage
            dvol=  gw(i)*gw(j)*det_jacb*thk
            rd  =  rd  - dvol*(phi_source*nd + 2.d0*lb*Gf/c0
     &          *  matmul(transpose(bd), matmul(bd, dd)))

            ! element matrices
            kdd =  kdd + dvol*((dphi_source)*dyadic(nd, nd, 4)
     &          +  2.d0*lb*Gf/c0*matmul(transpose(bd),bd))
            
            call energeticFunc(omega,domega,ddomega,phi)

            kuu =  kuu + dvol*matmul(matmul(transpose(b), omega*De), b)
          end do
        end do        
        ru = -matmul(kuu,uu) ! applies to hybrid formulation
        
        rr(1:8 ) = ru
        rr(9:12) = rd
        
        kk = 0.d0
        kk(1:8 , 1:8 ) = kuu
        kk(9:12, 9:12) = kdd
        
        rhs     = matmul(transpose(QQ),rr)
        amatrix = matmul(matmul(transpose(QQ),kk),QQ)
  
        return 
      end subroutine pfczm
      
!**********************************************************************************************************
!
      subroutine energeticFunc(omega,domega,ddomega,phi)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        implicit none
      
        real(rkind) :: omega, domega, ddomega, phi
        
        omega   =  (1.d0-phi)**2.d0      
        domega  =  -2.d0*(1.d0-phi)
        ddomega =  2.d0
        return
      end subroutine energeticFunc
      
!**********************************************************************************************************
!
      subroutine geometricFunc(dalpha,ddalpha,phi)
!
!**********************************************************************************************************
        use NumKind
        implicit none
        
        real(rkind) :: dalpha, phi, ddalpha
        
        dalpha  = 2.d0*phi
        ddalpha = 2.d0
        
        return 
      end subroutine geometricFunc  

!**********************************************************************************************************
!
      subroutine historyFunc(DH,DDH,strain,phi)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        implicit none
      
        real(rkind) :: DH, DDH, phi
        real(rkind) :: omega, domega, ddomega
        real(rkind) :: strain(3)

        call energeticFunc(omega,domega,ddomega,phi)

        DH = domega*dot_product(strain,matmul(De,strain))/2
        DDH = ddomega*dot_product(strain,matmul(De,strain))/2
        
        return
      end subroutine historyFunc

!**********************************************************************************************************
      subroutine UEL(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars,
     &               props, nprops, coords, mcrd, nnode, 
     &               u, du, v, a, jtype, time, dtime, kstep, 
     &               kinc, jelem, params, ndload, jdltyp, adlmag,
     &               predef, npredf, lflags, mlvarx, ddlmag, mdload,
     &               pnewdt, jprops,njprop,period)
!**********************************************************************************************************

        use NumKind
        use ModelParam
        implicit none

!**********************************************************************************************************
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! variables passed in
        integer (ikind), intent (in    ) :: ndofel, mlvarx, nrhs, 
     &    nsvars, nprops, mcrd, nnode, jtype, kstep, kinc, jelem, 
     &    ndload,npredf, mdload, njprop
     
        integer (ikind), intent (in    ) :: jdltyp(mdload,*), 
     &    lflags(*), jprops(njprop)
     
        real    (rkind), intent (in    ) :: props(nprops), 
     &    coords(mcrd,nnode), u(ndofel), du(mlvarx,*), v(ndofel), 
     &    a(ndofel), time(2), params(3), adlmag(mdload,*),
     &    ddlmag(mdload,*), predef(2,npredf,nnode), dtime, period
     
        ! variables can be updated
        real    (rkind), intent (in out) :: pnewdt
  
        ! variables to be updated (the update of energy(8) is optional)
        real    (rkind), intent (in out) :: rhs(mlvarx,nrhs), 
     &    amatrx(ndofel,ndofel), svars(nsvars), energy(8)
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !********************************************************************************************************
      !                   
      ! user coding to define rhs, amatrx, svars, energy and pnewdt (optional for the last two)
      !
        
        write (*, *) 'TIME = ', time(2)

        ! initialize parameters, etc.                     
        if (.not. bInitialized) call Initialize(props, nprops)
        
        ! right now only Q4 element is implemented
        call pfczm(rhs(:,1),amatrx,coords,u,svars)
        
        return

      end subroutine uel
