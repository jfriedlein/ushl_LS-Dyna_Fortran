c UEL-helper, e.g. get_initialNodalCoords      
      include "../UEL_helper_Fortran_LS-Dyna/UEL_helper.f"
c
!...
c
      subroutine ushl_e101(force,stiff,ndtot,istif,
     . x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,
     . fx1,fx2,fx3,fx4,
     . fy1,fy2,fy3,fy4,
     . fz1,fz2,fz3,fz4,
     . xdof,
     . dx1,dx2,dx3,dx4,dy1,dy2,dy3,dy4,dz1,dz2,dz3,dz4,
     . wx1,wx2,wx3,wx4,wy1,wy2,wy3,wy4,wz1,wz2,wz3,wz4,
     . dxdof,
     . thick,thck1,thck2,thck3,thck4,
     . hsv,ihsv,nhsv,
     . cm,lmc,
     . gl11,gl21,gl31,gl12,gl22,gl32,gl13,gl23,gl33,
     . cmtrx,lft,llt)
c
c Linear element     
c
      use UEL_helper
      include 'nlqparm'
      dimension force(nlq,ndtot),stiff(nlq,ndtot,ndtot)
      dimension x1(nlq),x2(nlq),x3(nlq),x4(nlq)
      dimension y1(nlq),y2(nlq),y3(nlq),y4(nlq)
      dimension z1(nlq),z2(nlq),z3(nlq),z4(nlq)
      dimension fx1(nlq),fx2(nlq),fx3(nlq),fx4(nlq)
      dimension fy1(nlq),fy2(nlq),fy3(nlq),fy4(nlq)
      dimension fz1(nlq),fz2(nlq),fz3(nlq),fz4(nlq)
      dimension xdof(nlq,8,*)
      dimension dx1(nlq),dx2(nlq),dx3(nlq),dx4(nlq)
      dimension dy1(nlq),dy2(nlq),dy3(nlq),dy4(nlq)
      dimension dz1(nlq),dz2(nlq),dz3(nlq),dz4(nlq)
      dimension wx1(nlq),wx2(nlq),wx3(nlq),wx4(nlq)
      dimension wy1(nlq),wy2(nlq),wy3(nlq),wy4(nlq)
      dimension wz1(nlq),wz2(nlq),wz3(nlq),wz4(nlq)
      dimension dxdof(nlq,8,*)
      dimension thick(nlq),thck1(nlq),thck2(nlq),thck3(nlq),thck4(nlq)
      dimension hsv(nlq,nhsv),ihsv(nlq,nhsv),cm(lmc)
      dimension gl11(nlq),gl21(nlq),gl31(nlq),
     .     gl12(nlq),gl22(nlq),gl32(nlq),
     .     gl13(nlq),gl23(nlq),gl33(nlq)
      dimension cmtrx(nlq,15,3)
c
      integer, parameter :: dimen = 2
      integer, parameter :: n_nodes = 4
      integer, parameter :: ndtot_loc = 8
c
      ! Current element ID
       integer :: ele
      ! Iterator for nodes
       integer :: iNode
      ! Undeformed coordinates of element's nodes
       real*8 :: xUndef(n_nodes, dimen)
      ! Undeformed coordinates of element's nodes
       real*8 :: xUndef8(8, 3)
      ! Deformed coordinates of element's nodes
       real*8 :: xDef(n_nodes, dimen)
      ! Total displacements of element's nodes
       real :: uDisp(n_nodes,dimen)
      ! Element force vector
       real*8 :: forceVector(ndtot_loc)
      ! Element stiffness matrix
       real*8 :: stiffMatrix(ndtot_loc,ndtot_loc)       
      ! List of all dofs (displacements and xdofs)
       real*8 :: pe(ndtot_loc)
      ! Index map: AceGen - LS-Dyna dof-indices
       integer :: index_xy_(ndtot_loc)     
c
      ! Determine total number of dofs per node (integer-division intended)
       ndofs_per_node = ndtot_loc / n_nodes
c
      index_xy_ = 
     &              (/
     &               ! x,  y,  z, xrot, yrot, zrot
     &                 1,  2, !3,    4,    5,    6, ! Node 1
     &                 7,  8, !9,   10,   11,   12, ! Node 2
     &                13, 14,!15,   16,   17,   18, ! Node 3
     &                19, 20 !21,   22,   23,   24  ! Node 4
     &               /)    
c
       ELEMENTS: do ele = lft, llt
         ! Fill matrix of deformed coordinates for current element 'ele'
          xDef = reshape(
     &              (/
     &                x1(ele), y1(ele), 
     &                x2(ele), y2(ele),
     &                x3(ele), y3(ele),
     &                x4(ele), y4(ele)
     &               /), (/n_nodes,dimen/), order=(/2,1/)
     &           )
         ! Retrieve matrix of undeformed coordinates for current element 'ele'
          xUndef8 = get_initialNodalCoords_1element_R102( ele )   
          xUndef = xUndef8(1:n_nodes,1:dimen)  
         ! Compute total displacements for current element 'ele'
          uDisp = xDef - xUndef          
         ! Combine all dofs in the list 'pe' (displacements and possibly xdofs)
          do iNode=1,n_nodes
            pe( 1+(iNode-1)*ndofs_per_node : iNode*ndofs_per_node )
     &        = (/ uDisp(iNode,:) /)
          enddo          
         ! Initialise output variables of UEL to zero, because of internal "AddIn"
          forceVector = 0.0
          stiffMatrix = 0.0
c
         ! Call UEL-routine for current element 'ele'
          call Q1X_2D( xUndef, pe,
     &                 cm, hsv(ele,1:nhsv),
     &                 nhsv, ndtot_loc, istif, ele,
     &                 forceVector, stiffMatrix, hsv(ele,1:nhsv) )

         ! Add internal force contribution of current element to 'force'
          force(ele,index_xy_) = force(ele,index_xy_)
     &                         + forceVector
c
         ! If stiffness matrix is requested (istiff==1), add it to 'stiff'
          if( istif == 1 ) then
            ! Check stiffMatrix for NaN
            ! and replace by 1e6 to avoid error "nan found in imasem_unsym"
            ! @todo-optimize Is 1e6 a good choice?
             do ii=1,ndtot_loc
               do jj=1,ndtot_loc
                  if (isNan( stiffMatrix(ii,jj) )) then
                     stiffMatrix(ii,jj) = 1e6
                  endif
               enddo
             enddo
            ! Add stiffness contribution of current element to 'stiff'
            ! @note Purposefully using "AddTo", to keep possible previous values of 'stiff' by LS-Dyna
             stiff(ele,index_xy_,index_xy_) = 
     &                                  stiff(ele,index_xy_,index_xy_)
     &                                  + stiffMatrix
c
          endif   
c                  
       end do ELEMENTS  
c
      return
      end
