c This subroutine loadud should replace the subroutine loadud in dyn21.f

      subroutine loadud(fnod,dt1,time,ires,x,d,v,a,ixs,
     . numels,ixb,numelb,idrflg,tfail,isf,p,npc,fval,iob,iadd64,numelh,
     . ixh,nhex_del,nbeam_del,nshell_del,hexarray,hextim,bemarray,
     . bemtim,shlarray,shltim,parm,numnp,fnodr,dr,vr,ndof,xmst,xmsr)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     input arrays
c
c     fnod - global nodal forces
c     fnodr- global nodal moment
c     dt1 - current time step size
c     time - current problem time
c     ires - restart flag, ( 0=solution phase     )
c                          (-n=input n parameters )
c                          ( 2=restart            )
c                          ( 3=write data into dump file  )
c                          ( 4=read data from restart file)
c
c     When data is read, DUMMY arrays are passed in the call.
c     Data should be read into a local common block which is
c     written into the restart database.
c
c
c     d - displacements
c     v - velocities
c     a - acceleations
c     dr- rotational displacements
c     vr- rotational velocities
c     xmst(numnp) =reciprocal of nodal translational masses in solution phase
c     xmsr(numnp) =reciprocal of nodal rotational masses in solution phase.
c                  this array is defined if and only if ndof=6
c
c     ixs - shell element connectivities (ixs(1,*)=part ID)
c                                        (ixs(2,*)=node 1)
c                                        (ixs(3,*)=node 2)
c                                        (ixs(4,*)=node 3)
c                                        (ixs(5,*)=node 4)
c
c     ixb - beam  element connectivities (ixb(1,*)=part ID)
c                                        (ixb(2,*)=node 1)
c                                        (ixb(3,*)=node 2)
c                                        (ixb(4,*)=orientation node)
c
c     ixh - shell element connectivities (ixh(1,*)=part ID)
c                                         or
c                                        (ixh(1,*)=0 implies element deleted)
c                                        (ixh(2,*)=node 1)
c                                        (ixh(3,*)=node 2)
c                                        (ixh(4,*)=node 3)
c                                        (ixh(5,*)=node 4)
c                                        (ixh(6,*)=node 5)
c                                        (ixh(7,*)=node 6)
c                                        (ixh(8,*)=node 7)
c                                        (ixh(9,*)=node 8)
c     numnp  - number of nodal points
c     numels - number of shell elements
c     numelb - number of beam  elements
c     numelh - number of solid elements
c     isf    - shell element failure flag (1=on)
c     tfail  - shell element failure time (eq.0:okay)
c                                         (ne.0:failure time)
c
c     idrflg - nonzero if dynamic relaxation phase
c     p      - load curve data pairs (abcissa,ordinate)
c     npc    - pointer into p.  (p(npc(lc)) points to the beginning
c              of load curve ID lc.  npoints=npc(lc+1)-npc(lc)=
c              number of points in the load curve.
c     fval   - fval(lc) is the value of load curve lc at t=time
c     iob    - i/o buffer
c     ndof   - number of degrees of freedom per node in solution phase,
c            - 0 in initialization phase
c
c
c     ELEMENT DELETION FROM THE USER LOADING SUBROUTINE
c
c     to delete elements the original file needs to have a
c     *DEFINE_ELEMENT_DEATH definition for one or more elements
c     of the type to be deleted.  The deletion time can be set
c     to a value that greatly exceeds the termination time for the
c     run.  With this keyword the necessary arrays are created that
c     allows elements of the same type to deleted during the run.
c
c     nhex_del if >0 element deletion option is active for solids
c     nbem_del if >0 element deletion option is active for beams
c     nshl_del if >0 element deletion option is active for shells
c     hexarray defines time to delete solid element,
c              the value should be >time
c     hextim   we check for solid element deletion when the current
c              time is greater to or equal to hextim
c     bemarray defines time to delete beam element,
c              the value should be >time
c     bemtim   we check for beam element deletion when the current
c              time is greater to or equal to bemim
c     shlarray defines time to delete shell element,
c              the value should be >time
c     shltim   we check for shell element deletion when the current
c              time is greater to or equal to shltim
c     parm     user loading parameter if ires<0
c
      use userinterface ! additionally necessary for use of getreal8ptr
c
      include 'iounits.inc'
      include 'bigprb.inc'
      include 'txtline.inc'
c
      parameter (NPARM=1000)
c     common/usrldv/parm(NPARM)
c
      integer*8 iadd64
      real*8 x
      real*8 d,dr
      dimension a(3,*),v(3,*),d(3,*),fnod(3,*),ixs(5,*),ixb(4,*),
     . x(3,*),tfail(*),p(*),npc(*),fval(*),iob(*),ixh(9,*),
     . hexarray(*),bemarray(*),shlarray(*),parm(*),fnodr(3,*),
     . vr(3,*),dr(3,*),xmst(*),xmsr(*)
c Custom declarations
       real*8 :: x3, x4, y3, y4, theta34, n34x, n34y, f_face, area
       integer*8 pid, pid_loaded, axi, undef
c
      ! array of initial/undeformed nodal coordinates   
       real*8, dimension(2,3) :: X_coords
       real*8, CONTIGUOUS, POINTER :: dm_rots(:) ! [POINTER_CONTIGUOUS]
       call getreal8ptr('dm_rots',dm_rots)
c
c     character*80 txts,mssg
c
      ! Retrieve the load curve ID, the PID that shall be loaded, 
      ! and whether the loading is axisymmetric or plane strain
      ! @note The numbering does not use the LCID or PID, but counts in
      !       the order of occurance. Therefore, a part defined as the
      !       second *PART will have a pid of 2 even though its actual ID
      !       is 10. The same for LCID.
      ! Keyword to activate user loading:
      ! *USER_LOADING
      ! $     LCID      PID      AXI    UNDEF
      !          2        2        1        0
       LCID = INT(PARM(1))
       pid_loaded = INT(PARM(2))
       axi = INT(PARM(3))
       undef = INT(PARM(4))
c
      if (ires.lt.0) then
        n=abs(ires)
c       write(iohsp,1030)
        call prludparm(0,parm,0,0)
        mssg='reading user loading subroutine'
        if (longs) then
          do 11 i=1,n,8
          call gttxsg (txts,lcount)
          read (txts,'(8e20.0)',err=400) (parm(j),j=i,min(i+3,n))
c         write(iohsp,1040) (j,parm(j),j=i,min(i+3,n))
          call prludparm(1,parm,i,min(i+3,n))
   11     continue
        else
          do 10 i=1,n,8
          call gttxsg (txts,lcount)
          read (txts,1020,err=400) (parm(j),j=i,min(i+7,n))
c         write(iohsp,1040) (j,parm(j),j=i,min(i+7,n))
          call prludparm(1,parm,i,min(i+7,n))
   10     continue
        endif
c       write(iohsp,1050)
        call prludparm(2,parm,0,0)
        return
      endif
c
      do 20 i=1,numels
        ! Retrieve the pid of the current element to check whether
        ! it belongs to the part that shall be loaded
         pid  = int(ixs(1,i))
        ! Check the pid
         if ( pid .eq. pid_loaded ) then
           ! Get node IDs
           ! @warning Here it is hardcoded that the pressurised surface is always spanned by node 3 and node 4 of the element. This must be enforced manually for the node numbering of the loaded elements.
            ixs4i=ixs(4,i) ! node 3 on the bottom
            ixs5i=ixs(5,i) ! node 4 on the bottom
           ! Get x and y coordinates of both nodes
            if ( undef==1 ) then
                ! Retrieve undeformed nodal coordinates
                ! (could be used for non-follower load or geometrically linear loading)
                 X_coords(:,:) = reshape(
     &                   dm_rots( 
     &                           (/
     &                              ((ixs4i-1)*3+1):(ixs4i*3),
     &                              ((ixs5i-1)*3+1):(ixs5i*3)
     &                            /)
     &                           ), (/2,3/), order=(/2,1/)
     &                         )
                x3=X_coords(1,1)
                y3=X_coords(1,2)
                x4=X_coords(2,1)
                y4=X_coords(2,2)
            else
               ! Retrieve deformed/current nodal coordinates
               ! e.g. for follower-load
                x3 =x(1,ixs4i)
                y3 =x(2,ixs4i)
                x4 =x(1,ixs5i)
                y4 =x(2,ixs5i)
            endif
c
           ! output y-coordinates to check correctness of nodes
            if ( abs(fval(LCID)) < 1e-8 ) then ! init step
                 write(*,*) "y3=",y3,"; y4=",y4
            endif
c
           ! Compute normal vector
           ! @note This must consider that when nodes 3 and 4 are interchanged,
           !       that the normal vector flips
             n34x = (y4-y3) / sqrt( (x3-x4)**2 + (y3-y4)**2 )
             n34y = -(x4-x3) / sqrt( (x3-x4)**2 + (y3-y4)**2 )
           ! Compute loaded area
            if ( axi==1 ) then
                ! @note Missing 2*pi as LS-Dyna uses force per radian for axisym
                 area = 1. !* 2* 4*atan(1.) ! 2* pi
     &                  * 0.5*(x3+x4) ! average radius
     &                  * sqrt( (x3-x4)**2 + (y3-y4)**2 ) ! edge length
            else ! plane strain
                 area = sqrt( (x3-x4)**2 + (y3-y4)**2 ) ! edge length
     &                  * 1. ! assuming unit thickness
            endif
           ! Pressure-induced effective force acting on edge node3-node4
           ! Computed from area and pressure
            f_face = area ! pressurised area of the current element
     &               * fval(LCID) ! pressure at current time
c
           ! Split up pressure-induced effective force on node 3 and 4
           ! (constant pressure distribution, split 50%-50% on both nodes, therefore factor "0.5")
            fnod(1,ixs4i)=fnod(1,ixs4i) + 0.5*f_face*n34x
            fnod(2,ixs4i)=fnod(2,ixs4i) + 0.5*f_face*n34y
            fnod(1,ixs5i)=fnod(1,ixs5i) + 0.5*f_face*n34x
            fnod(2,ixs5i)=fnod(2,ixs5i) + 0.5*f_face*n34y
c
      endif ! (pid)
   20 continue
c
c
      return
  400 call termin (txts,mssg,lcount,1)
 1020 format(8e10.0)
 1030 format(//' u s e r   d e f i n e d   l o a d i n g',
     .           '   p a r a m e t e r s',/)
 1040 format(
     1 5x,'   parameter number  ',i4,'=', e15.8)
 1050 format(//)
      end
