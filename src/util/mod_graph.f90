!***********************************************************************
      module mod_graph
!***********************************************************************
      use iso_fortran_env, only: r64=>real64
      implicit none
      contains
      !*****************************************************************
         recursive function is_connected(conn,v) result(res)
      !*****************************************************************
         implicit none
         integer :: &
            conn(:,:)
         integer, optional :: &
            v
         logical :: &
            res
         integer :: &
            u,&
            s
         !***
         res=.false.
         associate(first=>.not.present(v))
            if(first) then
               s=1
            else
               s=v
            endif
            associate(n=>size(conn,1))
               conn(s:,s)=-1
               do u=1,n
                  if(conn(s,u).eq.1) res=is_connected(conn,u)
               end do
               if(first) res=all([(conn(u,u),u=1,n)].eq.-1)
            end associate
         end associate
         !***
         end
      !*****************************************************************
         subroutine bipartite_matching(links,mapAB,weights,cost)
      !*****************************************************************
         implicit none
         integer :: &
            links(:,:)
         real(r64), optional :: &
            weights(:),&
            cost
         integer, allocatable, target :: &
            degsA(:)
         integer, allocatable :: &
            mapAB(:),&
            headsA(:,:),&
            mapBA(:),&
            ancmapB(:)
         real(r64), allocatable :: &
            weightsA(:,:),&
            potsA(:),&
            potsB(:),&
            min_dists(:),&
            distsA(:),&
            deltas(:)
         logical, allocatable :: &
            marksB(:),&
            marksA(:)
         integer :: &
            a,b,sink,&
            nA,nB,&
            ancA,ancB
         integer, pointer :: &
            d
         !***
         call init()
         do a=1,nA
            if(present(weights)) then
               call aug_dijkstra_map()
            else
               call aug_dfs_map()
            endif
            if(sink.eq.0) cycle
            b=sink
            do
               ancA=ancmapB(b)
               ancB=mapAB(ancA)
               mapAB(ancA)=b
               mapBA(b)=ancA
               b=ancB
               if(ancA.eq.a) exit
            end do
         end do
         !***
         contains
         !**************************************************************
            subroutine aug_dijkstra_map()
         !**************************************************************
            implicit none
            real(r64) :: sink_dist
            !***
            d=>degsA(a)
            associate(h => headsA(a,:d),&
                      c => weightsA(a,:d))
               marksA=spread(.false.,1,nA)
               marksA(a)=.true.
               distsA=spread(huge(0.d0),1,nA)
               distsA(a)=0.d0
               min_dists=spread(huge(0.d0),1,nB)
               min_dists(h)=max(0.d0,c-potsA(a)-potsB(h))
               marksB=spread(.false.,1,nB)
               ancmapB=spread(0,1,nB)
               ancmapB(h)=a
            end associate
            do
               sink_dist=minval(min_dists,mask=.not.marksB)
               sink=findloc(min_dists,sink_dist,1,&
                            mask=.not.marksB.and.&
                                 min_dists.ne.huge(0.d0))
               if(sink.eq.0) return
               marksB(sink)=.true.
               ancA=mapBA(sink)
               if(ancA.eq.0) exit
               marksA(ancA)=.true.
               distsA(ancA)=sink_dist
               d=>degsA(ancA)
               associate(h => headsA(ancA,:d),&
                         c => weightsA(ancA,:d))
                  deltas=max(0.d0,sink_dist+c-potsA(ancA)-potsB(h))
                  where(.not.marksB(h).and.min_dists(h).gt.deltas)
                     min_dists(h)=deltas
                     ancmapB(h)=ancA
                  end where
               end associate
            end do
            where(marksA) potsA=potsA+sink_dist-distsA
            where(marksB) potsB=potsB+min_dists-sink_dist
            if(present(cost)) cost=cost+sink_dist
            !***
            end
         !**************************************************************
            recursive subroutine aug_dfs_map(arg)
         !**************************************************************
            implicit none
            integer, optional :: &
               arg
            integer :: &
               ancA,&
               i
            !***
            if(present(arg)) then
               ancA=arg
            else
               ancA=a
               sink=0
               marksB=spread(.false.,1,nB)
               ancmapB=spread(0,1,nB)
            endif
            if(ancA.eq.0) then
               sink=b
               return
            endif
            d=>degsA(ancA)
            associate(h => headsA(ancA,:d))
               do i=1,d
                  b=h(i)
                  associate(mark => marksB(b),&
                            anc  => ancmapB(b),&
                            map  => mapBA(b))
                     if(mark) cycle
                     mark=.true.
                     anc=ancA
                     call aug_dfs_map(map)
                     if(sink.ne.0) return
                  end associate
               end do
            end associate
            !***
            end
         !**************************************************************
            subroutine init()
         !**************************************************************
            implicit none
            integer :: &
               i
            !***
            associate(nedges=>size(links,2))
               if(nedges.eq.0) then
                  nA=0
                  nB=0
               else
                  nA=maxval(links(1,:))
                  nB=maxval(links(2,:))
               endif
               mapAB=spread(0,1,nA)
               mapBA=spread(0,1,nB)
               marksB=spread(.false.,1,nB)
               ancmapB=spread(0,1,nB)
               if(present(weights)) then
                  potsA=spread(0.d0,1,nA)
                  potsB=spread(0.d0,1,nB)
                  marksA=spread(.false.,1,nA)
                  distsA=spread(0.d0,1,nA)
                  min_dists=spread(0.d0,1,nB)
                  deltas=spread(0.d0,1,nB)
               endif
               if(present(cost)) cost=0.d0
               allocate(degsA(nA))
               degsA=0
               do i=1,nedges
                  associate(a=>links(1,i))
                     degsA(a)=degsA(a)+1
                  end associate
               end do
               associate(mdegsA=>maxval(degsA))
                  allocate(headsA(nA,mdegsA))
                  if(present(weights)) &
                     allocate(weightsA(nA,mdegsA))
                  degsA=0
                  do i=1,nedges
                     associate(a => links(1,i),&
                               b => links(2,i))
                        d=>degsA(a)
                        d=d+1
                        headsA(a,d)=b
                        if(present(weights)) &
                           weightsA(a,d)=weights(i)
                     end associate
                  end do
               end associate
            end associate
            !***
            end
         !***
         end
      !***
      end
