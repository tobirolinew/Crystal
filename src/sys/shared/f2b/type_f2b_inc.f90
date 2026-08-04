      !**********************
       type site
      !*********************
         real(r64) :: &
            alpha=0.d0
         real(r64), allocatable :: &
            b_elst(:),&
            b_disp4(:),&
            b_disp6(:),&
            b_disp8(:),&
            b_disp10(:),&
            s_elst(:,:),&
            s_disp4(:,:),&
            s_disp6(:,:),&
            s_disp8(:,:),&
            s_disp10(:,:)
         integer, allocatable :: &
            lambda_elst(:,:),&
            lambda_disp4(:,:),&
            lambda_disp6(:,:),&
            lambda_disp8(:,:),&
            lambda_disp10(:,:)
       end type
      !*********************
       type site_pair
      !*********************
         real(r64) :: &
            beta=0.d0,&
            gamma=0.d0,&
            A12=0.d0,&
            delta1=0.d0,&
            delta3=0.d0,&
            delta4=0.d0,&
            delta6=0.d0,&
            delta8=0.d0,&
            delta10=0.d0,&
            C4=0.d0,&
            C6=0.d0,&
            C8=0.d0,&
            C10=0.d0
         real(r64), allocatable :: &
            b_exp(:),&
            s_exp(:,:)
         integer, allocatable :: &
            lambda_exp(:,:)
       end type
