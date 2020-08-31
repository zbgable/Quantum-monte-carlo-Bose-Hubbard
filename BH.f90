
  !*******************************************************************
  !

  !  codes for transverse Bose-Hubbard model.

  !  Author    : Ke-Hang Zhu
  !  Date      : July 16th, 2020
  !  Last Edit : Oct   16th, 2020
  !*******************************************************************


  
  
  
  !====================== START VARIABLE MODULE =======================
  MODULE my_vrbls
    IMPLICIT NONE

    !-- common parameters and variables ------------------------------
    ! THIS IS ALMOST PROJECT-INDEPENDENT 
    real(8)         ,parameter  :: tm32   = 1.d0/(2.d0**32)
    real(8)         ,parameter  :: eps    = 1.d-14      ! very small number
    logical                     :: prtflg               ! flag for write2file
    integer         ,parameter  :: Mxint  = 2147483647  ! maximum integer
    integer         ,parameter  :: Mnint  =-2147483647  ! minimum integer

    integer(8)                  :: Nsamp                ! # samples in unit 'NBlck'
    integer(8)                  :: Totsamp
    integer(8)                  :: Ntoss                ! # samples to be thrown away
    integer(8)                  :: itoss
    
    integer(8)                  :: i1=0,i2=0,i3=0,kw=0
    integer(8)                  :: N_measure
    integer(8)                  :: N_print
    integer(8)                  :: N_write

    character(10)               :: cnf_file  = 'cnf.dat'
    character(10)               :: stat_file = 'stat.dat'
    character(10)               :: print_file= 'print.txt'
    character(10)               :: green_file= 'green.txt'    ! save the data of green function

    !-- observe ----------------------------------------------------------
    ! THIS IS ALMOST PROJECT-INDEPENDENT
    integer         ,parameter  :: Max_block=4000
    integer         ,parameter  :: N_block=1000
    real            ,parameter  :: Max_cor=0.1
    logical                    :: flg_cor
    type obs                                            ! define obs type
        character(8)            :: nam                  ! the name of observables
        real(8)                 :: vnr                  ! current value (integer)
        real(8)                 :: val                  ! the accumulated value of observables
        real(8)                 :: cor                  ! correlation
        real(8)                 :: erb                  ! errobar
        real(8)                 :: blk(Max_block)       ! blocks for calculating the errorbars
        real(8)                 :: a                    !
        real(8)                 :: b                    ! obser=a*val+b
    end type


    integer         ,parameter  :: NObs_b=5             ! #basic observables
    integer         ,parameter  :: NObs_c=1             ! #composite observables
    integer         ,parameter  :: NObs=NObs_b+NObs_c
    type(Obs)                   :: Obser(NObs)
    integer(8)                  :: Npb                  ! # samples per block
    integer(8)                  :: Nbk                  ! the {Nbk}th block is the current one
    integer(8)                  :: Np                   ! the pointer

 
    !-----------------------------------------------------------------

    !-- parameters and variables -------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    real(8)         ,parameter  :: Pi=3.1415926536
    integer(1)      ,parameter  :: D = 2                ! dimensionality
    integer(8)                  :: Lx, Ly, Vol          ! the system size
    integer(8)                  :: Z                     ! partition sum

    real(8)         ,allocatable:: green(:)                 !green function
    real(8)                    :: G0
    real(8)       ,allocatable:: distance(:)                !distance in for different green function

    !calculatingg correlation time
    real(8)       ,allocatable:: gamma(:)                 !correlation function
    real(8)       ,allocatable:: M_corr(:)                 !used for calculating correlation function
    real(8)                 :: t0

    real(8)                     :: h
    real(8)                     :: t
    real(8)                     :: U
    real(8)                     :: mu
    real(8)                     :: beta
    real(8)                     :: delmu            !delmu=mu+U/2
    real(8)                     :: K                  !K=t*beta


    !-----------------------------------------------------------------

    !-- Lattice and State --------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 


    integer(1)                  :: nnb                  ! # neighboring vertices
    integer(1)                  :: nnb1                 ! nnb+1
    integer(1)                  :: nnb2                 ! nnb/2
    integer(8)      ,allocatable:: nnsite(:,:)          ! nearest neighbor site

    integer         ,allocatable:: backdir(:)           ! the backward direction

    integer(8)                  :: Ira,Masha
    logical                     :: GT_ZF                ! flag in G (T) or Z (F) space
    integer(8)                  :: delta_x, delta_y     !x (y) distance between ira and masha

    real(8)                     :: p,dE
    real(8)                     :: E
    integer(8)                  :: MP
    

    real(8)                     :: wormstep
    integer(8)                  :: step
    integer(8)                  :: test_i, test_j
    !-----------------------------------------------------------------


    !-- Random-number generator---------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    integer, parameter           :: mult=32781
    integer, parameter           :: mod2=2796203, mul2=125
    integer, parameter           :: len1=9689,    ifd1=471
    integer, parameter           :: len2=127,     ifd2=30
    integer, dimension(1:len1)   :: inxt1
    integer, dimension(1:len2)   :: inxt2
    integer, dimension(1:len1)   :: ir1
    integer, dimension(1:len2)   :: ir2
    integer                      :: ipnt1, ipnf1
    integer                      :: ipnt2, ipnf2
    integer, parameter           :: mxrn = 10000
    integer, dimension(1:mxrn)   :: irn(mxrn)

    integer                      :: Seed                   ! random-number seed
    integer                      :: nrannr                 ! random-number counter
    !-----------------------------------------------------------------

    !-- time-checking variables --------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    character( 8)           :: date
    character(10)           :: time
    character(5 )           :: zone
    integer, dimension(8)   :: tval
    double precision        :: t_prev, t_curr, t_elap
    integer                 :: h_prev, h_curr
    double precision        :: t_init, t_simu, t_meas, t_toss

    !-----------------------------------------------------------------
    ! quantum worm alg related
    type segment                                              ! define segment type
        integer(8)                 :: site                  !site
        integer(8)                 :: begin                !beginning point
        integer(8)                 :: length               !length
        integer(8)                 :: n_b                  !occupation number
        integer(8)                 :: pri_index
        integer(8)                 :: sub_index
        integer(8)                 :: kink_index
        integer(8)                 :: selfp
    end type

    type(segment)       ,allocatable:: seg(:)         ! set the maximum length of segments to be 500
    integer(8)           ,allocatable:: index(:)          ! nearest neighbor site
    integer(8)                      :: flg_used            !pointer for the used index and unused index
    integer(8)           ,allocatable:: supergrid(:)      !pointer for the used index and unused index
    integer(8)                      :: N_beta=1000000            !time slicing of imaginary axis
    real(8)                         :: omega_G=1.0/40000          !weight between G and Z space
    real(8)                         :: h_change            !diaginal part of the Hamiltonian
    real(8)                         :: N_kink ,n_k            !kink number
    real(8)                         :: nboson

    integer(8)                      :: multiple=20           !n time of Vol segments will be created
    real(8)                         :: p_delw,p_ts,p_crek,p_delk ! probabilty for calling the updates

    !-----------------------------------------------------------------
  END MODULE my_vrbls
  !======================== END VARIABLE MODULE =====================
  

  !======================== START MAIN PROGRAM ======================
  
  PROGRAM BH
    use my_vrbls
    implicit none
    integer :: isamp,r_cnf_stat,i_each,N_each,iblck

    read *, r_cnf_stat
    read *, Lx
    read *, t
    read *, U
    read *, mu
    read *, beta
    read *, p_delw
    read *, p_ts
    read *, p_crek
    read *, p_delk
    read *, Ntoss
    read *, Nsamp
    read *, N_print
    read *, N_write
    read *, N_each
    read *, N_measure
    read *, Seed

    print *, 'r_cnf_stat, Lx, t, U, mu, beta, p_delw ,p_ts ,p_crek ,p_delk, Ntoss, Nsamp, N_print, N_write, N_each, N_measure, Seed'
!    read  *,  r_cnf_stat, Lx, Ntoss, Nsamp, N_print, N_write, N_each, N_measure, t, h, beta, Seed
    print *,  r_cnf_stat, Lx, t, U, mu, beta, p_delw ,p_ts ,p_crek ,p_delk, Ntoss, Nsamp, N_print, N_write, N_each, N_measure, Seed

    
    !--- Initialization ----------------------------------------------    
    call set_time_elapse
    
    call initialize
    if(r_cnf_stat==0) then
        print*,"create a new simulation:"
        call init_stat                          !initiate set all observables to be 0
        call system("rm "//trim(print_file))
    else if(r_cnf_stat==1) then
        print*,"start a simulation on an old configuration:"
        call read_cnf()
        call init_stat
        call system("rm "//trim(print_file))
    else if(r_cnf_stat==2) then
        print*,"resume an old simulation:"
        call read_cnf()
        call init_stat
        call read_stat()
        call printing(6)
    else
        stop "r_cnf_stat should be 0, 1 or 2"
    end if

    call time_elapse                        !calculate the running time
    t_init = t_elap
    write(6,50) t_init
    50 format(/'        set up time:',f16.7,2x,'s')

    !print*, gamma
    !print*,"==================="
    !print*,"correlation time "
    !print*, t0

    !--- Thermialization ---------------------------------------------
    if(r_cnf_stat/=2) then
        prtflg=.False.
        t_toss=0.d0
        do isamp = 1, Ntoss             !Running MCMC NToss*N_each times
            do i_each = 1, N_each
                call markov
                call measure            !measure each time
            enddo

            i2=i2+1
            if(i2>=N_print) then
                call time_elapse;   t_toss=t_toss+t_elap        !calculate the running time
                itoss=itoss+N_print
                call writejjfile()
                call printing(6)
                open(1,file=trim(print_file),access="append")
                call printing(1)
                !open(1,file=trim(green_file),access="append")       !save green function 2020/7/31
                !call write_green
                close(1)
                call check
                i2=0
            end if
            i3=i3+1
            if(i3>=N_write) then
                print*, "write cnf. and stat."
                call write_cnf()
                call write_stat()
                print*, "write done"
                i3=0
            end if
        enddo

        call time_elapse                    !calculate the running time
        t_toss = t_toss+t_elap
        write(6,51) t_toss
        51 format(/'thermalization time:',f16.7,2x,'s')
        t_simu   = 0.d0;   t_meas = 0.d0
        wormstep=0; step=0
        !green=0         !set the green function to be 0 after thermalization
        !Z=1             !set the partition sum to be 0 after thermalization
        !G0=1

       
    else
        print*,"no need to thermalizate."
    end if
    print*, seg
    !--- Simulation --------------------------------------------------
    prtflg=.True.                  !flag for write2file
    do isamp = 1, Nsamp
        do i_each = 1, N_each
            call markov
            call measure
                
                    !if (mod(isamp*i_each,5000000)==0) then
                       ! print*,getN()/Vol
                   ! end if
                
!			call coll_data()
           
        end do
        
        i2=i2+1
        if(i2>=N_print) then
            call cal_ComObs()
            flg_cor=cal_cor_dev()    !if two sets of measurement has very big correlation, the flag would be false
            call time_elapse;    t_simu = t_simu+t_elap
            call writejjfile()
            call printing(6)
            open(1,file=trim(print_file),access="append")
            call printing(1)
            !open(1,file=trim(green_file),access="append")
            !call write_green                                    !save green function 2020/7/31
            close(1)
            call check
            print*, getN()
            i2=0
        end if
        i3=i3+1
        if(i3>=N_write) then
            print*,"write cnf. and stat."
            call write_cnf()
            call write_stat()
            print*,"write done"
            i3=0
        end if
    end do

    !call correlation
    !print*, gamma
    !print*,"==================="
    !print*,"correlation time "
    !print*, t0
    stop

CONTAINS

    !--- PROJECT-DEPENDENT ---------------------------------------------
  !============== START INITIALIZATION ======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE initialize
    implicit none

    !-- order should be changed --------------------------------------

    call allo_aux
    call def_latt

    call def_distance
    call set_RNG
    call init_cnf
    call set_Obs
    !Call set_segs
    call tst_and_prt

    return
  END SUBROUTINE initialize
 
  !==============Test and Print ======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE tst_and_prt
    implicit none

    !-- Test and Print -----------------------------------------------

    write(6,*)
    write(6,*) "Bose Hubbard model"
    write(6,40) Lx,Ly
    40 format('on ',i4,1x,"*",1x,i4,2x,'Square lattice')
    write(6,*)
    write(6,*) "H=-t*\sum_{<i,j>}b^+j bi+c.c.+U/2*sum_in_i(n_i-1)-mu*sum_i n_i"
    write(6,41) t,U,mu,beta
    41 format(',where t=',f12.8,' U=',f12.8,' mu=',f12.8,' beta=',f12.8)
    write(6,*)

    write(6,42) Nsamp,N_each
    42 format(' Will simulate      ',i10,'*',i10,'*',2x,'steps ')

    write(6,44) Ntoss,N_each,Vol
    44 format(' Throw away         ',i10,'*',i10,'*',2x,'steps')

    return
  END SUBROUTINE tst_and_prt


  !==============initialize configuration ======================================
  !! THIS IS PROJECT_DEPENDENT
  SUBROUTINE init_cnf
    implicit none
    integer(8) :: i,Vc,l

    !allocate(green(2*(Lx/2)**2))
    allocate(M_corr(Nsamp))
    allocate(gamma(Nsamp-1))

    allocate(seg(Vol*multiple))
    allocate(index(Vol*multiple))
    allocate(supergrid(Vol))

    supergrid=0
    index=0
!print*,index
    call set_segs


    l=0
    Ira=0;  Masha=0
    GT_ZF=.False.
    Z=0
    !green=0             !green function
    !G0=0                ! green (0)
    M_corr=0            !collected for claculating correction time
    N_kink=0.d0
    n_k=0.d0
    h_change=0.d0
    nboson=0.d0
    return
  END SUBROUTINE

!==============================================
!initialize the segments and index and pointer
    SUBROUTINE set_segs
      implicit none
      integer(8) :: ii,jj=1
!print*,seg

      do ii=1,Vol
!print*,'allocate done',jj
!jj=jj+1
        call init_segs_loop(ii)
        index(ii)=ii
        call griding(ii)
      end do

      do ii=Vol+1,multiple*Vol
        call init_segs(ii)
        index(ii)=ii
      end do

      flg_used=Vol

      return
    END SUBROUTINE


    SUBROUTINE init_segs_loop(i)
      implicit none
      integer(8) :: i

      seg(i)%site=i
      seg(i)%begin=0
      seg(i)%length=N_beta
      seg(i)%n_b=0
      seg(i)%pri_index=i
      seg(i)%sub_index=i
      seg(i)%kink_index=0
      seg(i)%selfp=i
      return
    END SUBROUTINE


    SUBROUTINE init_segs(i)
      implicit none
      integer(8) :: i

      seg(i)%site=0
      seg(i)%begin=0
      seg(i)%length=0
      seg(i)%n_b=0
      seg(i)%pri_index=0
      seg(i)%sub_index=0
      seg(i)%kink_index=0
      seg(i)%selfp=0
      return
    END SUBROUTINE


  !===============allocate auxiliary arrays ======================================
  !! THIS IS PROJECT_DEPENDENT
  SUBROUTINE allo_aux
    implicit none
    if(D==1) then
        Vol=Lx
    else if(D==2) then
        Ly=Lx;  Vol=Lx*Ly
    else
        stop "D should be 1 or 2"
    end if

    if(D==1) then
        allocate(nnsite(Vol,2))
        allocate(backdir(2))


    else if(D==2) then
        allocate(nnsite(Vol,4))
        allocate(backdir(4))
        allocate(distance(2*(Lx/2)**2))
        
    else
        stop "D should be 1 or 2"
    end if

    delmu=mu+U/2
    return
  END SUBROUTINE
  
  !==============define flipping probability =========================
  !! THIS IS PROJECT-DEPENDENT 

  !==============definite Lattice ====================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE def_latt
    implicit none
    integer(8)  :: i,j,ii
    if((Lx/2)*2/=Lx) then
        write(6,*) 'L be even?';   stop
    end if
    if(D==2) then                                       ! y
        nnb=4                                           ! ^
        do i=1,nnb                                      ! |  13 14 15 16
            backdir(i)=mod(i+1,4)+1                     ! |   9 10 11 12
        end do                                          ! |   5  6  7  8
        do j=1,Vol                                      ! |   1  2  3  4
            nnsite(j,1)=(j-1)/Lx*Lx+mod(j,Lx)+1         ! |---------------> x
            nnsite(j,2)=mod(j-Lx+Vol-1,Vol)+1           !
            nnsite(j,3)=(j-1)/Lx*Lx+mod(j-2+Lx,Lx)+1    !        4
            nnsite(j,4)=mod(j+Lx-1,Vol)+1               !        |
        end do                                          !   3 ---|--- 1
    else if(D==1) then                                   !        |
        nnb=2                                          !        2
        backdir(1)=2
        backdir(2)=1
        do j=1,Vol                                  ! periodic boundary condition
            if(j+1>Vol) then
                nnsite(j,1)=j+1-Vol
            else
                nnsite(j,1)=j+1
            end if
            if(j-1<1) then
                nnsite(j,2)=j-1+Vol
            else
                nnsite(j,2)=j-1
            end if
        end do
    else
        stop 'wrong! D must be 1 or 2'
    end if
    nnb1=nnb+1; nnb2=nnb/2
!    print*,j,nnsite(j,:)

!do ii=1, Vol
!print*, nnsite(ii,:)
!end do

    return
  END SUBROUTINE def_latt

SUBROUTINE def_distance
    implicit none
    integer(8)  :: i,j,k,l
    !real(8):: k,l
    distance=0
    do i=1,Lx+1
        do j=1,Lx+1
            k=i-(Lx/2+1)
            l=j-(Lx/2+1)
            if (k==0 .and. l==0) then
                G0=0
            else if (k**2+l**2>=Lx**2/4) then
                distance(k**2+l**2)=0
            else
                distance(k**2+l**2)=distance(k**2+l**2)+1
            end if
        end do
    end do
    
    return
  END SUBROUTINE def_distance
  !===================================================================

  !==========initialize state ======================
  !! THIS IS PROJECT-DEPENDENT

  SUBROUTINE init_stat           !initiate set all observables to be 0
    implicit none
    integer(8) :: Vc,Vcn,Rseg,Vseg
    integer :: dir

    kw=0
    wormstep=0.d0;  step=0

    delta_x=0
    delta_y=0
    t0=0

    return
  END SUBROUTINE init_stat

    FUNCTION getN()
    implicit none
    integer(8) :: i,ind
    real(8) :: leng,getN
        getN=0.d0
        do i=1, flg_used
            ind=index(i)

            leng=seg(ind)%length*1.0
            getN=getN+leng*seg(ind)%n_b/N_beta
        end do
    return
  END FUNCTION

  FUNCTION getE()
  implicit none
  integer(8) :: ii,ind
  real(8) :: leng,getE
      getE=0.d0
      do ii=1, flg_used
          ind=index(ii)

        leng=seg(ind)%length*1.0
        getE=getE+leng/N_beta*(U/2.0*seg(ind)%n_b**2-delmu*seg(ind)%n_b)
      end do
  return
END FUNCTION
  
  !==================== END INITIALIZATION =====================
  
  !====================== START MARKOV =========================
  
  
  !! THIS IS PROJECT-DEPENDENT
  SUBROUTINE markov
    implicit none
    integer :: istep
    do istep= 1, N_measure
        
        call Quaworm
       ! call check                 !check the program

    enddo
    return
  END SUBROUTINE markov


  !==============Worm ==========================================
  !! THIS IS PROJECT-DEPENDENT

  SUBROUTINE Quaworm           ! module for quantum Worm alg
    implicit none
    real(8) :: pt

    if (GT_ZF) then                      !if the configuration is in the G space
        pt=rn()
        if (pt<p_delw) then
           ! print*,'==========================','delete ira-masha update','=========================='

            Call Del_IM

            !if (.not. GT_ZF) then
              !   print*,'==========================','delete ira-masha succeeded','=========================='
            !end if
        else if (pt<p_delw+p_ts) then
           ! print*,'==========================','time shift update','=========================='

            Call Time_shift

        else if (pt<p_delw+p_ts+p_crek) then
            !print*,'==========================','create kink update','=========================='

            Call Cre_kink

        else if (pt<p_delw+p_ts+p_crek+p_delk) then
             !print*,'==========================','delete kink update','=========================='

            Call Del_kink
        end if

    else    !if the configuration is in the Z space

        ! print*,'==========================','create ira-masha update','=========================='
        Call Create_IM

    end if

    return
  END SUBROUTINE

!====================================================================================================================
    SUBROUTINE Create_IM           ! module for creating ira-masha pair
    implicit none
    integer(8) :: I,ind
    integer(8) :: nul=0 ! transfer fro Integer 4 to Interger 8
    integer(8) :: t_ira,t_masha,tt1,tt2
    integer(8) :: tm1,tm2,len
    integer(8) :: I_index, M_index,Subse
    integer(8) :: leng1,leng2,leng0
    integer(8) :: old_occ_num
    real(8) :: P_crew
    logical  :: is_loop

    I=dint(rn()*flg_used)+1                !randomly select one segment
    ind=index(I)
    Subse=seg(ind)%sub_index

    tm1=seg(ind)%begin
    len=seg(ind)%length
    tm2=seg(ind)%begin+len
    old_occ_num=seg(ind)%n_b
    is_loop=(len==N_beta)

    if (.not. is_loop) then                !if the segment is not a loop
            t_ira=dint(rn()*(len-1))+1+tm1  !select two time points
            t_masha=dint(rn()*(len-1))+1+tm1
            tt1=imatime_periodic(t_ira)      !satisfy the periodic boundary condition
            tt2=imatime_periodic(t_masha)
            
            ira=flg_used+1          !apply two pointers for ira and masha segment
            masha=flg_used+2
            I_index=index(ira)
            M_index=index(masha)

            if (t_ira>t_masha) then             !--masha======ira-
                P_crew=p_crew_noloop(t_ira,t_masha,old_occ_num,1,tm1,tm2)
                leng1=tm2-t_ira
                leng2=t_ira-t_masha
                leng0=t_masha-tm1
                if(leng0==0 .or.leng1==0 .or.leng2==0 ) return

                if (rn()<P_crew) then

                    Call cre_seg(seg(ind)%site,tt2,leng2,old_occ_num+1,index(I),I_index,nul,masha,M_index)
                    Call cre_seg(seg(ind)%site,tt1,leng1,old_occ_num,M_index,seg(ind)%sub_index,nul,ira,I_index)

                    seg(ind)%length=leng0         !update the original segment
                    seg(ind)%sub_index=M_index
                    seg(Subse)%pri_index=I_index
                    Call griding(ind)       ! supergriding
                    flg_used=flg_used+2     !set the pointer
                    Call update_energy_occupa_num(old_occ_num,1,seg(M_index)%length)
                    GT_ZF=.true.
                end if

            else if (t_ira<t_masha) then        !===ira-----masha==
                if (old_occ_num==0) return   ! if the occupation number <0,reject
                P_crew=p_crew_noloop(t_ira,t_masha,old_occ_num,-1,tm1,tm2)
                leng2=tm2-t_masha
                leng1=t_masha-t_ira
                leng0=t_ira-tm1
                if(leng0==0 .or.leng1==0 .or.leng2==0) return

                if (rn()<P_crew) then
                        Call cre_seg(seg(ind)%site,tt2,leng2,old_occ_num,I_index,seg(ind)%sub_index,nul,masha,M_index)
                        Call cre_seg(seg(ind)%site,tt1,leng1,old_occ_num-1,ind,M_index,nul,ira,I_index)
                        seg(ind)%length=leng0       !update the original segment
                        seg(ind)%sub_index=I_index
                        seg(Subse)%pri_index=M_index
                        ! supergriding
                        Call griding(ind)
                        flg_used=flg_used+2 !set the pointer
                        !measurement
                        Call update_energy_occupa_num(old_occ_num,-1,seg(I_index)%length)
                        GT_ZF=.true.
                end if
            end if
!=============================================!if the segment is a loop
    else                                !if the segment is a loop
        t_ira=dint(rn()*len)
        t_masha=dint(rn()*len)
        ira=flg_used+1

        masha=I
        I_index=index(ira)
        M_index=index(masha)
        if (t_ira>t_masha) then
            leng0=t_ira-t_masha
           if (rn()<1.0/2) then      ! choose at random counterclockwise or clockwise
                P_crew=p_crew_loop(leng0,old_occ_num,1)
                if (rn()<P_crew) then

                    Call cre_seg(seg(ind)%site,t_masha,leng0,old_occ_num+1,I_index,I_index,nul,masha,M_index)
                    Call cre_seg(seg(ind)%site,t_ira,N_beta-leng0,old_occ_num,M_index,M_index,nul,ira,I_index)
                    
                    flg_used=flg_used+1 !set the pointer

                    Call update_energy_occupa_num(old_occ_num,1,seg(M_index)%length)
                    GT_ZF=.true.
                end if
            else                    ! choose at random counterclockwise or clockwise
                if (old_occ_num==0) return
                    P_crew=p_crew_loop(N_beta-leng0,old_occ_num,-1)
                    if (rn()<P_crew) then

                        Call cre_seg(seg(ind)%site,t_masha,leng0,old_occ_num,I_index,I_index,nul,masha,M_index)
                        Call cre_seg(seg(ind)%site,t_ira,N_beta-leng0,old_occ_num-1,M_index,M_index,nul,ira,I_index)
                        flg_used=flg_used+1 !set the pointer
                        Call update_energy_occupa_num(old_occ_num,-1,seg(I_index)%length)
                        GT_ZF=.true.
                    end if

            end if

        else if (t_ira<t_masha) then
            leng0=t_masha-t_ira
            if (rn()<1.0/2) then      ! choose at random counterclockwise or clockwise
                P_crew=p_crew_loop(N_beta-leng0,old_occ_num,1)
              
                if (rn()<P_crew) then
                    Call cre_seg(seg(ind)%site,t_masha,N_beta-leng0,old_occ_num+1,I_index,I_index,nul,masha,M_index)
                    Call cre_seg(seg(ind)%site,t_ira,leng0,old_occ_num,M_index,M_index,nul,ira,I_index)
                    flg_used=flg_used+1 !set the pointer
                    Call update_energy_occupa_num(old_occ_num,1,seg(M_index)%length)
                    GT_ZF=.true.
                end if
            else                    ! choose at random counterclockwise or clockwise
                if (old_occ_num==0) return
                    P_crew=p_crew_loop(leng0,old_occ_num,-1)
                if (rn()<P_crew) then

                    Call cre_seg(seg(ind)%site,t_masha,N_beta-leng0,old_occ_num,I_index,I_index,nul,masha,M_index)
                    Call cre_seg(seg(ind)%site,t_ira,leng0,old_occ_num-1,M_index,M_index,nul,ira,I_index)
                    flg_used=flg_used+1 !set the pointer
                    Call update_energy_occupa_num(old_occ_num,-1,seg(I_index)%length)
                    GT_ZF=.true.
                end if
            end if
        end if
   ! if (GT_ZF)then
      !  print*,'create done'
  !  end if

    end if
    END SUBROUTINE
!====================================================================================================================
 SUBROUTINE Del_IM               ! module for deleting ira-masha pair
   implicit none
   integer(8) :: I,ind
   integer(8) :: nul=0
   integer(8) :: t1,t2
   integer(8) :: tm1,tm2,len,leng_tot
   integer(8) :: I_index, M_index,Prior,Subse
   integer(8) :: tira,tmasha
   real(8) :: P_delw_acc
   integer(8) :: old_occ_num

        I_index=index(ira)
        M_index=index(masha)
        tira=seg(I_index)%begin
        tmasha=seg(M_index)%begin
    
    if (seg(I_index)%sub_index==M_index .and. seg(I_index)%pri_index==M_index) then ! if Ira and masha are the only two segments composing the world line
        if (rn()<1.0/2) then     ! choose with equal prob to be n-1
            P_delw_acc=p_delw_loop(seg(M_index)%length,seg(M_index)%n_b,-1)

            if (rn()<P_delw_acc) then
                ! create a loop sgement
                ind=index(flg_used+1)
                flg_used=flg_used+1
                Call cre_seg(seg(I_index)%site,nul,N_beta,seg(I_index)%n_b,ind,ind,nul,flg_used,ind)

                Call update_energy_occupa_num(seg(M_index)%n_b,-1,seg(M_index)%length)
                
                Call delete_seg(ira)   !delete a segment ira
                Call delete_seg(masha) !delete a segment masha
                GT_ZF=.false.
                !print*, 'delete done'
            end if

         else           ! choose with equal prob to be n
           P_delw_acc=p_delw_loop(seg(I_index)%length,seg(I_index)%n_b,1)
            if (rn()<P_delw_acc) then
                ! create a loop sgement
                ind=index(flg_used+1)
                flg_used=flg_used+1
                Call cre_seg(seg(I_index)%site,nul,N_beta,seg(I_index)%n_b+1,ind,ind,nul,flg_used,ind)
                Call update_energy_occupa_num(seg(I_index)%n_b,1,seg(I_index)%length)

                Call delete_seg(ira)   !delete a segment ira
                Call delete_seg(masha) !delete a segment masha
                GT_ZF=.false.
            end if
         end if
!===============================================! if Ira and masha are not on a loop
     else if (seg(I_index)%sub_index==M_index) then

            Prior=seg(I_index)%pri_index
            Subse=seg(M_index)%sub_index
            old_occ_num=seg(I_index)%n_b
            leng_tot=seg(Prior)%length+seg(I_index)%length+seg(M_index)%length

            tira=get_starting_point(I_index)  !in case ira move through the beta
            P_delw_acc=p_delw_noloop(leng_tot,seg(I_index)%length,old_occ_num,+1)

            if (rn()<P_delw_acc) then
                seg(Prior)%length=leng_tot     ! update the Prior segment
                seg(Prior)%sub_index=Subse
                seg(Subse)%pri_index=Prior
                Call griding(Prior)  !supergrid
                Call update_energy_occupa_num(old_occ_num,1,seg(I_index)%length)

                Call delete_seg(ira)   !delete a segment ira
                Call delete_seg(masha) !delete a segment masha
                GT_ZF=.false.
            end if
     else if  (seg(M_index)%sub_index==I_index) then !t(ira)>t(masha)
            Prior=seg(M_index)%pri_index
            Subse=seg(I_index)%sub_index
            old_occ_num=seg(M_index)%n_b
            leng_tot=seg(Prior)%length+seg(I_index)%length+seg(M_index)%length

            tmasha=get_starting_point(M_index)  !in case masha move through the beta
            P_delw_acc=p_delw_noloop(leng_tot,seg(M_index)%length,old_occ_num,-1)

            if (rn()<P_delw_acc) then
                seg(Prior)%length=leng_tot   ! update the Prior segment
                seg(Prior)%sub_index=Subse
                seg(Subse)%pri_index=Prior
                Call griding(Prior)  !supergrid
                Call update_energy_occupa_num(old_occ_num,-1,seg(M_index)%length)

                Call delete_seg(ira)   !delete a segment ira
                Call delete_seg(masha) !delete a segment masha
                GT_ZF=.false.
            end if
     end if

if (.not. GT_ZF) then
    ira=0
    masha=0
end if
   return
 END SUBROUTINE

!====================================================================================================================
SUBROUTINE Time_shift               ! module for time shift ira
  implicit none
  integer(8) :: t1,tira
  integer(8) :: tmin,tmax,len,length_prior,length_cur
  integer(8) :: I_index,Prior
  integer(8) :: old_occ_num
  real(8)   :: P_tshift

    I_index=index(ira)
    Prior=seg(I_index)%pri_index
    tira=seg(I_index)%begin
    old_occ_num=seg(I_index)%n_b

    tmin=get_starting_point(Prior)
    tmax=seg(I_index)%begin+seg(I_index)%length
    len=tmax-tmin

    t1=dint(rn()*(len-1))+1         !pick up a point between tmax and tmin
    length_prior=t1
    length_cur=len-t1
    t1=t1+tmin

    P_tshift=prob_ts(old_occ_num,1,t1-tira)

    if (rn()<P_tshift) then !accept or not
        Call move_begin_end(I_index,Prior,t1,length_cur,length_prior)
        Call update_energy_occupa_num(old_occ_num,1,t1-tira)
    end if

  return
END SUBROUTINE
!====================================================================================================================
SUBROUTINE Cre_kink              ! module for creating a kink
  implicit none
  integer(8) :: dir,t_kink
  integer(8) :: site_ira,nn_of_ira,t_ira
  integer(8) :: tmin,tmax,len,old_occ_adja
  integer(8) :: I_index,Prior,adja
  real(8) :: P_crekink=0,ori,check_ind
  logical  :: backward,is_loop
    !print*,'maxtime &mintime error'
    I_index=index(ira)
    t_ira=seg(I_index)%begin     ! time of ira
    Prior=seg(I_index)%pri_index

    site_ira=seg(I_index)%site                  ! site of ira locates
    nn_of_ira=get_nearest_neighbor(site_ira)    ! randomly choose one of the n.n. sites of ira
    adja=adjacent(t_ira,nn_of_ira,I_index)      !find the adjacent segment of ira
    old_occ_adja=seg(adja)%n_b

    is_loop=(seg(adja)%length==N_beta)
    backward=(rn()<1.0/2)  ! determine with equal prob to shift ira backwards or forward
        ! ira is not the begin or the end of the adjacent segment
    if (t_ira==seg(adja)%begin .or. (t_ira==imatime_periodic(seg(adja)%begin+seg(adja)%length))) return
    if (backward) then      ! determine with equal prob to shift ira backwards
              !==========----adja-----        ---kink===new===---Ira---
              !===Prior======Ira---------     ===kink----------------
         Call suggest_to_create_kink_backward(I_index,Prior,adja,is_loop)
    else                   ! determine with equal prob to shift ira forward
        if (seg(adja)%n_b==0) return ! need to make sure that adja segment has occupation number >0, or ira would destroy it
                                    !=======adja==========     =adja=====Ira----kink====New====
                                    !===Prior======Ira------    ====Prior=======kink----------
         Call suggest_to_create_kink_forward(I_index,Prior,adja,is_loop)
    end if
  return
END SUBROUTINE

!====================================================================================================================
SUBROUTINE Del_kink              ! module for deleting a kink
  implicit none
  integer(8) :: site_ira,nnira,tira,I_index
  integer(8) :: tmin,tmax,tm1,tm2,t_kink
  integer(8) :: Prior, Subse, Connec,Pre_connec,Pre_pri
  integer(8) :: adja
  real(8) :: P_delkink,check_ii
  logical  :: backward,is_loop

    I_index=index(ira)
    site_ira=seg(I_index)%site      ! site of ira locates
    tira=seg(I_index)%begin     ! time of ira

    Prior=seg(I_index)%pri_index
    Subse=seg(I_index)%sub_index

    is_loop=(Prior==Subse)  !judge wgether this is a loop after deletion
    backward=rn()<1.0/2  ! determine with equal prob to shift ira backwards or forward

    if (backward) then      ! determine with equal prob to shift ira backwards
            Pre_pri=seg(Prior)%pri_index               !-----Pre_pri----==Prior==.ira.---
            Connec=seg(Prior)%kink_index               !===pre_connce===---Connec---------
           if (seg(Pre_pri)%n_b/=seg(I_index)%n_b) return
           if (Connec==0) return    ! take care of the special case where there is no kink on prior segment(masha)
           if (seg(Prior)%length>=seg(Connec)%length) return
                Pre_connec=seg(Connec)%pri_index
            !====================calculate the acceptance ratio================
            tmin=get_tmin(Connec,Pre_connec,I_index,Prior)
            tmax=tira
            t_kink=get_tkink(I_index,Prior,backward)
if (tmin>=tmax) then
    print*,'maxtime &mintime error',tmin,tmax
    check_ii=index(100000000)
end if
if (seg(Prior)%n_b==0) then
    print*,'occ error in prior seg'
    print*,seg(Prior)
    print*,seg(I_index)
    check_ii=index(100000000)
end if
            P_delkink=p_delk_acc(tmin,tmax,t_kink,seg(Pre_connec)%n_b,seg(Prior)%n_b,backward,I_index,Pre_connec)
            if (rn()<P_delkink) then
                n_k=n_k-1
                Call update_energy_occupa_num(seg(Prior)%n_b,-1,seg(Prior)%length)
                Call update_energy_occupa_num(seg(Connec)%n_b,1,seg(Prior)%length)
                Call update_del_kink(backward,Connec,Pre_pri,Pre_connec,Prior,I_index,is_loop)
                ira=seg(Connec)%selfp           !update the ira segment
            end if
        else                  ! determine with equal prob to shift ira forward
            if (seg(Prior)%n_b/=seg(Subse)%n_b) return !only then can we delete a kink
            Connec=seg(Subse)%kink_index            !===Prior===.ira.----===Subse===
                                                    !=======Pre_connec===----Connec----

                if (Connec==0) return      ! take care of the special case where there is no kink on prior segment
                Pre_connec=seg(Connec)%pri_index
                if (seg(I_index)%length>=seg(Pre_connec)%length) return
                adja=adjacent(tira,seg(Connec)%site,I_index)
                if (Pre_connec/=adja) return
                
                !==========================!calculate the acceptance ratio
                    tmax=get_tmax(Pre_Connec,Connec,I_index,Subse)
                    tmin=tira
if (tmin>=tmax) then
    print*,'maxtime &mintime error',tmin,tmax
    check_ii=index(100000000)
end if
                    t_kink=get_tkink(I_index,Subse,backward)
                    P_delkink= p_delk_acc(tmin,tmax,t_kink,seg(Pre_connec)%n_b,seg(Prior)%n_b,backward,I_index,Pre_connec)
                if (rn()<P_delkink) then
                    n_k=n_k-1  !measurement
                    Call update_energy_occupa_num(seg(I_index)%n_b,1,seg(I_index)%length)
                    Call update_energy_occupa_num(seg(Pre_Connec)%n_b,-1,seg(I_index)%length)
                    Call update_del_kink(backward,Connec,Prior,Pre_connec,Subse,I_index,is_loop)
                    ira=seg(Connec)%selfp           !update the ira segment
                end if
        end if

  return
END SUBROUTINE

FUNCTION get_tmin(adja,pre_adja,I_index,Prior)
  implicit none
  integer(8) :: get_tmin,t_adja,t_ira
  integer(8) :: adja,pre_adja,I_index,Prior

        t_ira=seg(I_index)%begin

            if (seg(adja)%begin<=t_ira .and. seg(adja)%begin+seg(adja)%length>=t_ira) then
                t_adja=seg(adja)%begin
            else if (seg(adja)%begin>t_ira) then
                t_adja=seg(adja)%begin-N_beta
            end if
            get_tmin=MAX(get_starting_point(Prior),t_adja-seg(pre_adja)%length)
  return
END FUNCTION

FUNCTION get_tmax(adja,after_adja,I_index,Subse)
  implicit none
  integer(8) :: get_tmax,t_adja,t_ira
  integer(8) :: adja,after_adja,I_index,Subse

        t_ira=seg(I_index)%begin

        if (seg(adja)%begin<=t_ira .and. seg(adja)%begin+seg(adja)%length>=t_ira) then
            t_adja=seg(adja)%begin
        else if (seg(adja)%begin>t_ira) then
            t_adja=seg(adja)%begin-N_beta
        end if
        get_tmax=MIN(t_adja+seg(adja)%length+seg(after_adja)%length,t_ira+seg(I_index)%length+seg(Subse)%length)

  return
END FUNCTION

FUNCTION get_tkink(I_index,kink_index,backward)
  implicit none
  integer(8) :: I_index,kink_index
  integer(8) :: get_tkink,t_ira
  logical :: backward

    t_ira=seg(I_index)%begin

    if (backward) then
        get_tkink=get_starting_point(kink_index)
    else
        get_tkink=t_ira+seg(I_index)%length
    end if
            
  return
END FUNCTION
!=====================================================================================================================
FUNCTION imatime_periodic(tt)              ! module for getting with periodic boundary condition
  implicit none
  integer(8) :: tt,imatime_periodic
            if (tt<0) then
                imatime_periodic=tt+N_beta
            else if (tt>=N_beta) then
                 imatime_periodic=tt-N_beta
            else
                imatime_periodic=tt
            end if
  return
END FUNCTION

!=====================================================================================================================
FUNCTION get_starting_point(ind_input)              ! module for getting with periodic boundary condition
  implicit none
  integer(8) ::Subse,ind_input,tt,get_starting_point
  integer(8) ::begin_ind, begin_subse

    if (seg(ind_input)%length==N_beta) then
        get_starting_point=0
    else
        Subse=seg(ind_input)%sub_index
        begin_ind=seg(ind_input)%begin
        begin_subse=seg(Subse)%begin

        if (begin_ind>=begin_subse) then
            get_starting_point=begin_ind-N_beta
        else
            get_starting_point=begin_ind
        end if
    end if
  return
END FUNCTION
!=====================================================================================================================
SUBROUTINE update_energy_occupa_num(old_occ_num,num_change,move_length)              ! module for getting with periodic boundary condition
  implicit none
  integer(8) :: new_occ_num,old_occ_num
  integer    :: num_change
  integer(8) :: move_length,ii
  real(8) ::delta_h0,leng

    delta_h0=delta_E(old_occ_num,num_change)
    leng=1.0*move_length
    leng=leng/N_beta

    h_change=h_change+delta_h0*leng!getE()!h_change+delta_h0*leng/N_beta
    nboson=nboson+leng*num_change!getN()!nboson+leng/N_beta*num_change

  return
END SUBROUTINE

!=====================================================================================================================
FUNCTION delta_E(old_occ_num,num_change)           ! module for deleting a segment
  implicit none
    integer(8) :: new_occ_num,old_occ_num
    integer    :: num_change
    real(8)    ::delta_E
    new_occ_num=old_occ_num+num_change

    delta_E=U/2.0*(new_occ_num**2-old_occ_num**2)-delmu*(new_occ_num-old_occ_num)
  return
END FUNCTION
!=====================================================================================================================
FUNCTION p_crew_noloop(t_ira,t_masha,old_occ_num,num_change,tm1,tm2)           ! module for deleting a segment
  implicit none
  integer(8) :: t_ira,t_masha,old_occ_num,tm1,tm2
  integer    :: num_change
  real(8) ::R,W_r,p_crew_noloop
  real(8) ::deltat,delta_h0
  real(8) ::max_min,leng1,inter_to_real

    inter_to_real=(tm2-tm1)*1.0
    max_min=beta*inter_to_real/N_beta
    delta_h0=delta_E(old_occ_num,num_change)

    if (t_masha>t_ira) then
        leng1=(t_masha-t_ira)*1.0
        deltat=ABS(beta*leng1/N_beta)
        W_r=old_occ_num*dexp(-deltat*delta_h0)
        !W_r=dexp(-deltat*delta_h0)
    else if (t_masha<t_ira) then
        leng1=(t_ira-t_masha)*1.0
        deltat=ABS(beta*leng1/N_beta)
        W_r=(old_occ_num+1)*dexp(-deltat*delta_h0)
        !W_r=dexp(-deltat*delta_h0)
    end if
    R=W_r*max_min**2/2.0*(flg_used)*p_delw*omega_G
    p_crew_noloop=MIN(1.0,R)
  return
END FUNCTION
!=====================================================================================================================

FUNCTION p_crew_loop(leng_IM,old_occ_num,num_change)           ! module for deleting a segment
  implicit none
  integer(8) ::leng_IM,old_occ_num
  integer :: num_change
  real(8) ::R,W_r,p_crew_loop
  real(8) ::deltat,inter_to_real
   real(8)   ::delta_h0

    inter_to_real=leng_IM*1.0
    deltat=beta*inter_to_real/N_beta
    delta_h0=delta_E(old_occ_num,num_change)

    if (num_change<0) then
        W_r=old_occ_num*dexp(-deltat*delta_h0)
        !W_r=dexp(-deltat*delta_h0)
    else
        W_r=(old_occ_num+1)*dexp(-deltat*delta_h0)
        !W_r=dexp(-deltat*delta_h0)
    end if

    R=W_r*beta*beta/2.0*(flg_used)*p_delw*omega_G
    p_crew_loop=MIN(1.0,R)
  return
END FUNCTION
!=====================================================================================================================
FUNCTION p_delw_loop(leng_IM,old_occ_num,num_change)
  implicit none
  integer(8) ::leng_IM,old_occ_num,N_fl
  integer :: num_change
  real(8) ::R,W_r,p_delw_loop,inter_to_real
  real(8) ::deltat,delta_h0

    inter_to_real=leng_IM*1.0
    deltat=beta*leng_IM*1.0/N_beta
     delta_h0=delta_E(old_occ_num,num_change)
    if (num_change>0) then
        W_r=1.0/(old_occ_num+1)*dexp(-deltat*delta_h0)
        !W_r=1.0*dexp(-deltat*(U/2*(2*n_m-1)-delmu))
    else
        W_r=1.0/old_occ_num*dexp(-deltat*delta_h0)
        !W_r=1.0*dexp(-deltat*(U/2*(-2*n_m+1)+delmu))
    end if
    N_fl=flg_used-1
    R=W_r/beta/beta/N_fl/p_delw/omega_G*2
    p_delw_loop=MIN(1.0,R)
  return
END FUNCTION

FUNCTION p_delw_noloop(leng_tot,leng_IM,old_occ_num,num_change)       ! module for deleting a segment
  implicit none
  integer(8) ::leng_IM,old_occ_num,N_fl,leng_tot
  integer :: num_change
  real(8) ::R,W_r,p_delw_noloop
  real(8) ::deltat,leng1
  real(8) ::max_min, delta_h0

    max_min=beta*(leng_tot-1.0)/N_beta
    deltat=beta*leng_IM*1.0/N_beta
     delta_h0=delta_E(old_occ_num,num_change)
    if (num_change>0) then
        W_r=1.0/(old_occ_num+1)*dexp(-deltat*delta_h0)
        !W_r=1.0*dexp(-deltat*(U/2*(2*n_m-1)-delmu))
    else
        W_r=1.0/old_occ_num*dexp(-deltat*delta_h0)
        !W_r=1.0*dexp(-deltat*(U/2*(-2*n_m+1)+delmu))
    end if

    N_fl=flg_used-2
    R=W_r/max_min/max_min/N_fl/p_delw/omega_G*2
    p_delw_noloop=MIN(1.0,R)
  return
END FUNCTION

!=====================================================================================================================
FUNCTION prob_ts(old_occ_num,num_change,move_length)          ! module for deleting a segment
  implicit none
  integer(8) :: move_length,old_occ_num
  integer :: num_change
  real(8) ::R,prob_ts
  real(8) ::deltat,leng1,delta_h0
    leng1=move_length*1.0
    deltat=beta*leng1/N_beta

    delta_h0=delta_E(old_occ_num,num_change)
    R=dexp(-deltat*(delta_h0))
    prob_ts=MIN(1.0,R)
  return
END FUNCTION
!=====================================================================================================================
FUNCTION p_crek_acc(tm1,tm2,tau1,n1,n2,backward,old_occ_adja)           ! module for calculating the probability of creating a kink
  implicit none
  integer(8) ::tm1,tm2,tau1,n1,n2,old_occ_adja
  real(8) ::R,W_r,p_crek_acc
  real(8) ::deltat,deltaE1,deltaE2
  real(8) ::max_min,inter_to_real,check_ii
  logical ::backward

    inter_to_real=(tm2-tm1)*1.0
    max_min=beta*inter_to_real/N_beta
    
    if (backward)then
        inter_to_real=1.0*(tm2-tau1)
        deltat=beta*inter_to_real/N_beta

        deltaE1=delta_E(n1,-1)
        deltaE2=delta_E(old_occ_adja,1)
        W_r=t*n2*dexp(-deltat*(deltaE1+deltaE2))
        !W_r=t*dexp(-deltat*(deltaE1+deltaE2))*sqrt(1.0*n2*n1)
        !W_r=t*n2*dexp(-deltat*(U/2.0*(-2*n1+2*n2)))
    else
        inter_to_real=1.0*(tau1-tm1)
        deltat=beta*inter_to_real/N_beta

        deltaE1=delta_E(n1-1,1)
        deltaE2=delta_E(old_occ_adja,-1)
        W_r=t*n2*dexp(-deltat*(deltaE1+deltaE2))
        !W_r=t*dexp(-deltat*(deltaE1+deltaE2))*sqrt(1.0*n2*n1)
        !W_r=t*n2*dexp(-deltat*(U/2.0*(+2*n1-2*n2)))
    end if

    R=W_r*max_min*P_delk/P_crek*2*D
    p_crek_acc=MIN(1.0,R)
  return
END FUNCTION
!=====================================================================================================================
FUNCTION p_delk_acc(tm1,tm2,tau1,n1,n2,backward,I_index,Pre_connec)           ! module for calculating the probability of deleting a kink
  implicit none
  integer(8) :: tm1,tm2,tau1,n1,n2,check_ii,I_index,Pre_connec
  real(8) ::R,W_r,p_delk_acc
  real(8) ::deltat,inter_to_real,deltaE1,deltaE2
  real(8) ::max_min
  logical ::backward

    inter_to_real=1.0*(tm2-tm1)
    max_min=beta*inter_to_real/N_beta

if (max_min<0) then
    print*,'maxtime &mintime error'
    check_ii=index(100000000)
end if

if (n2<0) then
    print*,'occ error'
    check_ii=index(100000000)
end if

    if (backward)then
        inter_to_real=1.0*(tm2-tau1)
        deltat=beta*inter_to_real/N_beta
        deltaE1=delta_E(seg(I_index)%n_b+1,-1)
        deltaE2=delta_E(seg(Pre_Connec)%n_b-1,1)
        !W_r=1.0/t/n2*dexp(-deltat*(U/2.0*(+2*n1-2*n2)))
        W_r=1.0/t/n2*dexp(-deltat*(deltaE1+deltaE2))
        !W_r=1.0/t*dexp(-deltat*(deltaE1+deltaE2))/sqrt(1.0*n2*n1)
    else
        inter_to_real=1.0*(tau1-tm1)
        deltat=beta*inter_to_real/N_beta
         !W_r=1.0/t/n2*dexp(-deltat*(U/2.0*(-2*n1+2*n2)))
        deltaE1=delta_E(seg(I_index)%n_b,1)
        deltaE2=delta_E(seg(Pre_Connec)%n_b,-1)
         W_r=1.0/t/n2*dexp(-deltat*(deltaE1+deltaE2))
        !W_r=1.0/t*dexp(-deltat*(deltaE1+deltaE2))/sqrt(1.0*n2*n1)
    end if

    R=W_r/max_min/P_delk*P_crek/2/D
    p_delk_acc=MIN(1.0,R)
!print*,'p_delk_acc',p_delk_acc
  return
END FUNCTION
!=====================================================================================================================
SUBROUTINE move_begin_end(cur_index,Prior_index,time_to_move,length_cur,length_pri)              ! module for moving starting point in imaginary axis
  implicit none                     !but only for segments without zero crossing
  integer(8) :: cur_index,Prior_index,time_to_move
  integer(8) :: length_pri,length_cur

    seg(Prior_index)%length=length_pri
    seg(cur_index)%begin=imatime_periodic(time_to_move)
    seg(cur_index)%length=length_cur
    
    Call griding(Prior_index)  !supergrid
    Call griding(cur_index)
  return
END SUBROUTINE

!=====================================================================================================================
SUBROUTINE cre_seg(s,l1,l2,n,i1,i2,i3,p,inx)              ! module for deleting a segment
  implicit none
  integer(8) :: s,l1,l2,n,i1,i2,i3,p,inx
  integer(8) :: cur
    cur=inx

    seg(cur)%site=s
    seg(cur)%begin=l1
    seg(cur)%length=l2
    seg(cur)%n_b=n
    seg(cur)%pri_index=i1
    seg(cur)%sub_index=i2
    seg(cur)%kink_index=i3
    seg(cur)%selfp=p

    Call griding(cur)
  return
END SUBROUTINE

SUBROUTINE delete_seg(pointee)              ! module for deleting a segment
  implicit none
  integer(8) :: pointee,err
  integer(8) :: temporary
    !delete a segment ira
    temporary=index(pointee)
    index(pointee)=index(flg_used)
        if (index(flg_used)==0)then
            print*,'index error'
print*,seg
            err=index(100000000)
        end if
        if (temporary==0)then
            print*,'index error'
                err=index(100000000)
        end if

    seg(index(flg_used))%selfp=seg(temporary)%selfp     !update the pointer to index array of the segemnt
    if (masha==flg_used) then
        masha=seg(temporary)%selfp
    end if

        if (seg(temporary)%selfp/=pointee) then
            print*, 'pointer error when deleting'
            err=index(100000000)
        end if
    index(flg_used)=temporary
    flg_used=flg_used-1

    Call null_seg(temporary)
  return
END SUBROUTINE

SUBROUTINE null_seg(inx)
  implicit none
  integer(8) :: cur,inx
    cur=inx

    seg(cur)%site=0
    seg(cur)%begin=0
    seg(cur)%length=0
    seg(cur)%n_b=0
    seg(cur)%pri_index=0
    seg(cur)%sub_index=0
    seg(cur)%kink_index=0
    seg(cur)%selfp=0
  return
END SUBROUTINE

!=====================================================================================================================
FUNCTION get_nearest_neighbor(site1)               ! module for time shift ira
  implicit none
    integer(8) :: dir,get_nearest_neighbor,site1

        dir=dint(nnb*rn())+1        ! randomly choose one of the n.n. sites of ira
        get_nearest_neighbor=nnsite(site1,dir)
  return
END FUNCTION

SUBROUTINE suggest_to_create_kink_backward(I_index,Prior,adja,is_loop)
implicit none
 integer(8) :: t_adja,nt_kink,t_kink,t_ira
 integer(8) :: tmin,tmax,len, length_new, length_ira, length_adja,old_occ_adja
 integer(8) :: I_index,Prior,adja
 integer(8) :: I_new, index_ira,index_new
 integer(8) :: nul=0
 real(8) :: P_crekink=0
 logical  :: is_loop,backward=.true.

    t_ira=seg(I_index)%begin
    old_occ_adja=seg(adja)%n_b
     if (is_loop) then
         tmin=get_starting_point(Prior)
         tmax=t_ira
         len=tmax-tmin
         t_kink=dint(rn()*(len-1))+1          !pick up a point betwen tmax and tmin
         t_kink=t_kink+tmin

         length_new=t_ira-t_kink
         length_ira=N_beta-length_new
     else
         if (seg(adja)%begin<=t_ira .and. seg(adja)%begin+seg(adja)%length>=t_ira) then
             t_adja=seg(adja)%begin
         else if (seg(adja)%begin>t_ira) then
             t_adja=seg(adja)%begin-N_beta
         end if
         tmin=MAX(get_starting_point(Prior),t_adja)
         tmax=t_ira
         len=tmax-tmin
         t_kink=dint(rn()*(len-1))+1
         t_kink=t_kink+tmin

         length_new=t_ira-t_kink
         length_adja=t_kink-t_adja
         length_ira=seg(adja)%length+t_adja-t_ira
     end if

     if (len<=1) return
     P_crekink=p_crek_acc(tmin,tmax,t_kink,seg(Prior)%n_b,seg(adja)%n_b+1,backward,old_occ_adja)
         if (rn()<P_crekink) then  !accept or not
             ira=flg_used+1      !update of ira
             I_new=flg_used+2    ! apply for new segment
             index_ira=index(ira)
             index_new=index(I_new)
             flg_used=flg_used+2

             nt_kink=imatime_periodic(t_kink)  !note t1 may be negative
             Call update_old_world_line(I_index,Prior,nt_kink,index_new,length_new,backward)
            Call update_nnsite_world_line(adja,t_ira,length_ira,length_new,t_kink,I_index,index_ira,&
            & index_new,I_new,is_loop,backward)

             n_k=n_k+1
             Call update_energy_occupa_num(seg(Prior)%n_b,-1,seg(index_new)%length)
             Call update_energy_occupa_num(old_occ_adja,1,seg(index_new)%length)
         end if
  return
END SUBROUTINE

SUBROUTINE suggest_to_create_kink_forward(I_index,Prior,adja,is_loop)
implicit none
  integer(8) :: I_index,Prior,adja
  integer(8) :: t_ira,t_adja,nt_kink,t_kink
  integer(8) :: tmin,tmax,len, length_new, length_ira, length_adja,old_occ_adja
  integer(8) :: I_new, index_ira,index_new,ii
  integer(8) :: nul=0
  real(8) :: P_crekink=0
  logical  :: is_loop,backward=.false.

        t_ira=seg(I_index)%begin
        old_occ_adja=seg(adja)%n_b
        if (is_loop) then
            tmax=t_ira+seg(I_index)%length
            tmin=t_ira
            len=tmax-tmin
            t_kink=dint(rn()*(len-1))+1
            t_kink=t_kink+tmin

            length_ira=t_kink-t_ira
            length_new=N_beta-length_ira
        else
            tmin=t_ira
            if (seg(adja)%begin<=t_ira .and. seg(adja)%begin+seg(adja)%length>=t_ira) then
                t_adja=seg(adja)%begin
            else if (seg(adja)%begin>t_ira) then
                t_adja=seg(adja)%begin-N_beta
            end if
            tmax=MIN(t_adja+seg(adja)%length,seg(I_index)%begin+seg(I_index)%length)
            len=tmax-tmin

            t_kink=dint(rn()*(len-1))+1
            t_kink=t_kink+tmin
            length_ira=t_kink-t_ira
            length_adja=get_length_back(t_ira,seg(adja)%begin)
            length_new=seg(adja)%length-length_ira-length_adja
        end if

        if (len<=1) return
        P_crekink=p_crek_acc(tmin,tmax,t_kink,seg(Prior)%n_b,seg(adja)%n_b,backward,old_occ_adja)
        if (rn()<P_crekink) then            !accecp or not
            
            ira=flg_used+1      !update of ira
            I_new=flg_used+2
            flg_used=flg_used+2
            index_ira=index(ira)
            index_new=index(I_new)
            nt_kink=imatime_periodic(t_kink)  !note t1 may be larger than N_beta

            Call update_old_world_line(I_index,Prior,nt_kink,index_new,length_ira,backward)
            Call update_nnsite_world_line(adja,t_ira,length_ira,length_new,t_kink,I_index,&
                & index_ira,index_new,I_new,is_loop,backward)
            n_k=n_k+1
            Call update_energy_occupa_num(seg(Prior)%n_b-1,1,seg(index_ira)%length)
            Call update_energy_occupa_num(old_occ_adja,-1,seg(index_ira)%length)
        end if
  return
END SUBROUTINE

SUBROUTINE update_old_world_line(segment1,segment2,time_to_move,ind_kink,length_mid,backward)
    implicit none
    integer(8) :: segment1,segment2,time_to_move,ind_kink,length1,length2,length_mid
    logical :: backward

        if (backward) then
            length1=seg(segment1)%length+length_mid
            length2=seg(segment2)%length-length_mid
        else
            length1=seg(segment1)%length-length_mid
            length2=seg(segment2)%length+length_mid
        end if

        Call move_begin_end(segment1,segment2,time_to_move,length1,length2)
        seg(segment1)%kink_index=ind_kink
  return
END SUBROUTINE

SUBROUTINE update_nnsite_world_line(adja,t_ira,length_ira,length_new,t_kink,index_kink,index_ira,index_new,I_new,is_loop,backward)
    implicit none
    integer(8) :: adja,index_ira
    integer(8) :: nul=0
    integer(8) :: index_new,I_new
    integer(8) :: t_ira,length_ira,length_new,t_kink,nt_kink,index_kink
    logical :: is_loop,backward

        nt_kink=imatime_periodic(t_kink)   !note t1 may be negative
        
    if (is_loop) then
        if (backward) then
            Call cre_seg(seg(adja)%site,t_ira,length_ira,seg(adja)%n_b,index_new,index_new,nul,ira,index_ira)
            Call cre_seg(seg(adja)%site,nt_kink,length_new,seg(adja)%n_b+1,index_ira,index_ira,index_kink,I_new,index_new)
            Call delete_seg(seg(adja)%selfp)
        else
            Call cre_seg(seg(adja)%site,t_ira,length_ira,seg(adja)%n_b-1,index_new,index_new,nul,ira,index_ira)
            Call cre_seg(seg(adja)%site,nt_kink,length_new,seg(adja)%n_b,index_ira,index_ira,index_kink,I_new,index_new)
            Call delete_seg(seg(adja)%selfp)
        end if
    else
        if (backward) then
            Call cre_seg(seg(adja)%site,t_ira,length_ira,seg(adja)%n_b,index_new,seg(adja)%sub_index,nul,ira,index_ira)
            Call cre_seg(seg(adja)%site,nt_kink,length_new,seg(adja)%n_b+1,adja,index_ira,index_kink,I_new,index_new)
            Call cut_adjacent_into_three(adja,index_new,index_ira,backward)
        else
            Call cre_seg(seg(adja)%site,t_ira,length_ira,seg(adja)%n_b-1,adja,index_new,nul,ira,index_ira)
            Call cre_seg(seg(adja)%site,nt_kink,length_new,seg(adja)%n_b,index_ira,seg(adja)%sub_index,index_kink,I_new,index_new)
            Call cut_adjacent_into_three(adja,index_new,index_ira,backward)
        end if
    end if

  return
END SUBROUTINE

SUBROUTINE cut_adjacent_into_three(adja,index_new,index_ira,backward)
    implicit none
    integer(8) :: adja,index_new,index_ira,sub_adja
    logical ::backward

    seg(adja)%length=seg(adja)%length-seg(index_new)%length-seg(index_ira)%length  !update of the adjacent segment
    sub_adja=seg(adja)%sub_index
    if (backward) then
         seg(adja)%sub_index=index_new
         seg(sub_adja)%pri_index=index_ira
    else
         seg(adja)%sub_index=index_ira
         seg(sub_adja)%pri_index=index_new
    end if

    Call griding(adja)
  return
END SUBROUTINE

!=====================================================================================================================
SUBROUTINE update_del_kink(backward,Connec,seg_index1,Pre_connec,seg_index2,I_index,is_loop)    ! module for deleting a segment
  implicit none
  logical :: backward
  integer(8) :: Connec,Pre_connec,I_index,seg_index1,seg_index2
  integer(8) :: Prior,Subse,Pre_pri,ii
  integer(8) :: nul=0
  integer(8) :: pointer_new,index_new,time_to_move,delta_length,length_connec,length_pre_connec
  logical :: is_loop
    
    !ira=seg(Connec)%selfp           !update the ira segment

    time_to_move=seg(I_index)%begin
    !print*, '------------------------------'
    !print*, 'segment of ira'
    !print*, seg(I_index)
    !print*,'backward??', backward
    if (backward) then

            delta_length=get_length_back(time_to_move,seg(Connec)%begin)
            !print*,'xxxxxxxxxxxxxx',time_to_move,seg(Connec)%begin
            !print*,'delta_length',delta_length
            !print*,seg(Connec)
            !print*,seg(Pre_Connec)
            length_connec=seg(Connec)%length-delta_length
            length_pre_connec=seg(Pre_connec)%length+delta_length

        if (length_connec<0 .or. length_connec>N_beta) then
            print*,'---------------------------'
            print*,'length error in move_begin_end'
            print*,'---------------------------'
             ii=index(10000000)
        end if
    
    else
            delta_length=get_length_forward(time_to_move,seg(Connec)%begin)
            length_connec=seg(Connec)%length+delta_length
            length_pre_connec=seg(Pre_connec)%length-delta_length
        
        if (length_connec<0 .or. length_connec>N_beta) then
            print*,'---------------------------'
            print*,'length error in move_begin_end'
            print*,'---------------------------'
             ii=index(10000000)
        end if
    end if
    
    seg(Connec)%kink_index=0   !update the prior segment of connec
    Call move_begin_end(Connec,Pre_connec,time_to_move,length_connec,length_pre_connec)

    if (is_loop) then
        Prior=seg_index2

        if (backward) then
            Call back_to_a_loop(I_index,Prior,seg(I_index)%n_b)
        else
            Call back_to_a_loop(I_index,Prior,seg(Prior)%n_b)
        end if

    else
        if (backward) then
            Pre_pri=seg_index1
            Prior=seg_index2

            Call merge_three_into_one(Pre_pri,Prior,I_index)
        else
            Subse=seg_index2
            Prior=seg_index1
            Call merge_three_into_one(Prior,I_index,Subse)
        end if

    end if

  return
END SUBROUTINE

FUNCTION get_length_back(time_to_move,cur_time)
    implicit none
    integer(8) :: time_to_move,cur_time,get_length_back

        if (time_to_move>=cur_time)then
           get_length_back=time_to_move-cur_time
        else
           get_length_back=N_beta+(time_to_move-cur_time)
        end if
    return
End FUNCTION

FUNCTION get_length_forward(time_to_move,cur_time)
        implicit none
        integer(8) :: time_to_move,cur_time,get_length_forward

            if (time_to_move<=cur_time)then
                get_length_forward=cur_time-time_to_move
            else
                get_length_forward=N_beta+(cur_time-time_to_move)
            end if
        return
End FUNCTION


SUBROUTINE back_to_a_loop(I_index,Prior,nb)
  implicit none
  integer(8) :: I_index,Prior,nb
  integer(8) :: nul=0
  integer(8) :: pointer_new,index_new

        pointer_new=flg_used+1
        index_new=index(pointer_new)
        flg_used=flg_used+1
        Call cre_seg(seg(I_index)%site,nul,N_beta,nb,index_new,index_new,nul,pointer_new,index_new)
        Call griding(index_new)
        Call delete_seg(seg(I_index)%selfp)     !delete the segment connected by the kink
        Call delete_seg(seg(Prior)%selfp)     !delete the segment connected by the kink
  return
END SUBROUTINE

SUBROUTINE merge_three_into_one(segment1,segment2,segment3)
  implicit none
  integer(8) :: segment1,segment2,segment3,sub3

        seg(segment1)%length=seg(segment1)%length+seg(segment2)%length+seg(segment3)%length
        seg(segment1)%sub_index=seg(segment3)%sub_index

        sub3=seg(segment3)%sub_index    !update the subsequent segment of segment3
        seg(sub3)%pri_index=segment1
        Call griding(segment1)

        Call delete_seg(seg(segment2)%selfp)
        Call delete_seg(seg(segment3)%selfp)
  return
END SUBROUTINE

!=====================================================================================================================
FUNCTION adjacent(tira,nnira,I_index)             ! module for searching for adjacent segment
  implicit none
  integer(8) :: tira,nnira,I_index,adjacent
  integer(8) :: grid_num,ii
  integer(8) :: ad_index,ad_begin,che
  logical :: flg_find=.false.

    ad_index=supergrid(nnira)

!if (seg(ad_index)%site==0) then
!print*,ad_index,tira
   ! ii=index(1000000)
!end if

    if (seg(ad_index)%begin==0) then
        ad_begin=0
    else if (seg(ad_index)%begin>0) then
        ad_begin=seg(ad_index)%begin-N_beta
    end if

    do while (ad_begin<=tira)
           if (ad_begin<=tira .and. ad_begin+seg(ad_index)%length >=tira) then
               adjacent=ad_index
                flg_find=.true.
               return
           else
               ad_begin=ad_begin+seg(ad_index)%length
               ad_index=seg(ad_index)%sub_index
           end if
       end do

if (.not. flg_find) then
        print*,'---------------------------'
        print*,'adjacent find error'
        print*,'---------------------------'
         ii=index(10000000)
end if
  return
END FUNCTION

!=====================================================================================================================
SUBROUTINE griding(indece)             ! module for searching for adjacent segment
  implicit none
  integer(8) :: indece
  integer(8) :: t1,t2
  integer(8) :: sit
    t1=seg(indece)%begin
    t2=seg(indece)%begin+seg(indece)%length
    sit=seg(indece)%site

    if (t1==0) then
        supergrid(sit)=indece
    else if (t1<=N_beta .and. t2>=N_beta) then
         supergrid(sit)=indece
    end if
    
  return
END SUBROUTINE

!=====================================================================================================================

! correlation time
SUBROUTINE correlation            !give the flag whether thw correlation of two data is too high
  implicit none
  integer(8) :: ii,jj,num_blk
  real(8) ::  M0,delta_t,Numerator,Denominator
    
  M0=0
  num_blk=Nsamp
  do ii=1,num_blk
    M0=M0+M_corr(ii)  !calculate the mean value of variable M
  end do
  M0=M0/num_blk
    
  do ii=1, num_blk-1
      Denominator = 0
      do jj=1, num_blk-ii
          Denominator = Denominator + (M_corr(jj) - M0) ** 2
      end do
        Numerator = 0
      do  jj=1, num_blk-ii
            Numerator = Numerator + (M_corr(jj) - M0) * (M_corr(jj+ii-1) - M0)
      end do
        gamma(ii)=Numerator / Denominator
   end do

    delta_t=N_each
    t0=delta_t/2
    do ii=1,num_blk-1
        t0 = t0+ gamma(ii)*delta_t
    end do

    return
END SUBROUTINE


  subroutine check()
    implicit none
    integer(8) :: ii,ind,leng,beginetime,index_kink
    integer(8) ::length_tot,ind_sub
    real(8) :: N_occ
    print*, "checking ..."
    if (flg_used>multiple*Vol) then
        print*, "index memory overflow"
    end if
if (GT_ZF) then
    if (seg(index(ira))%length==0) then
        print*,'ira length error'
        print*,seg(index(ira))
        ii=index(1000000000)
    end if
    if (seg(index(masha))%length==0) then
        print*,'masha length error'
        print*,seg(index(masha))
        ii=index(1000000000)
    end if

    if (seg(seg(index(ira))%pri_index)%n_b-seg(index(ira))%n_b/=1) then
        print*,'ira occupation number error'
        print*,seg(index(ira))
        ii=index(1000000000)
    end if

end if

    if (seg(index(flg_used+1))%length/=0 .or. seg(index(flg_used+1))%site/=0)then
        print*, "flag error"
        print*, seg(index(flg_used+1))
if (flg_used==Vol*multiple) then
    print*, "overflow"
end if
        ind_sub=index(10000000)
    end if

    if (GT_ZF) then
        if (seg(index(ira))%site==0) then
        print*, "delete ira error"
        ind_sub=index(10000000)
        end if
    end if

    do ii=1,flg_used
        ind=index(ii)
        leng=seg(ind)%length
        beginetime=seg(ind)%begin

        if(seg(ind)%begin==N_beta)then
            print*, "begin overflow"
            print*,seg(ind)
            ind_sub=index(10000000)
        end if

        if(seg(ind)%n_b<0)then
            print*, "occupation number overflow"
            print*,seg(ind)
            ind_sub=index(10000000)
        end if

        if (leng>N_beta .or. leng<=0) then
            print*, "length overflow"
            print*,seg(ind)
            ind_sub=index(10000000)
        end if
        if (beginetime>N_beta) then
            print*, "begin overflow> N"
        else if (beginetime<0) then
            print*, "begin overflow <0"
        end if
        !if (seg(ind)%pri_index/=seg(ind)%sub_index) then
            !print*, "loop error"
        !end if
        if (seg(index(seg(ind)%selfp))%selfp/=seg(ind)%selfp) then
            print*, "self pointer error"
        end if

        length_tot=seg(ind)%length
        ind_sub=seg(ind)%sub_index
        do while (ind/=ind_sub)
            length_tot=length_tot+seg(ind_sub)%length
            ind_sub=seg(ind_sub)%sub_index
        end do
        if (length_tot/=N_beta) then
            print*, "length error"
            print*, seg(ind)
            ind_sub=seg(ind)%sub_index
            print*,seg(ind_sub)
            do while (ind/=ind_sub)
                ind_sub=seg(ind_sub)%sub_index
                print*,seg(ind_sub)
            end do
            ind_sub=index(10000000)
        end if
        ind_sub=seg(ind)%sub_index
        if (MOD(seg(ind)%begin+seg(ind)%length,N_beta)/=seg(ind_sub)%begin) then
            print*, "begin+length error"
            print*, seg(ind)
            print*,'------------------------------------------------------------------------'
            print*, seg
            ind_sub=index(10000000)
        end if
        
        if (seg(ind)%length==N_beta .and. seg(ind)%begin/=0) then
            print*, "loop begin error"
            print*, seg(ind)
            print*,'------------------------------------------------------------------------'
            print*, seg
            ind_sub=index(10000000)
        end if
    
        if (seg(ind)%kink_index/=0) then
            if (seg(seg(ind)%kink_index)%length==N_beta) then
                  print*, "kink seg length error"
                print*,seg(ind)
                 ind_sub=index(10000000)
            end if
        end if

        if (seg(ind)%kink_index/=0) then
            index_kink=seg(ind)%kink_index
            if (seg(index_kink)%kink_index/=ind) then
                 print*, "kink index error"
                print*,seg(ind)
                print*,seg(index_kink)
                print*,seg(seg(index_kink)%kink_index)
                 print*,'------------------------------------------------------------------------'
                print*,seg
                ind_sub=index(10000000)
            end if
        end if
    end do
    !nboson=getN()
    N_occ=getN()
    if (dabs(nboson-N_occ)>1.d-7) then
        print*, "energy update error"
        print*,N_occ,nboson
do ii=1,flg_used
    print*, seg(index(ii))
end do
        ind_sub=index(10000000)
    end if


    !print*,nboson
    !print*,h_change
    !print*,n_k
    print*, "check done"
    return
  end subroutine
    
        
  !====================== END MARKOV =========================


  !================== START MEASUREMENT ======================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE measure
    implicit none


    step=step+1

!    i1=i1+1
!    if(i1>=kw) then
!        call sigma()
!        i1=0
!    end if

!	call check()

    if (GT_ZF) return
    !if (seg(ira)%site==seg(masha)%site) then
       ! N_kink=N_kink+1
    !end if

    Z=Z+1
    Obser(1)%vnr=nboson   !first observable would be the occu number
    Obser(2)%vnr=-n_k/beta+h_change
    Obser(3)%vnr=+n_k*(n_k-1)/beta**2+h_change**2-2*n_k/beta*h_change
    Obser(4)%vnr=n_k
    Obser(5)%vnr=0

    if(prtflg) then          !flag for write2file
        call coll_data()
    end if


!    i1=i1+1
!    if(i1>=kw) then
!        call sigma()
!        i1=0
!    end if


    wormstep=wormstep+step
    step=0

    return
  END SUBROUTINE measure
  



  SUBROUTINE sigma()   !write to file
    implicit none
    integer :: i,k,n,it

    open(1,file="Energy.txt",access="append")
   ! write(1,*) wormstep+step,E/Vol
    write(1,*) wormstep+step
    close(1)

    return
  END SUBROUTINE

  !================== END MEASUREMENT =======================
  
  !============== START WRITE TO FILES =====================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE writejjfile 
    implicit none
    integer     :: j
    real(8)     :: xx

    return
  END SUBROUTINE writejjfile

  SUBROUTINE printing(id)
    implicit none
    integer :: id,i
    real(8) :: temp
    write(id,*) "========================================"
    write(id,*) "Now	:"
    write(id,*) "Z=",Npb*(Nbk-1)+Np
    do i=1,NObs_b
        write(id,*) i,trim(Obser(i)%nam)//'	',Obser(i)%vnr*Obser(i)%a+Obser(i)%b
    end do
    if(prtflg) then
        write(id,*) "Average:"
        write(id,*) "ZZ=",Npb*(Nbk-1)
        do i=1,NObs
            write(id,*) i,'<'//trim(Obser(i)%nam)//'>	',Obser(i)%val*Obser(i)%a+Obser(i)%b,Obser(i)%erb*Obser(i)%a,Obser(i)%cor
        end do
        write(id,*) "flg_cor=",flg_cor,"Npb=",Npb,"Nbk=",Nbk,"Np=",Np
        write(id,*) "wormstep=",wormstep
      write(id,*) "simulation time :",t_simu
    else
        write(id,*) "therm step:",itoss*N_each
        write(id,*) "therm time:",t_toss
    end if
!    write(id,*) "measure time    :",t_meas
    write(id,*) "========================================"
    return
  END SUBROUTINE

  SUBROUTINE write_green
    implicit none
    integer :: i,j
    real(8) :: temp
    write(1,*) "========================================"
    write(1,*) "Now    :"
    write(1,*) "Z=",Z
    write(1,*) "========================================"
    write(1,*) 'distance','Green function '
    write(1,*) "========================================"
    write(1,*) 0,G0/Z

    do i=1,2*(Lx/2)**2
        if (distance(i)/=0) then
             write(1,*)  i,green(i)/Z
        end if
    end do

    if(prtflg) then
        write(1,*) "flg_cor=",flg_cor,"Npb=",Npb,"Nbk=",Nbk,"Np=",Np
        write(1,*) "wormstep=",wormstep
        write(1,*) "simulation time :",t_simu
    else
        write(1,*) "therm step:",itoss*N_each
        write(1,*) "therm time:",t_toss
    end if
!    write(id,*) "measure time    :",t_meas
    write(1,*) "========================================"
    return
  END SUBROUTINE


  SUBROUTINE write_cnf()
    implicit none
    integer(8) :: i,j
    open(1,file=trim(cnf_file),form='UNFORMATTED',access='sequential',status="replace")
    write(1) D,Lx,Ira,Masha,GT_ZF
    close(1)

!    print*,rtal

    return
  END SUBROUTINE

  SUBROUTINE read_cnf()
    implicit none
    integer(1) :: D1
    integer(8) :: Lx1, i, j
    open(1,file=trim(cnf_file),form='UNFORMATTED',access='sequential',status="old")
    read(1) D1,Lx1,Ira,Masha,GT_ZF
    if(D1/=D .or. Lx1/=Lx) stop "the cnf is not suitable for the mismatched Ly."
    close(1)
    return
  END SUBROUTINE

  SUBROUTINE write_stat()
    implicit none
    open(1,file=trim(stat_file),form='UNFORMATTED',access='sequential',status="replace")
    write(1) D,Lx,t,h,beta,N_measure,Seed
    write(1) i1,i2,prtflg,t_simu,t_meas,step,wormstep
    write(1) Npb,Nbk,Np
    write(1) Obser(1:NObs)
    write(1) ir1,ir2,ipnt1,ipnf1,ipnt2,ipnf2,irn,nrannr
    close(1)
    return
  END SUBROUTINE

  SUBROUTINE read_stat()
    implicit none
    integer(1) :: D1
    integer(8) :: Lx1,Ly1,N_measure1
    integer :: Seed1
    real(8) :: t1,h1,beta1
    open(1,file=trim(stat_file),form='UNFORMATTED',access='sequential',status="old")
    read(1) D1,Lx1,t1,h1,beta1,N_measure1,Seed1
    if(D1/=D) stop "the stat is not suitable for the mismatched Dimension."
    if(Lx1/=Lx) stop "the stat is not suitable for the mismatched Lx."
    if(N_measure1/=N_measure)  stop "the stat is not suitable for the mismatched N_measure."
    if(Seed1/=Seed) stop "the stat is not suitable for the mismatched Seed."
    read(1) i1,i2,prtflg,t_simu,t_meas,step,wormstep
    read(1) Npb,Nbk,Np
    read(1) Obser(1:NObs)

    read(1) ir1,ir2,ipnt1,ipnf1,ipnt2,ipnf2,irn,nrannr
    close(1)
    return
  END SUBROUTINE

  !====================== END WRITE TO FILES =============================================
  
  
  !================= Obs ========================================

  SUBROUTINE set_Obs
    implicit none
    integer :: i
    !call init_Obs(1,"energy  ",1.d0/Vol,0.d0)
    call init_Obs(1,"occupat ",1.d0/Vol,0.d0)
    call init_Obs(2,"energy  ",1.d0/Vol,0.d0)
    call init_Obs(3,"E^2     ",1.d0/Vol/Vol,0.d0)
    call init_Obs(4,"none    ",1.d0,0.d0)
    call init_Obs(5,"none    ",1.d0/Vol/Vol/Vol/Vol,0.d0)
    call init_Obs(NObs_b+1,"specheat",1.d0/Vol/Vol,0.d0)
    !call init_Obs(NObs_b+2,"none    ",1.d0/Vol/Vol,0.d0)
    !call init_Obs(NObs_b+3,"susceptibility  ",1.d0/Vol/Vol,0.d0)

    Npb=1
    Nbk=1
    Np=0
    return
  END SUBROUTINE
  
  SUBROUTINE init_Obs(i,nam0,a0,b0)
    implicit none
    integer :: i
    character(8) :: nam0
    real(8) :: a0,b0
    Obser(i)%nam=trim(nam0)
    Obser(i)%a=a0                   !a0 is the multiplier (physical quantity need to be divided by Volume)
    Obser(i)%b=b0                   !b0 is a constant added to the physical quantity
    Obser(i)%vnr=0.d0
    Obser(i)%val=0.d0
    Obser(i)%cor=1.d100
    Obser(i)%erb=1.d100
    Obser(i)%blk=0.d0
    return
  END SUBROUTINE
  
  SUBROUTINE coll_data()
    implicit none
    integer :: i,j
!   do i=1,NObs
!       Obser(i)%val=Obser(i)%val+Obser(i)%vnr
!   end do
    if(Np==Npb) then
        if(flg_cor) then                            !if the correlation between two consequtive data is small enough
            if(Nbk==Max_block) then
                call merge_blk()                    !combine data in two blocks and take the average in max_block
            end if
        else                                        !if the correlation between two consequtive data is big
            if(Nbk>=N_block .and. Nbk/2*2==Nbk) then !if the size of the block > standard block size and the component number is even
                call cal_ComObs()                   !calculeta composite variables
                if(cal_cor_dev()) then              !if the correlation within the data set is small now
                    flg_cor=.True.                  !flag the small correlation
                else
                    call merge_blk()
                end if
            end if
        end if
        Nbk=Nbk+1                                   !list the current position of block
        Np=0                                        !block size Np
    end if
    Np=Np+1
    do i=1,NObs_b
        Obser(i)%blk(Nbk)=Obser(i)%blk(Nbk)+Obser(i)%vnr !add the current value to the block summation
    end do
    return
  END SUBROUTINE
  
  SUBROUTINE merge_blk()
    implicit none
    integer :: i,j
    Nbk=Nbk/2
    Npb=2*Npb
    do i=1,NObs_b
        do j=1,Nbk
            Obser(i)%blk(j)=Obser(i)%blk(2*j-1)+Obser(i)%blk(2*j)
        end do
        Obser(i)%blk(Nbk+1:2*Nbk)=0.d0
    end do
    return
  END SUBROUTINE
  
  FUNCTION cal_cor_dev()            !give the flag whether thw correlation of two data is too high
    implicit none
    logical :: cal_cor_dev
    integer :: i,iblck,Nbk0
    real :: cor,dev,devp,devn
    cal_cor_dev=.True.
    if(Np/=Npb) then
        Nbk0=Nbk-1
    else
        Nbk0=Nbk
    end if
    do i=1,NObs
        Obser(i)%val=Sum(Obser(i)%blk(1:Nbk0))/Nbk0 !take the average of each block value
        cor=0.d0;   dev=0.d0;   devp=0.d0
        do iblck=1,Nbk0
            devn=Obser(i)%blk(iblck)-Obser(i)%val
            dev=dev+devn*devn
            cor=cor+devn*devp
            devp=devn
        end do
        Obser(i)%cor=cor/dev                !calculate correlation
        Obser(i)%val=Obser(i)%val/Npb       !calculate the average value
        Obser(i)%erb=dsqrt(dev/Nbk0/(Nbk0-1.d0))/Npb   !calculate the error bar
        if(dabs(Obser(i)%cor)>Max_cor) then     !if the correlation between two sets of data is larger than Max_cor
            cal_cor_dev=.False.                 !
        end if
    end do
    return
  END FUNCTION
  
  SUBROUTINE cal_ComObs()                           !calculate composite variables
    implicit none
    integer :: iObs(NObs),N
    real(8) ,external :: spec_heat   !define variable name
    iObs(1)=2;  iObs(2)=3;  N=2                  !select what kind of basic variable is needed
    call cal_ComObs_func(NObs_b+1,spec_heat,iObs,N)
    return
  END SUBROUTINE

  SUBROUTINE cal_ComObs_func(indx,func,iObs,N)
        implicit none
        integer :: indx,N,iblck,i
        integer :: iObs(N)
        real(8) :: a(N)
        real(8) ,external :: func
        do iblck=1,Nbk
            do i=1,N
                a(i)=Obser(iObs(i))%blk(iblck)/Npb
            end do
            Obser(indx)%blk(iblck)=func(a,N)*Npb
        end do
        return
  END SUBROUTINE


  !===============Shift register random number generator =============
  !  very long period sequential version
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE set_RNG
    implicit none
    integer :: i_r,k_r,k1_r
    integer :: iseed

    nrannr = mxrn
    iseed  = iabs(Seed)+1
    k_r    = 3**18+2*iseed
    k1_r   = 1313131*iseed
    k1_r   = k1_r-(k1_r/mod2)*mod2

    do i_r = 1, len1
      k_r  = k_r *mult
      k1_r = k1_r*mul2
      k1_r = k1_r-(k1_r/mod2)*mod2
      ir1(i_r) = k_r+k1_r*8193
    enddo

    do i_r = 1, len2
      k_r  = k_r *mult
      k1_r = k1_r*mul2
      k1_r = k1_r-(k1_r/mod2)*mod2
      ir2(i_r) = k_r+k1_r*4099
    enddo

    do i_r = 1, len1
      inxt1(i_r) = i_r+1
    enddo
    inxt1(len1) = 1
    ipnt1 = 1
    ipnf1 = ifd1+1

    do i_r = 1, len2
      inxt2(i_r) = i_r+1
    enddo
    inxt2(len2) = 1
    ipnt2 = 1
    ipnf2 = ifd2 + 1
    return
  END SUBROUTINE set_RNG 
  !===================================================================

 !===============Calculate next random number =======================
  !! THIS IS ALMOST PROJECT-INDEPENDENT 
  double precision function rn()
  !integer function rn()
    implicit none
    integer   :: i_r, l_r, k_r
    nrannr = nrannr +1
    if(nrannr>=mxrn) then
      nrannr = 1
      do i_r= 1, mxrn
        l_r = ieor(ir1(ipnt1),ir1(ipnf1))
        k_r = ieor(ir2(ipnt2),ir2(ipnf2))
        irn(i_r) = ieor(k_r,l_r)
        ir1(ipnt1)=l_r
        ipnt1 = inxt1(ipnt1)
        ipnf1 = inxt1(ipnf1)
        ir2(ipnt2) = k_r
        ipnt2 = inxt2(ipnt2)
        ipnf2 = inxt2(ipnf2)
      enddo
    endif 
    !rn = irn(nrannr)
    rn = irn(nrannr)*tm32+0.5d0
  end function rn
  !===================================================================


  !==============Trace elapsed time ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE set_time_elapse
    implicit none
    !-- read and calculate time (in seconds) -------------------------
    call date_and_time(date, time, zone, tval)
    t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
    h_curr = tval(5)
    t_prev = t_curr
    h_prev = h_curr
    return
  END SUBROUTINE set_time_elapse
    

  !==============Trace elapsed time ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE time_elapse
    implicit none
    
    !-- read and calculate time (in seconds) -------------------------
    call date_and_time(date, time, zone, tval)
    t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
    h_curr = tval(5)

    t_elap = t_curr-t_prev
    if(h_curr<h_prev) t_elap = t_elap+24*3600.d0
    t_prev = t_curr
    h_prev = h_curr 
    return
  END SUBROUTINE time_elapse
  !===================================================================

END PROGRAM

function spec_heat(a,N)
    implicit none
    integer :: N
    real(8) :: spec_heat,a(N)

    spec_heat=a(2)-a(1)**2
    return
end function

function binder(a,N)
    implicit none
    integer :: N
    real(8) :: binder,a(N)
    binder=a(2)/a(1)/a(1)
    return
end function

 !===================================================================



