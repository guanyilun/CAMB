    !---------------------------------------------------------------------------------------------------
    ! Recombination module for CAMB, using HyRec
    ! Note you will need to rename dtauda_ in history.c to exported_dtauda.
    ! To use with the python wrapper add -fPIC to the HYREC CCFLAGS (for gcc)
    ! YG: update to include baryon inhomogeneity
    !---------------------------------------------------------------------------------------------------

    module HyRecMod
    use HyRec
    use classes
    use results
    implicit none
    private

    type, extends(THyRec) :: THyRecMod
       real(dl) :: f1,f2,f3,d1,d2,d3
       class(THyRec), allocatable :: Rec1, Rec2, Rec3, Rec0
    contains
    procedure :: Init => THyRecMod_init
    procedure :: x_e => THyRecMod_xe
    procedure :: T_m => THyRecMod_tm !baryon temperature
    procedure, nopass :: SelfPointer => THyRecMod_SelfPointer
    end type THyRecMod

    class(CAMBdata), pointer :: CurrentState

    public THyRecMod
    contains


    function THyRecMod_tm(this,a)
    class(THyRecMod) :: this
    real(dl), intent(in) :: a
    real(dl) THyRecMod_tm
    THyRecMod_tm = this%Rec0%T_m(a)
    end function THyRecMod_tm

    function THyRecMod_xe(this,a)
    class(THyRecMod) :: this
    real(dl), intent(in) :: a
    real(dl) THyRecMod_xe, xe1, xe2, xe3
    xe1 = this%Rec1%x_e(a)
    xe2 = this%Rec2%x_e(a)
    xe3 = this%Rec3%x_e(a)
    ! sum up
    THyRecMod_xe = this%f1*this%d1*xe1 + this%f2*this%d2*xe2 + this%f3*this%d3*xe3
    end function THyRecMod_xe

    subroutine THyRecMod_init(this, State, WantTSpin, delta_in)
    implicit none
    class(THyRecMod), target :: this
    class(TCAMBdata), target :: State
    logical, intent(in), optional :: WantTSpin
    real(dl), intent(in), optional :: delta_in
    real(dl) :: d1,d2,d3,f1,f2,f3,b
    select type(State)
    class is (CAMBdata)
       ! given
       b = State%CP%recfast_b
       d1 = State%CP%recfast_d1
       d2 = State%CP%recfast_d2
       f2 = State%CP%recfast_f2
    class default
       call MpiStop('Wrong state type')
    end select
    ! derived
    d3 = (1+b-d1+d1*d2*f2-d2**2*f2)/(1-d1+d1*f2-d2*f2)
    f1 = (b-f2-b*f2+2*d2*f2-d2**2*f2)/(1+b-2*d1+d1**2-d1**2*f2+2*d1*d2*f2-d2**2*f2)
    f3 = -(1-d1+d1*f2-d2*f2)**2/(-1-b+2*d1-d1**2+d1**2*f2-2*d1*d2*f2+d2**2*f2)
    ! store
    this%d1 = d1
    this%d2 = d2
    this%d3 = d3
    this%f1 = f1
    this%f2 = f2
    this%f3 = f3

    allocate(THyRec::this%Rec0)
    allocate(THyRec::this%Rec1)
    allocate(THyRec::this%Rec2)
    allocate(THyRec::this%Rec3)

    ! initialize
    call this%Rec0%Init(State, WantTSpin)  ! original
    call this%Rec1%Init(State, WantTSpin, d1)
    call this%Rec2%Init(State, WantTSpin, d2)
    call this%Rec3%Init(State, WantTSpin, d3)
    end subroutine THyRecMod_init

    subroutine THyRecMod_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (THyRec), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine THyRecMod_SelfPointer

    end module HyRecMod
