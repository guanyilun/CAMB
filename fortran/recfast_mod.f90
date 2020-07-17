module RecfastMod
  use constants
  use classes
  use results
  use Recombination
  use MpiUtils, only : MpiStop
  implicit none
  private

  type, extends(TRecfast) :: TRecfastMod
     real(dl) :: f1,f2,f3,d1,d2,d3
     class(TRecfast), allocatable :: Rec1, Rec2, Rec3
   contains
     procedure :: Init => TRecfastMod_init
     procedure :: x_e => TRecfastMod_xe
     procedure, nopass :: SelfPointer => TRecfastMod_SelfPointer
  end type TRecfastMod

contains

  subroutine TRecfastMod_init(this, State, WantTSpin, delta)
    implicit none
    class(TRecfastMod), target :: this
    class(TCAMBdata), target :: State
    logical, intent(in), optional :: WantTSpin
    real(dl), intent(in), optional :: delta
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
    ! initialize three models
    allocate(TRecfast::this%Rec1)
    allocate(TRecfast::this%Rec2)
    allocate(TRecfast::this%Rec3)
    ! initialize
    call this%Rec1%Init(State, WantTSpin, d1)
    call this%Rec2%Init(State, WantTSpin, d2)
    call this%Rec3%Init(State, WantTSpin, d3)

  end subroutine TRecfastMod_init

  function TRecfastMod_xe(this, a)
    class(TRecfastMod) :: this
    real(dl), intent(in) :: a
    real(dl) TRecfastMod_xe, xe1, xe2, xe3
    ! get inidividual xe and sum them up
    xe1 = this%Rec1%x_e(a)
    xe2 = this%Rec2%x_e(a)
    xe3 = this%Rec3%x_e(a)
    ! sum up
    TRecfastMod_xe = this%f1*this%d1*xe1 + this%f2*this%d2*xe2 + this%f3*this%d3*xe3
  end function TRecfastMod_xe

  subroutine TRecfastMod_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TRecfastMod), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

  end subroutine TRecfastMod_SelfPointer


end module RecfastMod
