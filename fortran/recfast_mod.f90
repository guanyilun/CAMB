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
     procedure :: Validate => TRecfastMod_Validate
     procedure :: xe_Tm => TRecfastMod_xe_Tm !ionization fraction and baryon temperature
     procedure :: T_m => TRecfastMod_tm !baryon temperature
     procedure :: T_s => TRecfastMod_ts !Spin temperature
     procedure :: dDeltaxe_dtau => TRecfastMod_dDeltaxe_dtau
     procedure :: get_Saha_z => TRecfastMod_Get_Saha_z

     procedure :: Init => TRecfastMod_init
     procedure :: x_e => TRecfastMod_xe
     procedure :: ReadParams => TRecfastMod_ReadParams
     procedure, nopass :: SelfPointer => TRecfastMod_SelfPointer
  end type TRecfastMod

  public TRecfastMod

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

    ! allocate(TRecfast::this%Rec1)
    ! allocate(TRecfast::this%Rec2)
    ! allocate(TRecfast::this%Rec3)
    ! initialize
    call this%Rec1%Init(State, WantTSpin, d1)
    call this%Rec2%Init(State, WantTSpin, d2)
    call this%Rec3%Init(State, WantTSpin, d3)
  end subroutine TRecfastMod_init

  subroutine TRecfastMod_ReadParams(this, Ini)
    use IniObjects
    class(TRecfastMod) :: this
    class(TIniFile), intent(in) :: Ini
    ! real(dl) :: d1,d2,d3,f1,f2,f3,b

    ! initialize three models
    allocate(TRecfast::this%Rec1)
    allocate(TRecfast::this%Rec2)
    allocate(TRecfast::this%Rec3)

    ! allow individual routines to load parameters
    call this%Rec1%ReadParams(Ini)
    call this%Rec2%ReadParams(Ini)
    call this%Rec3%ReadParams(Ini)

  end subroutine TRecfastMod_ReadParams

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

  subroutine TRecfastMod_Validate(this, OK)
    class(TRecfastMod),intent(in) :: this
    logical, intent(inout) :: OK

    call this%Rec1%Validate(OK)
    call this%Rec2%Validate(OK)
    call this%Rec3%Validate(OK)

  end subroutine TRecfastMod_Validate

  function TRecfastMod_tm(this,a)
    class(TRecfastMod) :: this
    real(dl), intent(in) :: a
    real(dl) TRecfastMod_tm, xe1, xe2, xe3

    xe1 = this%Rec1%T_m(a)
    xe2 = this%Rec2%T_m(a)
    xe3 = this%Rec3%T_m(a)
    ! sum up
    TRecfastMod_tm = this%f1*this%d1*xe1 + this%f2*this%d2*xe2 + this%f3*this%d3*xe3
  end function TRecfastMod_tm

  function TRecfastMod_ts(this,a)
    class(TRecfastMod) :: this
    real(dl), intent(in) :: a
    real(dl) TRecfastMod_ts, xe1, xe2, xe3

    xe1 = this%Rec1%T_s(a)
    xe2 = this%Rec2%T_s(a)
    xe3 = this%Rec3%T_s(a)
    ! sum up
    TRecfastMod_ts = this%f1*this%d1*xe1 + this%f2*this%d2*xe2 + this%f3*this%d3*xe3
  end function TRecfastMod_ts

  subroutine TRecfastMod_xe_Tm(this, a, xe, Tm)
    class(TRecfastMod) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: xe, Tm

    call this%Rec1%xe_Tm(a, xe, Tm)
    call this%Rec2%xe_Tm(a, xe, Tm)
    call this%Rec3%xe_Tm(a, xe, Tm)
  end subroutine TRecfastMod_xe_Tm


  function TRecfastMod_dDeltaxe_dtau(this,a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb,adotoa)
    !d x_e/d tau assuming Helium all neutral and temperature perturbations negligible
    !it is not accurate for x_e of order 1
    class(TRecfastMod) :: this
    real(dl) TRecfastMod_dDeltaxe_dtau
    real(dl), intent(in):: a, Delta_xe, Delta_nH, Delta_Tm, hdot, kvb, adotoa
    real(dl) xe1, xe2, xe3
    xe1 = this%Rec1%dDeltaxe_dtau(a, Delta_xe, Delta_nH, Delta_Tm, hdot, kvb, adotoa)
    xe2 = this%Rec2%dDeltaxe_dtau(a, Delta_xe, Delta_nH, Delta_Tm, hdot, kvb, adotoa)
    xe3 = this%Rec3%dDeltaxe_dtau(a, Delta_xe, Delta_nH, Delta_Tm, hdot, kvb, adotoa)
    TRecfastMod_dDeltaxe_dtau = this%f1*this%d1*xe1 + this%f2*this%d2*xe2 + this%f3*this%d3*xe3
  end function TRecfastMod_dDeltaxe_dtau


  real(dl) function TRecfastMod_Get_Saha_z(this)
    class(TRecfastMod) :: this
    TRecfastMod_Get_Saha_z =  this%Calc%recombination_saha_z
  end function TRecfastMod_Get_Saha_z

  subroutine TRecfastMod_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TRecfastMod), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

  end subroutine TRecfastMod_SelfPointer

end module RecfastMod
