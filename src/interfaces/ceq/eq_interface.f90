module eq_interface
  use iso_c_binding
  use ceq_system
  implicit none

  integer(c_int) :: ns, ne 
  real(c_double), allocatable, dimension(:) :: wt
  type (sys_type), pointer :: sys


contains

subroutine eq_hp_init(ne_in, ns_in, wt_in, thermo, Ein) BIND(C)
  implicit none

  integer(c_int), value, intent (in) :: ne_in, ns_in
  real(c_double), dimension(ns_in), intent(in)    :: wt_in
  real(c_double), dimension(ns_in,15), intent(in) :: thermo
  real(c_double), dimension(ns_in,ne_in), intent(in) :: Ein

  integer, parameter :: ng=0, ncs=0
  integer :: iret, CS(ncs), lu_op, diag
  real(k_wp), dimension(ns_in,ng) :: Bg

  ns = ns_in
  ne = ne_in

  if(.not. allocated(wt)) allocate(wt(ns))
  wt = wt_in

  lu_op = 0 ! logical unit for output (used also by ceq_state)
  diag  = 1 ! level of diagnostic output

  call ceq_sys_init( ns, ne, ncs, ng, Ein, CS, Bg, thermo, lu_op, diag, sys, iret )

end subroutine eq_hp_init

subroutine eq_hp(p, T, y) BIND(C)
  implicit none
  real(c_double), value, intent(in) :: p
  real(c_double), intent(inout) :: T
  real(c_double), dimension(ns), intent(inout) :: y
  real(c_double), dimension(ns) :: xpre, xeq, yeq, ypre, zeq, zpre
  real(c_double) :: Teq, stats(20), HoR, mw
  integer :: iret

  ypre = y
  mw = sum(ypre)/sum(ypre/wt)
  xpre = ypre/wt*mw
  zpre = xpre/sum(xpre*wt)

  call ceq_state( sys, N=zpre, p_Pa=p, N_h=zpre, T_h=T, HoR_eq=HoR, T_eq=Teq, N_eq=zeq, info=iret )

  xeq = zeq/sum(zeq)
  yeq = xeq*wt/sum(xeq*wt)

  y = yeq
  T = Teq

  if( iret < 0 ) then!else
     print*,'Fail'
  end if
end subroutine eq_hp

subroutine eq_hp_free() BIND(C)
  !if(associated(sys)) call ceq_sys_rm(sys)
  if(associated(sys)) then
    if(sys%initialized) then
      !N.B. ceq_sys_rm only nullifies these,
      !     so we implement our own to avoid memory leaks
      if(sys%ncs >= 1) deallocate(sys%CS)
      if(sys%ng >= 1) deallocate(sys%Bg)
      deallocate(sys%sp_order)
      deallocate(sys%el_order)
      deallocate(sys%Ein)
      deallocate(sys%E)
      deallocate(sys%B)
      deallocate(sys%BR)
      deallocate(sys%A)
      deallocate(sys%thermo)
      deallocate(sys%vars)
    endif
    deallocate(sys)
  endif
  if(allocated(wt)) deallocate(wt)
end subroutine eq_hp_free

end module eq_interface
