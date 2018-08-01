
module qtree
implicit none
s



end module qtree

module MAC


end module MAC

module RK4


end module RK4

module galaxy
implicit none
contains
    subroutine initialize(v)
    real(8), intent(IN) :: v
    real(8), parameter ::  parsec = 3.085677581d16, R = parsec*50d3, m_sun = 2d30, M = m_sun*1d11, n=2000
    integer :: i
    real(8) :: grid_distance, x_grid(n), y_grid(n), z_grid(n)
    grid_distance = R/n
    x_grid = (/ (i,i=0,n)*grid_distance /)
    y_grid = (/ (i,i=0,n)*grid_distance /)
    z_grid = (/ (i,i=0,n)*grid_distance /)

    end subroutine intialize



end module galaxy

program gravitational_dynamics
    implicit none
    print *, 'hello World '
end program gravitational_dynamics
