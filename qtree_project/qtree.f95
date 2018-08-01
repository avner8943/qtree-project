
module qtree
    implicit none




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
        real(8), parameter ::  parsec = 3.085677581d16, R = parsec*50d3, m_sun = 2d30, M = m_sun*1d11, n=2d3, n2
        real(8), allocatable :: pos(:,3), velocity(:,3)
        integer :: i, j, k, l
        real(8) :: grid_distance, x_grid(n), y_grid(n), z_grid(n)
        grid_distance = R/n
        x_grid = (/ (i,i=0,n)*grid_distance-1d3 /)
        y_grid = (/ (i,i=0,n)*grid_distance-1d3 /)
        z_grid = (/ (i,i=0,n)*grid_distance-1d3 /)

        l=1
        do i=1,n
            do k=1,n
                do j=1,n
                    if (x_grid(i)**2+y_grid(k)**2+z_grid(j)**2<R**2) then
                        pos(l,:) = (x_grid(i), y_grid(k), z_grid(j))
                        l = l+1
                    end if
                end do
            end do
        end do

       n2 = floor(l/2.0_8)


    end subroutine intialize



end module galaxy

program gravitational_dynamics
    implicit none
    print *, 'hello World '
end program gravitational_dynamics
