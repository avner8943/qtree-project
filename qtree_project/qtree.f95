
module oc_tree
    implicit none

    contains

    recursive subroutine cube(x,y,z,pos)
    x_mid = sum(x)/2
    y_mid = sum(y)/2
    z_mid = sum(z)/2

    cube1 = newcube(x(1), x_mid, y(1), y_mid, z(1), z_mid)
    cube2 = newcube(x(1), x_mid, y(1), y_mid, z_mid, z(2))
    cube3 = newcube(x(1), x_mid, y_mid, y(2), z(1), z_mid)
    cube4 = newcube(x(1), x_mid, y_mid, y(2), z_mid, y(2))
    cube5 = newcube(x_mid, x(2), y(1), y_mid, z(1), z_mid)
    cube6 = newcube(x_mid, x(2), y(1), y_mid, z_mid, z(2))
    cube7 = newcube(x_mid, x(2), y_mid, y(2), z(1), z_mid)
    cube8 = newcube(x_mid, x(2), y_mid, y(2), z_mid, z(2))


    end subroutine cube

end module oc_tree

module MAC


end module MAC

module RK4


end module RK4

module galaxy
    implicit none
contains
    subroutine initialize(v,grid,pos,m, velocity)
        real(8), intent(IN) :: v
        integer, parameter :: n=21
        real(8), parameter ::  parsec = 3.085677581d16, R = parsec*50d3, m_sun = 2d30, M_tot = m_sun*1d11,  Pi = 3.1415927
        real(8), allocatable, intent(OUT) :: pos(:,:), velocity(:,:)
        real(8), intent(OUT) :: m, grid(n,3)
        real(8), allocatable :: temp(:,:)
        integer :: i, j, k, l, n2

        real(8) :: grid_distance, sigma,x1,x2,y1,y2,y3

        grid_distance = 2*R/(n-1._8)

        grid(:,1) = (/ (i,i=0,n-1) /)*grid_distance-R
        grid(:,2) = (/ (i,i=0,n-1) /)*grid_distance-R
        grid(:,3) = (/ (i,i=0,n-1) /)*grid_distance-R
        print *, R**2, grid(11,1)
        allocate(pos(int(1d4),3))

        l=1
        do i=1,n
            do k=1,n
                do j=1,n
                    if (grid(i,1)**2+grid(k,2)**2+grid(j,3)**2<R**2) then
                        pos(l,:) = (/ grid(i,1), grid(k,2), grid(j,3) /)
                        l = l+1
                    end if
                end do
            end do
        end do

        if (mod(l,2)==0) then
            call move_alloc(pos,temp)
            allocate (pos(l-2,3))
            pos(:,:) = temp(1:l-2,:)
            n2 = l-2
        else
            call move_alloc(pos,temp)
            allocate (pos(l-1,3))
            pos(:,:) = temp(1:l-1,:)
            n2 = l-1
        end if

       allocate(velocity(n2,3))
       sigma = v/sqrt(3._8)
       do i=1,n2
       call random_number(x1)
       call random_number(x2)
       y1 = sigma*(sqrt(-2*log(x1))*cos(2*pi*x2))
       y2 = sigma*(sqrt(-2*log(x1))*sin(2*pi*x2))
       velocity(i,1:2) = (/ y1,y2 /)
       call random_number(x1)
       call random_number(x2)
       y3 = sigma*(sqrt(-2*log(x1))*cos(2*pi*x2))
       velocity(i,3) = y3
       end do

       m = M_tot/n2

    end subroutine initialize



end module galaxy

program gravitational_dynamics
    use galaxy
    real(8) :: v=85d3, m, grid(21,3)
    real(8), allocatable ::  pos(:,:), velocity(:,:)
    integer :: i,n

    call initialize(v,grid,pos,m, velocity)
     n =size(velocity,1)
     print *, n, grid(:,1)
     open(unit=1 , file = 'test.csv')
     do i=1,n
     write (1, '(3es12.4)')velocity(i,1),velocity(i,2),velocity(i,3)
     end do
     close (unit=1)

     open(unit=2 , file = 'test1.csv')
     do i=1,n
     write (2, '(3es12.4)')pos(i,1),pos(i,2),pos(i,3)
     end do
     close (unit=2)


end program gravitational_dynamics
