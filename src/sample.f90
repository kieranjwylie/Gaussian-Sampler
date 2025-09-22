module sample

    implicit none

    public ! all routines used in main
    private :: gaussian
    
    type :: sampler_obj
        real(8), allocatable :: boundaries(:)
        real(8), allocatable :: height(:)
        real(8), allocatable :: area(:)
        integer              :: n
    end type sampler_obj

    type(sampler_obj) :: sampler 

    contains 

    
    function gaussian(a, b, x) result(res)
        implicit none

        real(8), intent(in) :: a, b, x
        real(8)             :: res

        res = a * exp(-x**2.0_8 / b)
    end function gaussian


    subroutine setup_sampler(xmin, xmax, a, b, n)
        implicit none

        real(8), intent(in) :: xmin, xmax, a, b ! a and b constant of gaussian
        integer, intent(in) :: n

        real(8) :: dx
        integer :: i

        allocate(sampler%boundaries(n+1))
        allocate(sampler%height(n))
        allocate(sampler%area(n))
        
        dx = (xmax - xmin) / float(n)

        do i=1, n + 1
            sampler%boundaries(i) = xmin + (i-1)*dx
        end do

        do i=1, n 
            sampler%height(i) = gaussian(a, b, max(sampler%boundaries(i), sampler%boundaries(i+1)))
            sampler%area(i) = dx * sampler%height(i)
        end do


    end subroutine setup_sampler

    ! Currently how ive setup the boxes is wrong
    subroutine ziggurat(xmin, xmax, n, a, b, res, its)
        implicit none

        real(8), intent(in)  :: xmin, xmax, a, b
        integer, intent(in)  :: n
        real(8), intent(out) :: res
        integer, intent(out) :: its

        real(8) :: sampled_point(2), gaus
        integer :: index, i

        its = 0
        sample_loop: do
            its = its + 1
            call random_number(sampled_point(1))
            call random_number(sampled_point(2))

            sampled_point(1) = sampled_point(1) * (xmax-xmin) + xmin

            index_loop: do i = 1, n 
                if (sampled_point(1) > sampler%boundaries(i) .and. sampled_point(1) < sampler%boundaries(i+1)) then
                    index = i
                    exit index_loop
                end if
            end do index_loop

            sampled_point(2) = sampled_point(2) * sampler%height(index)

            gaus = gaussian(a, b, sampled_point(1))
            if (sampled_point(2) < gaus) then
                res = sampled_point(2)
                
                exit sample_loop
            else
                cycle sample_loop
            end if
        end do sample_loop
        

    end subroutine ziggurat

    subroutine box_muller(a, b, res)

        real(8), intent(in) :: a,b
        real(8), intent(out) :: res

        real(8) :: z, u1, u2, r, pi 

        pi = 3.1415_8


        call random_number(u1)
        call random_number(u2)

        call random_number(z)

        ! Just for fun - can use either one
        if (z < 0.5_8) then
            res = (b/2.0_8)**0.5_8 * (-2.0_8*log(u1))**0.5_8 * cos(2.0_8*pi*u2)
        else
            res = (b/2.0_8)**0.5_8 * (-2.0_8*log(u1))**0.5_8 * sin(2.0_8*pi*u2)
        end if
    end subroutine box_muller

end module sample