! Samples from a gaussian of the form:
! p(x) = -a * exp(-x**2 / b)

program main

    use sample
    
    implicit none

    integer :: n, its, fn, i
    real(8) :: xmin, xmax, a, b, sampled

    
    ! Define parameters here
    n = 100000
    xmin = -100.0_8
    xmax = 100.0_8

    a = 0.3999_8
    b = 2.0_8

    call setup_sampler(xmin, xmax, a, b, n)

    !call ziggurat(xmin, xmax, n, a, b, sampled, its)

    open(newunit= fn, file="data.dat", status="new", action="write")
    
    do i=1, n
        call box_muller(a,b,sampled)
        write(fn, "(g0)") sampled
    end do

    write(*, "(a,g0)") "Sampled point: ", sampled
end program main