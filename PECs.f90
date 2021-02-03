module constants
    implicit none
    real(kind=8), parameter :: a_so_yb = 806.09d0
    integer, parameter :: ndat = 107, npec = 4, naso = 4
    integer, parameter :: nparam_pec = 4*npec, nparam_aso = 3*naso
    real(kind=8) :: cof_pec(nparam_pec), cof_aso(nparam_aso)
    real(kind=8) :: r(ndat)
    integer, parameter :: N1 = 3, N2=5, LWORK1 = max(1,3*N1-1), LWORK2 = max(1,3*N2-1)
    contains
    subroutine read_param
        implicit none
        integer :: i
        
        open(1,file='./fortran_pec_input.dat',status='old')
        open(2,file='./PECs_so.dat',status='old')
        
        read(1,*)cof_pec(:4)
        read(1,*)cof_pec(5:8)
        read(1,*)cof_pec(9:12)
        read(1,*)cof_pec(13:16)
        
        read(1,*)cof_aso(:4)
        read(1,*)cof_aso(5:8)
        read(1,*)cof_aso(9:12)

        do i=1,ndat
            read(2,*)r(i)
        enddo
        close(1)
        close(2)

    end subroutine read_param

end module constants

program pec
    use constants
    implicit none
    real(kind=8) :: v(npec), a_so(naso), w(1+n1+n2)
    integer :: i

    call read_param
    
    open(10,file='output_pecs_f90.dat',status='unknown')
    open(11,file='soc_function_f90.dat',status='unknown')
    open(12,file='model_pecs.dat',status='unknown')
    
    do i = 1, ndat
        call potential(r(i),cof_pec,v)
        call aso_func(r(i),cof_aso,a_so)
        call so_pecs(v,a_so,w)
        write(10,'(10f16.8)')r(i),w
        write(11,'(5f16.8)')r(i),a_so
        write(12,'(5f16.8)')r(i),v
    enddo
    close(10)
    close(11)
    close(12)
    
end program pec


subroutine potential(R,param,v)
    use constants, only: npec, nparam_pec
    implicit none
    real(kind=8), intent(in) :: R, param(nparam_pec)
    real(kind=8), intent(out) :: v(npec)
    real(kind=8) :: A(npec), B(npec), beta1(npec), beta2(npec)
    real(kind=8), parameter :: Eh = 2.1947463d5, a0 = 0.5291772d0
    real(kind=8), parameter :: c6pi = 4067.6*Eh*a0**6, c6sig = 4661.5*Eh*a0**6
    ! real(kind=8), parameter :: c6pi = 3915.3*Eh*a0**6, c6sig = 4486.9*Eh*a0**6
    real(kind=8), parameter, dimension(npec) :: c6 = [c6sig, c6pi, c6pi, c6sig]
    real(kind=8), parameter, dimension(npec) :: c8 = 80.0d0*a0**2*c6
    real(kind=8) :: damp
    integer :: i
    
    A = param(:4)*Eh
    B = param(5:8)*Eh
    beta1 = param(9:12)/a0
    beta2 = param(13:16)/a0
    v = A*exp(-beta1*R) - B*exp(-beta2*R)
    do i=1,npec
        v(i) = v(i) - damp(R,beta1(i),6)*c6(i)/R**6 - damp(R,beta1(i),8)*c8(i)/R**8
    enddo
    return
end subroutine potential

real(kind=8) function damp(R,beta,n)
    implicit none
    real(kind=8), intent(in) :: R, beta
    integer, intent(in) :: n
    integer :: i, fac
    real(kind=8) :: ss

    ss = 0.0d0
    do i=0,n
        ss = ss + (beta*R)**i/real(fac(i))
    enddo
    damp = 1.0d0 - exp(-beta*R)*ss

end function damp

integer function fac(n)
    implicit none
    integer, intent(in) :: n
    integer :: i
    fac = 1
    if (n==0) then
        fac=1
    else
        do i=1,n
            fac=fac*i
        enddo
    endif
end function fac

subroutine aso_func(R,param,aso)
    use constants, only: a_so_yb, naso, nparam_aso
    implicit none
    real(kind=8), intent(in) :: R, param(nparam_aso)
    real(kind=8), intent(out) :: aso(naso)
    real(kind=8) :: r0(naso), c(naso), alpha(naso)
    
    r0 = param(:4)
    c = param(5:8)
    alpha = param(9:12)
    ! aso = a_so_yb*(1.0d0 + c*exp(-alpha*(R-r0)**2))
    aso = a_so_yb + c*(1.0d0-tanh(alpha*(R-r0)))
    return    
end subroutine aso_func

subroutine so_pecs(v,aso,w)
    use constants, only: npec, naso, n1, n2, LWORK1, LWORK2
    implicit none
    real(kind=8), intent(in) :: v(npec), aso(naso)
    real(kind=8) :: w1(N1), work1(LWORK1), w2(N2), work2(LWORK2)
    real(kind=8) :: hso1(n1,n1),hso2(n2,n2)
    real(kind=8) :: ham1(n1,n1),ham2(n2,n2)
    real(kind=8) :: u1(n1), u2(n2)
    real(kind=8) :: a1,a2,a3,a4
    integer :: nst = 1 + n1 + n2, info, i
    real(kind=8), intent(out) :: w(nst)

    a1 = aso(1)
    a2 = aso(2)
    a3 = aso(3)
    a4 = aso(4)

    hso1(1,1) = a3/3.0d0
    hso1(1,2) = sqrt(2.0d0*a2*a3)/3.0d0
    hso1(1,3) = sqrt(2.0d0*a3*a4/3.0d0)
    hso1(2,2) = 2.0d0*a2/3.0d0
    hso1(2,3) = -sqrt(a2*a4/3.0d0)
    hso1(3,3) = 0.0d0

    do i=2,N1
        hso1(i,:i-1) = hso1(:i-1,i)
    enddo

    hso2(1,1) = -a3/3.0d0
    hso2(1,2) = sqrt(2.0d0*a2*a3)/3.0d0
    hso2(1,3) = 2.0d0*sqrt(2.0d0*a3*a4)/3.0d0
    hso2(1,4) = sqrt(a1*a3)/3.0d0
    hso2(1,5) = 0.0d0 
    hso2(2,2) = -2.0d0*a2/3.0d0
    hso2(2,3) = -sqrt(a2*a4)/3.0d0
    hso2(2,4) = 2.0d0*sqrt(2.0d0*a1*a2)/3.0d0
    hso2(2,5) = 0.0d0
    hso2(3,3) = 0.0d0
    hso2(3,4) = 0.0d0
    hso2(3,5) = sqrt(2.0d0*a3*a4/3.0d0)
    hso2(4,4) = 0.0
    hso2(4,5) = -sqrt(1.0d0*a1*a3/3.0d0)
    hso2(5,5) = -a3

    do i=2,N2
        hso2(i,:i-1) = hso2(:i-1,i)
    enddo

    u1 = [v(3),v(2),v(4)]
    u2 = [v(3),v(2),v(4),v(1),v(3)]
    ham1 = 0.0d0
    ham2 = 0.0d0
    forall(i=1:n1) ham1(i,i) = u1(i)
    forall(i=1:n2) ham2(i,i) = u2(i)
    ham1 = ham1 + hso1
    ham2 = ham2 + hso2
    call dsyev ('V', 'U', N1, ham1, N1, w1, WORK1, LWORK1, INFO)
    call dsyev ('V', 'U', N2, ham2, N2, w2, WORK2, LWORK2, INFO)
    w(1) = v(3)+aso(3)
    w(2:4) = w1
    w(5:9) = w2

    return
end subroutine so_pecs