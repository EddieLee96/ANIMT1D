program MT_ANISO_1D_Test
        implicit none

        real(8), parameter :: rkPi=3.141592653589793D0
        real(8), parameter :: rkMu=4.D0*rkPi*1.D-7
        integer, parameter :: n_layer=4
        complex(8), dimension(3) :: E_source
        real(8), dimension(0:n_layer,3,3) :: SigmaTensor
        real(8), dimension(0:n_layer) :: z_b
        integer :: ifreq, nfreq 
        real(8) :: minf, maxf, rfreq, w
        integer, parameter :: n_z=1
        real(8), dimension(n_z) :: z
        complex(8), dimension(n_z) :: Ex_z1, Ey_z1, Ez_z1, Hx_z1, Hy_z1, Hz_z1
        complex(8), dimension(n_z) :: Ex_z2, Ey_z2, Ez_z2, Hx_z2, Hy_z2, Hz_z2
        complex(8) :: Zxx, Zxy, Zyx, Zyy
        real(8) :: Rhoxy, Rhoyx, Phasexy, Phaseyx
        real(8), dimension(3,3) :: R1, R2, R3, R4, R5, R6, rPConTensor, rConTensor
        real(8) :: strike, dip, slant
        
        z_b(0)=-1000000.d0
        z_b(1)=0.d0
        z_b(2)=2000.d0
        z_b(3)=3680.d0
        z_b(4)=5680.d0

        !Air layer
        SigmaTensor(0,1,1) = 1.d0/10.d0**9 
        SigmaTensor(0,1,2) = 0.d0
        SigmaTensor(0,1,3) = 0.d0
        SigmaTensor(0,2,1) = SigmaTensor(1,1,2)
        SigmaTensor(0,2,2) = 1.d0/10.d0**9
        SigmaTensor(0,2,3) = 0.d0
        SigmaTensor(0,3,1) = SigmaTensor(1,1,3)
        SigmaTensor(0,3,2) = SigmaTensor(1,2,3)
        SigmaTensor(0,3,3) = 1.d0/10.d0**9

        !First layer
        SigmaTensor(1,1,1) = 1.d0/500.d0 
        SigmaTensor(1,1,2) = 0.d0
        SigmaTensor(1,1,3) = 0.d0
        SigmaTensor(1,2,1) = SigmaTensor(1,1,2)
        SigmaTensor(1,2,2) = 1.d0/500.d0
        SigmaTensor(1,2,3) = 0.d0
        SigmaTensor(1,3,1) = SigmaTensor(1,1,3)
        SigmaTensor(1,3,2) = SigmaTensor(1,2,3)
        SigmaTensor(1,3,3) = 1.d0/500.d0
 
        !Second layer
        strike=10.d0; dip=30.d0; slant=20.d0
        call RotZMat(R1, -strike)
        call RotXMat(R2, -dip)
        call RotZMat(R3, -slant)
        rPConTensor(:,:)=0.d0; rPConTensor(1,1)=1.d0/200.d0
        rPConTensor(2,2)=1.d0/1000.d0; rPConTensor(3,3)=1.d0/200.d0
        call RotZMat(R4, slant)
        call RotXMat(R5, dip)
        call RotZMat(R6, strike)
        rConTensor=matmul(matmul(matmul(matmul(matmul(matmul(R1,R2),R3),rPConTensor),R4),R5),R6)
        SigmaTensor(2,:,:)=rConTensor(:,:)
 
        !Third layer
        strike=10.d0; dip=20.d0; slant=30.d0
        call RotZMat(R1, -strike)
        call RotXMat(R2, -dip)
        call RotZMat(R3, -slant)
        rPConTensor(:,:)=0.d0; rPConTensor(1,1)=1.d0/1000.d0
        rPConTensor(2,2)=1.d0/200.d0; rPConTensor(3,3)=1.d0/1000.d0
        call RotZMat(R4, slant)
        call RotXMat(R5, dip)
        call RotZMat(R6, strike)
        rConTensor=matmul(matmul(matmul(matmul(matmul(matmul(R1,R2),R3),rPConTensor),R4),R5),R6)
        SigmaTensor(3,:,:)=rConTensor(:,:)
  
        !Fourth layer
        SigmaTensor(4,1,1) = 1.d0/500.d0 
        SigmaTensor(4,1,2) = 0.d0
        SigmaTensor(4,1,3) = 0.d0
        SigmaTensor(4,2,1) = 0.d0
        SigmaTensor(4,2,2) = 1.d0/500.d0
        SigmaTensor(4,2,3) = 0.d0
        SigmaTensor(4,3,1) = 0.d0
        SigmaTensor(4,3,2) = 0.d0
        SigmaTensor(4,3,3) = 1.d0/500.d0

        open(10, file='Ex_z_per1.dat')
        open(11, file='Ey_z_per1.dat')
        open(12, file='Ez_z_per1.dat')
        open(13, file='Hx_z_per1.dat')
        open(14, file='Hy_z_per1.dat')
        open(15, file='Hz_z_per1.dat')

        open(20, file='Ex_z_per2.dat')
        open(21, file='Ey_z_per2.dat')
        open(22, file='Ez_z_per2.dat')
        open(23, file='Hx_z_per2.dat')
        open(24, file='Hy_z_per2.dat')
        open(25, file='Hz_z_per2.dat')

        open(30, file='Rhoxy_per.dat')
        open(31, file='Rhoyx_per.dat')
        open(32, file='Phasexy_per.dat')
        open(33, file='Phaseyx_per.dat')

        nfreq=17; minf=0.01d0; maxf=100.d0
        z(1)=0.d0
        do ifreq=1, nfreq
            rfreq=10.d0**((dlog10(maxf)-dlog10(minf))/(dfloat(nfreq)-dfloat(1))*dfloat(ifreq-1)+dlog10(minf))
            w=2.d0*rkPi*rfreq
            E_source(1)=dcmplx(1.d0,0.d0); E_source(2)=dcmplx(0.d0,0.d0); E_source(3)=dcmplx(0.d0,0.d0)
            call MT_ANISO_1d(E_source,n_layer,SigmaTensor,z_b,w,n_z,z,Ex_z1,Ey_z1,Ez_z1,Hx_z1,Hy_z1,Hz_z1)
            write(10,*) 1.d0/rfreq, Ex_z1(1)
            write(11,*) 1.d0/rfreq, Ey_z1(1)
            write(12,*) 1.d0/rfreq, Ez_z1(1)
            write(13,*) 1.d0/rfreq, Hx_z1(1)
            write(14,*) 1.d0/rfreq, Hy_z1(1)
            write(15,*) 1.d0/rfreq, Hz_z1(1)

            E_source(1)=dcmplx(0.d0,0.d0); E_source(2)=dcmplx(1.d0,0.d0); E_source(3)=dcmplx(0.d0,0.d0)
            call MT_ANISO_1d(E_source,n_layer,SigmaTensor,z_b,w,n_z,z,Ex_z2,Ey_z2,Ez_z2,Hx_z2,Hy_z2,Hz_z2)
            write(20,*) 1.d0/rfreq, Ex_z2(1)
            write(21,*) 1.d0/rfreq, Ey_z2(1)
            write(22,*) 1.d0/rfreq, Ez_z2(1)
            write(23,*) 1.d0/rfreq, Hx_z2(1)
            write(24,*) 1.d0/rfreq, Hy_z2(1)
            write(25,*) 1.d0/rfreq, Hz_z2(1)

            Zxx=(Ex_z1(1)*Hy_z2(1)-Ex_z2(1)*Hy_z1(1))/(Hx_z1(1)*Hy_z2(1)-Hx_z2(1)*Hy_z1(1))
            Zxy=-(Ex_z1(1)*Hx_z2(1)-Ex_z2(1)*Hx_z1(1))/(Hx_z1(1)*Hy_z2(1)-Hx_z2(1)*Hy_z1(1))
            Zyx=(Ey_z1(1)*Hy_z2(1)-Ey_z2(1)*Hy_z1(1))/(Hx_z1(1)*Hy_z2(1)-Hx_z2(1)*Hy_z1(1))
            Zyy=-(Ey_z1(1)*Hx_z2(1)-Ey_z2(1)*Hx_z1(1))/(Hx_z1(1)*Hy_z2(1)-Hx_z2(1)*Hy_z1(1))
            Rhoxy=1.d0/(w*rkMu)*(dble(Zxy)**2+dimag(Zxy)**2)
            Rhoyx=1.d0/(w*rkMu)*(dble(Zyx)**2+dimag(Zyx)**2)
            Phasexy=datan(dimag(Zxy)/dble(Zxy))*180.d0/rkPi
            Phaseyx=datan(dimag(Zyx)/dble(Zyx))*180.d0/rkPi
            write(30,*) 1.d0/rfreq, Rhoxy
            write(31,*) 1.d0/rfreq, Rhoyx
            write(32,*) 1.d0/rfreq, Phasexy
            write(33,*) 1.d0/rfreq, Phaseyx
        end do
        close(10); close(11); close(12); close(13); close(14); close(15)
        close(20); close(21); close(22); close(23); close(24); close(25)
        close(30); close(31); close(32); close(33)
end program

subroutine RotZMat(R, theta)
        implicit none
        real(8), intent(out) :: R(3,3)
        real(8), intent(in) :: theta
        real(8), parameter :: rkPi=3.141592653589793D0

        R(1,1) = dcos(theta*rkPi/180.d0)
        R(1,2) = dsin(theta*rkPi/180.d0)
        R(1,3) = 0.d0
        R(2,1) = -dsin(theta*rkPi/180.d0)
        R(2,2) = dcos(theta*rkPi/180.d0)
        R(2,3) = 0.d0
        R(3,1) = 0.d0
        R(3,2) = 0.d0
        R(3,3) = 1.d0
end subroutine RotZMat


subroutine RotXMat(R, theta)
        implicit none
        real(8), intent(out) :: R(3,3)
        real(8), intent(in) :: theta
        real(8), parameter :: rkPi=3.141592653589793D0

        R(1,1) = 1.d0
        R(1,2) = 0.d0
        R(1,3) = 0.d0
        R(2,1) = 0.d0
        R(2,2) = dcos(theta*rkPi/180.d0)
        R(2,3) = dsin(theta*rkPi/180.d0)
        R(3,1) = 0.d0
        R(3,2) = -dsin(theta*rkPi/180.d0)
        R(3,3) = dcos(theta*rkPi/180.d0)
end subroutine RotXMat
