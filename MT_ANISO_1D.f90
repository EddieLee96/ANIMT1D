subroutine MT_ANISO_1d(E_source, n_layer, SigmaTensor, z_b, w, n_z, z, Ex_z, Ey_z, Ez_z, Hx_z, Hy_z, Hz_z)
        !------------------------------------------------ z0, E_source
        !0=air, SigmaTensor(0,:,:)
        !------------------------------------------------ z1
        !1, SigmaTensor(1,:,:)
        !------------------------------------------------ z2
        !...
        !------------------------------------------------ z_l
        !l, Ex=D_lp*exp(k(z-z_l+1))+D_ln*exp(-k(z-z_l))
        !   Ey=Q*Ex
        !------------------------------------------------ z_l+1
        !...
        !------------------------------------------------ z_nlayer
        !n_layer, SigmaTensor(n_layer,:,:)
        !infinite half space
        !...
        
        implicit none

        real(8), parameter :: rkPi=3.141592653589793D0
        real(8), parameter :: rkMu=4.D0*rkPi*1.D-7
        real(8), parameter :: rkEps=1.D-9/(36.D0*rkPi)
        complex(8), dimension(3),  intent(in) :: E_source !Ex, Ey, Ez
        integer, intent(in) :: n_layer !not include air layer
        real(8), intent(in) :: SigmaTensor(0:n_layer,3,3) !include air layer
        real(8), intent(in) :: z_b(0:n_layer) !include air layer, layer boundary z
        real(8), intent(in) :: w
        integer, intent(in) :: n_z !the number of desired z
        real(8), dimension(n_z), intent(in) :: z !the desired depths
        complex(8), dimension(n_z), intent(inout)  :: Ex_z, Ey_z, Ez_z, Hx_z, Hy_z, Hz_z !the desired EM fields

        complex(8), parameter :: ci=dcmplx(0.d0,1.d0)
        complex(8), dimension(4,4) :: Mtot, M0_z0, Mtmpinv, Mtmp1, Mtmp2
        complex(8), dimension(0:n_layer,4,4) :: D_up
        complex(8), dimension(4,1) :: D_nlayer, Ftmp, Dtmp
        integer :: i, j, n_z_layer

        Mtot(:,:)=dcmplx(0.d0,0.d0)
        Mtot(1,1)=dcmplx(1.d0,0.d0); Mtot(2,2)=dcmplx(1.d0,0.d0)
        Mtot(3,3)=dcmplx(1.d0,0.d0); Mtot(4,4)=dcmplx(1.d0,0.d0)
        D_up(n_layer,:,:)=dcmplx(0.d0,0.d0)
        D_up(n_layer,1,1)=dcmplx(1.d0,0.d0); D_up(n_layer,2,2)=dcmplx(1.d0,0.d0)
        D_up(n_layer,3,3)=dcmplx(1.d0,0.d0); D_up(n_layer,4,4)=dcmplx(1.d0,0.d0)
        do i=n_layer-1, 0, -1
            if (i.eq.n_layer-1) then
             call Matrix_M(z_b(i+1),w,SigmaTensor(i+1,:,:),z_b(i+1),0.d0,1,Mtmp1)
            else
              call Matrix_M(z_b(i+1),w,SigmaTensor(i+1,:,:),z_b(i+1),z_b(i+2),0,Mtmp1)
            end if        
            call Matrix_M(z_b(i+1),w,SigmaTensor(i,:,:),z_b(i),z_b(i+1),0,Mtmp2)
        print'(8ES12.4)', Mtmp2
            Mtmpinv=matinv4(Mtmp2)
            D_up(i,:,:)=matmul(matmul(Mtmpinv,Mtmp1),D_up(i+1,:,:))
        end do
        call Matrix_M(z_b(0),w,SigmaTensor(0,:,:),z_b(0),z_b(1),0,M0_z0)
        Mtot=matmul(M0_z0,D_up(0,:,:))
        D_nlayer(:,1)=dcmplx(0.d0,0.d0)
        D_nlayer(2,1)=(Mtot(2,4)*E_source(1)-Mtot(1,4)*E_source(2))/(Mtot(1,2)*Mtot(2,4)-Mtot(1,4)*Mtot(2,2))
        D_nlayer(4,1)=(-Mtot(2,2)*E_source(1)+Mtot(1,2)*E_source(2))/(Mtot(1,2)*Mtot(2,4)-Mtot(1,4)*Mtot(2,2))

        do i=1, n_z
!       do i=1, 1
            n_z_layer=-1
            do j=1, n_layer
                if (z(i).le.z_b(j)) then
                    n_z_layer=j-1
                    exit
                end if
            end do
            if (n_z_layer.eq.-1) n_z_layer=n_layer
            if (n_z_layer.eq.n_layer) then
              call Matrix_M(z(i),w,SigmaTensor(n_z_layer,:,:),z_b(n_z_layer),0.d0,1,Mtmp1)
            else
              call Matrix_M(z(i),w,SigmaTensor(n_z_layer,:,:),z_b(n_z_layer),z_b(n_z_layer+1),0,Mtmp1)
            end if
            Dtmp=matmul(D_up(n_z_layer,:,:),D_nlayer)
            Ftmp=matmul(Mtmp1,Dtmp)
            Ex_z(i)=Ftmp(1,1)
            Ey_z(i)=Ftmp(2,1)
            Ez_z(i)=-(dcmplx(SigmaTensor(n_z_layer,3,1),0.d0)*Ex_z(i)&
                      +dcmplx(SigmaTensor(n_z_layer,3,2),0.d0)*Ey_z(i))&
                     /(dcmplx(SigmaTensor(n_z_layer,3,3),w*rkEps))
            Hx_z(i)=Ftmp(3,1)
            Hy_z(i)=Ftmp(4,1)
            Hz_z(i)=dcmplx(0.d0,0.d0)
        end do
  contains
    function matinv4(A) result(B)
      !! Performs a direct calculation of the inverse of a 4×4 matrix.
      complex(8), intent(in) :: A(4,4)   !! Matrix
      complex(8)             :: B(4,4)   !! Inverse matrix
      complex(8)             :: detinv
      complex(8)             :: det

      ! Calculate the inverse determinant of the matrix
      det = &
      (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
      -A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
      +A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
      -A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))
        print*, 'det: ', det
        detinv = 1.d0/det

      ! Calculate the inverse of the matrix
      B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
      B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
      B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
      B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
      B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
      B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
      B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
      B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
      B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
      B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
      B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
      B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
      B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
      B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
      B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
      B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    end function
end subroutine


subroutine Matrix_M(z, w, SigmaTensor, z_b_upper, z_b_lower, l_znlayer,  M)
        implicit none
        real(8), intent(in) :: z, w
        real(8), dimension(3,3), intent(in) :: SigmaTensor
        real(8), intent(in) :: z_b_upper, z_b_lower
        integer, intent(in) :: l_znlayer
        complex(8), dimension(4,4), intent(inout) :: M
        
        complex(8), parameter :: ci=dcmplx(0.d0,1.d0)
        real(8), parameter :: rkPi=3.141592653589793D0
        real(8), parameter :: mu=4.D0*rkPi*1.D-7
        real(8), parameter :: rkEps=1.D-9/(36.D0*rkPi)
        complex(8) :: r
        complex(8) :: Axx, Axy, Ayy
        complex(8) :: k1, k2, Q1, Q2

        Axx=dcmplx(SigmaTensor(1,1),w*rkEps)-dcmplx(Sigmatensor(1,3)*SigmaTensor(3,1),0.d0)&
                                             /dcmplx(SigmaTensor(3,3),w*rkEps)
        Axy=dcmplx(SigmaTensor(1,2),0.d0)-dcmplx(Sigmatensor(1,3)*SigmaTensor(3,2),0.d0)&
                                          /dcmplx(SigmaTensor(3,3),w*rkEps)
        Ayy=dcmplx(SigmaTensor(2,2),w*rkEps)-dcmplx(Sigmatensor(2,3)*SigmaTensor(3,2),0.d0)&
                                             /dcmplx(SigmaTensor(3,3),w*rkEps)

        r=dcmplx(0.d0,1.d0/(w*mu))

        if (l_znlayer.eq.1) then
          if (Axy.eq.dcmplx(0.d0,0.d0)) then
            k1=cdsqrt(ci*dcmplx(w*mu,0.d0)*Axx)
            k2=cdsqrt(ci*dcmplx(w*mu,0.d0)*Ayy)
            M(:,:)=dcmplx(0.d0,0.d0)
            M(1,1)=dcmplx(0.d0,0.d0)
            M(1,2)=cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(2,3)=dcmplx(0.d0,0.d0)
            M(2,4)=cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(3,3)=dcmplx(0.d0,0.d0)
            M(3,4)=r*k2*cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(4,1)=dcmplx(0.d0,0.d0)
            M(4,2)=-r*k1*cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
!       print'(4ES)', k1, k2
!       print'(8ES12.4)', M
          else
            k1=cdsqrt(ci*dcmplx(w*mu/2.d0,0.d0)*(Axx+Ayy+cdsqrt((Axx-Ayy)**2+4.d0*Axy**2)))
            k2=cdsqrt(ci*dcmplx(w*mu/2.d0,0.d0)*(Axx+Ayy-cdsqrt((Axx-Ayy)**2+4.d0*Axy**2)))
            Q1=(k1**2-ci*dcmplx(w*mu,0.d0)*Axx)/(ci*dcmplx(w*mu,0.d0)*Axy)
            Q2=(k2**2-ci*dcmplx(w*mu,0.d0)*Axx)/(ci*dcmplx(w*mu,0.d0)*Axy)

            M(1,1)=dcmplx(0.d0,0.d0)
            M(1,2)=cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(1,3)=dcmplx(0.d0,0.d0)
            M(1,4)=cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(2,1)=dcmplx(0.d0,0.d0)
            M(2,2)=Q1*cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(2,3)=dcmplx(0.d0,0.d0)
            M(2,4)=Q2*cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(3,1)=dcmplx(0.d0,0.d0)
            M(3,2)=r*Q1*k1*cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(3,3)=dcmplx(0.d0,0.d0)
            M(3,4)=r*Q2*k2*cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(4,1)=dcmplx(0.d0,0.d0)
            M(4,2)=-r*k1*cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(4,3)=dcmplx(0.d0,0.d0)
            M(4,4)=-r*k2*cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
!       print'(4ES)', k1, k2
!       print'(4ES)', Q1, Q2
!       print'(8ES12.4)', M
          end if
        else if (l_znlayer.eq.0) then
          if (Axy.eq.dcmplx(0.d0,0.d0)) then
            k1=cdsqrt(ci*dcmplx(w*mu,0.d0)*Axx)
            k2=cdsqrt(ci*dcmplx(w*mu,0.d0)*Ayy)
            M(:,:)=dcmplx(0.d0,0.d0)
            M(1,1)=cdexp(k1*dcmplx(z-z_b_lower,0.d0))
            M(1,2)=cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(2,3)=cdexp(k2*dcmplx(z-z_b_lower,0.d0))
            M(2,4)=cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(3,3)=-r*k2*cdexp(k2*dcmplx(z-z_b_lower,0.d0))
            M(3,4)=r*k2*cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(4,1)=r*k1*cdexp(k1*dcmplx(z-z_b_lower,0.d0))
            M(4,2)=-r*k1*cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
!       print'(4ES)', k1, k2
!       print'(8ES12.4)', M
          else
            k1=cdsqrt(ci*dcmplx(w*mu/2.d0,0.d0)*(Axx+Ayy+cdsqrt((Axx-Ayy)**2+4.d0*Axy**2)))
            k2=cdsqrt(ci*dcmplx(w*mu/2.d0,0.d0)*(Axx+Ayy-cdsqrt((Axx-Ayy)**2+4.d0*Axy**2)))
            Q1=(k1**2-ci*dcmplx(w*mu,0.d0)*Axx)/(ci*dcmplx(w*mu,0.d0)*Axy)
            Q2=(k2**2-ci*dcmplx(w*mu,0.d0)*Axx)/(ci*dcmplx(w*mu,0.d0)*Axy)

            M(1,1)=cdexp(k1*dcmplx(z-z_b_lower,0.d0))
            M(1,2)=cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(1,3)=cdexp(k2*dcmplx(z-z_b_lower,0.d0))
            M(1,4)=cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(2,1)=Q1*cdexp(k1*dcmplx(z-z_b_lower,0.d0))
            M(2,2)=Q1*cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(2,3)=Q2*cdexp(k2*dcmplx(z-z_b_lower,0.d0))
            M(2,4)=Q2*cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(3,1)=-r*Q1*k1*cdexp(k1*dcmplx(z-z_b_lower,0.d0))
            M(3,2)=r*Q1*k1*cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(3,3)=-r*Q2*k2*cdexp(k2*dcmplx(z-z_b_lower,0.d0))
            M(3,4)=r*Q2*k2*cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
            M(4,1)=r*k1*cdexp(k1*dcmplx(z-z_b_lower,0.d0))
            M(4,2)=-r*k1*cdexp(-k1*dcmplx(z-z_b_upper,0.d0))
            M(4,3)=r*k2*cdexp(k2*dcmplx(z-z_b_lower,0.d0))
            M(4,4)=-r*k2*cdexp(-k2*dcmplx(z-z_b_upper,0.d0))
!       print'(4ES)', k1, k2
!       print'(4ES)', Q1, Q2
!       print'(8ES12.4)', M
          end if
        end if
end subroutine
