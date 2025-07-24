
module Hamiltonian

  use bacs
  use functions
  use ranpack
  implicit none

contains
  !---real 2DE Hamiltonian ---
  subroutine real_2DE(q, Nx, Ny, t, H)
    implicit none
    integer,intent(in)::q,Nx,Ny
    real(kind=double),intent(in)::t
    complex(kind=double),intent(out)::H(2*Nx*Ny,2*Nx*Ny)

    integer :: i, j, k
    integer :: spin
    integer :: m1, m2, n1, n2
    real(kind=double) :: flux 
 
    flux = 1.0d0/q
    do spin = 1, 2
        do m1 = 0, Nx-1
           m2 = mod(m1+1,Nx)
           do n1 = 0, Ny-1
              n2 = mod(n1+1,Ny)


              !--(m,n)~(m,n+1)--
              
              i = index(m1,n1,spin,Nx,Ny)
              j = index(m1,n2,spin,Nx,Ny)
              k = index(m2,n1,spin,Nx,Ny)
              
              H(j,i)=t*exp(2.0d0*pi*img*flux*(m1))
              H(k,i)=t*1.0d0

           end do
        end do
     end do
   end subroutine real_2DE
   !---end real 2DE hamiltonian---
    


  !--- real rashba Hamiltonian---
  subroutine realrashba(q, Nx, Ny, tr, H) !(p/q's q, system size, matrix)
    implicit none
    integer, intent(in) :: q, Nx, Ny
    real(kind=double), intent(in) :: tr
    complex(kind=double), intent(inout) :: H(2*Nx*Ny,2*Nx*Ny)
    complex(kind=double), allocatable :: Pmat(:,:,:)
    integer :: i, j, k
    integer :: spin1, spin2
    integer :: m1, m2, n1, n2
    real(kind=double) :: flux

    call Pauli_Matrix(Pmat)
    flux = 1.0d0/q
    
    do spin1 = 1, 2
       do spin2 = 1, 2
          do m1 = 0, Nx-1
             do n1 = 0, Ny-1
                m2 = mod(m1+1,Nx)
                n2 = mod(n1+1,Ny)
                
                i = index(m1,n1,spin1,Nx,Ny)
                j = index(m1,n2,spin2,Nx,Ny)
                k = index(m2,n1,spin2,Nx,Ny)
                
                H(k,i)=tr*img*Pmat(2,spin2,spin1)
                H(j,i)=-tr*img*Pmat(1,spin2,spin1)*exp(twopi*img*m1*flux)
                
             end do
          end do
       end do
    end do
    
  end subroutine realrashba

  !---end real rashba Hamiltonian---

  !---twisted boundary comdition real 2DE Hamiltonian ---
  subroutine tw_real_2DE(p, q, Nx, Ny, theta_x, theta_y, t, H)
    implicit none
    integer, intent(in) ::p, q, Nx, Ny
    real(kind=double), intent(in) :: t
    real(kind=double), intent(in) :: theta_x, theta_y
    complex(kind=double), intent(out) :: H(:,:)

    integer :: i, j, k
    integer :: spin
    integer :: m1, m2, n1, n2
    real(kind=double) :: flux

 
    flux = 1.0d0*p/q
    do spin = 1, 2
        do m1 = 0, Nx-1
           do n1 = 0, Ny-1
              m2 = mod(m1+1,Nx)
              n2 = mod(n1+1,Ny)

              !--(m,n)~(m,n+1)--
              
              i = index(m1,n1,spin,Nx,Ny)
              j = index(m1,n2,spin,Nx,Ny)
              k = index(m2,n1,spin,Nx,Ny)
              
              H(j,i)=t*exp(2.0d0*pi*img*flux*(m1))
              H(k,i)=t*1.0d0

              if((m1+1)/Nx == 1)then
                 H(k,i)=H(k,i)*exp(img*theta_x)
              end if
              
              if((n1+1)/Ny == 1)then
                 H(j,i)=H(j,i)*exp(img*theta_y)
              end if
              
           end do
        end do
     end do
   end subroutine tw_real_2DE
      !---end real 2DE hamiltonian---
    


   !---twisted boundary condition real rashba Hamiltonian---
   subroutine tw_realrashba(p, q, Nx, Ny, theta_x, theta_y, tr, H) !(p/q's q, system size, matrix)
    implicit none
    integer, intent(in) ::p, q, Nx, Ny
    real(kind=double), intent(in) :: tr
    real(kind=double), intent(in) :: theta_x, theta_y
    complex(kind=double), intent(inout) :: H(:,:)
    complex(kind=double), allocatable :: Pmat(:,:,:)
    integer :: i, j, k
    integer :: spin1, spin2
    integer :: m1, m2, n1, n2
    real(kind=double) :: flux

    call Pauli_Matrix(Pmat)
    flux = 1.0d0*p/q
      do spin1 = 1, 2
       do spin2 = 1, 2
          do m1 = 0, Nx-1
             do n1 = 0, Ny-1
                m2 = mod(m1+1,Nx)
                n2 = mod(n1+1,Ny)
                
                i = index(m1,n1,spin1,Nx,Ny)
                j = index(m1,n2,spin2,Nx,Ny)
                k = index(m2,n1,spin2,Nx,Ny)

                
                H(k,i)=tr*img*Pmat(2,spin2,spin1)
                H(j,i)=-tr*img*Pmat(1,spin2,spin1)*exp(2*pi*img*m1*flux)

                if((m1+1)/Nx == 1)then
                   H(k,i)=H(k,i)*exp(img*theta_x)
                end if
                
                if((n1+1)/Ny == 1)then
                   H(j,i)=H(j,i)*exp(img*theta_y)
                end if
                
             end do
          end do
       end do
    end do

  end subroutine tw_realrashba
  
  !---end real rashba Hamiltonian---

  !---bloch real 2DE Hamiltonian---
  subroutine BlochH(p, q, lx, ly, kx, ky, H)
    implicit none
    integer :: i
    integer, intent(in) :: p, q, lx, ly, kx, ky
    real(kind=double) :: flux
    complex(kind=double), intent(out) :: H(:,:,:,:)

    flux = 1.0d0*p/q
    
    do i =1,q
       H(i,i,lx,ly)=2.0d0*cos(ky-2.0d0*i*pi*flux)
    end do
    
    do i=1,q-1
       H(i,i+1,lx,ly)=1.0d0
    end do
    
    do i=2,q
       H(i,i-1,lx,ly)=1.0d0
    end do
    

    H(1,q,lx,ly)=exp(-img*kx)
    H(q,1,lx,ly)=exp(img*kx)
    
    
  end subroutine BlochH

  !---Bloch 2DE Hamiltonian ---

  !---Bloch 2DE Hamiltonian (w spin)
  subroutine Bloch_2DE(q,Nx,Ny,lx,ly,t,H)
    implicit none
    integer, intent(in) ::  q, Nx, Ny, lx, ly
    integer :: spin1, i, j
    integer :: m, n 
    real(kind=double), intent(in) :: t
    real(kind=double) :: flux
    real(kind=double)::kx
    real(kind=double)::ky
    complex(kind=double),intent(inout) :: H(2*q,2*q)
    
    flux=1.0d0/q
    kx = twopi/(Nx/q)*lx
    ky = twopi/Ny*ly
    
    do spin1 = 1, 2   !spin1=0...up, spin=1....down                     
       do m = 0, q-1
          i = indexB2(q,m,spin1)
          H(i,i)=t*exp(-img*ky)*exp(2*pi*img*flux*m)
       end do
      
       do m = 0, q-2
          i = indexB2(q,m+1,spin1)
          j = indexB2(q,m,spin1)
          H(i,j)=1.0d0*t
       end do
 
       H(indexB2(q,0,spin1),indexB2(q,q-1,spin1))=t*exp(-img*kx)
    end do

  end subroutine Bloch_2DE
  
  
  !---Bloch rashba Hamiltonian---

  subroutine Bloch_rashba(q,Nx,Ny,lx,ly,tr,H)
    implicit none
    integer, intent(in) :: q, Nx, Ny, lx, ly
    real(kind=double), intent(in) :: tr
    complex(kind=double), intent(inout) :: H(2*q,2*q)

    real(kind=double)::kx
    real(kind=double)::ky
    
    integer :: i, j, spin1, spin2
    integer :: m, n 
    real(kind=double)::flux
    complex(kind=double), allocatable :: Pmat(:,:,:)

    flux=1.0d0/q
    kx=twopi/(Nx/q)*lx
    ky=twopi/Ny*ly
    
    call Pauli_Matrix(Pmat)
    do spin1 = 1,2
       do spin2 = 1,2
          do m = 0, q-1
             i=indexB2(q,m,spin1)
             j=indexB2(q,m,spin2)
             H(i,j)=-tr*img*Pmat(1,spin1,spin2)*exp(-img*ky)*exp(2*pi*img*flux*m)
          end do
          
          do m = 0, q-2
             i = indexB2(q,m+1,spin1)
             j = indexB2(q,m,spin2)
             H(i,j)=tr*img*Pmat(2,spin1,spin2)
          end do
          
          H(indexB2(q,0,spin1),indexB2(q,q-1,spin2))=tr*img*Pmat(2,spin1,spin2)*exp(-img*kx)
       end do
    end do
  end subroutine Bloch_rashba

  !---2DE(bdg)---
  subroutine BdG_2DE(q,Nx,Ny,t,H)
    implicit none
    integer, intent(in) :: q, Nx, Ny
    real(kind=double), intent(in) :: t
    complex(kind=double), intent(inout) :: H(:,:)

    integer :: i,j,k
    integer :: spin, dagger
    integer :: m,n 
    integer :: m1, m2, n1, n2
    real(kind=double) :: flux

    flux = 1.0d0/q

    do dagger = 1, 2
       do spin = 1, 2
          do m1 = 0, Nx-1
             m2 = mod(m1+1,Nx)
             do n1 = 0, Ny-1
                n2 = mod(n1+1,Ny)
                i = index2(m1,n1,spin,dagger,Nx,Ny)
                j = index2(m2,n1,spin,dagger,Nx,Ny)
                k = index2(m1,n2,spin,dagger,Nx,Ny)

                H(k,i) = ((-1)**(dagger+1))*t*exp(((-1)**(dagger+1))*img*twopi*flux*m1)
                H(j,i) = ((-1)**(dagger+1))*t*1.0d0
             end do
          end do
       end do
    end do

  end subroutine BdG_2DE
  !---   ---

  !---Rashba(bdg)---
  subroutine BdG_Rashba(q,Nx,Ny,tr,H)
    implicit none
    integer,intent(in)::q,Nx,Ny
    real(kind=double),intent(in)::tr
    complex(kind=double),intent(inout)::H(:,:)
    complex(kind=double),allocatable::Pmat(:,:,:)
    integer :: i, j, k
    integer :: spin1, spin2
    integer :: dagger
    integer :: m1, m2, n1, n2
    real(kind=double) :: flux

    call Pauli_Matrix(Pmat)
    flux = 1.0d0/q
    do dagger = 1, 2
       do spin1 = 1, 2
          do spin2 = 1, 2
             do m1 = 0, Nx-1
                m2 = mod(m1+1,Nx)
                do n1 = 0, Ny-1
                   n2 = mod(n1+1,Ny)
                   
                   i = index2(m1,n1,spin1,dagger,Nx,Ny)
                   j = index2(m2,n1,spin2,dagger,Nx,Ny)
                   k = index2(m1,n2,spin2,dagger,Nx,Ny)

                   if(dagger == 1)then
                      H(j,i) = tr*img*Pmat(2,spin2,spin1)
                      H(k,i) = -tr*img*Pmat(1,spin2,spin1)*exp(img*twopi*flux*m1)
                   end if
                   if(dagger == 2)then
                      H(j,i) = tr*img*(Pmat(2,spin2,spin1))
                      H(k,i) = -tr*img*(Pmat(1,spin2,spin1))*exp(-twopi*img*flux*m1)
                   end if
                   
                end do
             end do
          end do
       end do
    end do

  end subroutine BdG_Rashba
  !--- ---

  subroutine BdG_SC(q,Nx,Ny,Delta_0,H)
    implicit none
    integer, intent(in) ::q, Nx, Ny
    real(kind=double), intent(in) :: Delta_0
    complex(kind=double), intent(inout) :: H(:,:)

    integer :: i, j, k
    integer :: m1, m2, n1, n2
    integer :: spin,spin1,spin2
    integer :: dagger
    integer :: sN 
    real(kind=double) :: flux
    complex(kind=double), allocatable :: Psi(:,:)
    complex(kind=double),allocatable::Pmat(:,:,:)
    character(50) :: filename,chq
    allocate(Psi(0:Nx-1,0:Ny-1))
    call Pauli_Matrix(Pmat)
    
    Psi = 0.0d0
    flux = 1.0d0/q
    sN = sqrt(Nx*Ny*2*flux)
    write(chq,'(i0)') q

    filename = "./result/vortex/vortex"//trim(chq)//".txt"
    open(100,file=filename)
    do m1 = 0, Nx-1
       do n1 = 0, Ny-1
          do k = -5, sN+5
             Psi(m1,n1) = Psi(m1,n1) + exp(twopi*img*sqrt(2.0d0*flux)*k*(n1))*exp(-2.0d0*pi*flux*(((m1)-1.0d0*k/sqrt(2.0d0*flux)))**2)
          end do
          write(100,*) m1, n1, abs(Psi(m1,n1))
       end do
       write(100,*)
    end do
    close(100)

    do spin1 = 1, 2
       do spin2 =1, 2
          do m1 = 0, Nx-1
             do n1 = 0, Ny-1
                i = index2(m1,n1,spin1,1,Nx,Ny)
                j = index2(m1,n1,spin2,2,Nx,Ny)
                H(i,j) = Delta_0*Psi(m1,n1)*Pmat(2,spin1,spin2)*img
             end do
          end do
       end do
    end do

    
  end subroutine BdG_SC
  !--- ---

  !---不純物---
  subroutine BdG_imp(Nx,Ny,iSeed,t,ls,Mu,H)                                                       
    implicit none
    integer, intent(in) :: Nx,Ny,iSeed
    real(kind=double), intent(in) :: ls,t,Mu
    complex(kind=double), intent(inout) :: H(:,:)

    integer :: i, j, k
    integer :: m1, n1
    integer :: spin
    integer :: dagger
    real(kind=double), allocatable :: V(:,:)

    allocate(V(0:Nx-1,0:Ny-1))
    
    call ran__dini(iSeed)
    open(1000,file="./result/Vimp.txt")
    do m1 = 0, Nx-1
       do n1 = 0, Ny-1
          V(m1,n1) = (ls*2)*ran__drnd()-ls
          write(1000,*) m1, n1, V(m1,n1)
       end do
       write(1000,*) 
    end do
    close(1000)
    
    do dagger = 1, 2
       do spin = 1, 2
          do m1 = 0, Nx-1
             do n1 = 0, Ny-1
                i = index2(m1,n1,spin,dagger,Nx,Ny)
                H(i,i)=((-1)**(dagger+1))*V(m1,n1)+eh(dagger)*4.0d0*t-eh(dagger)*Mu
             end do
          end do
       end do
    end do
    
  end subroutine BdG_imp
  !--- ---

  !---Bloch0---
  subroutine BdG_Bloch(q,Nx,Ny,t,tr,H)
    implicit none
    integer,intent(in)::q,Nx,Ny
    real(kind=double),intent(in)::t,tr
    complex(kind=double),intent(inout)::H(4*q,4*q,0:Nx/q-1,0:Ny-1)

    integer::i,j,k
    integer::m,m1,m2,n,n1,n2
    integer::spin,dagger
    integer::spin1,spin2,dagger1,dagger2
    integer::lx,ly
    real(kind=double)::kx,ky
    complex(kind=double), allocatable :: Pmat(:,:,:)
    complex(kind=double), allocatable :: H1(:,:),H2(:,:)
    real(kind=double)::flux
    real(kind=double)::del_kx
    real(kind=double)::del_ky

    allocate(H1(4*q,4*q),H2(4*q,4*q))

    flux=1.0d0/q

    H=0.0d0
    H1=0.0d0
    H2=0.0d0
    
    call Pauli_Matrix(Pmat)
    
    do lx = 0, Nx/q-1
       kx = twopi/(Nx/q)*lx
       do ly = 0, Ny-1
          ky = twopi/Ny*ly

          !---2DE---
          do dagger = 1, 2
             do spin = 1, 2
                do m = 0, q-1
                   i = indexB(q,m,spin,dagger)
                   H1(i,i)=eh(dagger)*t*exp(-eh(dagger)*img*ky)*exp(twopi*eh(dagger)*img*flux*(m))
                end do
                
                do m = 0,q-2
                   i = indexB(q,m+1,spin,dagger)
                   j = indexB(q,m,spin,dagger)
                   H1(i,j)=eh(dagger)*t
                end do
                
                H1(indexB(q,0,spin,dagger),indexB(q,q-1,spin,dagger)) = eh(dagger)*t*exp(-eh(dagger)*img*kx)
             end do
          end do
          !---***---
          !---rashba---
          do dagger = 1, 2
             do spin1 = 1, 2
                do spin2 = 1, 2
                   do m = 0, q-1
                      i = indexB(q,m,spin1,dagger)
                      j = indexB(q,m,spin2,dagger)

                      H2(i,j) = -tr*img*Pmat(1,spin1,spin2)*exp(-eh(dagger)*img*ky)*exp(twopi*eh(dagger)*img*flux*(m))
                   end do

                   do m = 0, q-2
                      i = indexB(q,m+1,spin1,dagger)
                      j = indexB(q,m,spin2,dagger)

                      H2(i,j) = tr*img*Pmat(2,spin1,spin2)
                   end do

                   H2(indexB(q,0,spin1,dagger),indexB(q,q-1,spin2,dagger))=tr*img*Pmat(2,spin1,spin2)*exp(-eh(dagger)*img*kx)
                end do
             end do
          end do
          !---***---

          do i = 1, 4*q
             do j = 1, 4*q
                H(i,j,lx,ly) = H1(i,j)+conjg(H1(j,i))+H2(i,j)+conjg(H2(j,i))
                
             end do
          end do

          do dagger = 1, 2
             do spin = 1, 2
                do m = 0, Nx-1
                   i = indexB(q,m,spin,dagger)
                   H(i,i,lx,ly)=H(i,i,lx,ly)+eh(dagger)*4.0d0*t
                end do
             end do
          end do
          
       end do
    end do

  end subroutine BdG_Bloch

  subroutine impurity(V,Nx,Ny,iSeed)
    implicit none
    integer, intent(in) :: Nx, Ny, iSeed
    integer :: m,n
    integer :: i,j,k
    complex(kind=double),intent(inout)::V(0:Nx-1,0:Ny-1)

    call ran__dini(iSeed)
    do m = 0, Nx-1
       do n = 0, Ny-1
          V(m,n) = 2*ran__drnd()-1
       end do
    end do
  end subroutine impurity

  subroutine Order_parameter(Psi,q,Nx,Ny)
    implicit none
    integer, intent(in) :: q, Nx, Ny
    complex(kind=double), intent(inout) :: Psi(0:Nx-1,0:Ny-1)

    integer :: i, j, k
    integer :: m1, m2, n1, n2
    integer :: spin
    integer :: dagger
    integer :: sN 
    real(kind=double) :: flux
    character(50) :: filename,chq

    Psi = 0.0d0
    flux = 1.0d0/q
    sN = sqrt(Nx*Ny*2*flux)
    write(chq,'(i0)') q

    filename = "./result/vortex"//trim(chq)//".txt"
    open(100,file=filename)
    do m1 = 0, Nx-1
       do n1 = 0, Ny-1
          do k = -5, sN+5
             Psi(m1,n1) = Psi(m1,n1) + exp(twopi*img*sqrt(2.0d0*flux)*k*(n1))*exp(-2.0d0*pi*flux*(((m1)-1.0d0*k/sqrt(2.0d0*flux)))**2)
          end do
          write(100,*) m1, n1, abs(Psi(m1,n1))
       end do
       write(100,*)
    end do
    close(100)

  end subroutine Order_Parameter


  !---tw-2DE(bdg)---
  subroutine BdG_2DE_tw(q,Nx,Ny,theta_x,theta_y,t,H)
    implicit none
    integer, intent(in) :: q, Nx, Ny
    real(kind=double), intent(in) :: t,theta_x,theta_y
    complex(kind=double), intent(inout) :: H(:,:)

    integer :: i,j,k
    integer :: spin, dagger
    integer :: m,n 
    integer :: m1, m2, n1, n2
    real(kind=double) :: flux

    flux = 1.0d0/q

    do dagger = 1, 2
       do spin = 1, 2
          do m1 = 0, Nx-1
             m2 = mod(m1+1,Nx)
             do n1 = 0, Ny-1
                n2 = mod(n1+1,Ny)
                i = index2(m1,n1,spin,dagger,Nx,Ny)
                j = index2(m2,n1,spin,dagger,Nx,Ny)
                k = index2(m1,n2,spin,dagger,Nx,Ny)

                H(k,i) = ((-1)**(dagger+1))*t*exp(((-1)**(dagger+1))*img*twopi*flux*m1)
                H(j,i) = ((-1)**(dagger+1))*t*1.0d0

                !---捻り境界条件---
                if(m2==0)then
                   H(j,i)=H(j,i)*exp(eh(dagger)*img*theta_x)
                end if
                if(n2==0)then
                   H(k,i)=H(k,i)*exp(eh(dagger)*img*theta_y)
                end if
                !------
                                
             end do
          end do
       end do
    end do

  end subroutine BdG_2DE_tw
  
  !---   ---

  !---tw-Rashba(bdg)---
  subroutine BdG_Rashba_tw(q,Nx,Ny,theta_x,theta_y,tr,H)
    implicit none
    integer,intent(in)::q,Nx,Ny
    real(kind=double),intent(in)::tr,theta_x,theta_y
    complex(kind=double),intent(inout)::H(:,:)
    complex(kind=double),allocatable::Pmat(:,:,:)
    integer :: i, j, k
    integer :: spin1, spin2
    integer :: dagger
    integer :: m1, m2, n1, n2
    real(kind=double) :: flux

    call Pauli_Matrix(Pmat)
    flux = 1.0d0/q
    do dagger = 1, 2
       do spin1 = 1, 2
          do spin2 = 1, 2
             do m1 = 0, Nx-1
                m2 = mod(m1+1,Nx)
                do n1 = 0, Ny-1
                   n2 = mod(n1+1,Ny)
                   
                   i = index2(m1,n1,spin1,dagger,Nx,Ny)
                   j = index2(m2,n1,spin2,dagger,Nx,Ny)
                   k = index2(m1,n2,spin2,dagger,Nx,Ny)

                   H(j,i) = tr*img*Pmat(2,spin2,spin1)
                   H(k,i) = -tr*img*Pmat(1,spin2,spin1)*exp(eh(dagger)*img*twopi*flux*m1)

                   if(m2==0)then
                      H(j,i)=H(j,i)*exp(eh(dagger)*img*theta_x)
                   end if
                   if(n2==0)then
                      H(k,i)=H(k,i)*exp(eh(dagger)*img*theta_y)
                   end if
                   
                end do
             end do
          end do
       end do
    end do

  end subroutine BdG_Rashba_tw
  !--- ---

  
end module Hamiltonian

