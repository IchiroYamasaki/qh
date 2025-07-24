program qh                                           
  use bacs
  use functions
  use ranpack
  use Hamiltonian
  implicit none

  !---変数--- 
  integer::i,j,k
  integer::m,n
  integer::m1,m2,n1,n2
  integer::spin,spin1,spin2
  integer::dagger,dagger1,dagger2
  integer::lx,ly
  integer::iSeed, mus, ds
  real(kind=double)::d_delta

  real(kind=double)::kx,ky
  real(kind=double)::theta_x,theta_y
  real(kind=double):: filling

  !---定数---
  integer,parameter::q = 18
  integer,parameter::Nx=q,Ny=q
  integer,parameter::qN=4*Nx*Ny
  real(kind=double),parameter::flux=1.0d0/q
  real(kind=double),parameter::t = 1.0d0
  real(kind=double),parameter::tr = 0.3d0
  real(kind=double),parameter::W = 0.2d0
  real(kind=double)::Delta_0
  real(kind=double)::Mu

  !---Hamiltonian
  complex(kind=double)::H(qN,qN),H1(qN,qN),H2(qN,qN)
  complex(kind=double)::H_qh(qN,qN)
  complex(kind=double)::H_soc(qN,qN)
  complex(kind=double)::H_sc(qN,qN)
  complex(kind=double)::H_imp(qN,qN)

  !---対角化---      
  complex(kind=double)::Evec(qN,qN),Evec1(qN,qN),Evec2(qN,qN)
  real(kind=double)::Eval(qN),Eval1(qN),Eval2(qN)
  !------
  !---kubo formura---
  complex(kind=double)::vx(qN,qN),vy(qN,qN)
  complex(kind=double)::vx_qh(qN,qN),vy_qh(qN,qN)
  complex(kind=double)::vx_soc(qN,qN),vy_soc(qN,qN)
  real(kind=double)::sigma_xy
  real(kind=double)::chern
  real(kind=double),parameter::eta=1.0d-6
  !------
  !---パウリ行列---
  complex(kind=double),allocatable::Pmat(:,:,:)
  !------

  !---コマンド因数---
  character(32)::arg1,arg2,arg3,arg4,arg5

  !---character---
  character(50) :: filename, ch_tr, ch_t, ch_Delta_0, ch_q, ch_Mu, ch_W,ch_iSeed,ch_ds

  call Pauli_Matrix(Pmat)

  write(ch_q,'(i0)') q
  write(ch_t,'(F3.1)') t
  write(ch_tr,'(F3.1)') tr
  write(ch_W,'(F3.1)') W
  
  call get_command_argument(1,arg1)
  call get_command_argument(2,arg2)
  call get_command_argument(3,arg3)

  read(arg1,*) iSeed
  read(arg2,*) ds
  read(arg3,*) d_delta

  write(ch_iSeed,'(i0)') iSeed
  write(ch_ds,'(i0)') ds
  
  call BdG_2DE(q,Nx,Ny,t,H_qh(:,:))
  call BdG_Rashba(q,Nx,Ny,tr,H_soc(:,:))
  call BdG_imp(Nx,Ny,iSeed,t,W,Mu,H_imp(:,:))

  do i = 1, qN
     do j = 1, qN
        H1(i,j)=H_qh(i,j)+conjg(H_qh(j,i))&
             +H_soc(i,j)+conjg(H_soc(j,i))&
             +H_imp(i,j)
     end do
  end do

  call diagonalize(qN,H1(:,:),Eval1(:),Evec1(:,:))

  filename = "./result/Eval_q"//trim(ch_q)//"t"//trim(ch_t)//"tr"//trim(ch_tr)//"W"//trim(ch_W)//"iSeed"//trim(ch_iSeed)//".txt"


  open(10,file=filename)
  do i = 1, qN
     write(10,*) i, Eval1(i)
  end do
  close(10)

  
  Mu = 0.0d0!(Eval1(qN/2+q/2+1)+Eval1(qN/2+q/2))/2
  filling = 1.0d0*(q/2)/q
  write(ch_Mu,'(F3.1)') Mu
  print*, Mu


  
  Delta_0 = d_delta*ds
  write(ch_Delta_0,'(F3.1)') Delta_0
  !---Hamiltonian---                                                       
  call BdG_SC(q,Nx,Ny,Delta_0,H_sc(:,:))
  call BdG_imp(Nx,Ny,iSeed,t,W,Mu,H_imp(:,:))

  do i = 1, qN
     do j = 1, qN
        H(i,j)=H_qh(i,j)+conjg(H_qh(j,i))&
             +H_soc(i,j)+conjg(H_soc(j,i))&
             +H_sc(i,j)+conjg(H_sc(j,i))&
             +H_imp(i,j)
     end do
  end do

  call diagonalize(qN,H(:,:),Eval(:),Evec(:,:))

  filename = "./result/q"//trim(ch_q)//"t"//trim(ch_t)//"tr"//trim(ch_tr)//"W"//trim(ch_W)//"iSeed"//trim(ch_iSeed)//"Delta"//trim(ch_Delta_0)//"Mu"//trim(ch_Mu)//".txt"

  open(20,file=filename)
  do i = 1, qN
     write(10,*) i, Eval(i)
  end do
  close(10)

  !-速度演算子-                                                                 
  !---vx,vy---                                                                  
  vx=0.0d0
  vy=0.0d0
  vx_qh=0.0d0
  vy_qh=0.0d0
  vx_soc=0.0d0
  vy_soc=0.0d0
  sigma_xy=0.0d0

  !---v_qh                                                                      
  do dagger = 1, 2
     do spin1 = 1, 2
        do n1 = 0, Ny-1
           n2 = mod(n1+1,Ny)

           i = index2(0,n1,spin1,dagger,Nx,Ny)
           j = index2(q-1,n1,spin1,dagger,Nx,Ny)
           vx_qh(i,j)=-img*t

           do m1 = 0, q-1
              i = index2(m1,n2,spin1,dagger,Nx,Ny)
              j = index2(m1,n1,spin1,dagger,Nx,Ny)

              vy_qh(i,j)=-img*t*exp(eh(dagger)*twopi*img*flux*m1)
           end do
        end do
     end do
  end do

  !------                                                                       
  !---v_soc---                                                                  
  do dagger = 1, 2
     do spin1 = 1, 2
        do spin2 = 1, 2
           do n1 = 0, Ny-1
              n2 = mod(n1+1,Ny)

              i = index2(0,n1,spin1,dagger,Nx,Ny)
              j = index2(q-1,n1,spin2,dagger,Nx,Ny)
              vx_soc(i,j)=eh(dagger)*tr*Pmat(2,spin1,spin2)

              do m1 = 0, q-1
                 i = index2(m1,n2,spin1,dagger,Nx,Ny)
                 j = index2(m1,n1,spin2,dagger,Nx,Ny)
                 vy_soc(i,j)=-eh(dagger)*tr*exp(eh(dagger)*img*twopi*flux*m1)*Pmat(1,spin1,spin2)
              end do
           end do
        end do
     end do
  end do


  do i = 1, qN
     do j = 1, qN
        vx(i,j)=vx_qh(i,j)+conjg(vx_qh(j,i))+vx_soc(i,j)+conjg(vx_soc(j,i))
        vy(i,j)=vy_qh(i,j)+conjg(vy_qh(j,i))+vy_soc(i,j)+conjg(vy_soc(j,i))
     end do
  end do
  !------                                                                
   !---kubo formula---                                                   
     do i = qN/2-q+1, qN/2
        do j = 1, qN
           if(i/=j)then
              sigma_xy=sigma_xy+twopi*2*AIMAG(&
                   ((dot_product(Evec(:,i),matmul(vx(:,:),Evec(:,j)))&
                   *dot_product(Evec(:,j),matmul(vy(:,:),Evec(:,i))))))&
                   /((Eval(i)-Eval(j))**2)/(Nx/q)/Ny
           end if
        end do
     end do
     !------
     print*, sigma_xy

     filename = "./result/iSeed"//trim(ch_iSeed)//"_ds"//trim(ch_ds)//"/chern_iSeed"//trim(ch_iSeed)//"_d"//trim(ch_ds)//".txt"
     open(iSeed,file=filename)
     write(iSeed,*) sigma_xy
     close(iSeed)
  
end program qh



