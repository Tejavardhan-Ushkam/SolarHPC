! ============================================================
! MODULE: integrators
! PROJECT: SolarHPC
! PURPOSE: Leapfrog (symplectic) and RK4 integrators for the
!          N-body ODE system. Both are OpenMP-parallelised.
!          Energy and position error diagnostics included.
!          Compiled FOURTH. Depends on: solar_constants,
!          omp_timer, nbody_force.
! ============================================================
module integrators
  use iso_c_binding
  implicit none

contains

  ! ============================================================
  ! LEAPFROG INTEGRATOR (FIXED: NO NESTED OPENMP)
  ! ============================================================
  subroutine leapfrog_step(n, pos, vel, mass, dt)
    implicit none
    integer, intent(in) :: n
    real(8), intent(inout) :: pos(n,3), vel(n,3)
    real(8), intent(in) :: mass(n), dt

    real(8) :: acc(n,3)
    real(8) :: acc_new(n,3)
    integer :: i, j

    ! --- Step 1: compute initial acceleration ---
    call compute_forces(n, pos, mass, acc)

    ! --- Step 2: half-step velocity update ---
    do i = 1, n
      do j = 1, 3
        vel(i,j) = vel(i,j) + 0.5d0 * dt * acc(i,j)
      end do
    end do

    ! --- Step 3: full-step position update ---
    do i = 1, n
      do j = 1, 3
        pos(i,j) = pos(i,j) + dt * vel(i,j)
      end do
    end do

    ! --- Step 4: compute new acceleration ---
    call compute_forces(n, pos, mass, acc_new)

    ! --- Step 5: second half velocity update ---
    do i = 1, n
      do j = 1, 3
        vel(i,j) = vel(i,j) + 0.5d0 * dt * acc_new(i,j)
      end do
    end do

  end subroutine leapfrog_step


  ! ============================================================
  ! RK4 INTEGRATOR (FIXED: NO NESTED OPENMP)
  ! ============================================================
  subroutine rk4_step(n, pos, vel, mass, dt)
    implicit none
    integer, intent(in) :: n
    real(8), intent(inout) :: pos(n,3), vel(n,3)
    real(8), intent(in) :: mass(n), dt

    real(8) :: k1v(n,3), k2v(n,3), k3v(n,3), k4v(n,3)
    real(8) :: k1x(n,3), k2x(n,3), k3x(n,3), k4x(n,3)
    real(8) :: temp_pos(n,3), temp_vel(n,3)
    integer :: i, j

    ! --- k1 ---
    call compute_forces(n, pos, mass, k1v)
    do i = 1, n
      do j = 1, 3
        k1x(i,j) = vel(i,j)
      end do
    end do

    ! --- k2 ---
    do i = 1, n
      do j = 1, 3
        temp_pos(i,j) = pos(i,j) + 0.5d0 * dt * k1x(i,j)
        temp_vel(i,j) = vel(i,j) + 0.5d0 * dt * k1v(i,j)
      end do
    end do
    call compute_forces(n, temp_pos, mass, k2v)
    do i = 1, n
      do j = 1, 3
        k2x(i,j) = temp_vel(i,j)
      end do
    end do

    ! --- k3 ---
    do i = 1, n
      do j = 1, 3
        temp_pos(i,j) = pos(i,j) + 0.5d0 * dt * k2x(i,j)
        temp_vel(i,j) = vel(i,j) + 0.5d0 * dt * k2v(i,j)
      end do
    end do
    call compute_forces(n, temp_pos, mass, k3v)
    do i = 1, n
      do j = 1, 3
        k3x(i,j) = temp_vel(i,j)
      end do
    end do

    ! --- k4 ---
    do i = 1, n
      do j = 1, 3
        temp_pos(i,j) = pos(i,j) + dt * k3x(i,j)
        temp_vel(i,j) = vel(i,j) + dt * k3v(i,j)
      end do
    end do
    call compute_forces(n, temp_pos, mass, k4v)
    do i = 1, n
      do j = 1, 3
        k4x(i,j) = temp_vel(i,j)
      end do
    end do

    ! --- Final update ---
    do i = 1, n
      do j = 1, 3
        pos(i,j) = pos(i,j) + dt/6.0d0 * (k1x(i,j) + 2.0d0*k2x(i,j) + 2.0d0*k3x(i,j) + k4x(i,j))
        vel(i,j) = vel(i,j) + dt/6.0d0 * (k1v(i,j) + 2.0d0*k2v(i,j) + 2.0d0*k3v(i,j) + k4v(i,j))
      end do
    end do

  end subroutine rk4_step


  ! ============================================================
  ! ENERGY COMPUTATION (CHECK UNITS CAREFULLY)
  ! ============================================================
  function compute_energy(n, pos, vel, mass) result(E)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: pos(n,3), vel(n,3), mass(n)
    real(8) :: E

    integer :: i, j
    real(8) :: kinetic, potential, rij, G
    real(8) :: dx, dy, dz

    G = 6.67430d-11  ! SI units

    kinetic = 0.0d0
    potential = 0.0d0

    ! --- Kinetic Energy ---
    do i = 1, n
      kinetic = kinetic + 0.5d0 * mass(i) * &
        (vel(i,1)**2 + vel(i,2)**2 + vel(i,3)**2)
    end do

    ! --- Potential Energy ---
    do i = 1, n
      do j = i+1, n
        dx = pos(i,1) - pos(j,1)
        dy = pos(i,2) - pos(j,2)
        dz = pos(i,3) - pos(j,3)
        rij = sqrt(dx*dx + dy*dy + dz*dz)

        if (rij > 0.0d0) then
          potential = potential - G * mass(i) * mass(j) / rij
        end if
      end do
    end do

    E = kinetic + potential

  end function compute_energy

end module integrators
