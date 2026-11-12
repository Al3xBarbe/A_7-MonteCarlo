

! ACTIVIDAD 7 - TECNICAS FISICAS III
! Ising 1-facilitado (1D). Monte Carlo convencional.
!
! Compilar: gfortran -Wall -Wextra -fcheck=all -O2 MonteCarlo.f90 -o MonteCarlo
! Ejecutar: ./MonteCarlo (macOS/Linux) o .\MonteCarlo.exe (Windows)

module simulaciones
   use, intrinsic :: iso_fortran_env, only: real64, int64
   implicit none
   private
   public :: metropolis

contains

!Subroutina para realizar la simulacion de Monte Carlo del modelo de Ising 1-facilitado en 1D
!Entran: n (numero de espines), temperatura, n_mcs (numero de pasos de Monte Carlo), 
!n_runs (numero de corridas independientes), archivo (nombre del archivo para guardar resultados)
   subroutine metropolis(n, temperatura, n_mcs, n_runs, archivo)
      implicit none

      ! Parametros que recibe del programa principal
      integer, intent(in) :: n, n_mcs, n_runs 
      real(real64), intent(in) :: temperatura
      character(len=*), intent(in) :: archivo


      ! Cadena de espines y variables del algoritmo
      integer, allocatable :: sigma(:)
      integer :: corrida, paso, intento, i, izq, der, n0
      integer :: k, unidad, signo
      integer(int64) :: aceptados
      real(real64) :: beta, r, p, c0, energia, c0_inicial


      ! Acumuladores para promediar las simulaciones
      real(real64), allocatable :: suma_c0(:), suma_c0_cuad(:)
      real(real64), allocatable :: suma_e(:), suma_e_cuad(:)
      real(real64) :: media_c0, error_c0, media_e, error_e
      real(real64) :: c0_teorico, e_equilibrio, fraccion_aceptada


      ! Comprobaciones iniciales
      !Como la cadena elegida es con el patron 1110, el numero de espines debe ser multiplo de 4 y mayor o igual a 4
      if (n < 4 .or. mod(n,4) /= 0) error stop 'N debe ser multiplo de 4 y >= 4'
      !Temperatura debe ser positiva
      if (temperatura <= 0.0_real64) error stop 'T debe ser positiva'
      !Al menos un paso de Monte Carlo
      if (n_mcs < 1) error stop 'n_mcs debe ser >= 1'
      !Al menos dos simulaciones independientes para poder calcular errores
      if (n_runs < 2) error stop 'n_runs debe ser >= 2 para calcular errores'


      beta = 1.0_real64 / temperatura  ! Unidades k_B=1
      c0_inicial = 0.25_real64         ! Patron inicial 1110 repetido
      e_equilibrio = -1.0_real64 / (1.0_real64 + exp(-beta))


      !Guardo memoria para los espines y los acumuladores de las observables
      allocate(sigma(n))
      allocate(suma_c0(0:n_mcs), suma_c0_cuad(0:n_mcs))
      allocate(suma_e(0:n_mcs), suma_e_cuad(0:n_mcs))
      !Inicializo los acumuladores a cero
      suma_c0 = 0.0_real64
      suma_c0_cuad = 0.0_real64
      suma_e = 0.0_real64
      suma_e_cuad = 0.0_real64
      aceptados = 0_int64


      call random_number(r)


      !Bucle para repetir la simulacion desde la misma configuracion inicial para promediar los resultados
      do corrida = 1, n_runs
        !Inicializo la cadena de espines con el patron 1110 repetido
         sigma = 1
         sigma(4:n:4) = 0              

         !Bucle para realizar los pasos de Monte Carlo
         do paso = 0, n_mcs

            ! Medimos observables en t=0,1,...,n_mcs
            !Calculo de la energía por espin
            energia = -real(sum(sigma), real64) / real(n, real64)
            !Calculo de la fracción de espines en estado 0 
            c0 = 0.0_real64
            do k = 1, n
               !Signo positivo para los espines pares y negativo para los impares
               signo = 1
               if (mod(k,2) /= 0) signo = -1
               c0 = c0 + real(signo*(1-sigma(k)), real64)
            end do
            !Normalizo la fracción de espines en estado 0 por el valor inicial
            c0 = (c0 / real(n, real64)) / c0_inicial

            !Actualizo los acumuladores para promediar las simulaciones
            suma_c0(paso) = suma_c0(paso) + c0
            suma_c0_cuad(paso) = suma_c0_cuad(paso) + c0**2
            suma_e(paso) = suma_e(paso) + energia
            suma_e_cuad(paso) = suma_e_cuad(paso) + energia**2

            !Salimos del bucle si ya hemos llegado al último paso de Monte Carlo
            if (paso == n_mcs) exit


            ! Un MCS equivale a N intentos de cambio
            do intento = 1, n
                !Generamos un número aleatorio para seleccionar un espín al azar
               call random_number(r)
               i = 1 + int(r * real(n, real64))


               ! Calculo de los vecinos teniendo en cuenta las condiciones de controno
               izq = modulo(i-2, n) + 1
               der = modulo(i, n) + 1
               n0 = (1-sigma(izq)) + (1-sigma(der))


               ! Metropolis con restriccion 1-facilitada
               ! p = (n0/2)*min(1,exp(-beta*DeltaE))
               ! Con E=-sum(sigma), equivale a:
               p = 0.5_real64 * real(n0,real64) * exp(-beta*real(sigma(i),real64))

               !Generamos otro número aleatorio para decidir si aceptamos el cambio de espín
               call random_number(r)
               if (r < p) then
                  sigma(i) = 1 - sigma(i)
                  aceptados = aceptados + 1_int64
               end if
            end do
         end do
      end do


      ! Guardar promedios y errores en la carpeta de ejecucion
      open(newunit=unidad, file=archivo, status='replace', action='write')
      ! Guardamos los encabezados de las columnas en el archivo CSV
      write(unidad,'(a)') &
         'mcs,c0_mc,error_c0,c0_teorico,energia_mc,error_energia,energia_equilibrio'

      !Bucle para calcular los promedios y errores de las observables y guardarlos en el archivo CSV
      do paso = 0, n_mcs
         media_c0 = suma_c0(paso)/real(n_runs,real64)
         media_e = suma_e(paso)/real(n_runs,real64)

         ! Error estandar de la media entre las simulaciones independientes
         error_c0 = sqrt(max(0.0_real64,(suma_c0_cuad(paso)/real(n_runs,real64)-media_c0**2) /real(n_runs-1,real64)))
            
         error_e = sqrt(max(0.0_real64, (suma_e_cuad(paso)/real(n_runs,real64)-media_e**2) /real(n_runs-1,real64)))

         !Calculo de la fracción de espines en estado 0 teórica
         c0_teorico = exp(-real(paso,real64)*exp(-beta))

         write(unidad,'(i0,6(",",es18.10))') paso, media_c0, error_c0, c0_teorico, media_e, error_e, e_equilibrio
      end do
      close(unidad)

      !Calculo de la fracción de cambios aceptados, nos da una idea de la eficiencia del algoritmo
      fraccion_aceptada = real(aceptados,real64) /  (real(n,real64)*real(n_mcs,real64)*real(n_runs,real64))

    !Imprimimos los resultados finales en la terminal para comprobar que la simulación se ha completado correctamente
      print *, 'Simulacion completada. Archivo: ', archivo
      print *, 'N = ', n, ' T = ', temperatura
      print *, 'MCS = ', n_mcs, ' Simulaciones = ', n_runs
      print *, 'Fraccion de cambios aceptados = ', fraccion_aceptada
      print *, 'Energia final media por espin = ', media_e
      print *, 'Energia teorica de equilibrio = ', e_equilibrio

      deallocate(sigma, suma_c0, suma_c0_cuad, suma_e, suma_e_cuad)
   end subroutine metropolis

end module simulaciones



program main

    use, intrinsic :: iso_fortran_env, only: real64
    use simulaciones, only: metropolis
    implicit none

    ! Parametros de la simulacion
    integer, parameter :: N = 256
    real(real64), parameter :: T = 0.5_real64
    integer, parameter :: N_MCS = 5000
    integer, parameter :: N_RUNS = 200

    ! Variables para construir el nombre del fichero
    character(len=20) :: str_N, str_T, str_MCS, str_RUNS
    character(len=150) :: archivo

    ! Convertir los parametros a texto
    write(str_N, '(i0)') N
    write(str_T, '(f10.3)') T
    write(str_MCS, '(i0)') N_MCS
    write(str_RUNS, '(i0)') N_RUNS

    ! Construir automaticamente el nombre del CSV
    archivo = 'Metropolis_N' // trim(str_N) // '_T' // trim(adjustl(str_T)) // '_MCS' // trim(str_MCS) // '_R' // trim(str_RUNS) // '.csv'

    ! Ejecutar la simulacion
    call metropolis(N, T, N_MCS, N_RUNS, trim(archivo))

end program main
