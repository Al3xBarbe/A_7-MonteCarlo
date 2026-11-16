!media_error(suma, suma2, n_runs, media, error)
!metropolis(n, temperatura, n_mcs, n_runs, media_c0, error_c0, media_e, error_e, aceptados, tcpu)
module simulaciones
   use, intrinsic :: iso_fortran_env, only: real64, int64
   implicit none

contains
   !Configuracion_1 = patron 1110 (C0(0) = 0.25)
   !Configuracion_2 = patron 0001 (C0(0) = 0.25)
   !Configuracion_3 = patron 1010 (C0(0) = 1.00)
   !Configuracion_4 = patron 0101 (C0(0) = 1.00)
   subroutine configuraciones(n, sigma, texto)
      integer, intent(in) :: n
      integer, allocatable, intent(out) :: sigma(:)
      character(len=100), intent(in) :: texto
      integer :: i
      allocate(sigma(n))
      if (texto == 'Configuracion_1') then
            sigma = 1
            sigma(4:n:4) = 0
      elif (texto == 'Configuracion_2') then
            sigma = 0
            sigma(4:n:4) = 1
      elif (texto == 'Configuracion_3') then
            sigma = 1
            sigma(2:n:2) = 0
      elif (texto == 'Configuracion_4') then
            sigma = 0
            sigma(2:n:2) = 1
      else
            print *, 'Error: configuracion desconocida'
            stop
      end if


   
      !---------------------------------------------------------------------------
   ! De acumuladores (suma, suma de cuadrados) a media y error estandar
   !---------------------------------------------------------------------------
   subroutine media_error(suma, suma2, n_runs, media, error)
      real(real64), intent(in) :: suma(0:), suma2(0:)
      integer, intent(in) :: n_runs
      real(real64), intent(out) :: media(0:), error(0:)
      real(real64) :: nr
      nr = real(n_runs, real64)
      media = suma / nr
      error = sqrt((suma2/nr - media**2) / real(n_runs-1, real64))
   end subroutine media_error

   !---------------------------------------------------------------------------
   ! Monte Carlo convencional (un intento de flip por paso, con rechazo)
   ! Salidas: medias y errores de C0 (normalizada) y de la energia por espin
   !---------------------------------------------------------------------------
   subroutine metropolis(n, temperatura, n_mcs, n_runs, media_c0, error_c0, &
                         media_e, error_e, aceptados, tcpu)
      integer, intent(in) :: n, n_mcs, n_runs
      real(real64), intent(in) :: temperatura
      real(real64), intent(out) :: media_c0(0:n_mcs), error_c0(0:n_mcs)
      real(real64), intent(out) :: media_e(0:n_mcs), error_e(0:n_mcs)
      integer(int64), intent(out) :: aceptados
      real(real64), intent(out) :: tcpu

      integer, allocatable :: sigma(:)
      integer :: corrida, paso, intento, i, izq, der, n0, k, signo
      real(real64) :: beta, r, p, c0, energia, c0_inicial, t_ini, t_fin
      real(real64) :: suma_c0(0:n_mcs), suma_c0_cuad(0:n_mcs)
      real(real64) :: suma_e(0:n_mcs), suma_e_cuad(0:n_mcs)

      beta = 1.0_real64 / temperatura
      c0_inicial = 0.25_real64           ! C0(0) para el patron 1110

      allocate(sigma(n))
      suma_c0 = 0.0_real64; suma_c0_cuad = 0.0_real64
      suma_e = 0.0_real64;  suma_e_cuad = 0.0_real64
      aceptados = 0_int64
      call cpu_time(t_ini)

      do corrida = 1, n_runs
         call configuraciones(n, sigma, 'Configuracion_1')
         do paso = 0, n_mcs
            ! Medidas en t = paso (MCS)
            energia = -real(sum(sigma), real64) / real(n, real64)
            c0 = 0.0_real64
            do k = 1, n
               signo = 1
               if (mod(k,2) /= 0) signo = -1
               c0 = c0 + real(signo*(1-sigma(k)), real64)
            end do
            c0 = (c0 / real(n, real64)) / c0_inicial

            suma_c0(paso) = suma_c0(paso) + c0
            suma_c0_cuad(paso) = suma_c0_cuad(paso) + c0**2
            suma_e(paso) = suma_e(paso) + energia
            suma_e_cuad(paso) = suma_e_cuad(paso) + energia**2

            if (paso == n_mcs) exit

            ! Un MCS = N intentos
            do intento = 1, n
               call random_number(r)
               i = 1 + int(r * real(n, real64))
               izq = modulo(i-2, n) + 1
               der = modulo(i, n) + 1
               n0 = (1-sigma(izq)) + (1-sigma(der))
               p = 0.5_real64 * real(n0, real64) * exp(-beta*real(sigma(i), real64))
               call random_number(r)
               if (r < p) then
                  sigma(i) = 1 - sigma(i)
                  aceptados = aceptados + 1_int64
               end if
            end do
         end do
      end do

      call cpu_time(t_fin)
      tcpu = t_fin - t_ini
      call media_error(suma_c0, suma_c0_cuad, n_runs, media_c0, error_c0)
      call media_error(suma_e, suma_e_cuad, n_runs, media_e, error_e)
      deallocate(sigma)
   end subroutine metropolis

   !---------------------------------------------------------------------------
   ! n-fold way (rejection-free). Misma dinamica que metropolis.
   ! Clases: clase = 1 + 3*sigma + n0  (6 clases), con p = 0.5*n0*exp(-beta*sigma)
   !---------------------------------------------------------------------------
   subroutine nfold(n, temperatura, n_mcs, n_runs, media_c0, error_c0, &
                    media_e, error_e, eventos, tcpu)
      integer, intent(in) :: n, n_mcs, n_runs
      real(real64), intent(in) :: temperatura
      real(real64), intent(out) :: media_c0(0:n_mcs), error_c0(0:n_mcs)
      real(real64), intent(out) :: media_e(0:n_mcs), error_e(0:n_mcs)
      integer(int64), intent(out) :: eventos
      real(real64), intent(out) :: tcpu

      integer, allocatable :: sigma(:), cls(:), loc(:,:), look(:), nncls(:)
      integer :: corrida, i, izq, der, j, k, c, siguiente
      real(real64) :: beta, eb, r, r2, r3, q, t, m, tnuevo, acum, objetivo
      real(real64) :: esum, csum, c0_inicial, w(6), p(6), t_ini, t_fin
      real(real64) :: c0v, ev
      real(real64) :: suma_c0(0:n_mcs), suma_c0_cuad(0:n_mcs)
      real(real64) :: suma_e(0:n_mcs), suma_e_cuad(0:n_mcs)

      beta = 1.0_real64 / temperatura
      eb = exp(-beta)
      p = [0.0_real64, 0.5_real64, 1.0_real64, 0.0_real64, 0.5_real64*eb, eb]

      allocate(sigma(n), cls(n), loc(n,6), look(n), nncls(6))
      suma_c0 = 0.0_real64; suma_c0_cuad = 0.0_real64
      suma_e = 0.0_real64;  suma_e_cuad = 0.0_real64
      eventos = 0_int64
      call cpu_time(t_ini)

      do corrida = 1, n_runs
         call configuraciones(n, sigma, 'Configuracion_1')
         nncls = 0
         esum = real(sum(sigma), real64)
         csum = 0.0_real64
         do i = 1, n
            call anadir(i)
            csum = csum + real(signo(i)*(1-sigma(i)), real64)
         end do

         t = 0.0_real64          ! tiempo en intentos
         siguiente = 0           ! siguiente instante de medida (MCS)

         do while (siguiente <= n_mcs)
            w = real(nncls, real64) * p
            q = sum(w) / real(n, real64)

            ! Tiempo (en intentos) hasta el siguiente flip aceptado
            if (q <= 0.0_real64) then
               tnuevo = huge(1.0_real64)            ! configuracion congelada
            else if (q >= 1.0_real64) then
               tnuevo = t + 1.0_real64
            else
               call random_number(r)
               r = 1.0_real64 - r                   ! r en (0,1]
               m = floor(log(r)/log(1.0_real64 - q)) + 1.0_real64
               tnuevo = t + m
            end if

            ! Medir con la configuracion actual todos los instantes que caen antes del salto
            do while (siguiente <= n_mcs)
               if (real(siguiente, real64)*real(n, real64) >= tnuevo) exit
               ev = -esum / real(n, real64)
               c0v = (csum / real(n, real64)) / c0_inicial
               suma_c0(siguiente) = suma_c0(siguiente) + c0v
               suma_c0_cuad(siguiente) = suma_c0_cuad(siguiente) + c0v**2
               suma_e(siguiente) = suma_e(siguiente) + ev
               suma_e_cuad(siguiente) = suma_e_cuad(siguiente) + ev**2
               siguiente = siguiente + 1
            end do
            if (siguiente > n_mcs) exit

            ! Elegir clase con probabilidad proporcional a n_i*p_i
            call random_number(r2)
            objetivo = r2 * sum(w)
            acum = 0.0_real64
            k = 6
            do c = 1, 6
               acum = acum + w(c)
               if (objetivo < acum) then
                  k = c
                  exit
               end if
            end do
            ! Elegir un espin de esa clase
            call random_number(r3)
            j = 1 + int(r3 * real(nncls(k), real64))
            if (j > nncls(k)) j = nncls(k)
            i = loc(j, k)
            izq = modulo(i-2, n) + 1
            der = modulo(i, n) + 1

            ! Flip y reclasificacion de los 3 espines afectados
            call quitar(i); call quitar(izq); call quitar(der)
            if (sigma(i) == 0) then
               esum = esum + 1.0_real64
               csum = csum - real(signo(i), real64)
            else
               esum = esum - 1.0_real64
               csum = csum + real(signo(i), real64)
            end if
            sigma(i) = 1 - sigma(i)
            call anadir(i); call anadir(izq); call anadir(der)

            eventos = eventos + 1_int64
            t = tnuevo
         end do
      end do

      call cpu_time(t_fin)
      tcpu = t_fin - t_ini
      call media_error(suma_c0, suma_c0_cuad, n_runs, media_c0, error_c0)
      call media_error(suma_e, suma_e_cuad, n_runs, media_e, error_e)
      deallocate(sigma, cls, loc, look, nncls)

   contains

      integer function signo(kk)       ! +1 en posiciones pares, -1 en impares
         integer, intent(in) :: kk
         signo = 1
         if (mod(kk,2) /= 0) signo = -1
      end function signo

      subroutine anadir(s)             ! clasifica el espin s y lo mete en su lista
         integer, intent(in) :: s
         integer :: n0, cc, l, rr
         l  = modulo(s-2, n) + 1
         rr = modulo(s, n) + 1
         n0 = (1-sigma(l)) + (1-sigma(rr))
         cc = 1 + 3*sigma(s) + n0
         cls(s) = cc
         nncls(cc) = nncls(cc) + 1
         loc(nncls(cc), cc) = s
         look(s) = nncls(cc)
      end subroutine anadir

      subroutine quitar(s)             ! saca s de su lista (intercambia con el ultimo)
         integer, intent(in) :: s
         integer :: cc, jj, ultimo
         cc = cls(s)
         jj = look(s)
         ultimo = loc(nncls(cc), cc)
         loc(jj, cc) = ultimo
         look(ultimo) = jj
         nncls(cc) = nncls(cc) - 1
      end subroutine quitar

   end subroutine nfold

   !---------------------------------------------------------------------------
   ! Evolucion completa: ejecuta metropolis y nfold con los mismos parametros,
   ! escribe un CSV con ambos y compara resultados y tiempos de CPU.
   !---------------------------------------------------------------------------
   subroutine evolucion_completa(n, temperatura, n_mcs, n_runs, archivo)
      integer, intent(in) :: n, n_mcs, n_runs
      real(real64), intent(in) :: temperatura
      character(len=*), intent(in) :: archivo

      real(real64), allocatable :: mc0_m(:), ec0_m(:), me_m(:), ee_m(:)
      real(real64), allocatable :: mc0_n(:), ec0_n(:), me_n(:), ee_n(:)
      integer(int64) :: acept, eventos
      real(real64) :: beta, eb, e_eq, t_m, t_n, dif, sig, max_sig_c0, max_sig_e
      real(real64) :: c0_teo
      integer :: paso, unidad

      call comprobar(n, temperatura, n_mcs, n_runs)
      allocate(mc0_m(0:n_mcs), ec0_m(0:n_mcs), me_m(0:n_mcs), ee_m(0:n_mcs))
      allocate(mc0_n(0:n_mcs), ec0_n(0:n_mcs), me_n(0:n_mcs), ee_n(0:n_mcs))

      beta = 1.0_real64 / temperatura
      eb = exp(-beta)
      e_eq = -1.0_real64 / (1.0_real64 + eb)

      call metropolis(n, temperatura, n_mcs, n_runs, mc0_m, ec0_m, me_m, ee_m, acept, t_m)
      call nfold(n, temperatura, n_mcs, n_runs, mc0_n, ec0_n, me_n, ee_n, eventos, t_n)

      open(newunit=unidad, file=archivo, status='replace', action='write')
      write(unidad,'(a)') 'mcs,c0_metropolis,err_c0_metropolis,c0_nfold,err_c0_nfold,' // &
         'c0_teorico,e_metropolis,err_e_metropolis,e_nfold,err_e_nfold,e_equilibrio'

      max_sig_c0 = 0.0_real64
      max_sig_e = 0.0_real64
      do paso = 0, n_mcs
         c0_teo = exp(-real(paso, real64)*eb)
         write(unidad,'(i0,10(",",es18.10))') paso, mc0_m(paso), ec0_m(paso), &
            mc0_n(paso), ec0_n(paso), c0_teo, me_m(paso), ee_m(paso), &
            me_n(paso), ee_n(paso), e_eq

         ! Diferencia entre metodos en unidades de su error combinado
         dif = abs(mc0_m(paso) - mc0_n(paso))
         sig = sqrt(ec0_m(paso)**2 + ec0_n(paso)**2)
         if (sig > 1.0e-12_real64) max_sig_c0 = max(max_sig_c0, dif/sig)
         dif = abs(me_m(paso) - me_n(paso))
         sig = sqrt(ee_m(paso)**2 + ee_n(paso)**2)
         if (sig > 1.0e-12_real64) max_sig_e = max(max_sig_e, dif/sig)
      end do
      close(unidad)

      print *, '=== Evolucion completa: N =', n, ' T =', temperatura, &
               ' MCS =', n_mcs, ' runs =', n_runs
      print *, 'Fraccion aceptada (metropolis)    =', &
               real(acept, real64)/(real(n, real64)*real(n_mcs, real64)*real(n_runs, real64))
      print *, 'Flips totales (n-fold)            =', eventos
      print *, 'CPU metropolis (s)                =', t_m
      print *, 'CPU n-fold (s)                    =', t_n
      if (t_n > 0.0_real64) print *, 'Aceleracion (CPU metro / CPU nfold) =', t_m/t_n
      print *, 'Max |dif|/error combinado, C0     =', max_sig_c0
      print *, 'Max |dif|/error combinado, E      =', max_sig_e
      print *, '(valores de unos pocos sigmas indican que ambos metodos coinciden)'
      print *, 'Archivo: ', archivo

      deallocate(mc0_m, ec0_m, me_m, ee_m, mc0_n, ec0_n, me_n, ee_n)
   end subroutine evolucion_completa

end module simulaciones


program main
   use, intrinsic :: iso_fortran_env, only: real64
   use simulaciones, only: evolucion_completa
   implicit none

   integer, parameter :: N = 256
   real(real64), parameter :: T = 0.5_real64
   integer, parameter :: N_MCS = 40
   integer, parameter :: N_RUNS = 200

   character(len=20) :: str_N, str_T, str_MCS, str_RUNS
   character(len=150) :: archivo

   write(str_N, '(i0)') N
   write(str_T, '(f10.3)') T
   write(str_MCS, '(i0)') N_MCS
   write(str_RUNS, '(i0)') N_RUNS

   archivo = 'Comparacion_N' // trim(str_N) // '_T' // trim(adjustl(str_T)) // &
             '_MCS' // trim(str_MCS) // '_R' // trim(str_RUNS) // '.csv'

   call evolucion_completa(N, T, N_MCS, N_RUNS, trim(archivo))
end program main