!===============================================================================
! ildg_lime -- lectura de configuraciones de norma en formato ILDG (contenedor LIME)
!
!  Formato:
!    Un archivo ILDG es una secuencia de registros LIME.  Cada registro tiene una
!    cabecera de 144 bytes en BIG-ENDIAN:
!
!       offset  tamano  contenido
!       ------  ------  ---------------------------------------------
!            0       4  numero magico  0x456789AB
!            4       2  version LIME
!            6       1  bits: MB = 0x80 (begin message), ME = 0x40 (end)
!            7       1  reservado
!           8       8  longitud de los datos (entero de 64 bits)
!           16     128  nombre del tipo de registro, terminado en NUL
!
!    Los datos siguen inmediatamente y se rellenan con ceros hasta el siguiente
!    multiplo de 8 bytes.
!
!    Registros relevantes:
!       ildg-format       XML con field / precision / lx ly lz lt
!       ildg-binary-data  el campo de norma
!       ildg-data-LFN     el nombre logico del archivo
!       scidac-checksum   XML con suma / sumb (CRC32 de SciDAC)
!       xlf-info          texto libre de tmLQCD (plaquette, beta, kappa, mu...)
!
!  Orden de los datos binarios (todo BIG-ENDIAN IEEE):
!       t (mas lento), z, y, x, mu = 1..4 (x,y,z,t), fila, columna, (re,im)
!    es decir, 4*3*3*2 numeros por sitio.
!
!  Convencion de almacenamiento en memoria (Fortran, primer indice el mas rapido):
!       u(a, b, mu, x, y, z, t)  con a = fila, b = columna, mu = 1..4 -> x,y,z,t
!       los indices de sitio van de 1 a L (no de 0 a L-1)
!
!  Antonio Garcia -- 2026
!===============================================================================
module ildg_lime
   use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, &
                                            real32, real64, error_unit
   implicit none
   private

   integer, parameter :: dp = real64

   !> Metadatos extraidos del archivo
   type, public :: ildg_metadata
      integer                       :: lx = 0, ly = 0, lz = 0, lt = 0
      integer                       :: prec = 0          !< 32 o 64 bits
      character(len=64)             :: field = ''        !< normalmente 'su3gauge'
      character(len=:), allocatable :: lfn               !< ildg-data-LFN
      character(len=:), allocatable :: format_xml
      character(len=:), allocatable :: xlf_info
      logical                       :: has_checksum = .false.
      integer(int64)                :: suma = 0, sumb = 0   !< del archivo
      integer(int64)                :: calc_suma = 0, calc_sumb = 0 !< recalculado
      integer(int64)                :: data_bytes = 0
   contains
      procedure :: volume => meta_volume
   end type ildg_metadata

   public :: ildg_list_records
   public :: ildg_read_metadata
   public :: ildg_read_gauge
   public :: ildg_plaquette
   public :: ildg_unitarity_error
   public :: hex8

   integer(int32), parameter :: LIME_MAGIC  = int(z'456789AB', int32)
   integer(int64), parameter :: MASK32      = int(z'FFFFFFFF', int64)
   integer,        parameter :: HDR_BYTES   = 144

   logical        :: crc_ready = .false.
   integer(int64) :: crc_table(0:255)

contains

!-------------------------------------------------------------------------------
   pure integer function meta_volume(self)
      class(ildg_metadata), intent(in) :: self
      meta_volume = self%lx*self%ly*self%lz*self%lt
   end function meta_volume

!===============================================================================
!  Conversion big-endian -> nativo
!===============================================================================
   pure function be_i16(b) result(v)
      integer(int8), intent(in) :: b(2)
      integer(int32) :: v
      v = ior(ishft(iand(int(b(1), int32), 255), 8), iand(int(b(2), int32), 255))
   end function be_i16

   pure function be_i32(b) result(v)
      integer(int8), intent(in) :: b(4)
      integer(int32) :: v
      integer :: i
      v = 0
      do i = 1, 4
         v = ior(ishft(v, 8), iand(int(b(i), int32), 255))
      end do
   end function be_i32

   pure function be_i64(b) result(v)
      integer(int8), intent(in) :: b(8)
      integer(int64) :: v
      integer :: i
      v = 0_int64
      do i = 1, 8
         v = ior(ishft(v, 8), iand(int(b(i), int64), 255_int64))
      end do
   end function be_i64

   !> Lee un real de 64 bits big-endian a partir del byte b(1)
   pure function be_r64(b) result(v)
      integer(int8), intent(in) :: b(8)
      real(real64) :: v
      v = transfer(b(8:1:-1), v)
   end function be_r64

   pure function be_r32(b) result(v)
      integer(int8), intent(in) :: b(4)
      real(real32) :: v
      v = transfer(b(4:1:-1), v)
   end function be_r32

!===============================================================================
!  Recorrido de registros LIME
!===============================================================================
   !> Abre el archivo, recorre las cabeceras y devuelve tipo/posicion/longitud
   subroutine scan_records(fname, ntypes, types, offsets, lengths, stat, msg)
      character(len=*),  intent(in)  :: fname
      integer,           intent(out) :: ntypes
      character(len=128), allocatable, intent(out) :: types(:)
      integer(int64),     allocatable, intent(out) :: offsets(:), lengths(:)
      integer,           intent(out) :: stat
      character(len=:), allocatable, intent(out) :: msg

      integer                        :: unit, ios, i, nul
      integer(int64)                 :: fsize, pos, dlen
      integer(int8)                  :: hdr(HDR_BYTES)
      character(len=128)             :: tname
      character(len=128), allocatable :: t_tmp(:)
      integer(int64),     allocatable :: o_tmp(:), l_tmp(:)
      integer, parameter             :: MAXREC = 256

      stat = 0
      msg = ''
      ntypes = 0
      allocate (t_tmp(MAXREC), o_tmp(MAXREC), l_tmp(MAXREC))

      inquire (file=fname, size=fsize)
      if (fsize <= 0) then
         stat = 1; msg = 'archivo vacio o inexistente: '//trim(fname); return
      end if

      open (newunit=unit, file=fname, access='stream', form='unformatted', &
            status='old', action='read', iostat=ios)
      if (ios /= 0) then
         stat = 1; msg = 'no se pudo abrir '//trim(fname); return
      end if

      pos = 1_int64
      do while (pos + HDR_BYTES - 1 <= fsize)
         read (unit, pos=pos, iostat=ios) hdr
         if (ios /= 0) exit

         if (be_i32(hdr(1:4)) /= LIME_MAGIC) then
            stat = 2
            msg = 'numero magico LIME incorrecto; el archivo no parece ILDG/LIME'
            close (unit); return
         end if

         if (be_i16(hdr(5:6)) /= 1) then
            write (error_unit, '(a,i0)') &
               'AVISO: version LIME inesperada: ', be_i16(hdr(5:6))
         end if

         dlen = be_i64(hdr(9:16))

         tname = ''
         nul = 128
         do i = 1, 128
            if (hdr(16 + i) == 0_int8) then
               nul = i - 1; exit
            end if
            tname(i:i) = achar(iand(int(hdr(16 + i), int32), 255))
         end do

         ntypes = ntypes + 1
         if (ntypes > MAXREC) then
            stat = 3; msg = 'demasiados registros LIME'; close (unit); return
         end if
         t_tmp(ntypes) = tname
         o_tmp(ntypes) = pos + HDR_BYTES          ! primer byte de datos (base 1)
         l_tmp(ntypes) = dlen

         ! avanza: datos + relleno hasta multiplo de 8
         pos = pos + HDR_BYTES + dlen
         if (mod(dlen, 8_int64) /= 0) pos = pos + (8_int64 - mod(dlen, 8_int64))
      end do

      close (unit)

      if (ntypes == 0) then
         stat = 4; msg = 'no se encontro ningun registro LIME'; return
      end if

      allocate (types(ntypes), offsets(ntypes), lengths(ntypes))
      types = t_tmp(1:ntypes)
      offsets = o_tmp(1:ntypes)
      lengths = l_tmp(1:ntypes)
   end subroutine scan_records

   !> Imprime la tabla de registros del archivo
   subroutine ildg_list_records(fname, stat, msg)
      character(len=*), intent(in)  :: fname
      integer,          intent(out) :: stat
      character(len=:), allocatable, intent(out) :: msg

      integer                        :: n, i
      character(len=128), allocatable :: types(:)
      integer(int64),     allocatable :: offs(:), lens(:)

      call scan_records(fname, n, types, offs, lens, stat, msg)
      if (stat /= 0) return

      write (*, '(a)') 'Registros LIME:'
      write (*, '(a)') '   #  tipo                        offset        bytes'
      write (*, '(a)') '   -----------------------------------------------------'
      do i = 1, n
         write (*, '(i4,2x,a24,i12,i13)') i, types(i) (1:24), offs(i) - 1, lens(i)
      end do
   end subroutine ildg_list_records

!===============================================================================
!  Lectura de un registro como texto
!===============================================================================
   subroutine read_text_record(unit, offset, length, text)
      integer,        intent(in)  :: unit
      integer(int64), intent(in)  :: offset, length
      character(len=:), allocatable, intent(out) :: text

      integer(int8), allocatable :: raw(:)
      integer :: i, n, ios

      n = int(min(length, 1000000_int64))
      allocate (raw(n))
      read (unit, pos=offset, iostat=ios) raw
      if (ios /= 0) then
         text = ''; return
      end if
      allocate (character(len=n) :: text)
      do i = 1, n
         if (raw(i) == 0_int8) then
            text = text(1:i - 1); return
         end if
         text(i:i) = achar(iand(int(raw(i), int32), 255))
      end do
   end subroutine read_text_record

!===============================================================================
!  Extraccion muy simple de <tag>valor</tag>
!===============================================================================
   function xml_value(text, tag) result(val)
      character(len=*), intent(in)  :: text, tag
      character(len=:), allocatable :: val
      integer :: i, j, k

      val = ''
      i = index(text, '<'//tag//'>')
      if (i == 0) return
      i = i + len(tag) + 2
      j = index(text(i:), '</'//tag//'>')
      if (j == 0) return
      val = adjustl(trim(text(i:i + j - 2)))
      ! recorta espacios/saltos de linea al final
      k = len(val)
      do while (k > 0)
         if (val(k:k) == ' ' .or. val(k:k) == achar(10) .or. val(k:k) == achar(13) &
             .or. val(k:k) == achar(9)) then
            k = k - 1
         else
            exit
         end if
      end do
      val = val(1:k)
   end function xml_value

!===============================================================================
!  Metadatos
!===============================================================================
   subroutine ildg_read_metadata(fname, meta, stat, msg)
      character(len=*),     intent(in)  :: fname
      type(ildg_metadata),  intent(out) :: meta
      integer,              intent(out) :: stat
      character(len=:), allocatable, intent(out) :: msg

      integer                        :: n, i, unit, ios
      character(len=128), allocatable :: types(:)
      integer(int64),     allocatable :: offs(:), lens(:)
      character(len=:),   allocatable :: txt, v
      logical                        :: found_fmt, found_bin

      call scan_records(fname, n, types, offs, lens, stat, msg)
      if (stat /= 0) return

      open (newunit=unit, file=fname, access='stream', form='unformatted', &
            status='old', action='read', iostat=ios)
      if (ios /= 0) then
         stat = 1; msg = 'no se pudo reabrir '//trim(fname); return
      end if

      found_fmt = .false.
      found_bin = .false.

      do i = 1, n
         select case (trim(types(i)))

         case ('ildg-format')
            call read_text_record(unit, offs(i), lens(i), txt)
            meta%format_xml = txt
            found_fmt = .true.
            v = xml_value(txt, 'field');     if (len(v) > 0) meta%field = v
            v = xml_value(txt, 'precision'); if (len(v) > 0) read (v, *, iostat=ios) meta%prec
            v = xml_value(txt, 'lx');        if (len(v) > 0) read (v, *, iostat=ios) meta%lx
            v = xml_value(txt, 'ly');        if (len(v) > 0) read (v, *, iostat=ios) meta%ly
            v = xml_value(txt, 'lz');        if (len(v) > 0) read (v, *, iostat=ios) meta%lz
            v = xml_value(txt, 'lt');        if (len(v) > 0) read (v, *, iostat=ios) meta%lt

         case ('ildg-binary-data')
            meta%data_bytes = lens(i)
            found_bin = .true.

         case ('ildg-data-LFN')
            call read_text_record(unit, offs(i), lens(i), txt)
            meta%lfn = trim(txt)

         case ('xlf-info')
            call read_text_record(unit, offs(i), lens(i), txt)
            meta%xlf_info = txt

         case ('scidac-checksum')
            call read_text_record(unit, offs(i), lens(i), txt)
            v = xml_value(txt, 'suma')
            if (len(v) > 0) then
               meta%suma = hex2int(v); meta%has_checksum = .true.
            end if
            v = xml_value(txt, 'sumb')
            if (len(v) > 0) meta%sumb = hex2int(v)

         end select
      end do

      close (unit)

      if (.not. found_fmt) then
         stat = 5; msg = 'falta el registro ildg-format'; return
      end if
      if (.not. found_bin) then
         stat = 6; msg = 'falta el registro ildg-binary-data'; return
      end if
      if (meta%prec /= 32 .and. meta%prec /= 64) then
         stat = 7; msg = 'precision no soportada en ildg-format'; return
      end if
      if (min(meta%lx, meta%ly, meta%lz, meta%lt) <= 0) then
         stat = 8; msg = 'dimensiones invalidas en ildg-format'; return
      end if
   end subroutine ildg_read_metadata

   !> Convierte una cadena hexadecimal a entero de 64 bits
   pure function hex2int(s) result(v)
      character(len=*), intent(in) :: s
      integer(int64) :: v
      integer :: i, d
      character :: c
      v = 0_int64
      do i = 1, len_trim(s)
         c = s(i:i)
         select case (c)
         case ('0':'9'); d = iachar(c) - iachar('0')
         case ('a':'f'); d = iachar(c) - iachar('a') + 10
         case ('A':'F'); d = iachar(c) - iachar('A') + 10
         case default;   cycle
         end select
         v = iand(ishft(v, 4) + int(d, int64), MASK32)
      end do
   end function hex2int

   !> Formatea un entero como 8 digitos hexadecimales
   pure function hex8(v) result(s)
      integer(int64), intent(in) :: v
      character(len=8) :: s
      integer :: i, d
      character(len=16), parameter :: digits = '0123456789abcdef'
      integer(int64) :: t
      t = iand(v, MASK32)
      do i = 8, 1, -1
         d = int(iand(t, 15_int64))
         s(i:i) = digits(d + 1:d + 1)
         t = ishft(t, -4)
      end do
   end function hex8

!===============================================================================
!  Lectura del campo de norma
!===============================================================================
   !> Lee ildg-binary-data en u(3,3,4,lx,ly,lz,lt) y recalcula el checksum SciDAC
   subroutine ildg_read_gauge(fname, meta, u, stat, msg)
      character(len=*),     intent(in)    :: fname
      type(ildg_metadata),  intent(inout) :: meta
      complex(dp), allocatable, intent(out) :: u(:, :, :, :, :, :, :)
      integer,              intent(out)   :: stat
      character(len=:), allocatable, intent(out) :: msg

      integer                        :: n, i, unit, ios, irec
      character(len=128), allocatable :: types(:)
      integer(int64),     allocatable :: offs(:), lens(:)
      integer(int8),      allocatable :: buf(:)
      integer(int64)                 :: expect, pos, rank
      integer                        :: pb, site_bytes, slice_bytes
      integer                        :: x, y, z, t, mu, a, b, off
      real(dp)                       :: re, im

      stat = 0; msg = ''

      call scan_records(fname, n, types, offs, lens, stat, msg)
      if (stat /= 0) return

      irec = 0
      do i = 1, n
         if (trim(types(i)) == 'ildg-binary-data') then
            irec = i; exit
         end if
      end do
      if (irec == 0) then
         stat = 6; msg = 'falta el registro ildg-binary-data'; return
      end if

      pb = meta%prec/8                       ! bytes por numero real
      site_bytes = 4*3*3*2*pb                ! bytes por sitio
      expect = int(meta%volume(), int64)*int(site_bytes, int64)

      if (lens(irec) /= expect) then
         stat = 9
         msg = 'el tamano de ildg-binary-data no coincide con lx*ly*lz*lt y la precision'
         return
      end if

      allocate (u(3, 3, 4, meta%lx, meta%ly, meta%lz, meta%lt), stat=ios)
      if (ios /= 0) then
         stat = 10; msg = 'memoria insuficiente para el campo de norma'; return
      end if

      slice_bytes = meta%lx*meta%ly*meta%lz*site_bytes
      allocate (buf(slice_bytes), stat=ios)
      if (ios /= 0) then
         stat = 10; msg = 'memoria insuficiente para el buffer de lectura'; return
      end if

      open (newunit=unit, file=fname, access='stream', form='unformatted', &
            status='old', action='read', iostat=ios)
      if (ios /= 0) then
         stat = 1; msg = 'no se pudo reabrir '//trim(fname); return
      end if

      call crc32_init()
      meta%calc_suma = 0_int64
      meta%calc_sumb = 0_int64
      rank = 0_int64
      pos = offs(irec)

      do t = 1, meta%lt
         read (unit, pos=pos, iostat=ios) buf
         if (ios /= 0) then
            stat = 11; msg = 'lectura incompleta del campo de norma'
            close (unit); return
         end if
         pos = pos + int(slice_bytes, int64)

         off = 0
         do z = 1, meta%lz
            do y = 1, meta%ly
               do x = 1, meta%lx

                  ! checksum SciDAC: CRC32 por sitio, rotado segun el rango
                  call scidac_accum(buf(off + 1:off + site_bytes), rank, &
                                    meta%calc_suma, meta%calc_sumb)
                  rank = rank + 1_int64

                  do mu = 1, 4
                     do a = 1, 3            ! fila
                        do b = 1, 3         ! columna
                           if (pb == 8) then
                              re = be_r64(buf(off + 1:off + 8)); off = off + 8
                              im = be_r64(buf(off + 1:off + 8)); off = off + 8
                           else
                              re = real(be_r32(buf(off + 1:off + 4)), dp); off = off + 4
                              im = real(be_r32(buf(off + 1:off + 4)), dp); off = off + 4
                           end if
                           u(a, b, mu, x, y, z, t) = cmplx(re, im, dp)
                        end do
                     end do
                  end do

               end do
            end do
         end do
      end do

      close (unit)
   end subroutine ildg_read_gauge

!===============================================================================
!  Observables de control
!===============================================================================
   !> Plaquette media: (1/(6 V Nc)) sum_x sum_{mu<nu} Re Tr U_{mu nu}(x)
   function ildg_plaquette(u) result(p)
      complex(dp), intent(in) :: u(:, :, :, :, :, :, :)
      real(dp) :: p

      integer  :: dims(4), c(4), cmu(4), cnu(4)
      integer  :: x, y, z, t, mu, nu, volume
      complex(dp) :: m(3, 3)
      real(dp) :: s

      dims = [size(u, 4), size(u, 5), size(u, 6), size(u, 7)]
      volume = dims(1)*dims(2)*dims(3)*dims(4)
      s = 0.0_dp

      !$omp parallel do collapse(4) private(x,y,z,t,mu,nu,c,cmu,cnu,m) reduction(+:s)
      do t = 1, dims(4)
         do z = 1, dims(3)
            do y = 1, dims(2)
               do x = 1, dims(1)
                  c = [x, y, z, t]
                  do mu = 1, 3
                     do nu = mu + 1, 4
                        cmu = c; cmu(mu) = mod(c(mu), dims(mu)) + 1
                        cnu = c; cnu(nu) = mod(c(nu), dims(nu)) + 1
                        m = matmul(u(:, :, mu, c(1), c(2), c(3), c(4)), &
                                   u(:, :, nu, cmu(1), cmu(2), cmu(3), cmu(4)))
                        m = matmul(m, conjg(transpose( &
                                   u(:, :, mu, cnu(1), cnu(2), cnu(3), cnu(4)))))
                        m = matmul(m, conjg(transpose( &
                                   u(:, :, nu, c(1), c(2), c(3), c(4)))))
                        s = s + real(m(1, 1) + m(2, 2) + m(3, 3), dp)
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      p = s/(6.0_dp*real(volume, dp)*3.0_dp)
   end function ildg_plaquette

   !> Maximo de |U^dagger U - 1| sobre todos los enlaces (control de lectura)
   function ildg_unitarity_error(u) result(e)
      complex(dp), intent(in) :: u(:, :, :, :, :, :, :)
      real(dp) :: e
      integer :: x, y, z, t, mu, a, b
      complex(dp) :: m(3, 3)
      real(dp) :: d

      e = 0.0_dp
      do t = 1, size(u, 7)
         do z = 1, size(u, 6)
            do y = 1, size(u, 5)
               do x = 1, size(u, 4)
                  do mu = 1, 4
                     m = matmul(conjg(transpose(u(:, :, mu, x, y, z, t))), &
                                u(:, :, mu, x, y, z, t))
                     do a = 1, 3
                        do b = 1, 3
                           d = abs(m(a, b) - merge((1.0_dp, 0.0_dp), &
                                                   (0.0_dp, 0.0_dp), a == b))
                           if (d > e) e = d
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end function ildg_unitarity_error

!===============================================================================
!  CRC32 (IEEE 802.3, igual que zlib) y checksum SciDAC
!===============================================================================
   subroutine crc32_init()
      integer(int64) :: c
      integer :: n, k
      if (crc_ready) return
      do n = 0, 255
         c = int(n, int64)
         do k = 1, 8
            if (iand(c, 1_int64) /= 0_int64) then
               c = ieor(int(z'EDB88320', int64), ishft(c, -1))
            else
               c = ishft(c, -1)
            end if
            c = iand(c, MASK32)
         end do
         crc_table(n) = c
      end do
      crc_ready = .true.
   end subroutine crc32_init

   pure function crc32(bytes) result(crc)
      integer(int8), intent(in) :: bytes(:)
      integer(int64) :: crc, c
      integer :: i, idx
      c = MASK32
      do i = 1, size(bytes)
         idx = int(iand(ieor(c, iand(int(bytes(i), int64), 255_int64)), 255_int64))
         c = iand(ieor(crc_table(idx), ishft(c, -8)), MASK32)
      end do
      crc = iand(ieor(c, MASK32), MASK32)
   end function crc32

   pure function rotl32(x, n) result(y)
      integer(int64), intent(in) :: x
      integer,        intent(in) :: n
      integer(int64) :: y
      integer :: m
      m = modulo(n, 32)
      if (m == 0) then
         y = iand(x, MASK32)
      else
         y = iand(ior(ishft(iand(x, MASK32), m), ishft(iand(x, MASK32), m - 32)), MASK32)
      end if
   end function rotl32

   !> Acumula el checksum SciDAC de un sitio con rango (indice lexicografico) rank
   pure subroutine scidac_accum(bytes, rank, suma, sumb)
      integer(int8),  intent(in)    :: bytes(:)
      integer(int64), intent(in)    :: rank
      integer(int64), intent(inout) :: suma, sumb
      integer(int64) :: w
      w = crc32(bytes)
      suma = ieor(suma, rotl32(w, int(modulo(rank, 29_int64))))
      sumb = ieor(sumb, rotl32(w, int(modulo(rank, 31_int64))))
   end subroutine scidac_accum

end module ildg_lime
