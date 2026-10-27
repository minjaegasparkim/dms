!____________________________________________________________________________________________________
!
! 9th module : minimal FITS writer (no external library needed)
!____________________________________________________________________________________________________
!
! Writes a primary HDU with a 2-d image or a 3-d cube in IEEE big-endian double precision
! (BITPIX = -64), following the FITS standard 4.0:
!   * header: 80-character ASCII cards in blocks of 2880 bytes, padded with blanks
!   * data  : big-endian, padded with zeros to a multiple of 2880 bytes
! The byte order is swapped by hand on little-endian machines, so no compiler-specific
! CONVERT= option is required.
!
! Usage:
!   call fits_header_init(hdr)
!   call fits_card_real(hdr, 'CDELT1', -0.01_r2, 'pixel size [deg]') ...
!   call fits_write_cube(filename, data(nx, ny, nz), hdr)      (nz = 1 gives a 2-d image)
!____________________________________________________________________________________________________

module fits_io

  use datatype

  implicit none

  private
  public :: fits_header, fits_header_init, fits_card_real, fits_card_int, fits_card_str, &
            fits_card_comment, fits_write_cube

  integer, parameter :: max_cards = 400
  integer, parameter :: ib1 = selected_int_kind(2), ib4 = selected_int_kind(9)

  type fits_header
    integer                                   :: n = 0
    character(len=80), dimension(max_cards)   :: card
  end type fits_header

contains

  subroutine fits_header_init(hdr)
    type(fits_header), intent(out) :: hdr
    hdr%n = 0
  end subroutine fits_header_init


  subroutine add_card(hdr, line)
    type(fits_header), intent(inout) :: hdr
    character(len=*), intent(in)     :: line
    if (hdr%n >= max_cards) return
    hdr%n = hdr%n + 1
    hdr%card(hdr%n) = line
  end subroutine add_card


  subroutine fits_card_real(hdr, key, val, comment)
    type(fits_header), intent(inout) :: hdr
    character(len=*), intent(in)     :: key, comment
    real(r2), intent(in)             :: val
    character(len=80)                :: line
    character(len=20)                :: sval
    write(sval, '(es20.12e3)') val
    line = ''
    line(1:8)   = key
    line(9:10)  = '= '
    line(11:30) = adjustr(sval)
    line(31:33) = ' / '
    line(34:80) = comment
    call add_card(hdr, line)
  end subroutine fits_card_real


  subroutine fits_card_int(hdr, key, val, comment)
    type(fits_header), intent(inout) :: hdr
    character(len=*), intent(in)     :: key, comment
    integer, intent(in)              :: val
    character(len=80)                :: line
    character(len=20)                :: sval
    write(sval, '(i20)') val
    line = ''
    line(1:8)   = key
    line(9:10)  = '= '
    line(11:30) = sval
    line(31:33) = ' / '
    line(34:80) = comment
    call add_card(hdr, line)
  end subroutine fits_card_int


  subroutine fits_card_str(hdr, key, val, comment)
    type(fits_header), intent(inout) :: hdr
    character(len=*), intent(in)     :: key, val, comment
    character(len=80)                :: line
    character(len=70)                :: sval
    integer                          :: nv, i
    character(len=140)               :: esc
    ! quotes inside the value are doubled (FITS standard 4.2.1); at most 66 characters
    esc = ''
    nv  = 0
    do i = 1, len_trim(val)
      if (nv >= 66) exit
      if (val(i:i) == "'") then
        if (nv >= 65) exit
        esc(nv+1:nv+2) = "''"
        nv = nv + 2
      else
        esc(nv+1:nv+1) = val(i:i)
        nv = nv + 1
      end if
    end do
    sval = "'"//esc(1:nv)//"'"
    if (nv < 8) sval = "'"//esc(1:nv)//repeat(' ', 8-nv)//"'"
    line = ''
    line(1:8)  = key
    line(9:10) = '= '
    line(11:)  = trim(sval)
    if (len_trim(line) < 31) then
      line(31:33) = ' / '
      line(34:80) = comment
    else if (len_trim(line) < 77) then
      line(len_trim(line)+1:) = ' / '//comment
    end if
    call add_card(hdr, line)
  end subroutine fits_card_str


  subroutine fits_card_comment(hdr, text)
    type(fits_header), intent(inout) :: hdr
    character(len=*), intent(in)     :: text
    character(len=80)                :: line
    line = 'COMMENT '//text
    call add_card(hdr, line)
  end subroutine fits_card_comment


  subroutine fits_write_cube(fname, dat, hdr_user)

    character(len=*), intent(in)              :: fname
    real(r2), dimension(:,:,:), intent(in)    :: dat           ! (nx, ny, nz)
    type(fits_header), intent(in)             :: hdr_user

    type(fits_header)                    :: hdr
    integer                              :: unit, i, nx, ny, nz, naxis
    integer(i8)                          :: nbytes, npad, k, ntot
    character(len=2880)                  :: block
    integer(kind=ib1), allocatable, dimension(:) :: bytes
    real(r2), allocatable, dimension(:)  :: flat
    logical                              :: little_endian
    integer(kind=ib1), dimension(8) :: b8

    nx = size(dat, 1)
    ny = size(dat, 2)
    nz = size(dat, 3)
    naxis = 3
    if (nz == 1) naxis = 2

    ! mandatory keywords first
    call fits_header_init(hdr)
    call add_card(hdr, card_logical('SIMPLE', .true., 'conforms to FITS standard'))
    call fits_card_int(hdr, 'BITPIX', -64, 'IEEE double precision')
    call fits_card_int(hdr, 'NAXIS', naxis, 'number of axes')
    call fits_card_int(hdr, 'NAXIS1', nx, 'x (detector xd)')
    call fits_card_int(hdr, 'NAXIS2', ny, 'y (detector yd)')
    if (naxis == 3) call fits_card_int(hdr, 'NAXIS3', nz, 'wavelength index or z')
    call add_card(hdr, card_logical('EXTEND', .true., 'extensions may be present'))
    do i = 1, hdr_user%n
      call add_card(hdr, hdr_user%card(i))
    end do
    call add_card(hdr, 'END')

    open(newunit=unit, file=fname, access='stream', form='unformatted', status='replace', &
         action='write')

    ! header blocks
    block = ''
    k = 0
    do i = 1, hdr%n
      block(k*80+1:k*80+80) = hdr%card(i)
      k = k + 1
      if (k == 36) then
        write(unit) block
        block = ''
        k = 0
      end if
    end do
    if (k > 0) write(unit) block

    ! data: column-major order of Fortran = FITS order (x fastest)
    ntot = int(nx, i8) * int(ny, i8) * int(nz, i8)
    allocate(flat(ntot))
    flat = reshape(dat, [ntot])

    little_endian = (transfer(1_ib4, b8(1)) == 1)

    allocate(bytes(8_i8 * ntot))
    do k = 1, ntot
      b8 = transfer(flat(k), b8)
      if (little_endian) b8 = b8(8:1:-1)
      bytes(8_i8*(k-1)+1 : 8_i8*k) = b8
    end do
    write(unit) bytes

    nbytes = 8_i8 * ntot
    npad   = mod(2880_i8 - mod(nbytes, 2880_i8), 2880_i8)
    if (npad > 0) then
      deallocate(bytes)
      allocate(bytes(npad))
      bytes = 0
      write(unit) bytes
    end if

    close(unit)
    deallocate(flat, bytes)

  end subroutine fits_write_cube


  function card_logical(key, val, comment) result(line)
    character(len=*), intent(in) :: key, comment
    logical, intent(in)          :: val
    character(len=80)            :: line
    line = ''
    line(1:8)  = key
    line(9:10) = '= '
    if (val) then
      line(30:30) = 'T'
    else
      line(30:30) = 'F'
    end if
    line(31:33) = ' / '
    line(34:80) = comment
  end function card_logical

end module fits_io
