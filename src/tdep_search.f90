!____________________________________________________________________________________________________
!
! Module tdep_search: files of temperature-dependent optical data
!____________________________________________________________________________________________________
!
! For a component with temperature-dependent optical constants the input gives a base name, e.g.
!     optical_data = 'kim2019/pure_ices/composite_cice'
! and DMS uses every file of the optical-data directory whose name is
!     <base>_<T>K   or   <base>_<T>K.<extension>
! with T an integer or decimal number (leading zeros allowed: _10K, _010K, _12.5K) and any
! extension (also none). The temperature is read from the name, and the sets are sorted by it.
! Two files with the same temperature (e.g. _010K.txt and _10K.dat) are an error. A lower-case
! 'k' is accepted as well (e.g. crys_100k.tab).
!
! Fortran has no standard routine for listing a directory; the directory is read with the POSIX
! routines opendir/readdir through a small C function (dir_listing.c, ISO_C_BINDING).
!____________________________________________________________________________________________________

module tdep_search

  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
  use datatype

  implicit none

  private
  public :: find_tdep_files, parse_tdep_name

  integer, parameter :: lname = 512

  interface
    function c_dir_open(path) bind(C, name="dms_dir_open") result(ierr)
      import :: c_int, c_char
      character(kind=c_char), dimension(*), intent(in) :: path
      integer(c_int) :: ierr
    end function c_dir_open
    function c_dir_next(name, maxlen) bind(C, name="dms_dir_next") result(n)
      import :: c_int, c_char
      character(kind=c_char), dimension(*), intent(inout) :: name
      integer(c_int), value :: maxlen
      integer(c_int) :: n
    end function c_dir_next
    subroutine c_dir_close() bind(C, name="dms_dir_close")
    end subroutine c_dir_close
  end interface

contains

  !-------------------------------------------------------------------------------------------------
  ! Is name = prefix // "_" // <T> // "K" [// "." // <extension>] ? If so, T [K].
  !-------------------------------------------------------------------------------------------------
  subroutine parse_tdep_name(name, prefix, ok, temp)

    character(len=*), intent(in) :: name, prefix
    logical, intent(out)         :: ok
    real(r2), intent(out)        :: temp

    integer :: lp, ln, i, i0, ndig, ndot, ios
    character(len=1) :: ch

    ok = .false.;  temp = -1.0_r2
    lp = len_trim(prefix);  ln = len_trim(name)
    if (lp < 1 .or. ln < lp + 3) return
    if (name(1:lp) /= prefix(1:lp)) return
    if (name(lp+1:lp+1) /= "_") return

    ! number: digits with at most one decimal point, at least one digit before the point
    i0 = lp + 2;  i = i0;  ndig = 0;  ndot = 0
    do while (i <= ln)
      ch = name(i:i)
      if (ch >= "0" .and. ch <= "9") then
        ndig = ndig + 1
      else if (ch == "." .and. ndot == 0 .and. ndig > 0) then
        ! a point is part of the number only if a digit follows
        if (i + 1 > ln) exit
        if (name(i+1:i+1) < "0" .or. name(i+1:i+1) > "9") exit
        ndot = 1
      else
        exit
      end if
      i = i + 1
    end do
    if (ndig == 0 .or. i > ln) return
    if (name(i:i) /= "K" .and. name(i:i) /= "k") return
    ! after the K: end of the name or an extension
    if (i < ln) then
      if (name(i+1:i+1) /= ".") return
    end if
    read(name(i0:i-1), *, iostat=ios) temp
    if (ios /= 0 .or. .not. (temp > 0.0_r2)) then
      temp = -1.0_r2;  return
    end if
    ok = .true.

  end subroutine parse_tdep_name

  !-------------------------------------------------------------------------------------------------
  ! All data sets of a base name: n files (relative to data_dir), sorted by temperature.
  ! errmsg /= "" on an error (directory not readable, two files with the same temperature).
  !-------------------------------------------------------------------------------------------------
  subroutine find_tdep_files(base, data_dir, n, temps, files, errmsg)

    character(len=*), intent(in)  :: base, data_dir
    integer, intent(out)          :: n
    real(r2), allocatable, dimension(:), intent(out)          :: temps
    character(len=lname), allocatable, dimension(:), intent(out) :: files
    character(len=*), intent(out) :: errmsg

    character(len=lname) :: dirpart, prefix, dirpath, entry
    character(kind=c_char), dimension(lname + 1) :: cbuf
    real(r2), allocatable, dimension(:)          :: t_tmp
    character(len=lname), allocatable, dimension(:) :: f_tmp
    integer  :: islash, k, j, nl, cap
    logical  :: ok
    real(r2) :: temp, tsw
    character(len=lname) :: fsw

    errmsg = "";  n = 0
    islash = index(base, "/", back=.true.)
    if (islash > 0) then
      dirpart = base(1:islash);  prefix = base(islash+1:)
    else
      dirpart = "";  prefix = base
    end if
    dirpath = trim(data_dir)//trim(dirpart)
    if (len_trim(dirpath) == 0) dirpath = "."
    if (len_trim(prefix) == 0) then
      errmsg = "empty base name for temperature-dependent optical data: "//trim(base)
      allocate(temps(0), files(0));  return
    end if

    if (c_dir_open(to_c(dirpath)) /= 0) then
      errmsg = "cannot read the directory "//trim(dirpath)//" (temperature-dependent data "//trim(base)//")"
      allocate(temps(0), files(0));  return
    end if

    cap = 64
    allocate(t_tmp(cap), f_tmp(cap))
    do
      nl = c_dir_next(cbuf, int(lname + 1, c_int))
      if (nl < 0) exit
      if (nl > lname .or. nl == 0) cycle
      entry = ""
      do k = 1, nl
        entry(k:k) = cbuf(k)
      end do
      if (entry(1:1) == ".") cycle                         ! hidden files, . and ..
      call parse_tdep_name(entry, prefix, ok, temp)
      if (.not. ok) cycle
      if (n == cap) call grow()
      n = n + 1
      t_tmp(n) = temp
      f_tmp(n) = trim(dirpart)//trim(entry)
    end do
    call c_dir_close()

    ! sort by temperature (insertion sort; a few hundred files at most)
    do k = 2, n
      tsw = t_tmp(k);  fsw = f_tmp(k);  j = k - 1
      do while (j >= 1)
        if (t_tmp(j) <= tsw) exit
        t_tmp(j+1) = t_tmp(j);  f_tmp(j+1) = f_tmp(j);  j = j - 1
      end do
      t_tmp(j+1) = tsw;  f_tmp(j+1) = fsw
    end do

    do k = 2, n
      if (abs(t_tmp(k) - t_tmp(k-1)) <= 1.0e-9_r2 * t_tmp(k)) then
        errmsg = "two files with the same temperature for "//trim(base)//": "//trim(f_tmp(k-1))// &
                 " and "//trim(f_tmp(k))
        exit
      end if
    end do

    allocate(temps(n), files(n))
    temps(:) = t_tmp(1:n);  files(:) = f_tmp(1:n)

  contains

    subroutine grow()
      real(r2), allocatable, dimension(:) :: t2
      character(len=lname), allocatable, dimension(:) :: f2
      allocate(t2(2*cap), f2(2*cap))
      t2(1:cap) = t_tmp;  f2(1:cap) = f_tmp
      call move_alloc(t2, t_tmp);  call move_alloc(f2, f_tmp)
      cap = 2 * cap
    end subroutine grow

  end subroutine find_tdep_files

  function to_c(str) result(cstr)
    character(len=*), intent(in) :: str
    character(kind=c_char), dimension(len_trim(str) + 1) :: cstr
    integer :: k
    do k = 1, len_trim(str)
      cstr(k) = str(k:k)
    end do
    cstr(len_trim(str) + 1) = c_null_char
  end function to_c

end module tdep_search
