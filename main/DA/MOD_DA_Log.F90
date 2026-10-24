#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Log
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Log file for data assimilation.
!
!    The log is always written by the master process to a single file in the
!    run directory (DEF_DA_LOG_FILE, default "DA.log"); there is no switch to
!    disable it. Three fixed message classes are provided:
!
!      DA_log        log file only   -- progress and per-step diagnostics
!      DA_log_key    log + stdout   -- run configuration, important events
!      DA_log_error  log + stdout   -- fatal errors; echoed on the calling rank
!
!    Every line carries the current model date and time so that the log can be
!    used for debugging. When the model date has not been set yet (very early
!    initialization), the wall-clock date is used instead.
!
! AUTHOR:
!   Zhilong Fan, 10/2026
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist,    only: DEF_DA_LOG_FILE
   USE MOD_TimeManager, only: julian2monthday
   IMPLICIT NONE
   SAVE

   PRIVATE
   PUBLIC :: init_DA_log
   PUBLIC :: end_DA_log
   PUBLIC :: DA_log_set_date
   PUBLIC :: DA_log
   PUBLIC :: DA_log_key
   PUBLIC :: DA_log_error

   integer :: unit_log = -1          ! log file unit
   logical :: log_open = .false.     ! true after the log file has been opened
   character(len=19) :: stamp = ''   ! current model date: 'YYYY-MM-DD HH:MM:SS'

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_log ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=256) :: logfile
   character(len=512) :: msg
   character(len=8)  :: cdate
   character(len=10) :: ctime
   integer :: ierr

!-----------------------------------------------------------------------------
      IF (.not. p_is_master) RETURN
      IF (log_open) RETURN

      logfile = trim(DEF_DA_LOG_FILE)
      IF (len_trim(logfile) == 0) logfile = 'DA.log'

      ! Append so that a restarted run keeps the previous segment in the file.
      OPEN (NEWUNIT=unit_log, FILE=trim(logfile), STATUS='unknown', &
            POSITION='append', ACTION='write', IOSTAT=ierr)
      IF (ierr /= 0) THEN
         WRITE (msg,'(A,A,A,I0)') '[CoLM-DA] ERROR: cannot open DA log file: ', &
            trim(logfile), ', iostat = ', ierr
         CALL CoLM_stop (msg)
      ENDIF
      log_open = .true.

      CALL DATE_AND_TIME (cdate, ctime)
      IF (len_trim(stamp) == 0) THEN
         stamp = cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8)//' '// &
                 ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
      ENDIF

      WRITE (unit_log,'(A)') ''
      WRITE (unit_log,'(A)') '[CoLM-DA] ===== DA log opened at '// &
         cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8)//' '// &
         ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)//' ====='
      FLUSH (unit_log)

   END SUBROUTINE init_DA_log

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_log ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (.not. p_is_master) RETURN
      IF (.not. log_open) RETURN

      WRITE (unit_log,'(A)') '[CoLM-DA] ===== DA log closed ====='
      CLOSE (unit_log)
      unit_log = -1
      log_open = .false.

   END SUBROUTINE end_DA_log

!-----------------------------------------------------------------------------

   SUBROUTINE DA_log_set_date (idate)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   integer, intent(in) :: idate(3)   ! calendar (year, julian day, seconds)

   integer :: month, mday, hour, minute, second, jday

!-----------------------------------------------------------------------------
      IF (.not. p_is_master) RETURN

      jday = idate(2)
      IF (jday < 1) jday = 1
      CALL julian2monthday (idate(1), jday, month, mday)

      hour   = idate(3) / 3600
      minute = MOD (idate(3), 3600) / 60
      second = MOD (idate(3), 60)
      WRITE (stamp,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         idate(1), month, mday, hour, minute, second

   END SUBROUTINE DA_log_set_date

!-----------------------------------------------------------------------------

   SUBROUTINE DA_log (msg, echo)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=*),           intent(in) :: msg
   logical,          optional, intent(in) :: echo   ! also write to stdout

   character(len=1024) :: line

!-----------------------------------------------------------------------------
      CALL init_DA_log ()

      line = '[CoLM-DA] '//stamp//' | '//trim(msg)
      IF (log_open) THEN
         WRITE (unit_log,'(A)') trim(line)
         FLUSH (unit_log)
      ENDIF
      IF (PRESENT(echo)) THEN
         IF (echo) WRITE (*,'(A)') trim(line)
      ENDIF

   END SUBROUTINE DA_log

!-----------------------------------------------------------------------------

   SUBROUTINE DA_log_key (msg)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=*), intent(in) :: msg

!-----------------------------------------------------------------------------
      ! Key messages are echoed once, by the master.
      IF (p_is_master) CALL DA_log (msg, echo = .true.)

   END SUBROUTINE DA_log_key

!-----------------------------------------------------------------------------

   SUBROUTINE DA_log_error (msg)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=*), intent(in) :: msg

   character(len=1024) :: line

!-----------------------------------------------------------------------------
      ! Errors may be detected on a worker; keep them visible on every rank
      ! that reports one, and record them in the log on the master.
      CALL init_DA_log ()

      IF (p_is_master) THEN
         line = '[CoLM-DA] '//stamp//' | '//trim(msg)
         IF (log_open) THEN
            WRITE (unit_log,'(A)') trim(line)
            FLUSH (unit_log)
         ENDIF
         WRITE (*,'(A)') trim(line)
      ELSE
         WRITE (*,'(A,I0,A)') '[CoLM-DA][rank ', p_iam_glb, '] '//trim(msg)
      ENDIF

   END SUBROUTINE DA_log_error

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Log
#endif
