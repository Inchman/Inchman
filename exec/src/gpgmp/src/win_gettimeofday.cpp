/*
gettimeofday

Imported into GPGMP by Aidan Lane on Mon Apr 16, 2012.

From
http://suacommunity.com/dictionary/gettimeofday-entry.php

Purpose
The gettimeofday() function obtain the current time, expressed
as seconds and microseconds since the Epoch, and stores it in
the timeval structure given in the parameters.

Discussion
There is no direct analog of the gettimeofday() in Windows;
developers porting codes that use this function must provide
a conversion. The typical conversion utilizes the Microsoft
Windows GetSystemTimeAsFileTime function, but there are some
conversion issues. One, the Windows routine returns the number
of 100 nanosecond intervals, not the number of microseconds.
Two, the Windows function defines the starting Epoch as
January 1, 1601 where the Unix gettimeofday() function defines
the starting Epoch as January 1, 1970. Below is an example of
wrapping GetSystemTimeAsFileTime that deals with both of these
issues, allowing developers to provide a replacement
gettimeofday() without requiring any other changes.
*/

#include "win_gettimeofday.h"

#include <time.h>
#include <windows.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif


int gettimeofday(struct timeval *tv, struct timezone *tz)
{
    // Define a structure to receive the current Windows filetime
    FILETIME ft;

    // Initialize the present time to 0 and the timezone to UTC
    unsigned __int64 tmpres = 0;
    static int tzflag = 0;

    if (NULL != tv)
    {
        GetSystemTimeAsFileTime(&ft);

        // The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
        // intervals since Jan 1, 1601 in a structure. Copy the high bits to 
        // the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        // Convert to microseconds by dividing by 10
        tmpres /= 10;

        // The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
        // in seconds from Jan 1 1601.
        tmpres -= DELTA_EPOCH_IN_MICROSECS;

        // Finally change microseconds to seconds and place in the seconds value. 
        // The modulus picks up the microseconds.
        tv->tv_sec = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }

    if (NULL != tz)
    {
        if (!tzflag)
        {
            _tzset();
            tzflag++;
        }

        // Adjust for the timezone west of Greenwich
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }

    return 0;
}
