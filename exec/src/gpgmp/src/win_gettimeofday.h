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

#ifndef __gpgmp_gettimeofday_h__
#define __gpgmp_gettimeofday_h__


struct timezone
{
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz);


#endif // !__gpgmp_gettimeofday_h__
