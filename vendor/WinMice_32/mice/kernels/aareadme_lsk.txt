

    LEAPSECONDS KERNEL VERSION NAIF0009

The file naif0009.tls is a unix-style text file. It is suitable for use
on all unix boxes, including PC/linux and MAC/OSX machines.

For PCs running Windows, use naif0009.tls.pc.

Use of one of these files is required for all SPICE computations
dealing with times on or after Jan 01, 2009 00:00:00 if you want
accurate conversions between UTC and Ephemeris Time (a.k.a. TDB).
Failure to use naif0009 for times after this epoch will result
in a one second time conversion error.

You may begin using naif0009 in place of naif0008 right now; there
is no need to wait until Jan 01, 2009. Time conversions for epochs
prior to Jan 01, 2009 will be the same using either naif0008 or
naif0009.


