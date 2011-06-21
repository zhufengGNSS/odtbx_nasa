
File:            aareadme.txt
Last Revised:    2010-03-03 (NJB)
Creation Date:   2006-08-04
Author:          Nat Bachman (JPL/NAIF)



Contents of path pub/naif/generic_kernels/pck
=============================================


Binary Earth PCKs
-----------------

   The following kernels provide Earth orientation data taking
   into account 

      Precession
      Nutation
      Nutation corrections
      UT1-TAI
      Polar motion

   See the comment areas of these kernels (using the SPICE utility
   COMMNT) for additional documentation.


   High accuracy EOP-based kernel, updated at least twice per week:

      earth_000101_yymmdd_yymmdd.bpc

   The first two dates in the file name are the file's coverage start
   and stop times. The third date in the file name is the epoch of the
   last datum in the source EOP file:  earth orientation from this date
   forward is predicted.


   Low accuracy, long term predict kernel. The extended predict region of this
   kernel---the time interval following the end of the predict region of the
   input EOP file---does not estimate UT1-TAI or polar motion. The dates in the
   file name are the file's coverage start and stop times:
 
      earth_070425_370426_predict.bpc 


   High accuracy, historical kernel.  The dates in the file name are the
   file's coverage start and stop times:

      earth_720101_070426.bpc


Binary Lunar PCK
----------------

   Latest version
   --------------

   DE-421 based kernel, providing orientation of Lunar Principal Axis (PA)
   reference frame:

      moon_pa_de421_1900-2050.bpc

   The data source for this kernel was the NAVIO file

      de421.nio

   The approximate time coverage is indicated by the dates in the
   file name.

   This kernel is to be used with the lunar frame kernel 

      moon_080317.tf
 
   This frame kernel is located in the path

      pub/naif/generic_kernels/fk/satellites

   See the cited frame kernel for a detailed discussion of lunar reference
   frames and associated SPICE kernels.


   Previous version
   ----------------

   DE-418 based kernel, providing orientation of Lunar Principal Axis (PA)
   reference frame:

      moon_pa_de418_1950-2050.bpc

   The data source for this kernel was the NAVIO file

      de418.nio

   The approximate time coverage is indicated by the dates in the
   file name.

   This kernel is to be used with the lunar frame kernel 

      moon_071218.tf
 
   This frame kernel is located in the path

      pub/naif/generic_kernels/fk/satellites

   See the cited frame kernel for a detailed discussion of lunar reference
   frames and associated SPICE kernels.
 


Generic text PCK
----------------

   This kernel contains orientation data for planets, natural satellites,
   the Sun, and selected asteroids:

      pck00009.tpc

   Most data are from article

     Seidelmann, P.K., Archinal, B.A., A'Hearn, M.F., 
     Conrad, A., Consolmagno, G.J., Hestroffer, D.,
     Hilton, J.L., Krasinsky, G.A., Neumann, G.,
     Oberst, J., Stooke, P., Tedesco, E.F., Tholen, D.J., 
     and Thomas, P.C. "Report of the IAU/IAG Working Group 
     on cartographic coordinates and rotational elements: 2006."

   This kernel is a text file; see the documentation in the file for
   details concerning the file's contents.


NORAD two-line element kernel
-----------------------------

   This kernel supports creation of SPK files containing 
   type 10 (NORAD two-line element) segments:
 
      geophysical.ker

   See the MKSPK User's Guide for further information.


GM/mass kernels
---------------
  
   The file 

      de-403-masses.tpc

   gives "GM" (gravitational constant times mass) values for the
   sun, planets, and planetary system barycenters.  These values
   are based on DE-403.

   The file 

      Gravity.tpc

   contains data of uncertain provenance.  Use at your own risk.
   NAIFers:  consider deleting this file?


EARTH_FIXED Alias Kernel
------------------------

   This is a kernel used to map the reference frame alias 
   EARTH_FIXED to a specific reference frame:

      earth_fixed.tf

   Using the alias EARTH_FIXED in code allows applications to 
   select the earth-fixed frame to use at run time via the 
   frame specification in this kernel.  The default association
   is with the frame 

      IAU_EARTH

   The file contains instructions for changing the association to
   ITRF93.



Getting help
============

   SPICE tutorials are available on the NAIF web site:

      http://naif.jpl.nasa.gov/naif/

   The following tutorials may be of particular interest:

      Intro to Kernel Files
      PCK (Planetary cartographic and physical constants)
      FK (Reference frame specifications)

   Also see the following Required Reading documents provided
   with the SPICE Toolkit:

      Frames      
      Kernel
      PCK

   Many SPICE kernels contain internal documentation. Text kernels may
   be viewed with a text editor or web browser. Comments in binary
   kernels may be viewed using the SPICE utility COMMNT. Coverage of
   binary PCKs may be viewed using the SPICE utility SPACIT. See the
   COMMNT and SPACIT User's Guides (provided with the SPICE Toolkit)
   for details.

