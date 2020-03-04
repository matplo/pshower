/* dn.i */
%module dn
%include "typemaps.i"
%include "cpointer.i"
%pointer_class(int, intp);
%pointer_class(long int, longintp);
%pointer_class(double, doublep);

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);
%array_class(long int, longintArray);

%array_class(doublep, doubleParrayClass);

 %{
 /* Put header files here or function declarations like below */
 #include "Parton_Shower_Lib.hh"
 #include "jet_fragmentation.hh"
 %}

%include "Parton_Shower_Lib.hh"
%include "jet_fragmentation.hh"
