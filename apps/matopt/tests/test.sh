#!/usr/bin/bash

{ cd  ./test1/ ; }
{ echo Test 1: ; }
{ python Monometallic_Nanocluster_Design.py ; }
{ cd .. ; }

{ cd  ./test2/ ; }
{ echo Test 2: ; }
{ python Bimetallic_Nanocluster_Design.py ; }
{ cd .. ; }

{ cd  ./test3/ ; }
{ echo Test 3: ; }
{ python Surface_Design.py ; }
{ cd .. ; }

{ cd  ./test4/ ; }
{ echo Test 4: ; }
{ python Bifunctional_Surface_Design.py ; }
{ cd .. ; }

{ cd  ./test5/ ; }
{ echo Test 5: ; }
{ python Metal_Oxide_Bulk_Design.py ; }
{ cd .. ; }