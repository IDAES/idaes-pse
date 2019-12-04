#!/usr/bin/bash

{ cd  ./test1/ ; }
{ echo Test 1: ; }
{ python3 Monometallic_Nanocluster_Design.py ; }
{ cd .. ; }

{ cd  ./test2/ ; }
{ echo Test 2: ; }
{ python3 Bimetallic_Nanocluster_Design.py ; }
{ cd .. ; }

{ cd  ./test3/ ; }
{ echo Test 3: ; }
{ python3 Surface_Design.py ; }
{ cd .. ; }

{ cd  ./test4/ ; }
{ echo Test 4: ; }
{ python3 Bifunctional_Surface_Design.py ; }
{ cd .. ; }

{ cd  ./test5/ ; }
{ echo Test 5: ; }
{ python3 Metal_Oxide_Bulk_Design.py ; }
{ cd .. ; }