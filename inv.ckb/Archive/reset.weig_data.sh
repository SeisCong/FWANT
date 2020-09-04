#!/bin/bash

Gd_list='inv_Gd_list'
Gd_list_bak="inv_Gd_list_`date +%Y%m%d`"

cp -f Gd_list_bak Gd_list_bak

awk '{print $1 " " $2 " " $3 " " 1 " " $4 " " $5 " " $6 " " $7}' $Gd_list_bak > $Gd_list
