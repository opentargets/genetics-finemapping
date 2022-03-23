#!/usr/bin/env bash

plink --bfile $1 --out $2 --make-bed --allow-extra-chr --extract range <(echo -e "$3")
