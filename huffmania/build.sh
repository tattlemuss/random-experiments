#!/bin/sh
rm a.out
cc main.c -g -Werror -lgcc
./a.out

