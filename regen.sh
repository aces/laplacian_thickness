#! /bin/bash

rm config.log
aclocal -I ~jason/share/aclocal/
autoheader
automake
autoconf

