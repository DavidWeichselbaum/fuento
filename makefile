CC = g++
CPPFLAGS = -std=c++0x
LIBS = -lcurl -lboost_regex -lboost_filesystem
PREFIX = /usr/local
BUILDDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
INSTALLDIR = /.fuento

ifeq ($(OPTIMIZE),no)
  CPPFLAGS += 
else
  CPPFLAGS += -O2 -w
endif


fuento:fuento.cpp
	$(CC) fuento.cpp $(LIBS) -o fuento $(CPPFLAGS)

.PHONY: install
install: fuento
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp -p $< $(DESTDIR)$(PREFIX)/bin/fuento
	#   This program is free software: you can redistribute it and/or modify
	#   it under the terms of the GNU General Public License as published by
	#   the Free Software Foundation, either version 3 of the License, or
	#   (at your option) any later version.
	#   This program is distributed in the hope that it will be useful,
	#   but WITHOUT ANY WARRANTY; without even the implied warranty of
	#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	#   GNU General Public License for more details.
	#   You should have received a copy of the GNU General Public License
	#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

.PHONY: uninstall
uninstall:
	rm  $(DESTDIR)$(PREFIX)/bin/fuento
	rm -r $(HOME)/.fuento
