@ECHO OFF
REM ------------------------------------------------------------------------------
REM This file is part of the OpenStructure project <www.openstructure.org>
REM
REM Copyright (C) 2008-2020 by the OpenStructure authors
REM
REM This library is free software; you can redistribute it and/or modify it under
REM the terms of the GNU Lesser General Public License as published by the Free
REM Software Foundation; either version 3.0 of the License, or (at your option)
REM any later version.
REM This library is distributed in the hope that it will be useful, but WITHOUT
REM ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
REM FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
REM details.
REM
REM You should have received a copy of the GNU Lesser General Public License
REM along with this library; if not, write to the Free Software Foundation, Inc.,
REM 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
REM ------------------------------------------------------------------------------
REM Windows startup script for a protein-centric user interface
REM Author: Juergen Haas

call "%~dp0\bin\dng.bat" %*
