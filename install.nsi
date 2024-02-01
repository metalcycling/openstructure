; NSIS script for DNG 
; Juergen Haas <juergen.haas@unibas.ch>
; Nov 2009
;=========================
  
; The name of the program
!define PRODUCT_NAME "openstructure"
!define PRODUCT_NAME_SMALL "dng"
!define PRODUCT_VERSION "0.9.0"
Name "${PRODUCT_NAME} v${PRODUCT_VERSION}"

; Installer filename
OutFile "dng_win32.exe"

;=========================

    LicenseText "If you accept the terms of the agreement, click I Agree to continue.\ 
	You must accept the agreement to install ${PRODUCT_NAME_SMALL} v${PRODUCT_VERSION}" 
    LicenseData "LICENSE.txt"
    Page license 
	
	DirText "Setup will install ${PRODUCT_NAME_SMALL} v${PRODUCT_VERSION} in the following folder.\ 
    To install in a different folder, click, browse and select another folder. Please do not choose\
	a path containing blank spaces as it will make ${PRODUCT_NAME} not work properly. Click Install to start the installation"
	InstallDir "C:\${PRODUCT_NAME_SMALL}"
    Page directory 
	
	;Page components
	
	ShowInstDetails Show
	ShowUninstDetails Show	
    Page instfiles
	

;=========================
;Sections
;=========================


Section "" ;No components page, name is not important

SetOutPath $INSTDIR
File "README"
Rename $INSTDIR\README $INSTDIR\README.txt

SetOutPath "$INSTDIR\bin"
File ".\stage\bin\*" 

SetOutPath "$INSTDIR\examples"
File ".\examples\entity\test_query.py"

WriteUninstaller "$INSTDIR\Uninstall.exe"

SectionEnd


 ; The uninstall section
Section "Uninstall"

;Exec "$INSTDIR\bin\unassoc.bat"

DeleteRegValue HKEY_CLASSES_ROOT ".pdb" ""
DeleteRegValue HKEY_CLASSES_ROOT ".pdb\DefaultIcon" ""
DeleteRegValue HKEY_CLASSES_ROOT ".pdb\Shell" ""
DeleteRegValue HKEY_CLASSES_ROOT ".pdb\Shell\Open" ""
DeleteRegValue HKEY_CLASSES_ROOT ".pdb\Shell\Open\Command" ""


DeleteRegKey HKCR ".pdb"

Delete "$INSTDIR\bin\*.*"
Delete "$INSTDIR\src\*.*"
Delete "$INSTDIR\examples\*.*"

RMDir "$INSTDIR"

RMDir "$INSTDIR\examples\entity"
RMDir "$INSTDIR\examples"

Delete "$SMPROGRAMS\${PRODUCT_NAME_SMALL}\dng.lnk"
Delete "$SMPROGRAMS\${PRODUCT_NAME_SMALL}\dng_cl.lnk"
Delete "$SMPROGRAMS\${PRODUCT_NAME_SMALL}\uninstall.lnk"
RMDir "$SMPROGRAMS\${PRODUCT_NAME_SMALL}"

Delete "$INSTDIR\README.txt"
Delete "$INSTDIR\README"


RMDir $INSTDIR

SectionEnd

Section "Shortcut in Start Menu"

SetOutPath "$INSTDIR\win32"

CreateDirectory "$SMPROGRAMS\${PRODUCT_NAME_SMALL}"
CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME_SMALL}\dng.lnk" "$INSTDIR\bin\dng.exe" "" "$INSTDIR\dng.ico"
CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME_SMALL}\dng_cl.lnk" "$INSTDIR\bin\dng_cl.exe"
CreateShortCut "$DESKTOP\dng.lnk" "$INSTDIR\bin\dng.exe" "" "$INSTDIR\win32\dng.ico"


CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME_SMALL}\uninstall.lnk" "$INSTDIR\Uninstall.exe" 
SectionEnd

Section "pdb file association"

;Exec "$INSTDIR\assoc.bat"

WriteRegStr HKEY_CLASSES_ROOT ".pdb" "" ""
WriteRegStr HKEY_CLASSES_ROOT ".pdb\DefaultIcon" "" "$INSTDIR\dng.ico"
WriteRegStr HKEY_CLASSES_ROOT ".pdb\Shell" "" ""
WriteRegStr HKEY_CLASSES_ROOT ".pdb\Shell\Open" "" ""
WriteRegStr HKEY_CLASSES_ROOT ".pdb\Shell\Open\Command" "" "$INSTDIR\dng_load.py $\"%1$\""

SectionEnd

;=========================
;Functions
;=========================

#Function .onInit
#FunctionEnd
Function un.onUninstSuccess

HideWindow
MessageBox MB_ICONINFORMATION|MB_OK "$(^Name) was successfully removed from your computer."

FunctionEnd

Function un.onInit

MessageBox MB_ICONQUESTION|MB_YESNO|MB_DEFBUTTON2 "Are you sure you want to completely remove $(^Name) and all of its components?" IDYES +2
Abort

FunctionEnd

 
