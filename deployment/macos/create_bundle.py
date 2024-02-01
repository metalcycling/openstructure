INFO_PLIST="""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple Computer//DTD PLIST 1.0//EN" 
    "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleDevelopmentRegion</key>
  <string>English</string>
  <key>CFBundleExecutable</key>
  <string>dng</string>
  <key>CFBundleName</key>
  <string>DNG</string>  
  <key>CFBundleDisplayName</key>
  <string>DNG</string>  
  <key>CFBundleHelpBookFolder</key>
  <string>OpenStructure Manual</string>
  <key>CFBundleHelpBookName</key>
  <string>_BUNDLE_ Help</string>
  <key>CFBundleIconFile</key>
  <string>icon.icns</string>
  <key>CFBundleIdentifier</key>
  <string>org.openstructure.OpenStructure</string>
  <key>CFBundleInfoDictionaryVersion</key>
  <string>6.0</string>
  <key>CFBundlePackageType</key>
  <string>APPL</string>
  <key>CFBundleSignature</key>
  <string>????</string>
  <key>CFBundleVersion</key>
  <string></string>
  <key>NSMainNibFile</key>
  <string>MainMenu</string>
  <key>NSPrincipalClass</key>
  <string>NSApplication</string>
</dict>
</plist>"""

import os
from shutil import copy
import stat
def _WritePkgInfo(bundle):
    pkg_info=open(os.path.join(bundle+'.app', 'Contents/PkgInfo'), 'w+')
    pkg_info.write('APPL????')
    pkg_info.close()

def _WriteInfoPList(bundle):
    info_plist=open(os.path.join(bundle+'.app', 'Contents/Info.plist'), 'w+')
    info_plist.write(INFO_PLIST)
    info_plist.close()

def _WriteScript(bundle):
  script="""#!/bin/sh
/usr/openstructure/bin/dng
  """
  bin_path=os.path.join('%s.app' % bundle, 'Contents/MacOS/dng')
  bin=open(bin_path, 'w+')
  bin.write(script)
  bin.close()
  os.system('chmod a+x %s' % bin_path)

def _CreateBundleSkeleton(bundle):
  bin_path=os.path.join('%s.app' % bundle, 'Contents/MacOS')
  if not os.path.exists(bin_path):
    os.makedirs(bin_path)
  _WritePkgInfo(bundle)
  _WriteInfoPList(bundle)
  _WriteScript(bundle)
  res_dir='%s.app/Contents/Resources' % bundle
  if not os.path.exists(res_dir):
    os.makedirs(res_dir)  
  copy('../../graphics/icon.icns', res_dir)

def CreateBundle(bundle_name):
  _CreateBundleSkeleton(bundle_name)

CreateBundle('DNG')