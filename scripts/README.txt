1. Setup the rpm build environment (See guides on the web)
2. Checkout openstructure source in a fodler called 'openstructure-1.1d'
3. Copy the compounds.chemlib into the folder
4. Remove '.git' and '.gitignore' from the folder 
5. COmpress the source with 'tar -cJvf openstructuere.tar.xz openstructure-1.1d'
6. Put the zipped source in '~/rpmbuild/SOURCES'
7. Put the desktop file in '~/rpmbuild/SOURCES'
8. Put the spec file in '~/rpmbuild/SPECS'
9. Go the the 'SPECS' directory and build with 'rpmbuild -ba openstructure.spec' 

