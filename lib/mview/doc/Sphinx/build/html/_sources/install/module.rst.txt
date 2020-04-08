Installation: Perl module
=========================

Perl module
^^^^^^^^^^^

1. Save the archive then uncompress and extract it (Linux, macOS, UNIX)::

        tar xvzf mview-VERSION.tar.gz

   or (Windows, using an archiver like 7-Zip, as here)::
   
        7z x mview-VERSION.zip
        
   This creates a sub-folder ``mview-VERSION`` containing all the files.
   
2. Change to this folder.

You can now use one of the following sets of instructions to do the install:

3. Run::

        perl Makefile.PL
        make install
        
   which attempts to install into the Perl distribution.

3. Or run something like::

        perl Makefile.PL INSTALL_BASE=/usr/local
        make install

   which attempts to install under the given folder. In this example you need
   write access to ``/usr/local`` and users will need ``/usr/local/bin`` on
   their ``PATH``.
   
3. Or, if you have a `local::lib <https://metacpan.org/pod/local::lib>`_
   setup, you can install MView there::

        perl Makefile.PL $PERL_MM_OPT
        make install

4. Finally, the unpacked archive can be deleted since the important components
   have been copied elsewhere.
